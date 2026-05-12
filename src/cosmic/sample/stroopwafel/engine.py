"""AdaptiveSampler: the main STROOPWAFEL engine.

Orchestrates the explore -> adapt -> refine -> weight calculation pipeline
using vectorized operations throughout. No Location objects are created.
"""
import os
import numpy as np
from scipy.stats import multivariate_normal

from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

from .mixture_model import GaussianMixture
from .result import STROOPWAFELResult
from .constants import MIN_ACTIVE_FRACTION


class AdaptiveSampler:
    """Adaptive importance sampler for binary population synthesis with COSMIC.

    Parameters
    ----------
    parameter_space : `ParameterSpace`
        Defines the sampling dimensions and their distributions.
    total_systems : `int`
        Total number of systems to simulate across all phases.
    batch_size : `int`
        Number of systems evolved per COSMIC call.
    BSEDict : `dict`
        COSMIC binary stellar evolution parameters.
    compute_derived : `callable`
        Function with signature
        ``(samples_physical, param_names) -> dict`` that computes
        derived quantities (e.g., ``'mass_2'``, ``'separation'``).
    reject_systems : `callable`
        Function with signature
        ``(samples_physical, derived, param_names) -> bool_mask``
        returning True for physically unacceptable systems.
    is_interesting : `callable`
        Function with signature ``(bpp) -> (n_hits, hit_bin_nums)``
        identifying systems of interest from COSMIC output.
    output_path : `str`, optional
        Directory for output files, by default ``'output'``
    nproc : `int`, optional
        Number of CPU cores for COSMIC, by default 1
    kappa : `float`, optional
        Gaussian width scaling factor, by default 1.0
    n_generations : `int`, optional
        Number of refinement generations, by default 1
    mc_only : `bool`, optional
        If True, only run exploration (standard Monte Carlo), by default
        False
    seed : `int` or None, optional
        Random seed for reproducibility, by default None
    """

    def __init__(self, parameter_space, total_systems, batch_size, BSEDict,
                 compute_derived, reject_systems, is_interesting,
                 output_path='output', nproc=1, kappa=1.0,
                 n_generations=1, mc_only=False, seed=None):
        self.param_space = parameter_space
        self.total_systems = total_systems
        self.batch_size = batch_size
        self.bse_dict = BSEDict
        self.compute_derived_fn = compute_derived
        self.reject_fn = reject_systems
        self.is_interesting_fn = is_interesting
        self.output_path = output_path
        self.nproc = nproc
        self.kappa = kappa
        self.n_generations = n_generations
        self.mc_only = mc_only
        self.rng = np.random.default_rng(seed)

        # State
        self.num_explored = 0
        self.num_hits = 0
        self.fraction_explored = 1.0
        self.finished = 0
        self.prior_fraction_rejected = 0.0
        self.mixture = None

        # Accumulators for all samples (in sampling space)
        self._all_samples = []
        self._all_is_hit = []
        self._all_generation = []
        self._all_gaussian_idx = []

    def run(self):
        """Run the full STROOPWAFEL pipeline.

        Returns
        -------
        `STROOPWAFELResult`
            Container holding all samples, weights, and metadata.
        """
        os.makedirs(self.output_path, exist_ok=True)

        self._explore()

        if not self.mc_only and self.num_hits > 0:
            self._adapt()
            self._refine()

        result = self._compute_weights()
        return result

    # ------------------------------------------------------------------
    # Exploration
    # ------------------------------------------------------------------
    def _explore(self):
        print("Exploration phase started")

        if not self.mc_only:
            self.prior_fraction_rejected = self._estimate_prior_rejection_rate()
            print(f"  Prior rejection rate: {self.prior_fraction_rejected:.4f}")

        while self._should_continue_exploring():
            n_oversample = int(2 * np.ceil(
                self.batch_size / max(1.0 - self.prior_fraction_rejected, MIN_ACTIVE_FRACTION)
            ))

            # Sample from prior
            samples, mask = self.param_space.sample(n_oversample, rng=self.rng)

            # Transform to physical space for rejection checking
            samples_phys = self.param_space.to_physical(samples)
            derived = self.compute_derived_fn(samples_phys, self.param_space.names)

            # Physical rejection
            phys_rejected = self.reject_fn(samples_phys, derived, self.param_space.names)

            # Combined mask: in bounds AND not physically rejected
            valid = mask & ~phys_rejected

            # Select up to batch_size valid samples
            valid_indices = np.where(valid)[0]
            self.rng.shuffle(valid_indices)
            selected = valid_indices[:self.batch_size]

            batch_samples = samples[selected]
            batch_samples_phys = samples_phys[selected]
            batch_derived = {k: v[selected] for k, v in derived.items()}

            # Evolve with COSMIC
            n_hits, hit_bin_nums, bpp, initC, kick_info = self._evolve_batch(
                batch_samples_phys, batch_derived
            )

            # Record which are hits
            is_hit = np.zeros(len(selected), dtype=bool)
            is_hit[hit_bin_nums] = True

            # Accumulate
            self._all_samples.append(batch_samples)
            self._all_is_hit.append(is_hit)
            self._all_generation.append(np.zeros(len(selected), dtype=int))
            self._all_gaussian_idx.append(np.full(len(selected), -1, dtype=int))

            self.num_hits += n_hits
            self.finished += len(selected)
            self.num_explored += len(selected)
            self._update_fraction_explored()

            self._print_progress()

        self.num_hits_exploratory = self.num_hits
        print(f"\nExploration done: {self.num_hits} hits / {self.num_explored} explored "
              f"(rate={self.num_hits/max(1, self.num_explored):.6f}, "
              f"f_expl={self.fraction_explored:.4f})")

    def _estimate_prior_rejection_rate(self, n_test=100000):
        """Estimate fraction of prior samples that are physically rejected.

        Parameters
        ----------
        n_test : `int`, optional
            Number of samples to use for the estimate, by default 100000

        Returns
        -------
        `float`
            Estimated fraction of samples rejected by bounds and physical
            criteria combined.
        """
        samples, mask = self.param_space.sample(n_test, rng=self.rng)
        rejected = n_test - np.sum(mask)

        valid_samples = samples[mask]
        if len(valid_samples) > 0:
            phys = self.param_space.to_physical(valid_samples)
            derived = self.compute_derived_fn(phys, self.param_space.names)
            phys_rejected = self.reject_fn(phys, derived, self.param_space.names)
            rejected += np.sum(phys_rejected)

        return rejected / n_test

    def _update_fraction_explored(self):
        if self.num_hits == 0 or self.num_explored == 0:
            return
        u = 1.0 / (self.fraction_explored * self.total_systems)
        r = self.num_hits / self.num_explored
        num = r * (np.sqrt(1.0 - r) - np.sqrt(u))
        den = np.sqrt(1.0 - r) * (np.sqrt(u * (1.0 - r)) + r)
        if den != 0:
            self.fraction_explored = 1 - num / den

    def _should_continue_exploring(self):
        if self.mc_only:
            return self.num_explored < self.total_systems
        return self.num_explored / self.total_systems < self.fraction_explored

    # ------------------------------------------------------------------
    # Adaptation
    # ------------------------------------------------------------------
    def _adapt(self):
        print("Adaptation phase started")

        # Gather all exploration hits in sampling space
        all_samples = np.vstack(self._all_samples)
        all_is_hit = np.concatenate(self._all_is_hit)
        hit_samples = all_samples[all_is_hit]

        average_density_one_dim = 1.0 / np.power(self.num_explored, 1.0 / self.param_space.ndim)

        self.mixture = GaussianMixture.from_hits(
            hit_samples, self.param_space, average_density_one_dim, kappa=self.kappa
        )

        print(f"  Created {self.mixture.n_components} Gaussian components")
        print("Adaptation phase finished")

    # ------------------------------------------------------------------
    # Refinement
    # ------------------------------------------------------------------
    def _refine(self):
        print("Refinement phase started")
        entropies = []

        for gen in range(self.n_generations):
            # Estimate rejection rate for current mixture
            dist_rejection_rate = self.mixture.compute_rejection_rate(
                self.param_space, self.compute_derived_fn, self.reject_fn,
                n_per_component=10000, rng=self.rng
            )

            n_per_gen = int((self.total_systems - self.num_explored) / self.n_generations)
            gen_samples_list = []
            gen_is_hit_list = []
            gen_finished = 0

            while gen_finished < n_per_gen and self.finished < self.total_systems:
                # Sample from mixture
                samples, mask, gauss_idx = self.mixture.sample(
                    self.batch_size, self.param_space,
                    consider_rejection=True, rng=self.rng
                )

                # Apply bounds and physical rejection
                valid_samples = samples[mask]
                valid_gauss_idx = gauss_idx[mask]

                if len(valid_samples) == 0:
                    continue

                phys = self.param_space.to_physical(valid_samples)
                derived = self.compute_derived_fn(phys, self.param_space.names)
                phys_rejected = self.reject_fn(phys, derived, self.param_space.names)

                keep = ~phys_rejected
                valid_samples = valid_samples[keep]
                valid_gauss_idx = valid_gauss_idx[keep]
                phys = phys[keep]
                derived = {k: v[keep] for k, v in derived.items()}

                # Trim to batch size with randomisation
                indices = np.arange(len(valid_samples))
                self.rng.shuffle(indices)
                n_take = min(len(valid_samples), self.batch_size)
                batch_samples = valid_samples[indices[:n_take]]
                batch_phys = phys[indices[:n_take]]
                batch_gauss_idx = valid_gauss_idx[indices[:n_take]]
                batch_derived = {k: v[indices[:n_take]] for k, v in derived.items()}

                # Evolve with COSMIC
                n_hits, hit_bin_nums, bpp, initC, kick_info = self._evolve_batch(
                    batch_phys, batch_derived
                )

                is_hit = np.zeros(n_take, dtype=bool)
                is_hit[hit_bin_nums] = True

                # Accumulate
                self._all_samples.append(batch_samples)
                self._all_is_hit.append(is_hit)
                self._all_generation.append(np.full(n_take, gen + 1, dtype=int))
                self._all_gaussian_idx.append(batch_gauss_idx)

                gen_samples_list.append(batch_samples)
                gen_is_hit_list.append(is_hit)

                self.num_hits += n_hits
                self.finished += n_take
                gen_finished += n_take
                self._print_progress()

            # EM update (if not the last generation)
            if gen < self.n_generations - 1 and len(gen_samples_list) > 0:
                gen_samples = np.vstack(gen_samples_list)
                gen_is_hit = np.concatenate(gen_is_hit_list)
                gen_priors = self.param_space.compute_prior(gen_samples)

                # Save current state in case we need to revert
                saved_mixture = GaussianMixture(
                    self.mixture.means.copy(),
                    self.mixture.covariances.copy(),
                    self.mixture.alphas.copy(),
                    self.mixture.rejection_rate
                )

                should_revert = self.mixture.update_em(
                    gen_samples, gen_is_hit, gen_priors,
                    self.prior_fraction_rejected,
                    tolerance=1e-10, entropies=entropies
                )

                if should_revert:
                    self.mixture = saved_mixture
                    print("  EM update reverted (insufficient entropy change)")

        n_refined = self.total_systems - self.num_explored
        if n_refined > 0:
            refine_hits = self.num_hits - self.num_hits_exploratory
            print(f"\nRefinement done: {refine_hits} hits / {n_refined} refined "
                  f"(rate={refine_hits/max(1, n_refined):.6f})")

    # ------------------------------------------------------------------
    # Weight calculation
    # ------------------------------------------------------------------
    def _compute_weights(self):
        """Compute importance sampling weights for all samples.

        Builds the final ``STROOPWAFELResult`` using the formula
        ``w(x) = π(x) / Q(x)`` where
        ``Q = f_e·π + (1 − f_e)·q`` is a mixture of the prior and the
        Gaussian proposal.

        Returns
        -------
        `STROOPWAFELResult`
            Container holding all samples, importance weights, and
            associated metadata.
        """
        all_samples = np.vstack(self._all_samples)
        all_is_hit = np.concatenate(self._all_is_hit)
        all_generation = np.concatenate(self._all_generation)
        all_gaussian_idx = np.concatenate(self._all_gaussian_idx)
        N = len(all_samples)

        # Prior probabilities
        pi_norm = 1.0 / max(1.0 - self.prior_fraction_rejected, MIN_ACTIVE_FRACTION)
        pi = self.param_space.compute_prior(all_samples) * pi_norm

        # Start with the exploration-phase contribution to denominator
        fraction_explored = self.num_explored / float(N)
        den = fraction_explored * pi

        # Add refinement-phase contributions from each generation's mixture
        if self.mixture is not None and not self.mc_only:
            q_norm = 1.0 / max(1.0 - self.mixture.rejection_rate, MIN_ACTIVE_FRACTION)
            # Evaluate the mixture PDF incrementally (memory-efficient)
            for k in range(self.mixture.n_components):
                xPDF_k = multivariate_normal.pdf(
                    all_samples, self.mixture.means[k], self.mixture.covariances[k],
                    allow_singular=True
                )
                den += (xPDF_k * self.mixture.alphas[k]
                        * (1 - fraction_explored) * q_norm) / self.n_generations

        weights = pi / den

        # Build result
        result = STROOPWAFELResult(
            samples=self.param_space.to_physical(all_samples),
            param_names=self.param_space.names,
            weights=weights,
            is_hit=all_is_hit,
            generation=all_generation,
            gaussian_idx=all_gaussian_idx,
            num_explored=self.num_explored,
            num_hits=self.num_hits,
            fraction_explored=self.fraction_explored,
        )

        print(f"\nTotal hits: {self.num_hits}")
        print(f"Sum of weights: {np.sum(weights):.4f}")
        print(f"Weighted hit rate: {result.hit_rate:.8f} +/- {result.hit_rate_uncertainty:.8f}")

        return result

    # ------------------------------------------------------------------
    # COSMIC interface
    # ------------------------------------------------------------------
    def _evolve_batch(self, samples_physical, derived):
        """Evolve a batch of binaries with COSMIC and identify hits.

        Parameters
        ----------
        samples_physical : `numpy.ndarray`
            (N, D) array of binary parameters in physical space.
        derived : `dict`
            Dictionary with keys ``'mass_2'``, ``'metallicity_1'``, etc.,
            each mapping to an (N,) array.

        Returns
        -------
        n_hits : `int`
            Number of systems classified as hits.
        hit_bin_nums : `numpy.ndarray`
            Integer array of 0-indexed positions within the batch that
            are hits.
        bpp : `pandas.DataFrame`
            COSMIC binary population parameters output.
        initC : `pandas.DataFrame`
            COSMIC initial conditions output.
        kick_info : `pandas.DataFrame`
            COSMIC natal kick information output.
        """
        n = len(samples_physical)
        idx = {name: i for i, name in enumerate(self.param_space.names)}

        batch_initial = InitialBinaryTable.InitialBinaries(
            m1=samples_physical[:, idx['mass_1']],
            m2=derived['mass_2'],
            porb=samples_physical[:, idx['porb']],
            ecc=samples_physical[:, idx['ecc']],
            tphysf=np.full(n, 13700.0),
            kstar1=np.full(n, 1),
            kstar2=np.full(n, 1),
            metallicity=derived['metallicity_1'],
        )

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=batch_initial,
            BSEDict=self.bse_dict,
            nproc=self.nproc,
        )

        # Apply user's hit identification
        n_hits, hit_bin_nums = self.is_interesting_fn(bpp)

        return n_hits, hit_bin_nums, bpp, initC, kick_info

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    def _print_progress(self):
        pct = 100 * self.finished / self.total_systems
        bar_len = 20
        filled = int(bar_len * self.finished // self.total_systems)
        bar = '|' * filled + '-' * (bar_len - filled)
        print(f'\r  progress |{bar}| {pct:.1f}% complete '
              f'({self.finished}/{self.total_systems})', end='\r')
