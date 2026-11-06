"""Vectorized Gaussian Mixture Model for adaptive importance sampling.

Stores all mixture components as numpy arrays and provides vectorized
sampling, PDF evaluation, and expectation-maximization (EM) updates.
"""
import numpy as np
from scipy.stats import multivariate_normal, entropy as scipy_entropy


class GaussianMixture:
    """A mixture of K multivariate Gaussians.

    Parameters
    ----------
    means : `numpy.ndarray`
        (K, D) array of component means.
    covariances : `numpy.ndarray`
        (K, D, D) array of covariance matrices.
    alphas : `numpy.ndarray`
        (K,) array of mixture weights (must sum to 1).
    rejection_rate : `float`, optional
        Fraction of samples that fall outside bounds, by default 0.0
    """

    def __init__(self, means, covariances, alphas, rejection_rate=0.0,
                 min_active_fraction=0.01, min_entropy_change=0.01):
        self.means = np.asarray(means)
        self.covariances = np.asarray(covariances)
        self.alphas = np.asarray(alphas, dtype=float)
        self.rejection_rate = rejection_rate
        self.min_active_fraction = min_active_fraction
        self.min_entropy_change = min_entropy_change

    @property
    def n_components(self):
        return len(self.means)

    @property
    def ndim(self):
        return self.means.shape[1]

    @classmethod
    def from_hits(cls, hit_samples, param_space, average_density_one_dim, kappa=1.0,
                  min_active_fraction=0.01, min_entropy_change=0.01):
        """Create a Gaussian mixture by placing one component at each hit.

        Parameters
        ----------
        hit_samples : `numpy.ndarray`
            (K, D) array of hit locations in sampling space.
        param_space : `ParameterSpace`
            Parameter space instance providing bounds and sigma computation.
        average_density_one_dim : `float`
            Characteristic inter-sample spacing,
            ``1 / num_explored ** (1 / D)``.
        kappa : `float`, optional
            Width scaling factor for the Gaussian covariances, by default 1.0
        min_active_fraction : `float`, optional
            Minimum value of (1 - rejection_rate) used in oversampling and
            normalisation calculations, by default 0.01
        min_entropy_change : `float`, optional
            Minimum change in normalised effective sample size (entropy) to
            avoid reverting to the previous mixture state, by default 0.01

        Returns
        -------
        `GaussianMixture`
            A new mixture with one component centred on each hit.
        """
        K = len(hit_samples)
        D = param_space.ndim

        # Compute sigma for each hit and dimension
        sigmas = param_space.compute_sigma(hit_samples, average_density_one_dim)

        # Build diagonal covariance matrices, scaled by kappa
        covariances = np.zeros((K, D, D))
        for k in range(K):
            covariances[k] = np.diag(sigmas[k] ** 2) * kappa ** 2

        # Equal mixture weights
        alphas = np.full(K, 1.0 / K)

        return cls(hit_samples.copy(), covariances, alphas,
                   min_active_fraction=min_active_fraction, min_entropy_change=min_entropy_change)

    def sample(self, n_total, param_space, consider_rejection=False, rng=None):
        """Sample from the mixture distribution.

        Parameters
        ----------
        n_total : `int`
            Total number of samples desired.
        param_space : `ParameterSpace`
            Parameter space used for bounds checking.
        consider_rejection : `bool`, optional
            If True, oversample to account for the current rejection
            rate, by default False
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None

        Returns
        -------
        samples : `numpy.ndarray`
            (M, D) array of samples in sampling space.
        mask : `numpy.ndarray`
            (M,) boolean array indicating which samples are in bounds.
        gaussian_indices : `numpy.ndarray`
            (M,) integer array indicating which component generated each
            sample.
        """
        rng = rng or np.random.default_rng()
        all_samples = []
        all_indices = []

        for k in range(self.n_components):
            n_k = int(np.ceil(n_total * self.alphas[k]))
            if consider_rejection and self.rejection_rate > 0.0:
                active = max(1.0 - self.rejection_rate, self.min_active_fraction)
                n_k = int(2 * np.ceil(n_k / active))
            if n_k <= 0:
                continue
            s = rng.multivariate_normal(self.means[k], self.covariances[k], size=n_k)
            all_samples.append(s)
            all_indices.append(np.full(n_k, k, dtype=int))

        if not all_samples:
            D = param_space.ndim
            return np.empty((0, D)), np.empty(0, dtype=bool), np.empty(0, dtype=int)

        samples = np.vstack(all_samples)
        gaussian_indices = np.concatenate(all_indices)
        mask = param_space.in_bounds(samples)

        return samples, mask, gaussian_indices

    def pdf(self, samples):
        """Evaluate the mixture PDF at the given samples.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N,) array of PDF values (weighted sum of components).
        """
        N = len(samples)
        result = np.zeros(N)
        for k in range(self.n_components):
            result += self.alphas[k] * multivariate_normal.pdf(
                samples, self.means[k], self.covariances[k], allow_singular=True
            )
        return result

    def component_pdfs(self, samples):
        """Evaluate each component's PDF at the given samples.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (K, N) array where element [k, n] is the PDF of component k
            evaluated at sample n.
        """
        K = self.n_components
        N = len(samples)
        xPDF = np.empty((K, N))
        for k in range(K):
            xPDF[k, :] = multivariate_normal.pdf(
                samples, self.means[k], self.covariances[k], allow_singular=True
            )
        return xPDF

    def compute_rejection_rate(self, param_space, reject_mask_fn,
                               n_per_component=10000, rng=None,
                               max_samples_per_call=1_000_000):
        """Estimate the rejection rate of the mixture.

        Samples from each component, transforms to physical space, applies
        rejection criteria, and computes the weighted rejection rate.

        Rather than calling ``reject_mask_fn`` once per component (which is
        slow when the rejection function has a large fixed cost per call,
        e.g. loading METISSE tracks in ``set_reff``), components are grouped
        into chunks of whole components holding at most
        ``max_samples_per_call`` samples, and ``reject_mask_fn`` is called
        once per chunk on all of its in-bounds samples.  Draws are made
        component by component in the same order regardless of chunking, so
        a fixed ``rng`` gives the same estimate for any chunk size.

        Parameters
        ----------
        param_space : `ParameterSpace`
            Parameter space for bounds checking and coordinate transforms.
        reject_mask_fn : `callable`
            Function ``(samples_physical) -> bool_mask`` returning True for
            physically rejected systems.
        n_per_component : `int`, optional
            Number of samples per component for the estimate, by default
            10000
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None
        max_samples_per_call : `int`, optional
            Maximum number of samples (before bounds rejection) passed to
            each ``reject_mask_fn`` call.  Larger values mean fewer calls but
            more memory (roughly ``max_samples_per_call * D`` floats plus
            whatever ``reject_mask_fn`` allocates).  Each chunk always holds
            at least one full component, so if ``n_per_component`` exceeds
            this value then one call is made per component.  By default
            1,000,000

        Returns
        -------
        `float`
            Estimated rejection rate (also stored as
            ``self.rejection_rate``).
        """
        rng = rng or np.random.default_rng()
        K = self.n_components
        n = n_per_component
        n_rejected = np.zeros(K)

        components_per_chunk = max(1, int(max_samples_per_call) // n)
        for start in range(0, K, components_per_chunk):
            ks = np.arange(start, min(start + components_per_chunk, K))
            s = np.concatenate([
                rng.multivariate_normal(self.means[k], self.covariances[k], size=n)
                for k in ks
            ])
            idx = np.repeat(np.arange(len(ks)), n)

            # Bounds rejection
            bounds_mask = param_space.in_bounds(s)
            rejected = np.bincount(idx[~bounds_mask], minlength=len(ks))

            # Physical rejection on in-bounds samples (all components at once)
            s_valid = s[bounds_mask]
            if len(s_valid) > 0:
                s_physical = param_space.to_physical(s_valid)
                phys_mask = np.asarray(reject_mask_fn(s_physical), dtype=bool)
                rejected += np.bincount(idx[bounds_mask][phys_mask], minlength=len(ks))

            n_rejected[ks] = rejected

        self.rejection_rate = float(np.sum(self.alphas * n_rejected / n))
        return self.rejection_rate

    def update_em(self, samples, is_hit, prior_probs, prior_fraction_rejected,
                  tolerance=1e-10, entropies=None):
        """Perform one (importance-weighted) expectation-maximization (EM) update.

        EM is the standard algorithm for fitting a mixture model.  It alternates
        an **E-step** -- computing each component's *responsibility* for every
        sample (the posterior probability that the sample was drawn from that
        component) -- with an **M-step** that re-estimates every component's
        weight, mean, and covariance as the responsibility-weighted moments of
        the samples.  Here the samples additionally carry importance weights
        (``prior_probs * is_hit / q``), so productive components (those near many
        hits) gain weight and recentre on where the hits actually are; components
        whose weight falls below ``tolerance`` are dropped.

        The update is only worth keeping if it improves the proposal.  This is
        measured by the normalised effective sample size ``exp(H(w)) / N`` (with
        ``H`` the Shannon entropy of the normalised weights): if it fails to
        increase by at least ``min_entropy_change`` over the previous generation,
        the method signals a revert (see Returns) -- STROOPWAFEL's indication
        that the mixture has stopped improving.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array of samples in sampling space.
        is_hit : `numpy.ndarray`
            (N,) boolean or integer array indicating hits.
        prior_probs : `numpy.ndarray`
            (N,) array of prior probabilities for each sample.
        prior_fraction_rejected : `float`
            Fraction of prior samples that are physically rejected.
        tolerance : `float`, optional
            Minimum mixture weight to keep a component, by default 1e-10
        entropies : `list`, optional
            List of previous entropy values (mutated in place for
            convergence tracking), by default None

        Returns
        -------
        `bool`
            True if the entropy check triggers reversion to the previous
            mixture state.
        """
        pi_norm = 1.0 / max(1.0 - prior_fraction_rejected, self.min_active_fraction)
        q_norm  = 1.0 / max(1.0 - self.rejection_rate, self.min_active_fraction)
        pi = prior_probs * pi_norm
        is_hit = np.asarray(is_hit, dtype=float)

        N = len(samples)
        K = self.n_components

        # Compute component PDFs: (K, N)
        xPDF = self.component_pdfs(samples).T  # -> (N, K)
        qPDF = xPDF * self.alphas * q_norm  # (N, K)

        # Responsibilities
        qPDF_sum = np.sum(qPDF, axis=1)  # (N,)
        rho = qPDF / qPDF_sum[:, None]  # (N, K)

        # Importance weights
        gaussian_weights = (pi * is_hit) / qPDF_sum  # (N,)
        weight_sum = np.sum(gaussian_weights)
        if weight_sum == 0:
            return False
        weights_normalized = (gaussian_weights / weight_sum)[:, None]  # (N, 1)

        # Update alphas
        new_alphas = np.sum(weights_normalized * rho, axis=0)  # (K,)

        # Remove insignificant components
        keep = new_alphas > tolerance
        if not np.any(keep):
            return False

        new_alphas = new_alphas[keep]
        rho = rho[:, keep]
        K_new = int(np.sum(keep))

        # Update means
        new_means = np.empty((K_new, self.ndim))
        for d in range(self.ndim):
            new_means[:, d] = np.sum(
                weights_normalized * samples[:, d:d+1] * rho, axis=0
            )
        new_means = new_means / new_alphas[:, None]

        # Update covariances
        old_covs = self.covariances[keep]
        new_covs = np.empty_like(old_covs)
        for k in range(K_new):
            diff = samples - new_means[k]                      # (N, D)
            outer = np.einsum('ni,nj->nij', diff, diff)        # (N, D, D)
            factor = weights_normalized[:, 0] * rho[:, k]     # (N,)
            new_covs[k] = np.sum(factor[:, None, None] * outer, axis=0) / new_alphas[k]

        # Normalised effective sample size: exp(H(w)) / N, where H is the
        # Shannon entropy of the normalised importance weights.  Ranges from
        # 1/N (all weight on one sample) to 1 (uniform weights).  We revert
        # if this metric changes by less than MIN_ENTROPY_CHANGE between
        # generations, which indicates the mixture has stopped improving.
        entropy_change = np.exp(scipy_entropy(weights_normalized[:, 0])) / N
        if entropies is not None:
            if len(entropies) >= 1 and entropy_change - entropies[-1] < self.min_entropy_change:
                return True  # Signal to revert
            entropies.append(entropy_change)

        # Apply updates
        self.means = new_means
        self.covariances = new_covs
        self.alphas = new_alphas

        return False
