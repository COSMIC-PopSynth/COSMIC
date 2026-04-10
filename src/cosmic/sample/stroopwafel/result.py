"""Container for STROOPWAFEL adaptive sampling results."""
import numpy as np


class STROOPWAFELResult:
    """Results from an AdaptiveSampler run.

    Parameters
    ----------
    samples : `numpy.ndarray`
        (N, D) array of all samples in physical space.
    param_names : `list` of `str`
        Parameter names defining the column ordering.
    weights : `numpy.ndarray`
        (N,) array of importance sampling weights.
    is_hit : `numpy.ndarray`
        (N,) boolean array indicating which samples are hits.
    generation : `numpy.ndarray`
        (N,) integer array (0 = exploration, 1+ = refinement).
    gaussian_idx : `numpy.ndarray`
        (N,) integer array (-1 = exploration, k = from Gaussian
        component k).
    num_explored : `int`
        Number of systems simulated during the exploration phase.
    num_hits : `int`
        Total number of hits found across all phases.
    fraction_explored : `float`
        Adaptive exploration fraction.
    bpp_frames : `list` of `pandas.DataFrame`, optional
        COSMIC binary population parameter tables, by default None
    initC_frames : `list` of `pandas.DataFrame`, optional
        COSMIC initial conditions tables, by default None
    kick_info_frames : `list` of `pandas.DataFrame`, optional
        COSMIC natal kick information tables, by default None
    """

    def __init__(self, samples, param_names, weights, is_hit, generation,
                 gaussian_idx, num_explored, num_hits, fraction_explored,
                 bpp_frames=None, initC_frames=None, kick_info_frames=None):
        self.samples = samples
        self.param_names = param_names
        self.weights = weights
        self.is_hit = is_hit
        self.generation = generation
        self.gaussian_idx = gaussian_idx
        self.num_explored = num_explored
        self.num_hits = num_hits
        self.fraction_explored = fraction_explored
        self.bpp_frames = bpp_frames or []
        self.initC_frames = initC_frames or []
        self.kick_info_frames = kick_info_frames or []

    @property
    def hit_rate(self):
        """Importance-weighted hit rate (sum of weights over hits / N).

        Returns
        -------
        `float`
            Weighted fraction of systems that are hits, or 0.0 if
            weights or hit flags are not available.
        """
        if self.weights is None or self.is_hit is None:
            return 0.0
        return np.sum(self.weights[self.is_hit]) / len(self.weights)

    @property
    def hit_rate_uncertainty(self):
        """Standard error on the importance-weighted hit rate.

        Returns
        -------
        `float`
            Standard error ``std(w_hits) / sqrt(N)``, or 0.0 if fewer
            than two hits are present.
        """
        if self.weights is None or self.is_hit is None:
            return 0.0
        w_hits = self.weights[self.is_hit]
        if len(w_hits) < 2:
            return 0.0
        return np.std(w_hits, ddof=1) / np.sqrt(len(self.weights))
