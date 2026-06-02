"""STROOPWAFEL adaptive importance sampling for COSMIC.

Provides a vectorized reimplementation of the STROOPWAFEL algorithm
for efficiently sampling rare binary stellar evolution outcomes.  The
three-phase pipeline (exploration → adaptation → refinement) places
Gaussian components at discovered hits, then importance-samples from
the resulting mixture to concentrate compute budget on interesting
regions of parameter space.

"""
from .engine import AdaptiveSampler
from .parameter_space import ParameterSpace, Parameter

__all__ = ['AdaptiveSampler', 'ParameterSpace', 'Parameter']
