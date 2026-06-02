.. _adaptive_checkpoint:

******************************
Saving and loading checkpoints
******************************

In this tutorial, we will cover how to save and load checkpoints during an adaptive importance sampling run in ``COSMIC``. This allows you to save the state of your simulation once the adaptation phase is complete, and then load it later to continue sampling from the same Gaussian mixture model (GMM) without having to redo the adaptation phase. This can be particularly useful if you want to run multiple sampling phases with the same adapted GMM (e.g. over multiple nodes on a cluster).

