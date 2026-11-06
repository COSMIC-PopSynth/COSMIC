.. _multiprocessing:

*************************
Multiprocessing in COSMIC
*************************

Both sampling and evolution in COSMIC can be parallelised across multiple CPU cores, which can
significantly reduce runtime when working with large populations. This guide explains how to use
multiprocessing for the independent sampler and for :meth:`~cosmic.evolve.Evolve.evolve`.

There are two ways to enable multiprocessing in COSMIC:

- Pass ``nproc`` — an integer specifying how many worker processes to use. COSMIC will create a
  pool internally, use it, and then shut it down automatically.
- Pass an existing ``pool`` — a :class:`multiprocessing.Pool` (or any pool-like object with a
  ``map`` and ``starmap`` interface) that you have already created. COSMIC will use it but will
  **not** close it, so you remain in control of when to close it.

.. note::

    Multiprocessing in Python requires that parallelised code runs inside a
    ``if __name__ == "__main__":`` guard when using the default ``fork``-free start methods
    (``spawn`` on Windows and macOS). All of the examples below define a ``main()`` function
    and call it under this guard. If you are running interactively in a Jupyter notebook, you
    may need to use a workaround such as ``multiprocessing.get_context("fork")`` or simply
    set ``nproc=1``.


How many cores should you use?
==============================

How many do you have?
---------------------

Before deciding on ``nproc``, it is worth knowing how many CPU cores are available on your
machine. You can check this with Python's :mod:`os` or :mod:`multiprocessing` modules:

.. code-block:: python

    import os
    import multiprocessing

    print(os.cpu_count())                  # logical cores (including hyperthreads)
    print(multiprocessing.cpu_count())     # same value, via the multiprocessing module

The number reported includes hyperthreaded logical cores. For CPU-bound work like binary
evolution the physical core count is often more relevant, though this varies by machine and
workload.

More cores is not always faster
-------------------------------

Spawning worker processes has a fixed overhead cost: Python needs to start each worker,
import modules, and copy data to and from the workers. For small populations this overhead can
easily exceed any time saved from parallelism, making multiprocessing *slower* than a single
process.

As a rough guide you will see speedups from multiprocessing once you have at least a few thousand binaries,
but the exact threshold depends on your machine and the details of your workload.

The best way to know whether more cores helps *your* specific workload is to benchmark it
directly. Here is a simple timing script you can adapt:

.. include:: ../_generated/default_bsedict.rst

.. code-block:: python

    import time
    import multiprocessing
    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve

    def benchmark(size=50000, n_cores=1):
        with multiprocessing.Pool(n_cores) as pool:
            t0 = time.perf_counter()
            ibt = InitialBinaryTable.sampler(
                'independent', [13, 14], [13, 14],
                binfrac_model=0.5, primary_model='kroupa01',
                ecc_model='sana12', porb_model='sana12',
                qmin=-1, SF_start=13700.0, SF_duration=0.0,
                met=0.02, size=size, pool=pool,
            )[0]
            print(f"    Sampling {size} binaries with nproc={n_cores} took {time.perf_counter() - t0:.1f}s")

            t0 = time.perf_counter()
            Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, pool=pool)
            elapsed = time.perf_counter() - t0
            print(f"    Evolving {size} binaries with nproc={n_cores} took {elapsed:.1f}s")


    if __name__ == "__main__":
        for size in [10, 1000]:
            print(f"Benchmarking with population size {size}:")
            for n in [1, 2, 4, 8, 16]:
                benchmark(size=size, n_cores=n)
            print()

Running this on your machine will show you the point of diminishing returns. You may find that
4 cores is nearly as fast as 16 for a small population,
while for a large population 16 cores is much faster than 4.

.. tip::

    For most research-quality simulations you'll likely just want to use the maximum number of cores available,
    but benchmarking can be useful if you want to know how long a run will take or if you are running on a shared cluster
    where you have to request a specific number of cores and may want to avoid requesting more than you need.


Parallelising the independent sampler
======================================

Using ``nproc``
---------------

The simplest way to parallelise sampling is to pass ``nproc`` directly to
:meth:`~cosmic.sample.initialbinarytable.InitialBinaryTable.sampler`. ``COSMIC`` splits the
requested sample into chunks, distributes them across ``nproc`` worker processes, and combines
the results.

.. code-block:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable


    def main():
        ibt, mass_singles, mass_binaries, n_singles, n_binaries = \
            InitialBinaryTable.sampler(
                'independent', [13, 14], [13, 14],
                binfrac_model=0.5, primary_model='kroupa01',
                ecc_model='sana12', porb_model='sana12',
                qmin=-1, SF_start=13700.0, SF_duration=0.0,
                met=0.02, size=100_000,
                nproc=4,
            )
        print(ibt)


    if __name__ == "__main__":
        main()

``COSMIC`` will only spawn as many chunks as make sense given the requested sample size — if the
sample is too small to benefit from parallelism it will silently cap the number of chunks, so
passing a large ``nproc`` on a small sample is harmless.

Using an existing pool
-----------------------

If you are already managing a :class:`multiprocessing.Pool` — for example because you want to
reuse it across multiple sampling calls, or with evolution as well — you can pass it via the ``pool`` argument instead.
When a ``pool`` is supplied, the ``nproc`` argument is ignored.

.. code-block:: python

    from multiprocessing import Pool
    from cosmic.sample.initialbinarytable import InitialBinaryTable


    def main():
        with Pool(4) as pool:
            ibt, mass_singles, mass_binaries, n_singles, n_binaries = \
                InitialBinaryTable.sampler(
                    'independent', [13, 14], [13, 14],
                    binfrac_model=0.5, primary_model='kroupa01',
                    ecc_model='sana12', porb_model='sana12',
                    qmin=-1, SF_start=13700.0, SF_duration=0.0,
                    met=0.02, size=100000,
                    pool=pool,
                )
            print(ibt)


    if __name__ == "__main__":
        main()


Parallelising evolution
========================

Using ``nproc``
---------------

:meth:`~cosmic.evolve.Evolve.evolve` accepts an ``nproc`` keyword argument. Each binary is
evolved independently, so the workload is split across ``nproc`` worker processes automatically.

.. include:: ../_generated/default_bsedict.rst

.. code-block:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve


    def main():
        ibt, mass_singles, mass_binaries, n_singles, n_binaries = \
            InitialBinaryTable.sampler(
                'independent', [13, 14], [13, 14],
                binfrac_model=0.5, primary_model='kroupa01',
                ecc_model='sana12', porb_model='sana12',
                qmin=-1, SF_start=13700.0, SF_duration=0.0,
                met=0.02, size=100,
            )

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=ibt,
            BSEDict=BSEDict,
            nproc=4,
        )
        print(bpp)


    if __name__ == "__main__":
        main()

Using an existing pool
-----------------------

If you want to reuse the same pool for evolution — or share one across sampling and evolution —
pass the pool object as the ``pool`` argument. Like the sampler, :meth:`~cosmic.evolve.Evolve.evolve`
will use the pool but leave it open when it finishes.

.. code-block:: python

    from multiprocessing import Pool
    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve


    def main():
        with Pool(4) as pool:
            ibt, mass_singles, mass_binaries, n_singles, n_binaries = \
                InitialBinaryTable.sampler(
                    'independent', [13, 14], [13, 14],
                    binfrac_model=0.5, primary_model='kroupa01',
                    ecc_model='sana12', porb_model='sana12',
                    qmin=-1, SF_start=13700.0, SF_duration=0.0,
                    met=0.02, size=100,
                    pool=pool,
                )

            bpp, bcm, initC, kick_info = Evolve.evolve(
                initialbinarytable=ibt,
                BSEDict=BSEDict,
                pool=pool,
            )

        print(bpp)


    if __name__ == "__main__":
        main()

Sharing a single pool avoids the overhead of spawning and tearing down worker processes twice,
which can be a meaningful saving when running many iterations.

.. tip::

    If you are running a parameter sweep and calling :meth:`~cosmic.evolve.Evolve.evolve` many
    times with different physics settings, creating one pool at the start and passing it to each
    call is more efficient than relying on ``nproc``.


Wrapping up
===========

Multiprocessing in COSMIC is basically just a matter of using ``nproc`` or passing a pool directly, but some
takeaways to keep in mind are:

- Use ``nproc`` for the simplest case — COSMIC handles pool creation and cleanup automatically.
- Use an existing ``pool`` when you want fine-grained control over worker processes or need to
  share the same pool across sampling and evolution.
- Both approaches work identically for the independent sampler and for
  :meth:`~cosmic.evolve.Evolve.evolve`.
- Always define a ``main()`` function and call it under ``if __name__ == "__main__":`` to
  avoid issues with Python's multiprocessing start methods on macOS and Windows.
- Benchmark before scaling up — more cores only helps once your population is large enough to
  outweigh the overhead of spawning workers.
