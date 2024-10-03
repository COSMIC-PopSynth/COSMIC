****************
Debugging COSMIC
****************

Sometimes it is useful to debug COSMIC using a debugger. This page provides a guide on how to set up debugging in Visual Studio Code.

VS Code background
------------------

If you don't have Visual Studio Code installed, you can download it from the `official website <https://code.visualstudio.com/>`_.
Additionally, if you are not familiar with debugging in Visual Studio Code, you may want to read the `general debugging documentation <https://code.visualstudio.com/docs/editor/debugging>`_ first.

Install extensions
------------------

First, you'll need a couple of extensions for VS code to make this work.
Open the extensions view by clicking on the square icon on the sidebar or pressing ``Ctrl+Shift+X``.
Search for the following extensions and install them:

- Modern Fortran
- Fortran Breakpoint Support
- C/C++

Debugging workflow explained
----------------------------

The debugging workflow uses a few different files to launch the debugging session. Let's go through them approximately in order:

``cosmic/src/test_bse.f``
    This is file that we're going to compile and run.
    This file is a simple test program that runs the BSE code on a single binary system. The input for the file
    is read from a file called ``binary.in`` and the output is written to a file called ``binary.data``.
    We'll talk about creating these files in the next section.

``cosmic/src/Makefile``
    This file sets up the compile instructions for the program. Note that we set the following flags:

    .. code-block:: makefile

        FFLAGS = -g -O0 -Wall -Wextra

    This tells the compiler to include debugging information in the executable and avoids optimisation.

``.vscode/tasks.json``
    This file contains a task that executes the compilation of the program by running ``make test`` in the
    ``cosmic/src`` directory and moves the executable to the debug directory.

``.vscode/launch.json``
    This is the file that launches the debugging sessions. If first re-compiles the test program and then
    runs it in the debugger.

Creating an input file
----------------------

To run the test program, you need to create an input file called ``binary.in``.
This file contains the input parameters for the BSE code. The format of this file can be inferred from the
``cosmic/src/test_bse.f`` file.

However, you don't need to create this file manually. In the ``debug`` directory, there is a file called
``create_binary_in.py`` which contains helper functions for creating this file. You can either manually pass
the parameters to the ``create_binary_in()`` function or use the ``convert_initC_row_to_binary_in`` to convert
a row from the ``initC`` file to a ``binary.in`` file.

Let's assume you've saved an ``initC`` file in the ``debug`` directory. You can convert the binary with
``bin_num = 42`` to a ``binary.in`` file with this code:

.. code-block:: python

    from create_binary_in import convert_initC_row_to_binary_in
    convert_initC_row_to_binary_in("initC.h5", 42)


Debugging in practice
---------------------

Now that you have the input file, you can start debugging. In VS Code, switch to the debug menu
(the bug and play button on the sidebar) and click on the green play button to start the debugging session.

Since you haven't set an breakpoints yet, the program will run until it finishes. VS Code has great documentation
about setting up `breakpoints <https://code.visualstudio.com/docs/editor/debugging#_breakpoints>`_ and more
details about `advanced breakpoints <https://code.visualstudio.com/docs/editor/debugging#_advanced-breakpoint-topics>`_
(e.g. where you can set conditions on when the breakpoint should stop). You can also use the panels to view
the value of any variable in the program and inspect the call stack.

Good luck debugging COSMIC!!