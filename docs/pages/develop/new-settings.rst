*********************
Creating new settings
*********************

*Last updated: 2026-02-17 by Tom Wagg*

This page guides you through adding an entirely new setting to COSMIC. If you just want to add a new option to an existing setting, you can read about how to do that `here <adding-options.html>`_.

.. admonition:: Summary checklist
    :class: tip
    
    Been here before and just making sure you're not missing anything? Here's a quick checklist:
    
    - ``src/cosmic/src/const_bse.h``: Initialise a new variable and add it to the relevant common block
    - ``src/cosmic/data/cosmic-settings.json``: Add the setting to the JSON file, with a description and default value
    - ``src/cosmic/evolve.py``
         - Add variable to ``INITIAL_CONDITIONS_BSE_COLUMNS``
         - Add variable in ``_evolve_single_system()`` list
    - ``src/cosmic/tests/data``: Update initC files and params.ini with new setting
         - ``src/cosmic/tests/data/initial_conditions_for_testing.hdf5``
         - ``src/cosmic/tests/data/kick_initial_conditions.h5``
         - ``src/cosmic/tests/data/Params.ini``
    - Actually use your new setting somewhere in the code! (I can't tell you where ;P)


Adding a new setting to ``COSMIC`` is a little more involved than just adding an options. Settings are stored in common blocks in the Fortran code so that they can be accessed anywhere in the code. Let's do it step by step.

Throughout this example, let's say we want to add a new setting called ``lbv_flag`` that sets the prescription for luminous blue variable mass loss.

Add new setting to a common block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we need to add the new setting to a common block in the Fortran code. This is done in ``src/cosmic/src/const_bse.h``. You can either add it to an existing common block or create a new one if it doesn't fit with any of the existing ones. For this example, let's say we want to add it to the ``WINDVARS`` common block since that's the most relevant.

In ``src/cosmic/src/const_bse.h``, we first declare the variable (in this case it's just an integer flag)

.. code-block:: fortran
    
           INTEGER lbv_flag

Then we add it to the common block declaration for ``WINDVARS``:

.. code-block:: fortran
    
           COMMON /WINDVARS/ neta,bwind,hewind,beta,xi,acc2,epsnov,
          &                  eddfac,gamma,lbv_flag


Allow user to specify the new setting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we need to allow the user to specify the new setting in a BSEDict or Params.ini. To enable this, we need to add it to ``src/cosmic/data/cosmic-settings.json``. This file is used to generate the default BSEDict and Params.ini file and also serves as a reference for all the settings available in COSMIC.

In ``src/cosmic/data/cosmic-settings.json``, we add a new entry for our setting:

.. code-block:: json
    
    "name": "lbv_flag",
    "description": "Flag for luminous blue variable wind mass loss prescription.",
    "type": "dropdown",
    "options-preface": "",
    "options": [
        {
            "name": 0,
            "description": "No LBV winds"
        },
        {
            "name": 1,
            "description": "LBV winds according to the prescription of <a href='https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract'>Hurley+2000</a>",
            "default": true
        },
        {
            "name": 2,
            "description": "LBV winds according to the prescription of <a href='https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract'>Belczynski+2008</a>"
        }
    ]

.. note::

    Notice that we set a ``type`` of ``dropdown`` since this setting is a flag that can only take a few discrete values. If your setting is a continuous variable, you would set the type to ``number`` instead. We also set a ``default`` value in the options, which is important for creating the default BSEDict and Params.ini files.


Pass the setting from evolve() to the Fortran code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we need to pass the setting from the Python code in ``evolve.py`` to the Fortran code. This is done in the ``_evolve_single_system()`` function in ``src/cosmic/evolve.py``. But first, we need to add the new setting to the list of columns in the initial binary table. This is done by adding it to the ``INITIAL_CONDITIONS_BSE_COLUMNS`` list at the top of ``src/cosmic/evolve.py``:

.. code-block:: python
    
    # dots are representing the existing columns
    INITIAL_CONDITIONS_BSE_COLUMNS = [..., 'lbv_flag']

Now we change the ``_evolve_single_system()`` function to pass the new setting to the Fortran code. In the list of flags that is passed to the Fortran code, we add our new setting:

.. code-block:: python
    
    _evolvebin.windvars.lbv_flag = f["lbv_flag"]

.. warning::

    Make sure you use the right common block here (``windvars`` in this case) otherwise the setting won't get properly passed but will also not throw an error which can be very confusing to debug! (I have made this mistake before and it was not fun to figure out what was going wrong)

Use the new setting in the Fortran code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This bit is up to you ;P In the case of our example setting, we would need to find the part of the Fortran code where luminous blue variable winds are calculated (in ``mlwind.f`` in this case) and modify the code to use the new setting to determine which prescription to use.

Update the test data
^^^^^^^^^^^^^^^^^^^^

.. note::

    This step will likely change when we improve the testing framework, but for now we need to manually update the test data to include our new setting. If it seems out of date by the time you read this, please make an issue to update this section with the new process for updating the tests.

Finally, we need to update the test data to include our new setting. This is important to ensure that the tests still pass and that the new setting is being properly tested. We need to update the initial conditions files and the Params.ini file in ``src/cosmic/tests/data``:
- ``src/cosmic/tests/data/initial_conditions_for_testing.hdf5``
- ``src/cosmic/tests/data/kick_initial_conditions.h5``
- ``src/cosmic/tests/data/Params.ini``

For the params.ini file, add a new entry for the new setting with a reasonable default value. For the initial conditions files, we need to add a new column for the new setting and fill it with reasonable values. You'll want to read in the initC files in Python, add the new column, fill it with values, and then write the files back out. Here's a quick helper function that does that for you

.. code-block:: python
    
    import pandas as pd
    
    def add_setting_to_test_files(setting_name, setting_value):
        initC = pd.read_hdf("src/cosmic/tests/data/initial_conditions_for_testing.hdf5", key="initC")
        initC[setting_name] = setting_value
        initC.to_hdf("src/cosmic/tests/data/initial_conditions_for_testing.hdf5", key="initC", mode="a")

        initC = pd.read_hdf("src/cosmic/tests/data/kick_initial_conditions.h5", key="initC")
        initC[setting_name] = setting_value
        initC.to_hdf("src/cosmic/tests/data/kick_initial_conditions.h5", key="initC", mode="a")

        with open("src/cosmic/tests/data/Params.ini", "a") as f:
            f.write(f"\n\n{setting_name} = {setting_value}\n\n")

    add_setting_to_test_files("lbv_flag", 1)

And that's it! You've added a new setting to COSMIC. Now you should absolutely run the tests to make sure everything is working properly! Better yet, add more tests that specifically test the new setting to make sure it's working as expected.