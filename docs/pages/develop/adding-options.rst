***********************************
Adding options to existing settings
***********************************

This page will guide you through the process of adding new options to existing settings in COSMIC.
This is a little simpler than adding a new setting, you can read more about how to do that `here <new-settings.html>`_.

As an example, we will add a new option for how to handle black hole supernova kicks.

Summary checklist
-----------------

Been here before and just making sure you're not missing anything? Here's a quick checklist:

- ``src/cosmic/src``: Add the new option to the relevant COSMIC code file and **test your changes**!
- ``docs/cosmic-settings.json``: Add the new option to the JSON file

Code changes
------------

First, we need to actually change the underlying COSMIC code to allow for this new option. Given that we're
going to change black hole kicks work, this will happen in the ``src/cosmic/src/kick.f`` file.

.. warning ::
    The exact file that you need to change may differ depending on what option you change!
    For example, if you were changing how mass transfer works, you'd need to change the ``src/cosmic/src/evolv2.f`` file instead.

    If you're not
    familiar with the codebase, you may want to search the ``src/cosmic/src`` directory for any mention of the
    setting to which you're adding an option.

In our case, we want to add to the ``bhflag`` setting, which is used to determine how black hole supernova kicks are handled.
At the time of writing, this setting is used as follows:

``src/cosmic/src/kick.f``

.. code-block:: fortranfixed
    :linenos:

        ....

          if(kw.eq.14.and.bhflag.eq.0)then
              vk2 = 0.d0
              vk = 0.d0
          elseif(kw.eq.14.and.bhflag.eq.1)then
              fallback = MIN(fallback,1.d0)
              vk = MAX((1.d0-fallback)*vk,0.d0)
              vk2 = vk*vk
          elseif(kw.eq.14.and.bhflag.eq.2)then
              vk = vk * mxns / m1n
              vk2 = vk*vk
          endif

        ....

There are currently 4 options implemented for the ``bhflag`` setting, which are 0, 1, 2, and 3
(Note that the code snippet above only shows 0, 1, and 2 because 3 means that the black hole kick is not changed
from the value for neutron stars). You can see conditions for each which check that ``kw`` is 14 (which ensures
this is a black hole) and then set the kick velocity according to the ``bhflag`` value.

Let's add a 5th option, which sets the black hole kick to **double** the value for neutron stars (this is
of course a rather contrived example, but it serves to illustrate the process).

``src/cosmic/src/kick.f``

.. code-block:: fortranfixed
    :linenos:
    :emphasize-lines: 13-15

        ....

          if(kw.eq.14.and.bhflag.eq.0)then
              vk2 = 0.d0
              vk = 0.d0
          elseif(kw.eq.14.and.bhflag.eq.1)then
              fallback = MIN(fallback,1.d0)
              vk = MAX((1.d0-fallback)*vk,0.d0)
              vk2 = vk*vk
          elseif(kw.eq.14.and.bhflag.eq.2)then
              vk = vk * mxns / m1n
              vk2 = vk*vk
          elseif(kw.eq.14.and.bhflag.eq.4)then
              vk = vk * 2.d0
              vk2 = vk*vk
          endif

        ....

Documentation changes
---------------------

Now we need to actual let COSMIC users that this new option exists. We'll need to update this in
a single place and this will propagate into the INI files and docs pages too (hurrah for structured data!)

Head over to the ``docs/cosmic-settings.json`` file and let's add a new option to ``bhflag``

``docs/pages/inifile.rst``

.. code-block:: json
    :linenos:
    :emphasize-lines: 25-28

        ....
            {
                "name": "bhflag",
                "description": "Sets the model for how SN kicks are applied to BHs, where bhflag != 0 allows for velocity kick at BH formation",
                "type": "dropdown",
                "options-preface": "",
                "options": [
                    {
                        "name": 0,
                        "description": "No BH kick"
                    },
                    {
                        "name": 1,
                        "description": "fallback-modulated kicks following <a class='reference external' href='https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract'>Fryer+2012</a>",
                        "default": true
                    },
                    {
                        "name": 2,
                        "description": "kicks decreased by ratio of BH mass to NS mass (1.44 Msun); conserves linear momentum"
                    },
                    {
                        "name": 3,
                        "description": "BH natal kicks are not decreased compared to NS kicks and are drawn from the same Maxwellian distribution with dispersion = <code>sigma</code> set above"
                    },
                    {
                        "name": 4,
                        "description": "A silly option that sets BH kicks to double the value for NSs"
                    }
                ]
            },
        ....

.. tip::

    If you're confused about how to format your addition to the JSON file, check out :ref:`json` to better understand this file.

This addition to the JSON file will then be used to update the Params.ini file and the config docs page with the new option.

And that's it! You've successfully added a new option to an existing setting in COSMIC, nice job!
