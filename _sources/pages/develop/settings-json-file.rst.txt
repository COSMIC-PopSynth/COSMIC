.. _json:

***************************************
``cosmic-settings.json`` file explained
***************************************

This file contains all of the settings information for every setting available in COSMIC.
We write all of this information in a single centralised JSON file and this can then be used to populate all other
parts of the COSMIC defaults and documentation.

In this guide let's take a look at the JSON schema for the file so it is clear how to add new settings or options.

File structure
==============

The high level structure of the file looks something like this:

.. code-block:: javascript

    {
        "category": "CATEGORY NAME",
        ...,
        "settings": [
            ...,
            {
                "name": "SETTING NAME",
                ...,
                "options": [
                    ...,
                    {
                        "name": "OPTION NAME",
                        ...,
                    },
                    ...
                ]
            }
            ...
        ]

where the ``...`` lines are extra parameters that we're not showing here (but do show below).
You can see there are a series of categories (e.g. Sampling, BSE),
settings (e.g. ``sampling_method``, ``alpha1``) and options (e.g. ``independent`` or ``multidim``) in the file.
I'm going to assume you'll not be creating new categories but it should be fairly similar to new settings or options - let's go into those two in more detail.

Schema for a setting
====================

Here is a list of the potential parameters in the setting object:

.. csv-table:: Schema for a setting
   :file: setting-schema.csv
   :header-rows: 1


In general any long strings (e.g. in the description or options-preface) you can use full HTML tags and they
will be rendered (e.g. bold some text, link a paper). See below if you're not familiar with this!
It may be simplest to just look at some existing settings if you're unsure.


Schema for an option
====================

And here's the same sort of thing for an option object, several of which comprise the ``setting["options"]`` array
in each setting.

.. csv-table:: Schema for an option
   :file: option-schema.csv
   :header-rows: 1


Quick HTML cheatsheet
=====================

As noted above, if you want to bold some text or link a paper in your setting/option descriptions 
then here's some quick HTML tags for doing that.
Once the options are converted into files these tags will get stripped out.
You can also add inline maths symbols following the last line
(**Notice the escaping of every backslash - i.e. two backslashes needed**).

Before (HTML code)
------------------

.. code-block:: HTML

    <b>This text is bold</b>
    <i>This text is italic</b>
    Here is <code>some code</code>
    This line is broken<br>by a new line
    <a href="https://google.com">This bit</a> is a link
    Here's a fancy equation \\( y = m x + c \\)

After (rendered results)
------------------------

These will render to approximately the following

    **This text is bold**

    *This text is italic*

    Here is ``some code``

    | This line is broken
    | by a new line

    `This bit <https://google.com>`_ is a link

    Here's a fancy equation :math:`y = m x + c`