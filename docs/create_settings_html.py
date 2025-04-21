"""
    create_settings_html.py

    This file is used to construct the HTML files for the config/settings page
    based on the information stored in cosmic-settings.json.

    Written by Tom Wagg
"""

import bs4
import json
import pandas as pd

# blame BS4 for me calling this soup
main_soup = """<!-- This file should be created using create_settings_html.py -->
<div class="container-fluid"></div>
"""

# a template for a settings group (we'll copy this and fill it in)
group_template = """<section class="card setting-card">
        <div class="card-header">
            <h2 class="setting-group-title"></h2>
            <p class="setting-group-description"></p>
        </div>
        <div class="card-body"></div>
    </section>"""

# as above but for a specific setting (this'll go in 'card-body' above)
settings_template = """<div class="setting">
                <div class="row align-items-center setting-chooser">
                    <div class="col-9">
                        <h3 class="name"><code></code></h3>
                        <p class="description"></p>
                        <p class="default"></p>
                    </div>
                    <div class="col-3"></div>
                </div>
                <div class="options-expander">Option details <i class="fa fa-chevron-down"></i></div>
                <div class="row options hide">
                    <p class="options-preface"></p>
                    <ul style="margin-left: 2rem; max-width: calc(100% - 2rem)"></ul>
                </div>
            </div>"""

# same as above, but for an option for a setting, this will go in the <ul> element
option_template = """<li><code class="docutils literal notranslate"><span class="pre opt-val"></span></code>: <span class="opt-desc"></span></li>"""

# read the settings file
with open("cosmic-settings.json") as f:
    settings = json.load(f)

# go through each major settings group
for group in settings:
    # start the overall main html file
    soup = bs4.BeautifulSoup(main_soup, 'html.parser')

    # start a new group, writing in its category name, label and description
    new_group = bs4.BeautifulSoup(group_template, 'html.parser')
    new_group.section["data-category"] = group["category"]
    title = new_group.select_one(".setting-group-title")
    title.string = group["category_label"]
    new_group.select_one(".setting-group-description").string = group["category_description"]

    # update the colours based on the choices in the JSON file
    new_group.select_one(".setting-card")["style"] = "border-color: " + group["docs-colour"] + ";"
    new_group.select_one(".card-header")["style"] = "background-color: " + group["docs-colour"] + ";"

    # proceed through the settings in this group
    for setting in group["settings"]:
        # if this is the start of a subsection (e.g. "Mass transfer")
        if "settings-section" in setting:
            # create a section and add a heading with the title
            settings_section = new_group.new_tag("div")
            settings_section["class"] = "settings-section"

            header = new_group.new_tag("h2")
            header["id"] = setting["settings-section"].lower().replace(" ", "-")
            header.string = setting["settings-section"]

            settings_section.append(header)
            settings_section.append(new_group.new_tag("hr"))

            # add an optional description if it's there (e.g. Common-envelope has this)
            if "settings-section-description" in setting:
                settings_section.append(bs4.BeautifulSoup(setting["settings-section-description"],
                                                          'html.parser'))

            new_group.select_one(".card-body").append(settings_section)

        # create the new settings with its details
        new_setting = bs4.BeautifulSoup(settings_template, 'html.parser')
        new_setting.select_one(".name").code.string = setting["name"]
        new_setting.select_one(".description").clear()
        new_setting.select_one(".description").append(bs4.BeautifulSoup(setting["description"],
                                                                        'html.parser'))

        # colour the sublinks the same as the border of the group
        new_setting.select_one(".options-expander")["style"] = "color: " + group["docs-colour"] + ";"

        # track the default value
        default = []

        if setting["type"] == "checkbox":
            insert_here = new_setting.select_one(".col-3")
            for i, option in enumerate(setting["options"]):
                # append to the default as necessary
                if "default" not in option:
                    option["default"] = False
                elif option["default"]:
                    default.append(option["name"])

                # create a checkbox with the option details
                new_option = new_setting.new_tag("input", type="checkbox", value=option["name"])
                if option["default"]:
                    new_option["checked"] = "true"
                new_option["name"] = f"{setting['name']}_opt{i}"
                new_label = new_setting.new_tag("label", for_=f"{setting['name']}_opt{i}")
                new_label.string = f" {option['name']}"

                new_container = new_setting.new_tag("div")
                new_container.append(new_option)
                new_container.append(new_label)

                insert_here.append(new_container)

            # create a secret hidden input which combines the checkbox inputs
            secret_input = new_setting.new_tag("input", type="string", value="")
            secret_input["class"] = "form-control hide-me"
            secret_input["value"] = str(default).replace("'", "")
            insert_here.append(secret_input)

        elif setting["type"] in ["string", "number"]:
            for option in setting["options"]:
                # create an input with the value of the default option
                if "default" in option and option["default"]:
                    default = option["name"]
                    new_option = new_setting.new_tag("input",
                                                     type="text" if setting["type"] == "string" else "number",
                                                     value=str(option["name"]))
                    new_option["class"] = "form-control"
                    new_setting.select_one(".col-3").append(new_option)
                    break

        elif setting["type"] == "dropdown":
            # create a dropdown input with all of the potential options
            new_select = new_setting.new_tag("select")
            new_select["class"] = "form-control"
            for i, option in enumerate(setting["options"]):
                new_option = new_setting.new_tag("option", value=option["name"])
                new_option.string = str(option["name"])
                if "default" in option and option["default"]:
                    new_option["selected"] = "true"
                    default = option["name"]
                new_select.append(new_option)
            new_setting.select_one(".col-3").append(new_select)

        # optionally there might be some extra explanation about the options
        if setting["options-preface"] != "":
            new_setting.select_one(".options-preface").clear()
            new_setting.select_one(".options-preface").append(bs4.BeautifulSoup(setting["options-preface"],
                                                                                'html.parser'))

        # add each of the options and its explanation to the list
        for i, option in enumerate(setting["options"]):
            new_option_expl = bs4.BeautifulSoup(option_template, 'html.parser')
            new_option_expl.select_one(".opt-val").string = str(option["name"])
            new_option_expl.select_one(".opt-desc").append(bs4.BeautifulSoup(option["description"], 'html.parser'))
            new_setting.select_one(".options").ul.append(new_option_expl)

        # convert the default options to a string and display it
        default_string = str(default)
        default_string = default_string.replace("'", "")
        new_setting.select_one(".default").string = f"Default: {default_string}"

        # special case: table of values for the critical mass ratios
        if setting["name"] == "qcflag":
            # read in the table and convert it to HTML
            qcrit_table_string = pd.read_csv("data/qcrit_table.csv").to_html(index=False, justify="center")

            # parse the html into some soup and add a couple of bootstrap classes to pretty it up
            qcrit_table = bs4.BeautifulSoup(qcrit_table_string, 'html.parser')
            qcrit_table.table["class"] = "table table-striped-columns table-hover table-responsive"
            qcrit_table.table["border"] = "0"
            qcrit_table.table["style"] = "text-align:center!important;"

            # add a caption to the table
            caption = new_setting.new_tag("caption")
            caption.string = "Comparison of qcrit Values (Donor Mass/Accretor Mass) For Each Donor Kstar Type Across Flag Options"

            # attach it to the main soup
            qcrit_table.table.append(caption)
            new_setting.select_one(".options").append(qcrit_table)

        # append this setting to the group
        new_group.select_one(".card-body").append(new_setting)

    # append this group to the soup (guess who's a poet and didn't even know it)
    soup.select_one(".container-fluid").append(new_group)

    # write the soup out to an HTML file for this category
    with open(f"pages/config/config_insert_{group['category']}.html", "w") as f:
        f.write(str(soup))
