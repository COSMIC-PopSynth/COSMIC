import bs4
import json
import pandas as pd

main_soup = """<!-- This file should be created using create_settings_html.py -->
<div class="container-fluid"></div>
"""

group_template = """<section class="card setting-card">
        <div class="card-header">
            <h2 class="setting-group-title"></h2>
            <p class="setting-group-description"></p>
        </div>
        <div class="card-body"></div>
    </section>"""

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

option_template = """<li><code class="docutils literal notranslate"><span class="pre opt-val"></span></code>: <span class="opt-desc"></span></li>"""

with open("cosmic-settings.json") as f:
    settings = json.load(f)

for group in settings:
    soup = bs4.BeautifulSoup(main_soup, 'html.parser')

    print(f"Creating group {group['category']}")
    new_group = bs4.BeautifulSoup(group_template, 'html.parser')
    new_group.section["data-category"] = group["category"]
    title = new_group.select_one(".setting-group-title")
    title.string = group["category_label"]
    new_group.select_one(".setting-group-description").string = group["category_description"]
    new_group.select_one(".setting-card")["style"] = "border-color: " + group["docs-colour"] + ";"
    new_group.select_one(".card-header")["style"] = "background-color: " + group["docs-colour"] + ";"

    for setting in group["settings"]:
        if "settings-section" in setting:
            settings_section = new_group.new_tag("div")
            settings_section["class"] = "settings-section"

            header = new_group.new_tag("h2")
            header["id"] = setting["settings-section"].lower().replace(" ", "-")
            header.string = setting["settings-section"]

            settings_section.append(header)
            settings_section.append(new_group.new_tag("hr"))

            if "settings-section-description" in setting:
                settings_section.append(bs4.BeautifulSoup(setting["settings-section-description"], 'html.parser'))

            new_group.select_one(".card-body").append(settings_section)

        print(f"Creating setting {setting['name']}")
        new_setting = bs4.BeautifulSoup(settings_template, 'html.parser')
        new_setting.select_one(".name").code.string = setting["name"]
        new_setting.select_one(".description").string = ""
        new_setting.select_one(".description").append(bs4.BeautifulSoup(setting["description"], 'html.parser'))
        new_setting.select_one(".options-expander")["style"] = "color: " + group["docs-colour"] + ";"

        default = []

        if setting["type"] == "checkbox":
            insert_here = new_setting.select_one(".col-3")
            for i, option in enumerate(setting["options"]):
                if "default" not in option:
                    option["default"] = False
                elif option["default"]:
                    default.append(option["name"])
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

            secret_input = new_setting.new_tag("input", type="string", value="")
            secret_input["class"] = "form-control hide-me"
            secret_input["value"] = str(default).replace("'", "")
            insert_here.append(secret_input)

        elif setting["type"] in ["string", "number"]:
            for option in setting["options"]:
                if "default" in option and option["default"]:
                    default = option["name"]
                    new_option = new_setting.new_tag("input",
                                                     type="text" if setting["type"] == "string" else "number",
                                                     value=str(option["name"]))
                    new_option["class"] = "form-control"
                    new_setting.select_one(".col-3").append(new_option)
                    break

        elif setting["type"] == "dropdown":
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

        new_setting.select_one(".options-preface").string = ""
        new_setting.select_one(".options-preface").append(bs4.BeautifulSoup(setting["options-preface"], 'html.parser'))

        for i, option in enumerate(setting["options"]):
            print(f"Creating option {option['name']}")
            new_option_expl = bs4.BeautifulSoup(option_template, 'html.parser')
            new_option_expl.select_one(".opt-val").string = str(option["name"])
            new_option_expl.select_one(".opt-desc").append(bs4.BeautifulSoup(option["description"], 'html.parser'))
            new_setting.select_one(".options").ul.append(new_option_expl)

        default_string = str(default)
        default_string = default_string.replace("'", "")
        new_setting.select_one(".default").string = f"Default: {default_string}"

        if setting["name"] == "qcrit":
            qcrit_table_string = pd.read_csv("data/qcrit_table.csv").to_html(index=False, justify="center")
            qcrit_table = bs4.BeautifulSoup(qcrit_table_string, 'html.parser')
            qcrit_table.table["class"] = "table table-striped-columns table-hover table-responsive"
            qcrit_table.table["border"] = "0"
            qcrit_table.table["style"] = "text-align:center!important;"

            caption = new_setting.new_tag("caption")
            caption.string = "Comparison of qcrit Values (Donor Mass/Accretor Mass) For Each Donor Kstar Type Across Flag Options"

            qcrit_table.table.append(caption)
            new_setting.select_one(".options").append(qcrit_table)

        new_group.select_one(".card-body").append(new_setting)

    soup.select_one(".container-fluid").append(new_group)

    with open(f"pages/config/config_insert_{group['category']}.html", "w") as f:
        f.write(str(soup))

