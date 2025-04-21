import json
import re

# regular expression for html tags
html_tag = r"<[^>]*>"

# convertor for html tags to ini file decorators
replacer = {
    "code": "``",
    "b": "**",
    "i": "*"
}


def replace_tag(m):
    """Replace HTML tags with some decoration in the ini file

    Parameters
    ----------
    m : `re.Match`
        A match from a regular expression (containing an HTML tag)

    Returns
    -------
    r : `str`
        Replacement string
    """
    # grab the tag out of the regular expression match
    match = m.group()
    full_tag = match.replace("/", "").replace("<", "").replace(">", "")
    tag = full_tag.split(" ")[0]

    # by default delete the tag, but replace if it's in the dictionary
    replacement = ""
    if tag in replacer:
        replacement = replacer[tag]

    return replacement


def construct_ini_from_json(config, include_comments=True):
    """Construct an INI file from the JSON

    Parameters
    ----------
    config : `dict`
        Dictionary built from the JSON file
    include_comments : `bool`, optional
        Whether to include comments in the INI file, by default True

    Returns
    -------
    ini_file : `str`
        INI file for writing to a file
    """
    ini_file = "; COSMIC INI file\n"

    for category in config:
        ini_file += f"\n[{category['category']}]\n"

        for setting in category["settings"]:
            if "settings-section" in setting:
                subheading = setting["settings-section"]
                ini_file += '\n;;;;' + ';' * len(subheading) + ';;;;\n'
                ini_file += f';;; {subheading} ;;;\n'
                ini_file += '\n;;;;' + ';' * len(subheading) + ';;;;\n'

            comments = ""
            if "settings-section-description" in setting:
                comments += f'; {setting["settings-section-description"]}\n'

            comments += f'\n; {setting["name"]}\n'
            comments += f'; {setting["description"]}\n'
            if "options-preface" in setting:
                comments += f'; {setting["options-preface"]}\n'

            comments += "; Options: \n"

            default = []
            for option in setting["options"]:
                comments += f';    {option["name"]} - {option["description"]}\n'

                if "default" in option and option["default"]:
                    if setting["type"] == "checkbox":
                        default.append(str(option["name"]))
                    else:
                        default = str(option["name"])

            default_string = f"[{', '.join(default)}]" if isinstance(default, list) else f"{default}"
            comments += f"; Default: {default_string}\n"

            if include_comments:
                ini_file += comments

            ini_file += f'{setting["name"]} = {default_string}\n'

    return ini_file


def main():
    # read the config file
    with open("cosmic-settings.json") as f:
        config = json.load(f)

    # convert it to an INI file, replace HTML tags and save the output
    with open("../examples/Params.ini", "w") as f:
        f.write(re.sub(pattern=html_tag,
                       repl=replace_tag,
                       string=construct_ini_from_json(config, include_comments=True)))


if __name__ == "__main__":
    main()
