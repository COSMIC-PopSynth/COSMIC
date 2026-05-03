import json
from pathlib import Path

def get_default_BSE_settings(to_python=False):
    """Get a copy of the default BSE settings from the COSMIC settings JSON file"""
    json_path = (
        Path(__file__).resolve().parent.parent
        / "src"
        / "cosmic"
        / "data"
        / "cosmic-settings.json"
    )

    with json_path.open(encoding="utf-8") as f:
        cosmic_settings = json.load(f)

    defaults = {}
    # loop through settings to find bse
    for cat in cosmic_settings:
        if cat["category"] != "bse":
            continue
        
        # go through each setting, finding the default option
        for setting in cat["settings"]:
            for option in setting["options"]:
                if option.get("default", False):
                    defaults[setting["name"]] = option["name"]

    if to_python:
        # ensure array settings are converted from strings to lists
        for setting in ["qcrit_array", "natal_kick_array", "fprimc_array", "alpha1", "acc_lim"]:
            # this one requires special handling because of the fractions
            if setting == "fprimc_array":
                parts = defaults[setting].strip("[]").split(",")
                defaults[setting] = [float(p.split("/")[0]) / float(p.split("/")[1]) for p in parts]
            else:
                defaults[setting] = json.loads(defaults[setting])

    return defaults


def generate_rst_bsedict(MAX_LINE_LENGTH=80):
    """Generate the default BSE settings dictionary for use in documentation"""
    defaults = get_default_BSE_settings()

    indent_str = " " * 8

    with open("_generated/default_bsedict.rst", "w") as f:
        lines = [
            ".. ipython:: python",
            "",
            "    BSEDict = {",
        ]

        first = True
        current_line = indent_str
        for k, v in defaults.items():
            entry = f'"{k}": {v}'
            if not first:
                entry = ", " + entry
            first = False

            # check if adding this entry would exceed line length
            if len(current_line) + len(entry) > MAX_LINE_LENGTH:
                lines.append(current_line.rstrip() + ",")
                current_line = indent_str + entry.lstrip(", ")
            else:
                current_line += entry

        if current_line.strip():
            lines.append(current_line.rstrip())
        lines.append("    }")

        f.write("\n".join(lines))

if __name__ == "__main__":
    generate_rst_bsedict()