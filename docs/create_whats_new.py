"""
    create_whats_new.py
    
    This file constructs the "What's New" section of the documentation, which highlights new features and settings in the latest version of COSMIC. It ingests the changelog and converts it to an rst file for the docs.

    Each heading that starts with ## is converted to a new section in the "What's New" page, and each bullet point is converted to a new feature.
    The script also links any headings of the form "## x.x.x" to the corresponding release in GitHub

    Written by Tom Wagg
"""

import re
from pathlib import Path
from os.path import exists

def parse_changelog(changelog_path):
    with open(changelog_path, 'r') as f:
        changelog = f.read()

    # Split the changelog into sections based on "## " headings
    sections = re.split(r'##\s+', changelog)[1:]  # Skip the first split which is before the first heading

    parsed_sections = []
    for section in sections:
        lines = section.strip().split('\n')
        heading = lines[0].strip()
        features = [line for line in lines[1:]]

        # if the first character of the heading is not a number, skip
        if not re.match(r'^\d', heading):
            continue

        # for each feature, replace markdown links of the form [text](link) with rst links of the form `text <link>`_
        features = [re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'`\1 <\2>`_', feature) for feature in features]

        parsed_sections.append((heading, features))

    return parsed_sections

def generate_rst(parsed_sections, output_path):
    lines = ["What's New in COSMIC\n====================\n"]

    for heading, features in parsed_sections:
        lines.append(f"v{heading}\n{'-' * (len(heading) + 1)}\n")

        # Check if the heading is a version number (e.g., "4.1.0")
        if re.match(r'^\d+\.\d+\.\d+$', heading):
            github_link = f"https://github.com/COSMIC-PopSynth/COSMIC/releases/tag/v{heading}"
            lines.append(f"`GitHub release <{github_link}>`_\n")

        for feature in features:
            lines.append(feature)
        lines.append("")

    content = "\n".join(lines)
    if exists(output_path):
        with open(output_path) as f:
            if f.read() == content:
                print(f"  No changes to {output_path}, skipping write.")
                return

    print(f"  Writing 'What's New' content to {output_path}...")
    with open(output_path, "w") as f:
        f.write(content)

if __name__ == "__main__":
    changelog_path = Path("../changelog.md")
    output_path = Path("_generated/whats_new.rst")

    parsed_sections = parse_changelog(changelog_path)
    generate_rst(parsed_sections, output_path)
