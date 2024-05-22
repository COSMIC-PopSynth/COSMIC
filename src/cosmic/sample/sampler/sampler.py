# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWpy.
#
# GWpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""Fetch registration for database queries
"""

import re

from six import string_types

from astropy.io.registry import IORegistryError
from astropy.table import Table

_SAMPLERS = {}

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"


def register_sampler(data_format, data_class, function, force=False, usage=None):
    """Register a new method to InitialBinaryTable.sampler() for a given format

    Parameters
    ----------
    data_format : `str`
        name of the format to be registered

    data_class : `type`
        the class that the sampler returns

    function : `callable`
        the method to call from :meth:`InitialBinaryTable.sampler`

    force : `bool`, optional
        overwrite existing registration for ``data_format`` if found,
        default: `False`
    """
    key = (data_format, data_class)
    if key not in _SAMPLERS or force:
        _SAMPLERS[key] = (function, usage)
    else:
        raise IORegistryError(
            "Fetcher for format '{0}' and class '{1}' "
            "has already been "
            "defined".format(data_format, data_class)
        )
    _update__doc__(data_class)


def get_sampler(data_format, data_class):
    """Return the :meth:`~InitialBinaryTable.sampler` function for the given format

    Parameters
    ----------
    data_format : `str`
        name of the format

    data_class : `type`
        the class that the sampler returns

    Raises
    ------
    astropy.io.registry.IORegistryError
        if not registration is found matching ``data_format``
    """
    try:
        return _SAMPLERS[(data_format, data_class)][0]
    except KeyError:
        formats = "\n".join(_SAMPLERS.keys())
        raise IORegistryError(
            "No sampler definer for format %r. "
            "The available formats are:\n%r" % (data_format, formats)
        )


def _update__doc__(data_class):
    header = "The available named formats are:"
    sampler = data_class.sampler

    # if __doc__ isn't a string, bail-out now
    if not isinstance(sampler.__doc__, string_types):
        return

    # remove the old format list
    lines = sampler.__doc__.splitlines()
    try:
        pos = [i for i, line in enumerate(lines) if header in line][0]
    except IndexError:
        pass
    else:
        lines = lines[:pos]

    # work out the indentation
    matches = [re.search(r"(\S)", line) for line in lines[1:]]
    indent = min(match.start() for match in matches if match)

    # now re-write the format list
    formats = []
    for fmt, cls in sorted(_SAMPLERS, key=lambda x: x[0]):
        if cls is not data_class:
            continue
        usage = _SAMPLERS[(fmt, cls)][1]
        formats.append((fmt, "``sampler(%r, %s)``" % (fmt, usage)))
    format_str = Table(rows=formats, names=["Format", "Basic usage"]).pformat(
        max_lines=-1, max_width=80, align=(">", "<")
    )
    format_str[1] = format_str[1].replace("-", "=")
    format_str.insert(0, format_str[1])
    format_str.append(format_str[0])

    lines.extend([" " * indent + line for line in [header, ""] + format_str])
    # and overwrite the docstring
    try:
        sampler.__doc__ = "\n".join(lines)
    except AttributeError:
        sampler.__func__.__doc__ = "\n".join(lines)
