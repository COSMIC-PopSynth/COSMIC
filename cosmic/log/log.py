# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of hveto.
#
# hveto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hveto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hveto.  If not, see <http://www.gnu.org/licenses/>.

"""Utilities for logging output from Omicron in python
"""

import logging

import datetime

COLORS = dict((c, 30 + i) for i, c in enumerate(
    ['black', 'red', 'green', 'yellow',
     'blue', 'magenta', 'cyan', 'white']))
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"
LEVEL_COLORS = {
    'WARNING': COLORS['yellow'],
    'INFO': COLORS['white'],
    'DEBUG': COLORS['blue'],
    'CRITICAL': COLORS['red'],
    'ERROR': COLORS['red'],
}


class ColoredFormatter(logging.Formatter):
    """A `~logging.Formatter` that supports coloured output
    """
    def __init__(self, msg, use_color=True, **kwargs):
        logging.Formatter.__init__(self, msg, **kwargs)
        self.use_color = use_color

    def format(self, record):
        record.gpstime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        levelname = record.levelname
        if self.use_color and levelname in LEVEL_COLORS:
            record.levelname = color_text(levelname, LEVEL_COLORS[levelname])
        return logging.Formatter.format(self, record)


class Logger(logging.Logger):
    """`~logging.Logger` with a nice format
    """
    FORMAT = ('[{bold}%(name)s{reset} %(gpstime)d] %(levelname)+19s: '
              '%(message)s'.format(bold=BOLD_SEQ, reset=RESET_SEQ))
    def __init__(self, name, level=logging.DEBUG):
        try:
            super(Logger, self).__init__(name, level=level)
        except TypeError:
            logging.Logger.__init__(self, name, level=level)
        colorformatter = ColoredFormatter(self.FORMAT)
        console = logging.StreamHandler()
        console.setFormatter(colorformatter)
        self.addHandler(console)


def color_text(text, color):
    if not isinstance(color, int):
        color = COLORS[color]
    return COLOR_SEQ % color + str(text) + RESET_SEQ
