#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Katie Breivik (2017 - 2019)
#
# This file is part of COSMIC
#
# COSMIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COSMIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with COSMIC.  If not, see <http://www.gnu.org/licenses/>

# load tables
from astropy.table import (Column, Table)
from .initialbinarytable import InitialBinaryTable

# attach unified I/O
from . import sampler
