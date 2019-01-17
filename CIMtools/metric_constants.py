# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#


# length
m = 1
sm = .01
mm = .001
dm = .1
km = 1000
inch = 2.54 * sm
foot = ft = 12 * inch
yard = 3 * foot
mile = 5280 * foot
nautical_mile = 1.852 * km

# squares
a = 100
ha = 100 * a
acre = 4046.86

# volumes
liter = litre = L = 1 * dm ** 3
pint = .56826125 * L
gallon = 8 * pint
bushel = bsh = bu = 8 * gallon
barrel = 36 * gallon

us_pint_liq = .473176473 * L
us_pint_dry = .5506104713575 * L
us_gallon_liq = 8 * us_pint_liq
us_gallon_dry = 8 * us_pint_dry
us_barrel_liq = 31.5 * us_gallon_liq
us_barrel_dry = 105 * us_gallon_dry / 4

beer_barrel = 31 * us_gallon_liq
oil_barrel = bbls = 42 * us_gallon_liq

# masses
kg = 1
g = .001
centner = 100
tonne = t = 1000
ounce = oz = 28.349523125 * g
troy_ounce = ozt = 31.1034768 * g
pound = lb = 16 * oz

# temperatures


class _Celsius:
    def __rmul__(self, other):
        return 273.15 + other

    def __rtruediv__(self, other):
        return other - 273.15


class _Fahrenheit:
    def __rmul__(self, other):
        return 273.15 + (other - 32) / 1.8

    def __rtruediv__(self, other):
        return (other - 273.15) * 1.8 + 32


K = 1
celsius = C = _Celsius()
fahrenheit = F = _Fahrenheit()

# times
second = s = 1
minute = M = 60
hour = h = 60 * M
day = 23 * hour + 56 * minute + 4 * second
year = 365.26 * day

# forces
gravity = 9.80665
newton = N = 1
kgf = gravity
lbf = gravity * lb

# pressures
pascal = Pa = 1
bar = 10 ** 5
at = kgf / sm ** 2
atm = 1.033 * at
m_water = .1 * at
mm_hg = 133.322
psi = lbf / inch ** 2
