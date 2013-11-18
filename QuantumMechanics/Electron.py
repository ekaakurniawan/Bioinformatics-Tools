# Copyright (C) 2013 by Eka A. Kurniawan
# eka.a.kurniawan(ta)gmail(tod)com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

class Electron:
    # Electron mass in kilogram (kg)
    MASS = 9.10938291e-31
    # Elementary charge in coulomb (C)
    ELEMENTARY_CHARGE = 1.602176565e-19

    # mass is in kg
    def __init__(self, mass = MASS, e = ELEMENTARY_CHARGE):
        # Electrom mass in kilogram (kg)
        self.mass = mass
        # Elementary charge in coulomb (C)
        self.e = e

