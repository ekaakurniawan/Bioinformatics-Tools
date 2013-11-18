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

# References:
#  - Miller, D.A.B., 2013. QMSE01: Quantum Mechanics for Scientists and
#    Engineers. [Online] Available at:
#    https://class.stanford.edu/courses/Engineering/QMSE01/Quantum_Mechanics_for_Scientists_and_Engineers/about
#    [Accessed 17 November 2013]

# Example:
"""
1. Calculating actual electric field in an infinitely high well with 1nm
   thickness knowing that the dimensionless electric field is 4 units.

from Electron import Electron
from InfinitePotentialWell import InfinitePotentialWell
electron = Electron()
# Well thickness (L): 1.0e-9 meter
# Energy eigenstate (n) = 1
well = InfinitePotentialWell(L = 1.0e-9, particle = electron, n = 1)
# Electric dimensionless field unit (f): 4
f = 4
# Get actual field (E) in volt per meter (V/m)
E = well.actual_field(f)
print '{0:.2e}'.format(E) # 1.50e+09 V/m
"""

from math import pi

class InfinitePotentialWell:
    # Planck's constant (h) in joule-second (J s) or m^2 kg / s
    H = 6.62606957e-34
    # Reduced Planck's constant or Dirac's constant or
    # Planck's constant over 2pi (h-bar) in joule-second (J s) or m^2 kg / s
    H_BAR = 1.054571726e-34
    # pi square
    PI_SQ = pi * pi

    def __init__(self, L = 1.0e-9, particle = None, n = 1, \
                 h = H, h_bar = H_BAR):
        # Well thickness in meter (m)
        self.L = L
        self.particle = particle
        # Energy eigenstate (integer starts from 1, dimensionless)
        self.n = n
        # Planck's constant (h) in joule-second (J s) or m^2 kg / s
        self.h = h
        # Reduced Planck's constant or Dirac's constant or
        # Planck's constant over 2pi (h-bar) in joule-second (J s) or m^2 kg / s
        self.h_bar = h_bar

    # Return eigenenergy or eigenvalue
    def eigenenergy(self, n = None):
        # Energy eigenstate (integer starts from 1, dimensionless)
        if n == None:
            n = self.n

        return ((self.h_bar ** 2) * (n ** 2) * (self.PI_SQ)) / \
               (2 * self.particle.mass * (self.L ** 2))

    # Unit energy of potential charge from one side to another (E_0)
    def unit_energy(self):
        return self.eigenenergy(n = 1) / (self.particle.e * self.L)

    # Dimensionless field unit (f)
    def dimensionless_field_unit(self, actual_field = None):
        # Actual field (E)
        return actual_field / self.unit_energy()

    # Actual field (E)
    def actual_field(self, dimensionless_field_unit = None):
        # Dimensionless field unit (f)
        if dimensionless_field_unit is not None:
            return dimensionless_field_unit * self.unit_energy()

