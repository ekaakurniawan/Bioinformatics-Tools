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

# Likelihood-Ratio Test

# References:
#  - http://en.wikipedia.org/wiki/Likelihood-ratio_test

import numpy as np

# Number of pairs of two ungapped DNA sequences x and y
#         A   C   G   T
n_xy = [[25,  6,  7,  9],  # A
        [ 8, 10, 10,  6],  # C
        [ 7,  5, 30,  8],  # G
        [ 6,  9, 10, 15]]  # T

n =   np.sum([np.sum(nij) for nij in n_xy])
n_x = [np.sum(ni) for ni in n_xy]
n_y = [np.sum(nj) for nj in zip(*n_xy)]

sum_x = np.sum([ni * np.log(float(ni) / n) for ni in n_x])
sum_y = np.sum([nj * np.log(float(nj) / n) for nj in n_y])
sum_xy = np.sum([np.sum([nj * np.log(float(nj) / n) for nj in ni]) for ni in n_xy])

# Log-likelihood
S = (sum_x + sum_y) - sum_xy
# Test Statistic (equivalent to a chi-squared distribution)
D = -2 * S

print "Log-likelihood = %10.3f" % S
print "Test Statistic = %10.3f" % D

