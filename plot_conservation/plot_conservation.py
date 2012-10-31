# Copyright (C) 2012 by Eka A. Kurniawan
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

# Tested on:
#  - Python 2.7.3
#  - NumPy 1.6.2
#  - MatPlotLib 1.1.1

import numpy as np
import matplotlib.pyplot as plot

files = [['H1N1 - Avian - protein_conservation.txt',    'H1N1 - Avian'],
         ['H1N1 - Human - protein 1a_conservation.txt', 'H1N1 - Human 1'],
         ['H1N1 - Human - protein 1b_conservation.txt', 'H1N1 - Human 2'],
         ['H1N1 - Human - protein 2a_conservation.txt', 'H1N1 - Human 3'],
         ['H1N1 - Human - protein 2b_conservation.txt', 'H1N1 - Human 4'],
         ['H1N1 - Human - protein 3a_conservation.txt', 'H1N1 - Human 5'],
         ['H1N1 - Human - protein 3b_conservation.txt', 'H1N1 - Human 6'],
         ['H1N1 - Swine - protein_conservation.txt',    'H1N1 - Swine'],
         ['H3N2 - Avian - protein_conservation.txt',    'H3N2 - Avian'],
         ['H3N2 - Human - protein 1_conservation.txt',  'H3N2 - Human 1'],
         ['H3N2 - Human - protein 2_conservation.txt',  'H3N2 - Human 2'],
         ['H3N2 - Human - protein 3_conservation.txt',  'H3N2 - Human 3'],
         ['H3N2 - Swine - protein_conservation.txt',    'H3N2 - Swine'],
         ['H5N1 - Avian - protein_conservation.txt',    'H5N1 - Avian'],
         ['H5N1 - Human - protein_conservation.txt',    'H5N1 - Human'],
         ['H5N1 - Swine - protein_conservation.txt',    'H5N1 - Swine']]

conservations = []
totalFile = len(files)

for file in files:
    inFile = open(file[0], 'r')
    conservations.append(np.array(inFile.read().split(',')[:-1], \
                                  dtype = np.float))
    inFile.close()

plot.boxplot([np.asarray(cs) for cs in conservations])
plot.title('Conservation Box Plot of Different Viruses')
plot.ylabel('Score (0 to 11)')
plot.xticks(np.arange(totalFile + 1), [''] + [file[1] for file in files], \
            rotation = -90)
plot.show()
