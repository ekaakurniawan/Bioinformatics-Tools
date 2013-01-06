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


# -----------------------------------------------------------------------------
# Pseudo-Random Number Generator using Linear Feedback Shift Register (LFSR)
# -----------------------------------------------------------------------------
#
# This method carries the behaviour of counter which the value will be resetted
# back after certain number of iterrations. The longer the bit register is the
# lasting the randomness generated. Following calculations make sure that the
# random number will continue for more than 15 days of running using 200MHz
# clock frequency. Based on the result, we need at least 48-bit LFSR for such
# requirement.
#
# total15 = (200000000 x 60 x 60 x 24 x 15) + 1
#                |                      |
#                |                      +--------> days
#                +-------------------------------> clock frequency (Hz)
# total15 =  2.5920e+14
# log2(total15) = 47.881 = 48 bits
#
# Note:  - The formula is intended for FPGA implementation or any other devices
#          performing one LFSR operation per clock.
#        - The implementation follows Wikipedia version (using XOR) instead of
#          Xilinx version (using XNOR).
#
# References:
#  - http://en.wikipedia.org/wiki/Linear_feedback_shift_register
#  - http://www.xilinx.com/support/documentation/application_notes/xapp052.pdf
#
# Tested on:
#  - Python 2.7.3
#  - MatPlotLib 1.1.1
# -----------------------------------------------------------------------------

import matplotlib.pyplot as plot

# Selections: 4, 8, 16 or 48 bits
bits = 16

lfsr_seed = 0x1
lfsr = lfsr_seed
lfsrs = []
# 4 bits
if bits == 4:
    while True:
        bit = ((lfsr >> 0) ^ (lfsr >> 1)) & 1
        lfsr = (lfsr >> 1) | (bit << 3)
        lfsrs.append(lfsr)
        if lfsr == lfsr_seed:
            break
# 8 bits
if bits == 8:
    while True:
        bit = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 4)) & 1
        lfsr = (lfsr >> 1) | (bit << 7)
        lfsrs.append(lfsr)
        if lfsr == lfsr_seed:
            break
# 16 bits
if bits == 16:
    while True:
        bit = ((lfsr >> 0) ^ (lfsr >> 1) ^ (lfsr >> 3) ^ (lfsr >> 12)) & 1
        lfsr = (lfsr >> 1) | (bit << 15)
        lfsrs.append(lfsr)
        if lfsr == lfsr_seed:
            break
# 48 bits
if bits == 48:
    while True:
        bit = ((lfsr >> 0) ^ (lfsr >> 1) ^ (lfsr >> 27) ^ (lfsr >> 28)) & 1
        lfsr = (lfsr >> 1) | (bit << 47)
        lfsrs.append(lfsr)
        if lfsr == lfsr_seed:
            break
print lfsrs
print len(lfsrs)

x_axis = list(xrange(len(lfsrs)))
plot.title('Random Number Generator using %s-bit LFSR' % bits)
plot.ylabel('LFSR Result')
plot.xlabel('Iteration')
plot.scatter(x_axis, lfsrs, s = 1, color = 'blue')
plot.show()
