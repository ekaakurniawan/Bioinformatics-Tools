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

# Genetic Algorithm Implementation for Optimizing Component Selection
# Suitable for integer constrain problem
# 
# References:
#  - http://en.wikipedia.org/wiki/Genetic_algorithm
#  - http://www.mathworks.com/videos/optimal-component-selection-using-the-mixed-integer-genetic-algorithm-68956.html

# Tested on:
#  - Python 2.7.3
#  - NumPy 1.6.2
#  - MatPlotLib 1.1.1

# TODO:
#  - Selection: Use random method to choose between elite and non-elite for
#               breeding (elite has higher chance)
#  - Breeding: Modify the script to produce two offsprings instead of one
#  - Migration: Adding migration event (fraction is 0.2 and interval is every
#               20 generations)

DEBUG = 1

import sys
from numpy import *
import matplotlib.pyplot as plot
import copy, time

def display(x0, y0, x, y):
    plot.scatter(x0, y0, s=20, color='blue')
    plot.scatter(x, y, s=5, color='red')
    plot.show()

def voltageCurve(Tdata, x, Res, ThVal, ThBeta):
    y = zeros(8)
    y[0] = Res[x[0]]
    y[1] = Res[x[1]]
    y[2] = Res[x[2]]
    y[3] = Res[x[3]]
    y[4] = ThVal[x[4]]
    y[5] = ThBeta[x[4]]
    y[6] = ThVal[x[5]]
    y[7] = ThBeta[x[5]]
    return tempCompCurve(y, Tdata)

def tempCompCurve(x, Tdata):
    # Imput voltage
    Vin = 1.1

    # Thermistor Calculations
    # x = R1, R2, R3, R4, RTH1(T_25degc), Beta1, RTH2(T_25degc), Beta2
    # RTH(T) = RTH(T_25degc) / exp(Beta * (T - T_25) / (T * T_25))
    T_25 = 298.15
    T_off = 273.15
    Beta1 = x[5]
    Beta2 = x[7]
    RTH1 = x[4] / exp(Beta1 * (((Tdata + T_off) - T_25) / ((Tdata + T_off) * T_25)))
    RTH2 = x[6] / exp(Beta2 * (((Tdata + T_off) - T_25) / ((Tdata + T_off) * T_25)))

    # Calculate equivalent circuits for parallel R's and RTH's
    R1_eq = (x[0] * RTH1) / (x[0] + RTH1)
    R3_eq = (x[2] * RTH2) / (x[3] + RTH2)

    # Calculate voltage at point B
    return Vin * ((R3_eq + x[3]) / (R1_eq + x[1] + R3_eq + x[3]))

def createComunity_random(pop_num, pop_size, gene_len, lb, ub):
    comunity = random.rand(pop_num, pop_size, gene_len)
    for x in xrange(pop_num):
        for y in xrange(pop_size):
            comunity[x][y][0:4] = lb[0] + (comunity[x][y][0:4] * (ub[0] - lb[0])).astype(int)
            comunity[x][y][4:6] = lb[4] + (comunity[x][y][4:6] * (ub[4] - lb[4])).astype(int)
    return comunity

def evaluateFitness(comunity, pop_num, pop_size, elite_cnt, \
                    Vdata, Tdata, Res, ThVal, ThBeta):
    # Calculate score
    score = zeros((pop_num, pop_size))
    V_len = len(Vdata)
    for x in xrange(pop_num):
        for y in xrange(pop_size):    
            Vnew = voltageCurve(Tdata, comunity[x][y], Res, ThVal, ThBeta)
            score[x][y] = sum(((Vdata - Vnew) ** 2) / V_len)
    # Get minimum
    elite_set = zeros((pop_num, elite_cnt, 6))
    elite_score = zeros((pop_num, elite_cnt))
    for x in xrange(pop_num):
        for indx, y in enumerate(score[x].argsort()[:3]):
            elite_set[x][indx] = comunity[x][y]
            elite_score[x][indx] = score[x][y]
    return elite_set, elite_score

def crossover_onePoint(P1, P2, gene_len, xover_frc):
    if random.rand() <= xover_frc:
        point = int(random.rand() * gene_len)
        P1[point:] = P2[point:]
    return P1

def mutation_uniform(C, gene_len, lb, ub):
    pt = int(random.rand() * gene_len)
    C[pt] = lb[pt] + (random.rand() * (ub[pt] - lb[pt])).astype(int)
    return C

def breed(elite_set, comunity, pop_num, pop_size, gene_len, elite_cnt, \
          xover_frc, lb, ub):
    for x in xrange(pop_num):
        for y in xrange(pop_size):
            parents = (random.rand(2) * elite_cnt).astype(int)
            child = crossover_onePoint(elite_set[x][parents[0]], \
                                       elite_set[x][parents[1]], \
                                       gene_len, xover_frc)
            comunity[x][y] = mutation_uniform(child, gene_len, lb, ub)
    return comunity

def main(args):
    Res = array([ 300,
                  330,
                  360,
                  390,
                  430,
                  470,
                  510,
                  560,
                  620,
                  680,
                  750,
                  820,
                  910,
                 1000,
                 1100,
                 1200,
                 1300,
                 1500,
                 1600,
                 1800,
                 2000,
                 2200,
                 2400,
                 2700,
                 3000,
                 3300,
                 3600,
                 3900,
                 4300,
                 4700,
                 5100,
                 5600,
                 6200,
                 6800,
                 7500,
                 8200,
                 9100,
                10000,
                11000,
                12000,
                13000,
                15000,
                16000,
                18000,
                20000,
                22000,
                24000,
                27000,
                30000,
                33000,
                36000,
                39000,
                43000,
                47000,
                51000,
                56000,
                62000,
                68000,
                75000,
                82000,
                91000,
               100000,
               110000,
               120000,
               130000,
               150000,
               160000,
               180000,
               200000,
               220000])
    ThVal = array([   50,
                     220,
                    1000,
                    2200,
                    3300,
                    4700,
                   10000,
                   22000,
                   33000])
    ThBeta = array([   2750,
                       3680,
                       3560,
                       3620,
                       3528,
                       3930,
                       3960,
                       4090,
                       3740])

    Tdata = r_[-40:86:5]
    Vdata = 1.026E-1 + -1.125E-4 * Tdata + 1.125E-5 * Tdata ** 2

    # Setting lower and upper bounds
    lb = array([0, 0, 0, 0, 0, 0])
    ub = array([70, 70, 70, 70, 9, 9])

    # GA Options:
    #  - Creation Function
    #  - Crossover Function                 : One Point
    #  - Crossover Fraction                 : xover_frc
    #  - Elite Count                        : elite_cnt
    #  - Fitness Limit -> Halt
    #  - Generations -> Halt                : gen_num
    #  - Initial Population                 : pop_num
    #  - Migration Direction
    #  - Migration Fraction
    #  - Migration Interval
    #  - Muation Function                   : Uniform
    #  - Output Function
    #  - Population Size                    : pop_size
    #  - Selection Function
    #  - Stall Generation Limit -> Halt
    #  - Time Limit -> Halt

    pop_num = 5
    pop_size = 100
    gen_num = 200
    elite_cnt = 3
    elite_set = zeros((pop_num, elite_cnt, 6))
    elite_score = zeros((pop_num, elite_cnt))
    xover_frc = 0.80

    topScore = float('inf')
    topCoeff = array([])
    # Initial populations
    comunity = createComunity_random(pop_num, pop_size, 6, lb, ub)

    print "Number of generations    = %s" % gen_num
    tic = time.time()
    for gen in xrange(gen_num):
        if DEBUG:
            print comunity[4][99] #bar

        # Evaluate fitness of each individual and return the best-fit set
        elite_set, elite_score = evaluateFitness(comunity, pop_num, pop_size, \
                                                 elite_cnt, Vdata, Tdata, Res, \
                                                 ThVal, ThBeta)

        # Display the minimum
        min_indx = unravel_index(elite_score.argmin(), elite_score.shape)
        topCurrentScore = elite_score[min_indx]
        if topCurrentScore < topScore:
            topScore = topCurrentScore
            topCoeff = copy.deepcopy(elite_set[min_indx])
            if DEBUG:
                Vnew = voltageCurve(Tdata, topCoeff, Res, ThVal, ThBeta)
                display(Tdata, Vdata, Tdata, Vnew)

        # Breed the best-fit individuals 
        comunity = breed(elite_set, comunity, pop_num, pop_size, 6, elite_cnt, \
                         xover_frc, lb, ub)

    # Final result
    print "Elapsed time             = %s" % (time.time() - tic)
    print "Final Score              = %s" % topScore
    print "Final Coefficient        = %s" % topCoeff
    Vnew = voltageCurve(Tdata, topCoeff, Res, ThVal, ThBeta)
    display(Tdata, Vdata, Tdata, Vnew)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
