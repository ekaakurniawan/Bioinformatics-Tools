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

# Determination of the Starting Exon using Markov Model

DEBUG = False
PLOT = True

import re
import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import fsolve
from math import sqrt
import time

def normalize_matrix(m):
    m_new = []
    for i in m:
        m_new.append(np.divide(i, float(np.sum(i))))
    return m_new

def print_trans_matrix(m):
    for i in m:
        line = ""
        for j in i:
            line += ("%6d" % j)
        print line

def print_trans_prob_matrix(m):
    for i in m:
        line = ""
        for j in i:
            line += ("%8.4f" % j)
        print line

# Curate the sequences if the stop codons are within the sequence, remove them!
def fit_plus_model_file(fileName, a_plus_ij):
    with open(fileName, 'r') as p_file:
        for line in p_file:
            line = line.split()[0]
            a_plus_ij, seq_len = fit_plus_model_line(line, a_plus_ij)
    return a_plus_ij

def fit_plus_model_line(line, a_plus_ij):
    seq_len = 0
    # Check starting exon (ATG)
    st_exon = line[:3]
    if st_exon[:2] == "AT":
        a_plus_ij[0][1] += 1
    if st_exon[1:] == "TG":
        a_plus_ij[1][2] += 1

    seq_len += 3

    # Check start transition
    st_trans = line[2:4]
    if   st_trans == "GA":
        a_plus_ij[2][3] += 1
    elif st_trans == "GC":
        a_plus_ij[2][4] += 1
    elif st_trans == "GG":
        a_plus_ij[2][5] += 1
    elif st_trans == "GT":
        a_plus_ij[2][6] += 1

    # Check codons
    is_stop_codon = False
    is_complete_codon = True
    i = -1
    j = -1
    for codon in re.findall('..?.?', line[3:]):
        codon_len = len(codon)
        if codon_len < 3:
            is_complete_codon = False
            seq_len += codon_len
        else:
            seq_len += 3
        # Check backedges
        if i != -1:
            if   codon[0] == "A":
                j = 3
                a_plus_ij[i][j] += 1
            elif codon[0] == "C":
                j = 4
                a_plus_ij[i][j] += 1
            elif codon[0] == "G":
                j = 5
                a_plus_ij[i][j] += 1
            elif codon[0] == "T":
                j = 6
                a_plus_ij[i][j] += 1
        if codon_len == 1: break

        # Check stop codons (TAA, TAG, TGA)
        if (codon_len == 2) and (codon in ["TA", "TG"]):
            if   codon == "TA":
                a_plus_ij[6][11] += 1
            elif codon == "TG":
                a_plus_ij[6][12] += 1
            break
        if (codon_len == 3) and (codon in ["TAA", "TAG", "TGA"]):
            if   codon == "TAA":
                a_plus_ij[6][11] += 1
                a_plus_ij[11][17] += 1
                a_plus_ij[17][19] += 1
            elif codon == "TAG":
                a_plus_ij[6][11] += 1
                a_plus_ij[11][18] += 1
                a_plus_ij[18][19] += 1
            elif codon == "TGA":
                a_plus_ij[6][12] += 1
                a_plus_ij[12][17] += 1
                a_plus_ij[17][19] += 1
            is_stop_codon = True
            break

        # Check normal codons
        if   codon[0] == "A":
            i = 3
        elif codon[0] == "C":
            i = 4
        elif codon[0] == "G":
            i = 5
        elif codon[0] == "T":
            i = 6

        if   codon[1] == "A":
            if i == 6:
                j = 11
                a_plus_ij[i][j] += 1
                i = 11
            else:
                j = 7
                a_plus_ij[i][j] += 1
                i = 7
        elif codon[1] == "C":
            j = 8
            a_plus_ij[i][j] += 1
            i = 8
        elif codon[1] == "G":
            if i == 6:
                j = 12
                a_plus_ij[i][j] += 1
                i = 12
            else:
                j = 9
                a_plus_ij[i][j] += 1
                i = 9
        elif codon[1] == "T":
            j = 10
            a_plus_ij[i][j] += 1
            i = 10
        if codon_len == 2: break

        if   codon[2] == "A":
            j = 13
            a_plus_ij[i][j] += 1
            i = 13
        elif codon[2] == "C":
            j = 14
            a_plus_ij[i][j] += 1
            i = 14
        elif codon[2] == "G":
            j = 15
            a_plus_ij[i][j] += 1
            i = 15
        elif codon[2] == "T":
            j = 16
            a_plus_ij[i][j] += 1
            i = 16
    # Check end transition
    if is_stop_codon or (not is_complete_codon):
        pass
    else:
        a_plus_ij[i][19] += 1

    return a_plus_ij, seq_len

# Curate the sequence if the length is not a multiple of three.
# Use the same length for training and testing.
def fit_minus_model_file(fileName, a_minus_ij):
    with open(fileName, 'r') as p_file:
        for line in p_file:
            line = line.split()[0]
            fit_minus_model_line(line, a_minus_ij)
    return a_minus_ij

def fit_minus_model_line(line, a_minus_ij):
    i = -1
    j = -1
    # Check start transition
    st_trans = line[0]
    if   st_trans == "A":
        a_minus_ij[0][1] += 1
        i = 1
    elif st_trans == "C":
        a_minus_ij[0][2] += 1
        i = 2
    elif st_trans == "T":
        a_minus_ij[0][3] += 1
        i = 3
    elif st_trans == "G":
        a_minus_ij[0][4] += 1
        i = 4

    # Check codons
    for codon in re.findall('..?.?', line[3:]):
        if len(codon) < 3: break

        for dna in codon:
            if   dna == "A":
                a_minus_ij[i][1] += 1
                i = 1
            elif dna == "C":
                a_minus_ij[i][2] += 1
                i = 2
            elif dna == "T":
                a_minus_ij[i][3] += 1
                i = 3
            elif dna == "G":
                a_minus_ij[i][4] += 1
                i = 4

    a_minus_ij[i][5] += 1

    return a_minus_ij

def calc_log_odds_ratio(fileName, n_plus, a_plus_ij, n_minus, a_minus_ij):
    lor = []
    lor_n = []
    a_plus_ij_log = np.log(a_plus_ij)
    a_minus_ij_log = np.log(a_minus_ij)
    with open(fileName, 'r') as p_file:
        for line in p_file:
            line = line.split()[0]
            plus_m = [[0 for j in xrange(n_plus)] for i in xrange(n_plus)]
            minus_ma = [[0 for j in xrange(n_minus)] for i in xrange(n_minus)]
            # Calculate probability of sequence given both '+' and '-' models
            plus_m, plus_seq_len = fit_plus_model_line(line, plus_m)
            plus_m = np.multiply(plus_m, a_plus_ij_log)
            p_plus = np.sum(plus_m) / plus_seq_len
            minus_ma = np.multiply(fit_minus_model_line(line, minus_ma), \
                                   a_minus_ij_log)
            p_minus = np.sum(minus_ma) / len(line)
            # Calculate length-normalized log-oods ratio
            s = p_plus - p_minus
            lor.append(s)
    return lor

def plot_histogram(lor_plus_n, lor_minus_n, title):
    fig = plot.figure()
    ax = fig.add_subplot(1, 1, 1)
    # Plot histogram
    n, bins, patches = ax.hist([lor_plus_n, lor_minus_n], 40, \
                               color = ['b','r'], alpha = 0.75)
    # Get bin certers
    bincenters = 0.5 * (bins[1:] + bins[:-1])
    # Get curve fitting for '+' data and plot line graph
    y = np.polyfit(np.array(bincenters), np.array(n[0]), 20)
    f1 = np.poly1d(y)
    p1 = ax.plot(bincenters, f1(bincenters), 'b-o', label = 'plus data', \
                 linewidth = 1)
    # Get curve fitting for '-' data and plot line graph
    y = np.polyfit(np.array(bincenters), np.array(n[1]), 20)
    f2 = np.poly1d(y)
    p2 = ax.plot(bincenters, f2(bincenters), 'r-*', label = 'minus data', \
                 linewidth = 1)
    # Show details
    plot.legend()
    ax.grid(True)
    ax.set_title(title)
    ax.set_xlabel('Normalized Log-Odds Ratio')
    ax.set_ylabel('Number of Sequences')
    # Plot
    if PLOT:
        plot.show()

    return fsolve(lambda x : f1(x) - f2(x), 0.0)

def eval_plus(lor_plus_n, threshold):
    TP = len((np.array(lor_plus_n) >= threshold).nonzero()[0])
    FN = len((np.array(lor_plus_n) < threshold).nonzero()[0])
    return TP, FN

def eval_minus(lor_minus_n, threshold):
    TN = len((np.array(lor_minus_n) <= threshold).nonzero()[0])
    FP = len((np.array(lor_minus_n) > threshold).nonzero()[0])
    return TN, FP

def measure_performance(TP, FN, TN, FP):
    # Accuracy
    accuracy = (float(TP + TN) / (TP + TN + FP + FN))

    # Sensitivity
    sensitivity = (float(TP) / (TP + FN))

    # Specificity
    specificity = (float(TN) / (TN + FP))

    # Correlation Coefficient
    cc = (((TP * TN) - (FN * FP)) / \
          sqrt((TP + FN) * (TN + FP) * (TP + FP) * (TN + FN)))

    # Precision
    precision = (float(TP) / (TP + FP))

    # Recall
    recall = (float(TP) / (TP + FN))

    # F-measure
    f_measure = ((2 * precision * recall) / (precision + recall))

    return accuracy, sensitivity, specificity, cc, precision, recall, f_measure

def print_performance(train_TP, train_FN, train_TN, train_FP, \
                      test_TP, test_FN, test_TN, test_FP, \
                      train_accuracies, train_sensitivity, \
                      train_specificity, train_cc, \
                      train_precision, train_recall, train_f_measure, \
                      test_accuracies, test_sensitivity, \
                      test_specificity, test_cc, \
                      test_precision, test_recall, test_f_measure):

    print "            TP      FN      TN      FP"
    print " Train: %6d  %6d  %6d  %6d" % (train_TP, train_FN, train_TN, train_FP)
    print " Test : %6d  %6d  %6d  %6d" % (test_TP, test_FN, test_TN, test_FP)
    print "          Accuracies   Sensitivity   Specificity   Correlation   " \
          "  Precision        Recall     F-measure"
    print "                                                   Coefficient"
    print " Train: %12.2f  %12.2f  %12.2f  %12.2f  %12.2f  %12.2f  %12.2f" % \
        (train_accuracies, train_sensitivity, train_specificity, train_cc, \
         train_precision, train_recall, train_f_measure)
    print " Test : %12.2f  %12.2f  %12.2f  %12.2f  %12.2f  %12.2f  %12.2f" % \
        (test_accuracies, test_sensitivity, test_specificity, test_cc, \
         test_precision, test_recall, test_f_measure)

# ----------------------------------------------------------------- Training ---
tic = time.time() * 1000
# Total '+' states
n_plus = 20
# Transition matrix '+' model from i state to j state
a_plus_ij = [[1 for j in xrange(n_plus)] for i in xrange(n_plus)]
# Total '-' states
n_minus = 6
# Transition matrix '-' model from i state to j state
a_minus_ij = [[1 for j in xrange(n_minus)] for i in xrange(n_minus)]

# Train '+' model
a_plus_ij = fit_plus_model_file("./Data/train_plus.txt", a_plus_ij)
if DEBUG:
    print_trans_matrix(a_plus_ij)
a_plus_ij = normalize_matrix(a_plus_ij)
if DEBUG:
    print_trans_prob_matrix(a_plus_ij)

# Train '-' model
a_minus_ij = fit_minus_model_file("./Data/train_minus.txt", a_minus_ij)
if DEBUG:
    print_trans_matrix(a_minus_ij)
a_minus_ij = normalize_matrix(a_minus_ij)
if DEBUG:
    print_trans_prob_matrix(a_minus_ij)

# Calculate log-odds ratio '+' data
train_lor_plus_n = calc_log_odds_ratio("./Data/train_plus.txt", \
                                       n_plus, a_plus_ij, \
                                       n_minus, a_minus_ij)

# Calculate log-odds ratio '-' data
train_lor_minus_n = calc_log_odds_ratio("./Data/train_minus.txt", \
                                        n_plus, a_plus_ij, \
                                        n_minus, a_minus_ij)

# Plot training data using length-normalized log-odds ratio
# The plot function also returns threshold
threshold = plot_histogram(train_lor_plus_n, train_lor_minus_n, \
                           'Histogram of Training Data')
print "Threshold: %s" % threshold

# Get true positive and false negative
train_TP, train_FN = eval_plus(train_lor_plus_n, threshold)
# Get true negative and false positive
train_TN, train_FP = eval_minus(train_lor_minus_n, threshold)
if DEBUG:
    print len(train_lor_plus_n)
    print train_TP, train_FN
    print len(train_lor_minus_n)
    print train_TN, train_FP

# Train Data Performance
train_accuracies, train_sensitivity, train_specificity, train_cc, \
    train_precision, train_recall, train_f_measure  = \
    measure_performance(train_TP, train_FN, train_TN, train_FP)
if DEBUG:
    print train_accuracies, train_sensitivity, train_specificity, train_cc

toc = time.time() * 1000
print "Training Time: %4d milliseconds" % (toc - tic)

# ------------------------------------------------------------------ Testing ---
tic = time.time() * 1000
# Calculate log-odds ratio '+' data
test_lor_plus_n = calc_log_odds_ratio("./Data/test_plus.txt", \
                                      n_plus, a_plus_ij, \
                                      n_minus, a_minus_ij)

# Calculate log-odds ratio '-' data
test_lor_minus_n = calc_log_odds_ratio("./Data/test_minus.txt", \
                                       n_plus, a_plus_ij, \
                                       n_minus, a_minus_ij)

# Plot training data using length-normalized log-odds ratio
# The plot function also returns threshold
plot_histogram(test_lor_plus_n, test_lor_minus_n, 'Histogram of Testing Data')

# Get true positive and false negative
test_TP, test_FN = eval_plus(test_lor_plus_n, threshold)
# Get true negative and false positive
test_TN, test_FP = eval_minus(test_lor_minus_n, threshold)
if DEBUG:
    print len(test_lor_plus_n)
    print test_TP, test_FN
    print len(test_lor_minus_n)
    print test_TN, test_FP

# Test Data Performance
test_accuracies, test_sensitivity, test_specificity, test_cc, \
    test_precision, test_recall, test_f_measure = \
    measure_performance(test_TP, test_FN, test_TN, test_FP)
if DEBUG:
    print test_accuracies, test_sensitivity, test_specificity, test_cc

# Print performance
print_performance(train_TP, train_FN, train_TN, train_FP, \
                  test_TP, test_FN, test_TN, test_FP, \
                  train_accuracies, train_sensitivity, \
                  train_specificity, train_cc, \
                  train_precision, train_recall, train_f_measure, \
                  test_accuracies, test_sensitivity, \
                  test_specificity, test_cc, \
                  test_precision, test_recall, test_f_measure)

toc = time.time() * 1000
print "Testing Time: %4d milliseconds" % (toc - tic)

# ---------------------------------------------------------------------- ROC ---
tic = time.time() * 1000

min_lor = np.min(train_lor_plus_n + train_lor_minus_n + \
                 test_lor_plus_n + test_lor_minus_n)
max_lor = np.max(train_lor_plus_n + train_lor_minus_n + \
                 test_lor_plus_n + test_lor_minus_n)

train_1_specificities = []
train_sensitivities = []
test_1_specificities = []
test_sensitivities = []
for threshold in np.linspace(min_lor, max_lor, 200):
    # Train
    TP, FN = eval_plus(train_lor_plus_n, threshold)
    TN, FP = eval_minus(train_lor_minus_n, threshold)
    train_1_specificities.append(1 - (float(TN) / (TN + FP)))
    train_sensitivities.append((float(TP) / (TP + FN)))
    # Test Data
    TP, FN = eval_plus(test_lor_plus_n, threshold)
    TN, FP = eval_minus(test_lor_minus_n, threshold)
    test_1_specificities.append(1 - (float(TN) / (TN + FP)))
    test_sensitivities.append((float(TP) / (TP + FN)))

plot.plot(train_1_specificities, train_sensitivities, 'o-', linewidth = 2, \
          label = 'Training Data')
plot.plot(test_1_specificities, test_sensitivities, '*-', linewidth = 2, \
          label = 'Testing Data')
plot.title('ROC Space of Training and Testing Data')
plot.xlabel('FPR (1 - specificity)')
plot.ylabel('TPR (sensitivity)')
plot.legend()
plot.grid()
if PLOT:
    plot.show()

toc = time.time() * 1000
print "ROC Calculation Time: %4d milliseconds" % (toc - tic)

