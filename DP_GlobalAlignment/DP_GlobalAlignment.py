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

# Global Alignment - Dynamic Programming

import copy

# Seq1 = "AAGC"
# Seq2 = "ATAGGC"

# Seq1 = "CCAA"
# Seq2 = "CA"

# Seq1 = "TCAAATCAACCAAGATGGAAGCAAAACTGTTTGTAC"
# Seq2 = "ATGAAGGCAATACTATTAGTCTTGCTATATACATTC"

Seq1 = "CTTGAT"
Seq2 = "GCAT"

V = "0" + Seq1
lenV = len(V)
W = "0" + Seq2
lenW = len(W)

S = [[0 for j in xrange(lenW)] for i in xrange(lenV)]
D = copy.deepcopy(S)
P = copy.deepcopy(S)

mu = 1
sigma = 2

t = 0

def getMaxScore(i, j):
    if V[i] == W[j]:
        m = S[i-1][j-1] + 1
    else:
        m = S[i-1][j-1] - mu

    scores = [S[i-1][j] - sigma,
              S[i][j-1] - sigma,
              m]
    maxScore = max(scores)
    
    direction = 0
    totalPath = 0
    if maxScore == scores[0]:
        direction += 1
        totalPath += P[i-1][j]
    if maxScore == scores[1]:
        direction += 2
        totalPath += P[i][j-1]
    if maxScore == scores[2]:
        direction += 4
        totalPath += P[i-1][j-1]

    return maxScore, direction, totalPath

def updateSequence(seq, p, indx, i, j):
    if seq == 0:
        R[p][seq][i] = indx + 1
    if seq == 1:
        R[p][seq][j] = indx + 1

def traceBack(i, j, p, indx, d):
    global t
    while (i is not 0) or (j is not 0):
        if d == 0: d = D[i][j]
        
        if d == 1:
            updateSequence(0, p, indx, i, j)
            i -= 1
            indx += 1
        if d == 2:
            updateSequence(1, p, indx, i, j)
            j -= 1
            indx += 1
        if d == 4:
            updateSequence(0, p, indx, i, j)
            i -= 1
            updateSequence(1, p, indx, i, j)
            j -= 1
            indx += 1
        if d == 3:
            traceBack(i, j, p, indx, 1)
            t += 1
            R[t] = copy.deepcopy(R[p])
            p = t
            traceBack(i, j, p, indx, 2)
            return
        if d == 5:
            traceBack(i, j, p, indx, 1)
            t += 1
            R[t] = copy.deepcopy(R[p])
            p = t
            traceBack(i, j, p, indx, 4)
            return
        if d == 6:
            traceBack(i, j, p, indx, 2)
            t += 1
            R[t] = copy.deepcopy(R[p])
            p = t
            traceBack(i, j, p, indx, 4)
            return
        if d == 7:
            traceBack(i, j, p, indx, 1)
            t += 1
            R[t] = copy.deepcopy(R[p])
            p = t
            traceBack(i, j, p, indx, 2)
            t += 1
            R[t] = copy.deepcopy(R[p])
            p = t
            traceBack(i, j, p, indx, 4)
            return
        
        d = 0

def swapSequence():
    for p in xrange(len(R)):
        maxIndx = max(R[p][0] + R[p][1])
        for indx in xrange(1, max([lenV, lenW])):
            if R[p][0][indx] is not 0:
               R[p][0][indx] = maxIndx - R[p][0][indx] + 1
            if R[p][1][indx] is not 0:
               R[p][1][indx] = maxIndx - R[p][1][indx] + 1

def displayDirection():
    for i in xrange(0, lenV):
        line = ""
        for j in range(0, lenW):
            if D[i][j] == 1: line += "  |, "
            if D[i][j] == 2: line += "  _, "
            if D[i][j] == 3: line += " _|, "
            if D[i][j] == 4: line += "  \, "
            if D[i][j] == 5: line += " \|, "
            if D[i][j] == 6: line += " _\, "
            if D[i][j] == 7: line += "_\|, "
        print line

def displaySequence():
    for p in xrange(len(R)):
        Vseq = ""
        Wseq = ""
        for indx in xrange(max(R[p][0] + R[p][1])):
            indx += 1
            if indx in R[p][0]:
                Vseq += V[R[p][0].index(indx)]
            else :
                Vseq += "-"
            if indx in R[p][1]:
                Wseq += W[R[p][1].index(indx)]
            else :
                Wseq += "-"
        print Vseq
        print Wseq
        print ""

for i in range(lenV):
    S[i][0] = i * -sigma
    D[i][0] = 1
    P[i][0] = 1

for j in range(lenW):
    S[0][j] = j * -sigma
    D[0][j] = 2
    P[0][j] = 1

for i in xrange(1, lenV):
    for j in range(1, lenW):
        S[i][j], D[i][j], P[i][j] = getMaxScore(i, j)

R = [[[0 for k in xrange(max([lenV, lenW]))] for j in xrange(2)] for i in xrange(P[lenV-1][lenW-1])]
traceBack(lenV-1, lenW-1, 0, 0, 0)
swapSequence()

print S #bar
print D #bar
displayDirection() #bar
print P #bar
print R #bar
displaySequence() #bar

