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

# Block Alignment using Look Up Table - Divide and Conquer
# Sub-quadratic Time

# References:
#  - Neil C. Jones, Pavel A. Pevzner. An Introduction to Bioinformatics Algorithms. Cambridge: The MIT Press, 2004.

# TODO:
#  - Adding backtracking

import math

# Seq1 = "CACCC"
# Seq2 = "CATC"

# Seq1 = "CAC"
# Seq2 = "CATC"

# Seq1 = "CTTGAT"
# Seq2 = "GCAT"

# Seq1 = "TCAAATCAACCAAGATGGAAGCAAAACTGTTTGTAC"
# Seq2 = "ATGAAGGCAATACTATTAGTCTTGCTATATACATTC"

# Seq1 = "MEAKLFVLFCTFTVLKADTICVGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCSLNGIAPLQLGKCNVAGWLLGNPECDLLLTANSWSYIIETSNSENGTCYPGEFIDYEELREQLSSVSSFEKFEIFPKANSWPNHETTKGVTAACSYSGASSFYRNLLWITKKGTSYPKLSKSYTNNKGKEVLVLWGVHHPPTTSEQQSLYQNTDAYVSVGSSKYNRRFTPEIAARPKVRGQAGRMNYYWTLLDQGDTITFEATGNLIAPWYAFALNKGSDSGIITSDAPVHNCDTRCQTPHGALNSSLPFQNVHPITIGECPKYVKSTKLRMATGLRNVPSIQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNNLERRIENLNKKVDDGFLDVWTYNAELLVLLENERTLDFHDSNVRNLYEKVRSQLRNNAKELGNGCFEFYHKCDDECMESVKNGTYDYPKYSEESKLNREEIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI"
# Seq2 = "MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGTSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"

# Seq1 = "TTAAG"
# Seq2 = "AAGT"

Seq1 = "TCAAATCAACCAAGAT"
Seq2 = "ATGAAGGCAATACTAT"


BASES = "ATGC"

def getScoreLCS(i, j):
    if V[i] == W[j]:
        m = S[i-1][j-1] + 1
        return max([S[i-1][j],
                    S[i][j-1],
                    m])
    else:
        return max([S[i-1][j],
                    S[i][j-1]])

def convertIndxToBase(indx, card):
    asgn = [0 for x in xrange(card)]
    for x in xrange(card, 0, -1):
        asgn[x - 1] = indx % 4
        indx /= 4
    return "".join([BASES[x] for x in asgn])

def getLUTScore():
    for y in xrange(1, t + 1):
        for x in xrange(1, t + 1):
            S[y][x] = getScoreLCS(y, x)
    return S[y][x]

n = max(len(Seq1), len(Seq2))
t = int(math.ceil(math.log(n, 2) / 2))

lut_n = 4 ** t
lut = [[0 for j in xrange(lut_n)] for i in xrange(lut_n)]
S = [[0 for y in xrange(t + 1)] for x in xrange(t + 1)]
for i in xrange(lut_n):
    for j in xrange(lut_n):
        V = "0" + convertIndxToBase(i, t)
        W = "0" + convertIndxToBase(j, t)
        lut[i][j] = getLUTScore()

print "LUT:\n%s" % lut #bar

def getBlockScore(i_block, j_block):
    i = (i_block - 1) * t
    V = Seq1[i: i + t]
    j = (j_block - 1) * t
    W = Seq2[j: j + t]

    m = S[i_block-1][j_block-1] + lut[convertBaseToIndex(V)][convertBaseToIndex(W)]
    return max([S[i_block-1][j_block],
                S[i_block][j_block-1],
                m])

def convertBaseToIndex(base):
    indx = 0
    for x in xrange(len(base)):
        indx *= 4
        indx += BASES.index(base[x])
    return indx

a = n / t
b = n / t
S = [[0 for y in xrange(b + 1)] for x in xrange(a + 1)]
for i_block in xrange(1, b + 1):
    for j_block in xrange(1, a + 1):
        S[i_block][j_block] = getBlockScore(i_block, j_block)

print "Final Score:\n%s" % S #bar
