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

import copy

# Seq1 = "CACCC"
# Seq2 = "CATC"

# Seq1 = "CAC"
# Seq2 = "CATC"

# Seq1 = "CTTGAT"
# Seq2 = "GCAT"

# Seq1 = "TCAAATCAACCAAGATGGAAGCAAAACTGTTTGTAC"
# Seq2 = "ATGAAGGCAATACTATTAGTCTTGCTATATACATTC"

Seq1 = "MEAKLFVLFCTFTVLKADTICVGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCSLNGIAPLQLGKCNVAGWLLGNPECDLLLTANSWSYIIETSNSENGTCYPGEFIDYEELREQLSSVSSFEKFEIFPKANSWPNHETTKGVTAACSYSGASSFYRNLLWITKKGTSYPKLSKSYTNNKGKEVLVLWGVHHPPTTSEQQSLYQNTDAYVSVGSSKYNRRFTPEIAARPKVRGQAGRMNYYWTLLDQGDTITFEATGNLIAPWYAFALNKGSDSGIITSDAPVHNCDTRCQTPHGALNSSLPFQNVHPITIGECPKYVKSTKLRMATGLRNVPSIQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNNLERRIENLNKKVDDGFLDVWTYNAELLVLLENERTLDFHDSNVRNLYEKVRSQLRNNAKELGNGCFEFYHKCDDECMESVKNGTYDYPKYSEESKLNREEIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI"
Seq2 = "MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGTSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"
# Sample output:
# MEAK--LFVLFCTF-TVLKADTICVGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCSLNGIAPLQLGKCNVAGW-LLGNPECD-LLL-TANSWSYII-ETSNSENGTCYPGEFIDYEELREQLSSVSSFEKFEIFPKANSWPNHETTKGVTAACSYSGA-SSFYRNLLWIT--KKGTSYPKLSKSYTNN-KGKEVLVLWGVHHPP-TTSE-QQSLYQNTDAYVSVG-SSKYNRRFTPEIAA-RPKVRGQAGRMNYYWTLL-DQGDTITFEATGNLIAPWYAFAL--NKGSDSGIITSDAPVHNCD--TRCQTPHGALN-SSLPFQNVHPITIGECPKYVKSTKLRMATGLRNVPSIQSRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAIDGITNKVNSVIEKMNTQFTAVGKEFNN-LE-RRIENLNKKVDDGFLDVWTYNAELLVLLENERTLDFHDSNVRNLYEKVRSQLRNNAKELGNGCFEFYHKCDDE-CMESVKNGTYDYPKYSEESKLNREEIDGVKLESMGVYQILAIYSTVASSLVLL-VSLGAISFWMCSNGSLQCRICI
# M--KAILVVLLYTFAT-ANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWI-LGNPECESL--STASSWSY-IVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKS-FYKNL--IWLVKKGNSYPKLSKSYIN-DKGKEVLVLWGIHHP-S-TSADQQSLYQNADAYVFVGT-SRYSKKFKPEIA-IRPKVRDQEGRMNYYWTL-VEPGDKITFEATGNLVVPRYAFAMERNAG--SGIIISDTPVH--DCNTTCQTPKGAINTS-LPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFN-HLEKR-IENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCD-NTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVL-VVSLGAISFWMCSNGSLQCRICI

mu = 1
sigma = 1

def getMaxScore(i, j):
    if V[i] == W[j]:
        m = S[i-1][j-1] + 1
    else:
        m = S[i-1][j-1] - mu

    scores = [S[i-1][j] - sigma,
              S[i][j-1] - sigma,
              m]
    return max(scores)

def calculatePrefix(source, sink, i):
    global V, W, S

    V = "0" + Seq1[source[0]:i + 1]
    W = "0" + Seq2[source[1]:sink[1]]

    lenV = len(V)
    lenW = len(W)

    S = [[0 for j in xrange(lenW)] for i in xrange(lenV)]
    for a in range(lenV): S[a][0] = a * -sigma
    for b in range(lenW): S[0][b] = b * -sigma

    for a in xrange(1, lenV):
        for b in range(1, lenW):
            S[a][b] = getMaxScore(a, b)
    return S[lenV - 1][1:lenW + 1]

def calculateSuffix(source, sink, i):
    global V, W, S

    V = "0" + Seq1[i:sink[0]][::-1]
    W = "0" + Seq2[source[1]:sink[1]][::-1]

    lenV = len(V)
    lenW = len(W)

    S = [[0 for j in xrange(lenW)] for i in xrange(lenV)]
    for a in range(lenV): S[a][0] = a * -sigma
    for b in range(lenW): S[0][b] = b * -sigma

    for a in xrange(1, lenV):
        for b in range(1, lenW):
            S[a][b] = getMaxScore(a, b)
    return S[lenV - 1][1:lenW + 1][::-1]

def getPath(source, sink):
    end = False
    if (sink[0] - source[0]) <= 2:
        if D[source[0]] == None:
            mid_i = source[0]
        elif D[source[0] + 1] == None:
            mid_i = source[0] + 1
        else:
            return
        end = True
    else:
        mid_i = source[0] + ((sink[0] - source[0]) / 2)

    prefix = calculatePrefix(source, sink, mid_i)
    suffix = calculateSuffix(source, sink, mid_i)

    sumScore = [prefix[b] + suffix[b] for b in xrange(sink[1] - source[1])]
    maxScore = max(sumScore)
    mid_k = source[1] + sumScore.index(maxScore)
    D[mid_i] = maxScore
    K[mid_i] = mid_k

    if end:
        return

    getPath(source, [mid_i + 1, mid_k + 1])
    getPath([mid_i, mid_k], sink)

def generateSequence():
    indx = 0
    k_indx = 0
    for i in xrange(0, n):
        if i in K[k_indx:]:
            total = sum([1 for j in K[k_indx:] if j == i])
            if total > 1:
                R[0] += [indx + j + 1 for j in xrange(total)]

                startIndx = k_indx + K[k_indx:].index(i)
                maxVal = max(D[startIndx:startIndx+total])
                R[1] += [indx + D[startIndx:startIndx+total].index(maxVal) + 1]
                
                indx += total
                k_indx += total
            else:
                R[0] += [indx + 1]
                R[1] += [indx + 1]
                indx += 1
                k_indx += 1
        else:
            R[1] += [indx + 1]
            indx += 1

def displaySequence():
    V = "0" + Seq1
    W = "0" + Seq2
    Vseq = ""
    Wseq = ""
    for indx in xrange(max(R[0] + R[1])):
        indx += 1
        if indx in R[0]:
            Vseq += V[R[0].index(indx)]
        else :
            Vseq += "-"
        if indx in R[1]:
            Wseq += W[R[1].index(indx)]
        else :
            Wseq += "-"
    print Vseq
    print Wseq
    print ""

m = len(Seq1)
n = len(Seq2)

S = []
V = ""
W = ""
D = [None for i in xrange(m)]
K = copy.deepcopy(D)
R = [[0], [0]]
getPath([0,0], [m,n])
generateSequence()
print R #bar
displaySequence()
