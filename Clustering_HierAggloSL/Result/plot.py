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

# To run:
# python plot.py dataFile1.data out90000.txt

import sys
import string, re
import matplotlib.pyplot as plot

if len(sys.argv) != 3:
    print "Use: python plot.py points_file clusters_file"
    print "Example: python plot.py dataFile1.data out90000.txt"
    sys.exit()

inFile = open(sys.argv[1], 'r')
points_x = []
points_y = []
for line in inFile:
    point = re.findall(r'[0-9]+.[0-9]+', line)
    points_x.append(float(point[0]))
    points_y.append(float(point[1]))
inFile.close()

#(0, 0, 255),      //Blue
#(255, 0, 0),      //Red
#(0, 255, 0),      //Green
#(255, 255, 0),    //Yellow
#(255, 0, 255),    //Magenta
#(255, 128, 128),  //Pink
#(128, 128, 128),  //Gray
#(128, 0, 0),      //Brown
#(255, 128, 0),    //Orange
colors = ['#0000ff', '#ff0000', '#00ff00', \
          '#ffff00', '#ff00ff', '#000000', \
          '#808080', '#800000', '#ff8000']
inFile = open(sys.argv[2], 'r')
for i in xrange(7):
    inFile.readline()
clusters = re.findall(r'[0-9]+', inFile.readline())
inFile.close()
ttl_clusters = len(clusters)

# Convert cluster ID to color ID. Maximum 9 clusters.
cluster_val = []
for i in xrange(ttl_clusters):
    if not clusters[i] in cluster_val:
        cluster_val.append(clusters[i])

data = [[[] for i in xrange(3)] for i in xrange(len(cluster_val))]
for i in xrange(ttl_clusters):
    cluster_idx = cluster_val.index(clusters[i])
    data[cluster_idx][0].append(points_x[i])
    data[cluster_idx][1].append(points_y[i])
    data[cluster_idx][2].append(colors[cluster_idx])

for i in xrange(len(data)):
    plot.scatter(data[i][0], data[i][1], s = 3, color = data[i][2])
plot.show()
