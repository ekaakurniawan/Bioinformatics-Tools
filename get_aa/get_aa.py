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

import sys
import string, re

if len(sys.argv) == 1:
    print "no input file"
    sys.exit()
else:
    inFileName = sys.argv[1]
    fileName = inFileName.split(".")[0]

inFile = open(inFileName, 'r')
for line in inFile:
    result = re.findall(r'[a-zA-Z]', line)
    if result:
        break
inFile.close()

consensus_aa = ""
outFileName = fileName + "_aa.fa"
outFile = open(outFileName, 'w')
outFile.write(">" + fileName.replace(" ", "_") + "\n" + \
              re.sub("(.{60})", "\\1\n", consensus_aa.join(result)))
outFile.close()
