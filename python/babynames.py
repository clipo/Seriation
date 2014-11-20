__author__ = 'carllipo'

import os

for year in range(1910,2013,1):
    cmd = "python /Volumes/Macintosh\ HD/Users/carllipo/PycharmProjects/Seriation/python/IDSS.py --inputfile=/Volumes/Macintosh\ HD/Users/carllipo/PycharmProjects/Seriation/testdata/namesbyyear/" + str(year)+".txt --continuity=1 --noheader=0 --graphs=1  --debug=1"
    print cmd
    os.system(cmd)
