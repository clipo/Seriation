__author__ = 'carllipo'

import os

## local path to project
path = "/Users/clipo/PycharmProjects/Seriation/"

for year in range(1911,2013,1):
    print "working on year:  ", year
    cmd = "python " +  path + "python/IDSS.py " \
          "--xyfile="+ path + "/testdata/namesbyyear/usstatesXY.txt " \
          "--inputfile="+path + "/testdata/namesbyyear/" + str(year)+".txt " \
          "--shapefile=1 --continuity=1 --noheader=0  --debug=1"
    print cmd
    os.system(cmd)
