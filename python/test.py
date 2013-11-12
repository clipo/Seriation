__author__ = 'carllipo'

''''
An example of using the IDSS module within another python program.
Right now, the seriation module doesnt return anything useful other than True (done). But I suppose it could.
'''''

import IDSS

seriation = IDSS()

args={}
args['inputfile'] ="../testdata/pfg.txt"
args['xyfile']="../testdata/pfgXY.txt"
args['largestonly']=1
args['screen']=1
args['mst']=1

results=seriation.seriate(args)