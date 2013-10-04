__author__ = 'carllipo'


import numpy as np
import networkx as nx
import pprint
import Dumper
from pylab import *
import matplotlib.pyplot as plt
import time
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic
import random
import scipy
import argparse
import datetime
import curses
import pprint
import csv

DEBUG = False
debug = 0
filterflag              = 0 ## do you want to try to output all of the solutions (filtered for non trivial)
largestonly             = 0 #  only output the largest set of solutions
individualfileoutput    = 0 ## create files for all the indivdual networks
bootstrapCI             = 0 ## flag for the CI bootstrap
bootstrapSignificance   = 99.5
man                     = 0
help                    = 0
inputfile = ''
threshold               = 0.5
noscreen                = 0     ## flag for screen output
excel                   = 0       ## flag for excel file output (not implemented yet)
xyfile                  = ""
mst                     = 0 ## minimum spanning tree
pairwiseFile            = ""
stats                   = 0 ## output stats, histograms of counts, etc
nosum                   = 0
allSolutions            = 0
memusage                = 0 ## memory usage flag


# start the clock to track how long this run takes
start = datetime.datetime.now()

# start prettyprint (python Dumper)
pp = pprint.PrettyPrinter(indent=4)

def main():
    parser = argparse.ArgumentParser(description='Conduct seriation analysis')
    parser.add_argument('--debug')
    parser.add_argument('--bootstrapCI')
    parser.add_argument('--bootstrapSignificance', type=int)
    parser.add_argument('--filtered')
    parser.add_argument('--largestonly')
    parser.add_argument('--individualfileoutput')
    parser.add_argument('--excel')
    parser.add_argument('--threshold')
    parser.add_argument('--noscreen')
    parser.add_argument('--xyfile')
    parser.add_argument('--pairwiseFile')
    parser.add_argument('--mst')
    parser.add_argument('--stats')
    parser.add_argument('--nosum')
    parser.add_argument('--allSolutions')
    parser.add_argument('--memusage')
    args = parser.parse_args()
    #pp.pprint(args)
    scr=curses.initscr()
    # Frame the interface area at fixed VT100 size
    scr.border(0)
    scr.addstr(12, 25, "Python curses in action!")
    #scr.refresh()
    #scr.getch()
    #screen.box()
    scr.addstr(14,21,"Step number:                     ")
    scr.addstr(5,5,"Bsdfsdfsdfsdfsdfsdfsdfsdfsdf      ")
    scr.refresh()
    #scr.getch()
    time.sleep(20)
if __name__ == "__main__":
    main()