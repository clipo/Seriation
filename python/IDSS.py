#!/usr/bin/env python
# Copyright (c) 2013.  Carl P. Lipo <clipo@csulb.edu>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.
__author__ = 'carllipo'

## NOTE the library base for this program is the anaconda distribution of python. There are additional
## libraries needed:
##      svglib
##      reportlab
##      curses       for windows use:http://www.lfd.uci.edu/~gohlke/pythonlibs/#curses
##      svgwrite

import csv
from datetime import datetime
import operator
import argparse
import sys
import logging as logger
import itertools
import math
import random as rnd
import curses  # for windows use:http://www.lfd.uci.edu/~gohlke/pythonlibs/#curses
from itertools import chain
import traceback
import collections
import operator
import time
from datetime import datetime
import os
import re
import uuid
import pickle
import multiprocessing
import numpy as np
import scipy as sp
import scipy.stats
import networkx as nx
import networkx.algorithms.isomorphism as iso
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import xlsxwriter
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic
import MST
import shapefile
if os.name != 'nt':
    import memory
from frequencySeriationMaker import frequencySeriationMaker
import seriationEvaluation
from occurrenceSeriationMaker import occurrenceSeriationMaker



class IDSS():
    color = ["b", "r", "m", "y", "k", "w", (0.976, 0.333, 0.518), (0.643, 0.416, 0.894),
             (0.863, 0.66, 0.447), (0.824, 0.412, 0.118)]

    def __init__(self):
        self.pairGraph = nx.Graph(is_directed=False)
        self.nodeSizeFactor = 500 # factor to scale the size of the nodes when drawing
        self.solutionCount=0
        self.args={}
        self.inputFile = ""
        self.solutionCount=0
        self.outputDirectory = ""
        if os.name != "nt":
            self.mem = memory.Memory()
        self.start = time.time()
        self.assemblageSize = {}
        self.assemblageFrequencies = {}
        self.assemblages = {}
        self.countOfAssemblages = 0
        self.assemblageValues = {}
        self.labels = []
        self.numberOfClasses = 0
        self.maxSeriationSize = 0
        self.xAssemblage = {}
        self.yAssemblage = {}
        self.xyAssemblages = []
        self.distanceBetweenAssemblages = {}
        self.largestX = 0
        self.largestY = 0
        self.distanceBetweenAssemblages = {}
        self.validComparisonsHash = {}
        self.typeFrequencyLowerCI = {}
        self.typeFrequencyUpperCI = {}
        self.typeFrequencyMeanCI = {}
        self.pairwise = {}
        self.pairwiseError = {}
        self.sumOfDifferencesBetweenPairs = {}
        self.typeNames=[]
        self.FalseList=[None,0,False,"None","0","False"]
        logger.debug("Start time:  %s ", self.start)
        self.scr = None
        self.args={}
        self.totalAssemblageSize = 0 ## total for all assemblages
        self.solutionsChecked = 0 ## total # of seriations evaluated to find the final subset

    def saveGraph(self, graph, filename):
        nx.write_gml(graph, filename)

    def all_pairs(self, lst):
        return list((itertools.permutations(lst, 2)))

    def all_tuples(self, lst):
        tuples = list(itertools.combinations(lst, 3))
        useable_tuples = []
        for e in tuples:
            useable_tuples.append(e)
        return useable_tuples

    def openFile(self, filename):
        try:
            logger.debug("trying to open: %s ", filename)
            file = open(filename, 'rU')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", filename, e)
            sys.exit('file %s does not open: %s') % ( filename, e)

        reader = csv.reader(file, delimiter='\t', quotechar='|')
        values = []
        rowcount=0
        for row in reader:
            row = map(str, row)
            if rowcount==0 and self.args['noheader'] <> 1:
                rowcount=1
                row.pop(0)
                self.typeNames=row
            else:
                if len(row) > 1:
                    label = row[0]
                    self.labels.append(label)
                    row.pop(0)
                    row = map(float, row)
                    self.numberOfClasses = len(row)
                    freq = []
                    rowtotal = sum(row)
                    for r in row:
                        if self.args['occurrence'] not in self.FalseList:
                            if float(r)>0:
                                values.append(1)
                                freq.append(1)
                            else:
                                values.append(0)
                                freq.append(0)
                        else:
                            freq.append(float(float(r) / float(rowtotal)))
                            values.append(float(r))
                    self.assemblages[label] = freq
                    self.assemblageFrequencies[label] = freq
                    self.assemblageValues[label] = values
                    self.assemblageSize[label] = rowtotal
                    self.totalAssemblageSize += rowtotal
                    self.countOfAssemblages += 1
        self.maxSeriationSize = self.countOfAssemblages
        self.nodeSizeFactor *= self.countOfAssemblages
        return True

    def aggregateIdenticalAssemblages(self):
        '''function combines assemblages that are the same.
        It takes no parameters and uses arrays/dicts for the class
        '''
        tempList={}
        self.occurrenceSeriationList={}

        for e in self.assemblages:
            combined=""
            for val in self.assemblages[e]:
                combined += str(val)
            #print "combined:", combined
            self.occurrenceSeriationList[e]=combined
            val = tempList.get(combined, 0)
            if val==0:
                newval = []
                newval.append(e)
                tempList[combined]=newval
                #print "combined: ", combined, " value: ", tempList[combined]
            else:
                tempList[combined].append(e)
                #print "now combined: ", combined, "now value: ", tempList[combined]
        self.assemblages={}
        self.assemblageFrequencies = {}
        self.labels=[]

        for t in tempList:
            newval = str("/".join(tempList[t]))
            #print "t:", t
            newarray = []
            for s in t:
                #print "s:", s
                newarray.append(float(s))
            #print "new array:", newarray
            #print "NewKey:", newval, "Value:", t
            self.assemblages[newval]=newarray
            self.assemblageFrequencies[newval]=newarray
            self.labels.append(newval)
            self.assemblageSize[newval]=int(sum(newarray))

    def preCalculateSumOfDifferencesBetweenPairs(self):
        logger.debug("Precalculate differences between pairs")
        pairs = self.all_pairs(self.assemblages)
        for pair in pairs:
            diff = self.calculateSumOfDifferences(pair[0], pair[1])
            key1 = pair[0] + "*" + pair[1]
            key2 = pair[1] + "*" + pair[0]
            self.sumOfDifferencesBetweenPairs[key1] = diff
            self.sumOfDifferencesBetweenPairs[key2] = diff
            if diff == 0:
                logger.info("Potential problem:  %s and %s are identical in their frequencies. ", pair[0],pair[1])

    def preCalculateComparisons(self):
        logger.debug("Precalculating the comparisons between all pairs of assemblages...")
        pairs = self.all_pairs(self.assemblages)
        for pair in pairs:
            self.pairGraph.add_node(pair[0], name=pair[0])
            self.pairGraph.add_node(pair[1], name=pair[1])
            columns = range(len(self.assemblages[pair[0]]))
            ass1 = self.assemblages[pair[0]]
            ass2 = self.assemblages[pair[1]]
            comparison = ""
            for i in columns:
                val1 = ass1[i]
                val2 = ass2[i]
                logger.debug("\t\tComparing Assemblage: %s  and    Assemblage: %s  ########", pair[0], pair[1])
                logger.debug("\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,
                             i, i, i, i, i, i)

                if self.args['bootstrapCI'] not in self.FalseList:
                    upperCI_test = self.typeFrequencyUpperCI[pair[0]][i]
                    lowerCI_test = self.typeFrequencyLowerCI[pair[0]][i]
                    upperCI_end = self.typeFrequencyUpperCI[pair[1]][i]
                    lowerCI_end = self.typeFrequencyLowerCI[pair[1]][i]

                    if upperCI_test < lowerCI_end:
                        comparison += "D"
                    elif lowerCI_test > upperCI_end:
                        comparison += "U"
                    else:
                        comparison += "M"
                else:
                    if val1 < val2:
                        comparison += "D"
                    if val1 > val2:
                        comparison += "U"
                    if val1 == val2:
                        comparison += "M"
                logger.debug("Type %d: - comparison is: %s ", i, comparison[i])

            logger.debug("Comparison for %s and %s is: %s ", pair[0], pair[1], comparison)
            self.pairGraph.add_edge(pair[0], pair[1], weight=comparison)

    def openPairwiseFile(self, filename):
        logger.debug("Opening pairwise file %", filename)
        try:
            pw = open(filename, 'r')
        except csv.Error as e:
            sys.exit('pairwise file %s does not open: %s') % ( filename, e)

        reader = csv.reader(pw, delimiter='\t', quotechar='|')
        for row in reader:
            pair = row[0] + "#" + row[1]
            self.pairwise[pair] = row[2]
            self.pairwiseError[pair] = row[3]
        return True

    def openXYFile(self, filename):
        logger.debug("Opening XY file %s", filename)
        ## open the xy file
        try:
            xyf = open(filename, 'r')
        except csv.Error as e:
            sys.exit('file %s does not open: %s') % ( filename, e)

        reader = csv.reader(xyf, delimiter='\t', quotechar='|')

        for row in reader:
            label = row[0]
            self.xyAssemblages.append(label)
            self.yAssemblage[label] = float(row[1])
            self.xAssemblage[label] = float(row[2])

        assemblagePairs = self.all_pairs(self.xyAssemblages)
        ## Go through all of the combinations
        for combo in assemblagePairs:
            pairname = combo[0] + "*" + combo[1]
            xdistance = float(self.xAssemblage[combo[0]]) - float(self.xAssemblage[combo[1]])
            xxdistance = xdistance * xdistance
            ydistance = float(self.yAssemblage[combo[0]]) - float(self.yAssemblage[combo[1]])
            yydistance = ydistance * ydistance
            distance = math.sqrt(xxdistance + yydistance)
            self.distanceBetweenAssemblages[pairname] = distance
        largestXname = max(self.xAssemblage.iterkeys(), key=(lambda key: self.xAssemblage[key]))
        largestYname = max(self.yAssemblage.iterkeys(), key=(lambda key: self.yAssemblage[key]))
        self.largestX = self.xAssemblage[largestXname]
        self.largestY = self.yAssemblage[largestYname]
        return True


    #############################################  THRESHOLD DETERMINATION ####################################
    ## first compare each assemblage and see if the threshold is exceeded.
    ## By this I mean the maximum difference between frequencies of any Type %ds greater than what is specified.
    ## When threshold = 0, all combinations are used. Later the "ends" of solutions are not evaluated if the
    ## difference between the last assemblage and any other free assemblage is > the threshold.
    ## This arbitrary setting is to keep from arbitrary solutions being stuck on that come from the "ends"
    ## of solutions.
    ##
    ## Precalculate all of the max differences between types in assemblage pairs.

    def thresholdDetermination(self, threshold):
        assemblageComparison = {}
        ##  get all the combinations of 2
        pairs = self.all_pairs(self.assemblages)

        ## Go through all of the combinations
        for combo in pairs:
            logger.debug("comparing combination of %s and %s ", combo[0], combo[1])
            pairname1 = combo[0] + "*" + combo[1]
            pairname2 = combo[1] + "*" + combo[0]
            maxDifference = 0
            assemblage1 = self.assemblages[combo[0]]
            assemblage2 = self.assemblages[combo[1]]
            i = -1
            columns = len(self.assemblages[combo[0]])
            logger.debug("Number of columns: %d", columns)
            ## calculate the maximum differences between the pairs of assemblages (across all types)
            for i in (0, columns - 1):
                logger.debug("i: %d ", i)
                ass1 = float(assemblage1[i])
                ass2 = float(assemblage2[i])
                diff = abs(ass1 - ass2)
                logger.debug("assemblage1: %f assemblage2: %f diff: %f", ass1, ass2, diff)
                if diff > maxDifference:
                    maxDifference = diff
            assemblageComparison[pairname1] = maxDifference
            assemblageComparison[pairname2] = maxDifference

        ############## pre calculate the valid pairs of assemblages and stuff into hash ############################
        for assemblage1 in self.assemblages:
            cAssemblages = []
            for assemblage2 in self.assemblages:
                if not assemblage1 == assemblage2:
                    testpair = assemblage1 + "*" + assemblage2
                    logger.debug("Pairs: %s and %s", assemblage1, assemblage2)
                    logger.debug("Comp value of pairs: %s:  %f and threshold is: %f", testpair,
                                 assemblageComparison[testpair], threshold)
                    if assemblageComparison[testpair] <= threshold:
                        logger.debug("Appending %s to the list of valid comparisons for %s ", assemblage1, assemblage2)
                        cAssemblages.append(assemblage2)
            self.validComparisonsHash[assemblage1] = cAssemblages
        return True

    def confidence_interval(self, data, confidence=0.05):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * sp.stats.t._ppf((1 + confidence) / 2., n - 1)
        return m, m - h, m + h

    ########################################### BOOTSTRAP CI SECTION ####################################
    def bootstrapCICalculation(self, bootsize=1000, confidenceInterval=0.05):

        if self.args['screen']:
            self.scr.addstr(1, 40, "STEP: Bootstrap CIs...        ")
            self.scr.refresh()

        ## for each assemblage
        logger.debug("Calculating bootstrap confidence intervals")
        # for each assemblage

        for currentLabel in sorted(self.assemblages.iterkeys()):
            assemblage = self.assemblages[currentLabel]
            types = len(self.assemblages[currentLabel])
            currentAssemblageSize = self.assemblageSize[currentLabel]

            ## create an array of arrays - one for each type
            arrayOfStats = []
            for c in range(0, types):
                array = []
                arrayOfStats.append([])

            ## size of bootstrapping (how many assemblages to create)
            loop = bootsize
            for counter in range(0, bootsize):

                assemsize = currentAssemblageSize
                # clear and set the array
                cumulate = []
                for d in range(0, types):
                    cumulate.append(0.0)

                index = 0
                count = 0
                ## now count through the classes and set up the frequencies
                for typeFrequency in assemblage:
                    index += typeFrequency
                    cumulate[count] = index  ## this is the index of the frequency for this class
                    ## index should be total # of types at end
                    count += 1

                ## set new_assemblage
                new_assemblage = []
                for c in range(0, types):
                    new_assemblage.append(0.0)

                for sherd in range(0, int(currentAssemblageSize)):
                    rand = random()             ## random number from 0-1
                    classVar = 0
                    typeIndex = types - 1
                    found = 0
                    total = sum(cumulate)
                    for t in reversed(cumulate):
                        if rand <= t:
                            found = typeIndex
                        typeIndex -= 1
                    new_assemblage[found] += 1

                ## count new assemblage frequencies
                counter = 0
                new_assemblage_freq = []
                newassemblage_size=sum(new_assemblage)
                for g in new_assemblage:
                    new_assemblage_freq.append(float(g / float(newassemblage_size)))
                    arrayOfStats[counter].append(float(g / float(newassemblage_size)))
                    counter += 1
                    ## this should result in a new assemblage of the same size

            lowerCI = []
            upperCI = []
            meanCI = []
            for freq in arrayOfStats:
                upper = 0.0
                lower = 0.0
                mean = 0.0
                if sum(freq) > 0.0:
                    mean, lower, upper = self.confidence_interval(freq, confidence=float(confidenceInterval))
                else:
                    mean = lower = upper = 0
                if math.isnan(lower) is True:
                    lower = 0.0
                if math.isnan(upper) is True:
                    upper = 0.0
                if math.isnan(mean) is True:
                    mean = 0.0
                lowerCI.append(lower)
                upperCI.append(upper)
                meanCI.append(mean)

            self.typeFrequencyLowerCI[currentLabel] = lowerCI
            self.typeFrequencyUpperCI[currentLabel] = upperCI
            self.typeFrequencyMeanCI[currentLabel] = meanCI

        return True

    ########################################### FIND ALL THE VALID TRIPLES  ####################################
    ########################################### #################################### ###########################
    def findAllValidTriples(self):
        triples = []
        error = 0
        numberOfTriplets = 0

        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Find valid triples....      ")
            self.scr.refresh()
        permutations = self.all_tuples(self.assemblages)

        for permu in permutations:
            if self.args['screen'] not in self.FalseList:
                c = self.scr.getch()
                if c == ord('q'):
                    curses.endwin()
                    curses.resetty()
                    sys.exit("Quitting as requested.\n\r")
            logger.debug("Triple test: %s * %s * %s", permu[0], permu[1], permu[2])

            comparison12 = ""
            comparison23 = ""
            error = 0
            columns = len(self.assemblages[permu[0]])
            logger.debug("Columns: %d", columns)
            difscore = 0
            difscore2 = 0
            comparison12 = ""
            comparison23 = ""
            for i in range(0, columns):
                ass1 = self.assemblages[permu[0]][i]
                ass2 = self.assemblages[permu[1]][i]
                ass3 = self.assemblages[permu[2]][i]
                logger.debug("ass1: %f ass2: %f ass3: %f", ass1, ass2, ass3)
                if self.args['bootstrapCI'] not in self.FalseList:
                    low1 = self.typeFrequencyLowerCI[permu[0]][i]
                    low2 = self.typeFrequencyLowerCI[permu[1]][i]
                    low3 = self.typeFrequencyLowerCI[permu[2]][i]
                    high1 = self.typeFrequencyUpperCI[permu[0]][i]
                    high2 = self.typeFrequencyUpperCI[permu[1]][i]
                    high3 = self.typeFrequencyUpperCI[permu[2]][i]
                    #mean1 = self.typeFrequencyMeanCI[permu[0]][i]
                    #mean2 = self.typeFrequencyMeanCI[permu[1]][i]
                    #mean3 = self.typeFrequencyMeanCI[permu[2]][i]

                    # compare 1 and 2
                    if high1 < low2:
                        comparison12 = "D"
                    elif low1 > high2:
                        comparison12 = "U"
                    else:
                        comparison12 = "M"

                    ## compare 2 and 3
                    if high2 < low3:
                        comparison23 = "U"
                    elif low2 > high3:
                        comparison23 = "D"
                    else:
                        comparison23 = "M"

                    ## check for problems
                    if comparison12 == "D" and comparison23 == "D":
                        comparison12 += "X"
                        comparison23 += "X"
                    elif comparison12 == "U" and comparison23 == "U":
                        error += 1

                else:     ### not bootstrap
                    if ass1 < ass2 < ass3:
                        comparison12 += "U"
                        comparison23 += "U"
                    elif ass1 < ass2 > ass3:
                        comparison12 += "X"
                        comparison23 += "X"
                    elif ass1 < ass2 == ass3:
                        comparison12 += "U"
                        comparison23 += "M"
                    elif ass1 > ass2 < ass3:
                        comparison12 += "D"
                        comparison23 += "U"
                        error += 1
                    elif ass1 > ass2 > ass3:
                        comparison12 += "D"
                        comparison23 += "D"
                    elif ass1 > ass2 == ass3:
                        comparison12 += "D"
                        comparison23 += "M"
                    elif ass1 == ass2 == ass3:
                        comparison12 += "M"
                        comparison23 += "M"
                    elif ass1 == ass2 > ass3:
                        comparison12 += "M"
                        comparison23 += "D"
                    elif ass1 == ass2 < ass3:
                        comparison12 += "M"
                        comparison23 += "U"
                    else:
                        logger.debug(
                            "\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f" % ass1, ass2,
                            ass3)
                        print "\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f \n\r" % ass1, ass2, ass3
                        print "I must quit. Debugging required.\n\r"
                        sys.exit()

                    logger.debug("Comparison12: %s Comparison23: %s", comparison12, comparison23)

            comparison = comparison12 + comparison23
            test = re.compile('DU').search(comparison)
            if test not in self.FalseList:
                error += 1

            if error == 0:
                # uses NetworkX
                graphID = uuid.uuid4().urn
                net = nx.Graph(name=str(numberOfTriplets), GraphID=graphID, End1=permu[0], End2=permu[2],
                               Middle=permu[1], is_directed=False)
                net.add_node(permu[0], name=permu[0], site="end", end=1, connectedTo=permu[1])
                net.add_node(permu[1], name=permu[1], site="middle", end=0, connectedTo="middle")
                net.add_node(permu[2], name=permu[2], site="end", end=1, connectedTo=permu[1])
                net.add_edge(permu[0], permu[1], weight=comparison12, GraphID=graphID, end=1)
                net.add_edge(permu[2], permu[1], weight=comparison23, GraphID=graphID, end=1)
                logger.debug("VALID TRIPLE SOLUTION: %s * %s * %s ", permu[0], permu[1], permu[2])
                logger.debug("VALID TRIPLE SOLUTION: %s  <--->   %s", comparison12, comparison23)
                logger.debug("VALID TRIPLE SOLUTION: %s ", net.adjacency_list())
                path = nx.shortest_path(net, source=permu[0], target=permu[2])
                logger.debug("VALID TRIPLE SOLUTION: Ends are: %s and %s", permu[0], permu[2])
                logger.debug("VALID TRIPLE SOLUTION: Shortest Path: %s ", path)
                triples.append(net)
                numberOfTriplets += 1
                logger.debug("Current number of triplets: %d", numberOfTriplets)
        return triples

    def filter_list(self, full_list, excludes):
        ''' This function filters two lists to find uniques
        :param a: A list to be compared
        :param b: A list to be compared
        :returns: a list of the things that are unique
        '''
        s = set(excludes)
        return (x for x in full_list if x not in s)

    def iso(self, G1, glist):
        """Quick and dirty nonisomorphism checker used to check isomorphisms."""
        for G2 in glist:
            if isomorphic(G1, G2):
                return True
        return False

    def MST(self, sGraph, filename):

        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family']='sans-serif'
        plt.figure(filename, figsize=(8, 8))
        M = nx.minimum_spanning_tree(sGraph)

        os.environ["PATH"] += ":/usr/local/bin:"
        pos = nx.graphviz_layout(M)
        #pos=nx.graphviz_layout(M,prog="twopi",root=self.args['graphroot'])
        edgewidth = []
        weights = nx.get_edge_attributes(M, 'weight')
        for w in weights:
            edgewidth.append(weights[w])
        maxValue = max(edgewidth)
        widths = []
        for w in edgewidth:
            widths.append(((maxValue - w) + 1) * 5)
        assemblageSizes = []
        sizes = nx.get_node_attributes(M, 'size')
        for s in sizes:
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(M, pos, alpha=0.3, width=widths)
        sizes = nx.get_node_attributes(M, 'size')
        nx.draw_networkx_nodes(M, pos, node_size=assemblageSizes, node_color='w', alpha=0.4)
        nx.draw_networkx_edges(M, pos, alpha=0.4, node_size=0, width=1, edge_color='k')
        nx.draw_networkx_labels(M, pos, fontsize=10,font_family='sans-serif', fontname='Helvetica')
        font = {'fontname': 'Helvetica',
                'color': 'k',
                'fontweight': 'bold',
                'font_family':'sans-serif',
                'fontsize': 10}
        plt.axis('off')
        plt.savefig(filename, dpi=75)
        self.saveGraph(sGraph, filename + ".gml")
        if self.args['shapefile'] is not None and self.args['xyfile'] is not None:
            self.createShapefile(M, self.outputDirectory + filename + ".shp")

    def iso(self, G1, glist):
        """Quick and dirty nonisomorphism checker used to check isomorphisms."""
        for G2 in glist:
            if isomorphic(G1, G2):
                return True
        return False

    def sumGraphsByWeight(self, filteredarray):
        sumGraph = nx.Graph(is_directed=False)
        # First add all the nodes to the sumgraph
        for node in self.assemblages:
            xCoordinate = 0
            yCoordinate = 0
            name = node
            if self.args['xyfile'] is not None:
                xCoordinate = self.xAssemblage[name]
                yCoordinate = self.yAssemblage[name]
            sumGraph.add_node(name, name=name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                              size=self.assemblageSize[name]/self.totalAssemblageSize*self.nodeSizeFactor)

        ## first find global max weight and global min weight
        globalMaxWeight=0
        globalMinWeight=1000
        for g in filteredarray:
            for e in g.edges():
                fromAssemblage = e[0]
                toAssemblage = e[1]
                currentWeight = self.sumOfDifferencesBetweenPairs[fromAssemblage+"*"+toAssemblage]
                if currentWeight>globalMaxWeight:
                    globalMaxWeight=currentWeight
                if currentWeight<globalMinWeight:
                    globalMinWeight=currentWeight

        ## Now create the summary graph by going through the edges of all the graphs
        for g in filteredarray:
            ## go through all the edges for each graph
            for e in g.edges_iter():
                d = g.get_edge_data(*e)
                # get the pair of data
                fromAssemblage = e[0]
                toAssemblage = e[1]
                currentWeight = self.sumOfDifferencesBetweenPairs[fromAssemblage+"*"+toAssemblage]
                normalizedWeight = ((globalMaxWeight-currentWeight)/((globalMaxWeight-globalMinWeight)+1))+1
                sumGraph.add_path([fromAssemblage, toAssemblage], sumDiffWeight=currentWeight, weight=normalizedWeight, inverseweight=(1/normalizedWeight))

        return sumGraph

    def sumGraphsByCount(self, filteredarray):
        sumGraph = nx.Graph(is_directed=False)
        ## go through all the graphs
        for g in filteredarray:
            ## go through all the edges for each graph
            for node in g.nodes(data=True):
                xCoordinate = 0
                yCoordinate = 0
                name = node[0]
                if self.args['xyfile'] is not None:
                    xCoordinate = self.xAssemblage[name]
                    yCoordinate = self.yAssemblage[name]
                sumGraph.add_node(name, name=name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                                  size=self.assemblageSize[name]/self.totalAssemblageSize*self.nodeSizeFactor)

            maxWeight = 0
            for e in g.edges_iter():
                d = g.get_edge_data(*e)
                fromAssemblage = e[0]
                toAssemblage = e[1]
                exists = False
                currentWeight = 1
                for e in sumGraph.edges():
                    dd = sumGraph.get_edge_data(*e)
                    if fromAssemblage in e and toAssemblage in e:   ## if exists
                        exists = True
                    currentWeight = 1
                    if exists is True:
                        currentWeight += 1

                if currentWeight > maxWeight:
                    maxWeight = currentWeight
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight)

            for e in sumGraph.edges_iter():
                d = sumGraph.get_edge_data(*e)
                currentWeight = int(d['weight'])
                inverseWeight = (maxWeight + 1) - currentWeight
                fromAssemblage = e[0]
                toAssemblage = e[1]
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight, inverseweight=inverseWeight)

        return sumGraph


    def finalGoodbye(self):
        if self.args['screen'] is not None:
            curses.endwin()
            curses.resetty()
            curses.nl()
            curses.echo()
            os.system("reset")

    #################################################### set up all the output files ####################################################
    def setupOutput(self):
        outputFile = self.outputDirectory + self.inputFile[0:-4] + ".vna"
        OUTMSTFILE = OUTMSTDISTANCEFILE = ""
        try:
            OUTFILE = open(outputFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outputFile, e
            sys.exit(msg)

        outpairsFile = self.outputDirectory + self.inputFile[0:-4] + "-pairs.vna"
        try:
            OUTPAIRSFILE = open(outpairsFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outpairsFile, e
            sys.exit(msg)

        outmstFile = self.outputDirectory + self.inputFile[0:-4] + "-mst.vna"
        outmst2File = self.outputDirectory + self.inputFile[0:-4] + "-mst-distance.vna"

        if self.args['mst'] not in self.FalseList:
            try:
                OUTMSTFILE = open(outmstFile, 'w')
                OUTMSTDISTANCEFILE = open(outmst2File, 'w')
            except csv.Error as e:
                msg = "Can't open file %s to write: %s" % outputFile, e
                sys.exit(msg)

        return OUTFILE, OUTPAIRSFILE, OUTMSTFILE, OUTMSTDISTANCEFILE

    def createShapefile(self, graph, shapefilename):
        w = shapefile.Writer(shapefile.POLYLINE)  # 3= polylines
        #turn geometry/attribute autoBalanace on
        #w.autoBalance = 1
        # add a field, any field
        w.field('FIELD','C','1')
        w.autoBalance = 1
        xCoordinates = nx.get_node_attributes(graph, "xCoordinate")
        yCoordinates = nx.get_node_attributes(graph, "yCoordinate")
        num=0
        lineparts=[]
        for e in graph.edges_iter():
            num += 1
            d = graph.get_edge_data(*e)
            node1 = e[0]
            node2 = e[1]
            x1 = float(xCoordinates[node1])
            y1 = float(yCoordinates[node1])
            x2 = float(xCoordinates[node2])
            y2 = float(yCoordinates[node2])
            lineparts.append([[x1,y1],[x2,y2]])
            #w.poly(parts=[[[x1, y1], [x2, y2]]])
        w.line(parts=lineparts)
        w.record('')
        # create the PRJ file
        #prj = open("%s.prj" % shapefilename[0:-4], "w")
        #epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
        #prj.write(epsg)
        #prj.close()
        w.save(shapefilename)
        w = shapefile.Writer(shapefile.POINT)
        w.field('Label','C','40')
        for n,d in graph.nodes_iter(data=True):
            w.point(float(d['xCoordinate']),float(d['yCoordinate']))
            w.record(d['name'],'Point')
        w.save(shapefilename[0:-4]+"-points.shp")

    def iso(self, G1, glist):
        """Quick and dirty nonisomorphism checker used to check isomorphisms."""
        for G2 in glist:
            if isomorphic.is_isomorphic(G1, G2):
                return True
        return False

    def createAtlasOfSolutions(self, filteredarray, type):
        plt.rcParams['font.family']='sans-serif'
        plt.figure(self.inputFile[0:-4] + "-" + str(type) + "-atlas.png", figsize=(8, 8))
        num = 0
        for g in filteredarray:
            t = 0
            for n in g.nodes():
                g.node[n]['label'] = g.node[n]['name'] ##+" ("+str(num)+")"
                #print "label: ", g.node[n]['label']
                t += 1
            num += 1
        UU = nx.Graph(is_directed=False)
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist = self.iso_filter_graphs(filteredarray) # list of nonisomorphic graphs
        atlasGraph = nx.disjoint_union_all(nlist)
        os.environ["PATH"] += ":/usr/local/bin:"
        pos = nx.graphviz_layout(atlasGraph, prog="twopi")
        #labels=nx.draw_networkx_labels(filteredarray,pos)
        C = nx.connected_component_subgraphs(atlasGraph)
        for g in C:
            c = [random()] * nx.number_of_nodes(g)
            #nodes, names = zip(*nx.get_node_attributes(g, 'label').items())
            names = nx.get_node_attributes(g, 'label')
            #print names
            nx.draw(g,
                    pos,
                    node_size=40,
                    node_color=c,
                    vmin=0.0,
                    vmax=1.0,
                    alpha=.2,
                    font_family='sans-serif',
                    fontname='Helvetica',
                    font_size=7,
                    with_labels=True,
                    labels=names
            )
        atlasFile = self.outputDirectory + self.inputFile[0:-4] + "-" + str(type) + "-atlas.png"
        plt.savefig(atlasFile, dpi=250)
        #plt.show() # display

    def outputExcel(self, filteredarray, filename, type):
        csv.register_dialect('excel_tab', delimiter='\t',lineterminator="\n")
        textFileName = filename +"-"+type+"-seriations.txt"
        f=open(textFileName, 'wb')
        writer=csv.writer(f,'excel_tab')
        excelFileName= filename + "-" + type + ".xlsx"
        workbook = xlsxwriter.Workbook(excelFileName)
        worksheet = workbook.add_worksheet()
        row = 0
        worksheet.write(row, 0, 'Seriation_Number')
        worksheet.write(row, 1, 'Assemblage')
        outputRow =[]
        outputRow.append('Seriation_Number')
        outputRow.append('Assemblage')
        if self.args['noheader'] in (1,True,"yes"):
            for type in range(2, self.numberOfClasses + 2):
                typename = "Type_" + str(type - 1)
                worksheet.write(row, type, typename)
                outputRow.append(typename)
        else:
            colcount=2
            for typename in self.typeNames:
                worksheet.write(row, colcount, typename)
                outputRow.append(typename)
                colcount +=1
        writer.writerow(outputRow)
        sernum = 0
        for g in filteredarray:
            col = 0
            sernum += 1
            row += 1
            for node in nx.shortest_path(g, g.graph['End1'], g.graph['End2']):
                outputRow=[]
                worksheet.write(row, 0, sernum)
                outputRow.append(sernum)
                worksheet.write(row, 1, node)
                outputRow.append(node)
                col = 2
                for a in self.assemblageFrequencies[node]:
                    val = int(self.assemblageSize[node] * a)
                    worksheet.write(row, col, val)
                    outputRow.append(val)
                    col += 1
                writer.writerow(outputRow)
                #print outputRow

            row += 1
            writer.writerow('')

        workbook.close()
        return excelFileName,textFileName

    def createAtlas(self, filteredarray):
        # remove isolated nodes, only connected graphs are left
        U = nx.Graph(is_directed=False) # graph for union of all graphs in atlas
        for G in filteredarray:
            U = nx.disjoint_union(U, G)
        # list of graphs of all connected components
        C = nx.connected_component_subgraphs(U)

        UU = nx.Graph(is_directed=False)
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist = [] # list of nonisomorphic graphs
        for G in C:
            # check against all nonisomorphic graphs so far
            if not self.iso(G, nlist):
                nlist.append(G)
                UU = nx.disjoint_union(UU, G) # union the nonisomorphic graphs
        return UU

    def calculateSumOfDifferences(self, assemblage1, assemblage2):
        diff = 0
        for type in range(0, self.numberOfClasses):
            diff += pow((float(self.assemblageFrequencies[assemblage1][type]) - float(
                self.assemblageFrequencies[assemblage2][type])),2)
        return pow(diff,0.5)


    def iso_filter_graphs(self, list):
        newlist = []
        compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
        for g in list:
            addList = True
            for h in newlist:
                #print "g: ",g.nodes(),"--- h: ",h.nodes(), "----> ", compare(g.nodes(),h.nodes())
                if compare(g.nodes(), h.nodes()) == True:
                    addList = False
            if addList == True:
                newlist.append(g)

        return newlist

    def continunityMaximizationSeriation(self):
        graphList = []
        numGraphs = 0

        ## special case for the first time through
        ## set up all the initial pairs of two closest assemblages.
        for ass in self.assemblages:
            numGraphs += 1
            g = nx.Graph(startAssemblage=ass, End1=ass, is_directed=False)
            g.add_node(ass, name=ass, size=self.assemblageSize[ass], xCoordinate=self.xAssemblage[ass],
                       yCoordinate=self.yAssemblage[ass])
            minMatch = 10
            newNeighbor = ""
            ## now find the smallest neighbor from the rest of the assemblages.
            for potentialNeighbor in self.assemblages:
                if potentialNeighbor is not ass:
                    diff = self.calculateSumOfDifferences(potentialNeighbor, ass)
                    if diff < minMatch:
                        minMatch = diff
                        newNeighbor = potentialNeighbor
            g.add_node(newNeighbor, name=newNeighbor, xCoordinate=self.xAssemblage[newNeighbor],
                       yCoordinate=self.xAssemblage[newNeighbor], size=self.assemblageSize[newNeighbor])
            if minMatch == 0:  ## prevent divide by zero errors
                minMatch = 10
            g.add_path([ass, newNeighbor], weight=minMatch, inverseweight=(1 / minMatch ))
            g.graph['End2'] = newNeighbor
            graphList.append(g)   ## create a starting graph for each of assemblage put into an array
        ## filter list so that we just have non-isomorphic graphs (1<->2 and not 2<->1).
        filtered_graph_list = self.iso_filter_graphs(graphList)
        ## Now go through list looking at each one and increasing as long as I can. Add graphs when there are equivalent solutions
        for current_graph in filtered_graph_list:
            globalMinMatch = 10
            endMinMatch = {"End1": 10, "End2": 10}
            currentMinimumMatch = {}
            matchEnd = ""
            matchEndAssemblage = {}   ## this will contain the end assemblages and the differences
            #print "Now on GRAPH: ", current_graph

            ## go through this the # of times of the assemblages -2 (since the network already has 2 nodes)
            for count in range(0, len(self.assemblages) - 2):
                ## examine both ends to see which is the smallest summed difference.
                matchEndAssemblage = {}
                matchEnd = ""
                currentMinimumMatch = {}
                endMinMatch = {"End1": 10, "End2": 10}
                for assEnd in ("End1", "End2"):
                    ## set the current end assemblages
                    endAssemblage = current_graph.graph[assEnd]
                    ## go through the other assemblages (not already in the graph.
                    for assemblage in self.assemblages:
                        if assemblage not in current_graph.nodes():
                            #print "assemblage: ", assemblage, " is not in : ", current_graph.nodes()
                            diff = self.calculateSumOfDifferences(endAssemblage, assemblage)
                            if diff < endMinMatch[assEnd]:
                                endMinMatch[assEnd] = diff
                                currentMinimumMatch[assEnd] = assemblage
                                matchEndAssemblage[assEnd] = endAssemblage
                                matchEnd = assEnd

                ## at this point we should have the minimum distance match for each end.
                ## we then need to compare each end to find which one is the smallest
                ## three possibilities -- end1, end2 and both (i.e., the diff is the same)
                smallestMatchEnd = []
                assemblagesMatchedToEnd = []
                #print "endminmatch-end1: ", endMinMatch['End1'], "  endminmatch-end2: ", endMinMatch['End2']
                if endMinMatch['End1'] < endMinMatch['End2']:
                    globalMinMatch = endMinMatch['End1']
                    smallestMatchEnd.append('End1')
                    assemblagesMatchedToEnd.append(currentMinimumMatch['End1'])
                    #print "new match is to end1:  ", currentMinimumMatch['End1']

                elif endMinMatch['End2'] < endMinMatch['End1']:
                    globalMinMatch = endMinMatch['End2']
                    smallestMatchEnd.append('End2')
                    assemblagesMatchedToEnd.append(currentMinimumMatch['End2'])
                    #print "new match is to end2:  ", currentMinimumMatch['End1']
                elif endMinMatch['End1'] < 10 and endMinMatch['End2'] < 10:
                    #print endMinMatch['End2'], "<--", endMinMatch['End1']
                    #print "matchEnd: ", matchEnd
                    #print currentMinimumMatch['End1'], "---", currentMinimumMatch['End2']
                    globalMinMatch = endMinMatch['End1']
                    smallestMatchEnd.append('End1')
                    smallestMatchEnd.append('End2')
                    assemblagesMatchedToEnd.append(currentMinimumMatch[matchEnd])

                ## find out if there are others that have the same minimum value
                for b in self.assemblages:
                    if b not in current_graph.nodes() and b not in assemblagesMatchedToEnd:
                        diff = self.calculateSumOfDifferences(b, endAssemblage)
                        if diff == globalMinMatch:
                            ## add this as a matched equivalent assemblage. We will then deal with more than one match
                            assemblagesMatchedToEnd.append(b)
                loop = 1
                firstOne = True

                original_network = current_graph.copy()
                for match in assemblagesMatchedToEnd:
                    for endAss in smallestMatchEnd:
                        # for the first time we need to simply add it to the right end but after this we copy...
                        if firstOne == True:
                            firstOne = False
                            current_graph.add_node(match, name=match, xCoordinate=self.xAssemblage[match],
                                                   yCoordinate=self.xAssemblage[match],
                                                   size=self.assemblageSize[match])
                            if globalMinMatch == 0:
                                globalMinMatch = 10
                            current_graph.add_path([match, matchEndAssemblage[endAss]], weight=globalMinMatch,
                                                   inverseweight=(1 / globalMinMatch ))
                            current_graph.graph[endAss] = match

                        ## if there are more than one we need to copy first before adding node
                        else:
                            loop += 1
                            #print "Loop: ", loop
                            new_network = original_network.copy()
                            new_network.add_node(match, name=match, xCoordinate=self.xAssemblage[match],
                                                 yCoordinate=self.xAssemblage[match], size=self.assemblageSize[match])
                            if globalMinMatch == 0:
                                globalMinMatch = 10
                            new_network.add_path([matchEndAssemblage[endAss], match], weight=globalMinMatch,
                                                 inverseweight=(1 / globalMinMatch ))
                            new_network.graph[endAss] = match
                            filtered_graph_list.append(new_network)
                            numGraphs += 1
                            #print "Number of graphs: ", numGraphs, " -- ", len(filtered_graph_list)

        return filtered_graph_list

    ## Output to file and to the screen
    def graphOutput(self, sumGraph, sumgraphfilename):

        ## Now make the graphic for set of graphs
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family']='sans-serif'
        newfilename = self.outputDirectory + sumgraphfilename
        gmlfilename = self.outputDirectory + sumgraphfilename + ".gml"
        self.saveGraph(sumGraph, gmlfilename)
        if self.args['shapefile'] is not None and self.args['xyfile'] is not None:
            self.createShapefile(sumGraph, self.outputDirectory + newfilename[0:-4] + ".shp")
        plt.figure(newfilename, figsize=(8, 8))
        os.environ["PATH"] += ":/usr/local/bin:"
        pos = nx.graphviz_layout(sumGraph)
        edgewidth = []

        ### Note the weights here are biased where the *small* differences are the largest (since its max value - diff)
        weights = nx.get_edge_attributes(sumGraph, 'weight')
        for w in weights:
            edgewidth.append(weights[w])
        maxValue = max(edgewidth)
        widths = []
        for w in edgewidth:
            widths.append(((maxValue - w) + 1) * 5)

        assemblageSizes = []
        sizes = nx.get_node_attributes(sumGraph, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s]/self.totalAssemblageSize*self.nodeSizeFactor)
        nx.draw_networkx_edges(sumGraph, pos, alpha=0.3, width=widths)
        sizes = nx.get_node_attributes(sumGraph, 'size')
        nx.draw_networkx_nodes(sumGraph, pos, node_size=assemblageSizes, font_family='sans-serif',node_color='w', alpha=0.4)
        nx.draw_networkx_edges(sumGraph, pos, alpha=0.4, node_size=0, width=1, edge_color='k')
        nx.draw_networkx_labels(sumGraph, pos, font_family='sans-serif', fontname='Helvetica', fontsize=10)
        font = {'fontname': 'Helvetica',
                'font_family':'sans-serif',
                'color': 'k',
                'fontweight': 'bold',
                'fontsize': 10}
        plt.axis('off')
        plt.savefig(newfilename, dpi=75)
        self.saveGraph(sumGraph, newfilename + ".gml")


    ## Output to file and to the screen
    def sumGraphOutput(self, sumGraph, sumgraphfilename):

        nodeList = sumGraph.nodes()
        for a in self.assemblages:
            if a not in nodeList:
                sumGraph.add_node(a, name=a, xCoordinate=self.xAssemblage[a], yCoordinate=self.yAssemblage[a],
                                  size=self.assemblageSize[a]/self.totalAssemblageSize*self.nodeSizeFactor)
        sumgraphOutputFile = self.outputDirectory + sumgraphfilename + ".vna"
        SUMGRAPH = open(sumgraphOutputFile, 'w')
        SUMGRAPH.write("*Node data\n")
        SUMGRAPH.write("ID AssemblageSize X Y Easting Northing\n")

        for node in sumGraph.nodes(data=True):
            nodeName = node[0]
            x = 0
            y = 0
            northing = 0
            easting = 0
            if self.args['xyfile'] is not None:
                x = float(self.xAssemblage[nodeName]) / 1000000.0
                y = (float(self.largestY) - float(self.yAssemblage[nodeName])) / 100000.0
                easting = self.xAssemblage[nodeName]
                northing = self.yAssemblage[nodeName]
            msg = nodeName + " " + str(self.assemblageSize[nodeName]) + " " + str(x) + " " + str(y) + " " + str(
                easting) + " " + str(northing) + "\n"
            SUMGRAPH.write(msg)
        SUMGRAPH.write("*Tie data\nFrom To Edge Weight InverseWeight\n")
        edgeCount = 0
        for e in sumGraph.edges_iter():
            d = sumGraph.get_edge_data(*e)
            edgeCount += 1
            text = e[0] + " " + e[1] + " " + str(edgeCount) + " " + str(d['weight']) + " " + str(
                d['inverseweight']) + "\n"
            SUMGRAPH.write(text)

        ## Now make the graphic for the sumgraph
        newfilename = self.outputDirectory + sumgraphfilename + "-weight.png"
        self.saveGraph(sumGraph, sumgraphfilename + ".gml")
        plt.figure(newfilename, figsize=(8, 8))
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family']='sans-serif'
        os.environ["PATH"] += ":/usr/local/bin:"
        pos = nx.graphviz_layout(sumGraph, prog="neato")

        edgewidth = []
        weights = nx.get_edge_attributes(sumGraph, 'weight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths = []
        for w in edgewidth:
            widths.append(((maxValue - w) + 1) * 5)

        assemblageSizes = []
        sizes = nx.get_node_attributes(sumGraph, 'size')
        #print sizes
        for s in sizes:
            assemblageSizes.append(sizes[s]/self.totalAssemblageSize*self.nodeSizeFactor)
        nx.draw_networkx_edges(sumGraph, pos, alpha=0.3, width=widths)
        sizes = nx.get_node_attributes(sumGraph, 'size')
        nx.draw_networkx_nodes(sumGraph, pos, node_size=assemblageSizes, node_color='w', alpha=0.4)
        nx.draw_networkx_edges(sumGraph, pos, alpha=0.4, node_size=0, font_family='sans-serif',width=1, edge_color='k')
        nx.draw_networkx_labels(sumGraph, pos,font_family='sans-serif',fontname='Helvetica', fontsize=10)
        font = {'fontname': 'Helvetica',
                'color': 'k',
                'fontweight': 'bold',
                'font_family':'sans-serif',
                'fontsize': 10}
        plt.axis('off')
        plt.savefig(newfilename, dpi=75)

        if self.args['shapefile'] is not None and self.args['xyfile'] is not None:
            self.createShapefile(sumGraph, self.outputDirectory + sumgraphfilename + "-weight.shp")


    #################################################### OUTPUT SECTION ####################################################
    def output(self, filteredArray, OUTFILE, OUTPAIRSFILE, OUTMSTFILE, OUTMSTDISTANCEFILE, maxEdges):
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(13, 1, "Now printing output file... ")
            self.scr.addstr(1, 40, "STEP: Output files...         ")
            self.scr.refresh()

        OUTFILE.write("*Node data\n")
        OUTFILE.write("ID AssemblageSize X Y Easting Northing\n")
        OUTPAIRSFILE.write("*Node data\n")
        OUTPAIRSFILE.write("ID AssemblageSize X Y Easting Northing\n")
        count = 0
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Printing list of nodes....     ")
            self.scr.refresh()
        ## note this assumes the use of UTM coordinates (northing and easting)
        for l in self.assemblages:
            x = 0
            y = 0
            northing = 0
            easting = 0
            if self.args['xyfile'] not in self.FalseList and self.args['occurrence'] in self.FalseList:
                x = float(self.xAssemblage[l]) / 1000000.0
                y = (float(self.largestY) - float(self.yAssemblage[l])) / 100000.0
                easting = self.xAssemblage[l]
                northing = self.yAssemblage[l]

            msg = l + " " + str(self.assemblageSize[l]) + " " + str(x) + " " + str(y) + " " + str(easting) + " " + str(
                northing) + "\n"
            OUTFILE.write(msg)
            OUTPAIRSFILE.write(msg)
            if self.args['mst'] not in self.FalseList:
                OUTMSTFILE.write(msg)
                OUTMSTDISTANCEFILE.write(msg)

        OUTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
        OUTPAIRSFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

        if self.args['mst'] not in self.FalseList:
            OUTMSTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
            OUTMSTDISTANCEFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Printing list of nodes attributes... ")
            self.scr.refresh()
        for l in self.assemblages:
            easting = 0
            northing = 0
            x = 0
            y = 0
            if self.args['xyfile'] not in self.FalseList:
                x = float(self.xAssemblage[l]) / 1000000
                y = (float(self.largestY) - float(self.yAssemblage[l])) / 100000
                easting = self.xAssemblage[l]
                northing = self.yAssemblage[l]
            msg = l + " " + str(self.assemblageSize[l]) + " " + str(x) + " " + str(y) + " " + str(easting) + " " + str(
                northing) + "\n"
            OUTFILE.write(msg)
            OUTPAIRSFILE.write(msg)
            if self.args['mst'] not in self.FalseList:
                OUTMSTFILE.write(msg)
                OUTMSTDISTANCEFILE.write(msg)

        ## This prints out counts of the edges as they appear in ALL of the solutions
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Going through and counting pairs...     ")
            self.scr.refresh()
        OUTPAIRSFILE.write("*Tie data\nFrom To Edge Count\n")
        if self.args['mst'] not in self.FalseList:
            OUTMSTFILE.write("*Tie data\nFrom To Edge End Weight ID\n")
            OUTMSTDISTANCEFILE.write("*Tie data\nFrom To Edge End Weight ID\n")

        ## first count up all of the edges by going through the solutions and each edge
        ## put the edge count in a hash of edges
        edgeHash = {}

        for network in filteredArray:
            for e in network.edges_iter():
                pairname = e[0] + "*" + e[1]
                edgeHash[pairname] = 0

            for e in network.edges_iter():
                pairname = e[0] + "*" + e[1]
                edgeHash[pairname] += 1

        ## now go through the edgeHash and print out the edges
        ## do this is sorted order of the counts. For fun.
        if self.args['screen'] is not None:
            self.scr.addstr(1, 40, "STEP: Doing the pair output...                ")
            self.scr.refresh()

        sorted_pairs = sorted(edgeHash.iteritems(), key=operator.itemgetter(1))

        for key, value in sorted_pairs:
            ass1, ass2 = key.split("*")
            msg = ass1 + "\t" + ass2 + "\t" + " 1 " + str(value) + "\n"
            OUTPAIRSFILE.write(msg)

        OUTFILE.write("*Tie data\nFrom To Edge Weight Network End pValue pError meanSolutionDistance\n")
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Eliminating duplicates...     ")
            self.scr.addstr(1, 40, "STEP: Printing edges...     ")
            self.scr.refresh()

        uniqueArray = set(filteredArray)
        distanceHash = {}
        seriationHash = {}
        ## only print unique ones...
        pairwise = {}
        pairwiseError = {}
        for network in filteredArray:
            if self.args['screen'] not in self.FalseList:
                self.scr.addstr(14, 1, "Now on solution: ")
                self.scr.addstr(14, 18, str(network.graph["GraphID"]))
                #print "now on solution: ", network["GraphID"],"\n"
            if self.args['largestonly'] not in self.FalseList and len(network.edges()) == maxEdges - 1:
                edgeCount = len(network.edges())
                groupDistance = 0
                meanDistance = 0.0
                eCount = 0
                if self.args['xyfile'] not in self.FalseList:
                    for e in network.edges_iter():
                        pairname = e[0] + "*" + e[1]
                        groupDistance += self.distanceBetweenAssemblages[pairname]
                        eCount += 1
                    meanDistance = groupDistance / eCount         ##use the average distance as the metric
                else:
                    meanDistance = "0"

                ## initialize edges
                for e in network.edges_iter():
                    pairname = e[0] + "#" + e[1]
                    pairwise[pairname] = 0
                    pairwiseError[pairname] = 0

                for e in network.edges_iter():
                    pVal = 0.0
                    pErr = 0.0
                    if self.args['pairwisefile'] is not None:
                        pairname = e[0] + "#" + e[1]
                        pVal = pairwise[pairname]
                        pErr = pairwiseError[pairname]
                    else:
                        pVal = 0.0
                        pErr = 0.0
                    text = e[0] + " " + e[1] + " 1 " + str(edgeCount) + " " + str(network.graph["GraphID"]) + " " + \
                           str(pVal) + " " + str(pErr) + " " + str(meanDistance) + "\n"
                    OUTFILE.write(text)

                network.graph["meanDistance"] = meanDistance
                distanceHash[text] = meanDistance

            else:  ## not just the largest, but ALL seriation solutions
                edgeCount = len(network.edges())
                groupDistance = 0
                meanDistance = 0.0
                eCount = 0
                if self.args['xyfile'] is not None:
                    for e in network.edges_iter():
                        pairname = e[0] + "*" + e[1]
                        groupDistance += self.distanceBetweenAssemblages[pairname]
                        eCount += 1
                    meanDistance = groupDistance / eCount         ##use the average distance as the metric
                else:
                    meanDistance = "0"

                ## initialize edges
                for e in network.edges_iter():
                    pairname = e[0] + "#" + e[1]
                    pairwise[pairname] = 0
                    pairwiseError[pairname] = 0

                for e in network.edges_iter():
                    pVal = 0.0
                    pErr = 0.0
                    if self.args['pairwisefile'] is not None:
                        pairname = e[0] + "#" + e[1]
                        pVal = pairwise[pairname]
                        pErr = pairwiseError[pairname]
                    else:
                        pVal = 0.0
                        pErr = 0.0
                    text = e[0] + " " + e[1] + " 1 " + str(edgeCount) + " " + str(network.graph["GraphID"]) + " " + \
                           str(pVal) + " " + str(pErr) + " " + str(meanDistance) + "\n"
                    OUTFILE.write(text)

                network.graph["meanDistance"] = meanDistance
                distanceHash[text] = meanDistance

    # from a "summed" graph, create a "min max" solution - but use weights not counts
    def createMinMaxGraphByWeight(self, **kwargs):
        ## first need to find the pairs with the maximum occurrence, then we work down from there until all of the
        ## nodes are included
        ## the weight
        weight = kwargs.get('weight', "weight")
        input_graph = kwargs.get('input_graph')

        output_graph = nx.Graph(is_directed=False)

        ## first add all of the nodes
        for name in self.assemblages:
            output_graph.add_node(name, name=name, label=name, xCoordinate=self.xAssemblage[name],
                    yCoordinate=self.yAssemblage[name], size=self.assemblageSize[name]/self.totalAssemblageSize*self.nodeSizeFactor)

        pairsHash={}

        for e in input_graph.edges_iter():
            d = input_graph.get_edge_data(*e)
            fromAssemblage = e[0]
            toAssemblage = e[1]
            key = fromAssemblage+"*"+toAssemblage
            value = input_graph[fromAssemblage][toAssemblage]['weight']
            #print "Original value: ",self.sumOfDifferencesBetweenPairs[key], " New Value: ", value
            #value = self.sumOfDifferencesBetweenPairs[key]
            pairsHash[key]=value

        for key, value in sorted(pairsHash.iteritems(), key=operator.itemgetter(1), reverse=True ):
            ass1, ass2 = key.split("*")
            #print ass1, "-", ass2, "---",value
            edgesToAdd={}
            if nx.has_path(output_graph, ass1, ass2) == False:
                edgesToAdd[key]=value
                 ## check to see if any other pairs NOT already represented that have the same value
                for p in pairsHash:
                    if pairsHash[p] == value:
                        k1,k2 = p.split("*")
                        if nx.has_path(output_graph, k1,k2) == False:
                            edgesToAdd[p]=pairsHash[p]
                ## now add all of the edges that are the same value if they dont already exist as paths
                for newEdge in edgesToAdd:
                    a1,a2 = newEdge.split("*")
                    weight = self.sumOfDifferencesBetweenPairs[newEdge]
                    if weight in self.FalseList:
                        weight=0.000000000001
                    output_graph.add_path([a1, a2], weight=weight, inverseweight=(1/weight))
        return output_graph

    # from a "summed" graph, create a "min max" solution -- using Counts
    def createMinMaxGraphByCount(self, **kwargs):
        ## first need to find the pairs with the maximum occurrence, then we work down from there until all of the
        ## nodes are included
        ## the weight
        weight = kwargs.get('weight', "weight")
        input_graph = kwargs.get('input_graph')
        maxWeight = 0
        pairsHash = {}
        output_graph = nx.Graph(is_directed=False)

        for e in input_graph.edges_iter():
            d = input_graph.get_edge_data(*e)
            fromAssemblage = e[0]
            toAssemblage = e[1]
            if weight == "weight":
                currentWeight = int(d['weight'])
            else:
                currentWeight = int(d['weight'])
            pairsHash[fromAssemblage + "*" + toAssemblage] = currentWeight
            label = fromAssemblage + "*" + toAssemblage

        matchOnThisLevel = False
        currentValue = 0
        for key, value in sorted(pairsHash.iteritems(), key=operator.itemgetter(1), reverse=True):
            #print key, "->", value
            if value==0:
                value=.0000000000001
            if currentValue == 0:
                currentValue = value
            elif value < currentValue:
                matchOnThisLevel = False  ## we need to match all the connections with equivalent weights (otherwise we
                ## would stop after the nodes are included the first time which would be arbitrary)
                ## here we set the flag to false.
            ass1, ass2 = key.split("*")
            #print ass1, "-", ass2, "---",value
            if ass1 not in output_graph.nodes():
                output_graph.add_node(ass1, name=ass1, xCoordinate=self.xAssemblage[ass1],
                                      yCoordinate=self.yAssemblage[ass1],
                                      size=self.assemblageSize[ass1]/self.totalAssemblageSize*self.nodeSizeFactor)
            if ass2 not in output_graph.nodes():
                output_graph.add_node(ass2, name=ass2, xCoordinate=self.xAssemblage[ass2],
                                      yCoordinate=self.yAssemblage[ass2],
                                      size=self.assemblageSize[ass2]/self.totalAssemblageSize*self.nodeSizeFactor)
            if nx.has_path(output_graph, ass1, ass2) == False or matchOnThisLevel == True:
                matchOnThisLevel = True   ## setting this true allows us to match the condition that at least one match was
                ## made at this level

                output_graph.add_path([ass1, ass2], weight=value, inverseweight=(1/value ))

        return output_graph

    def filterSolutions(self, end_solutions, all_solutions):
        ################################################# FILTERING  ####################################
        # now do some weeding. Basically start with the last network ( largest),
        #  and work backwards to smaller and smaller solutions. Ignore any
        # network that is already represented larger since these are trivial (e.g., A->B->C->D already covers
        # A->C->D.) That should then leave all the unique maximum solutions (or so it seems)
        ################################################# FILTERING  ####################################

        filteredarray = []
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Filter to get uniques... ")
        logger.debug("--- Filtering solutions so we only end up with the unique ones.")
        logger.debug("--- Start with %d solutions.", len(end_solutions))
        filteredarray.append(end_solutions[-1])
        newOne = 0
        for tnetwork in reversed(all_solutions):
            exists = 0
            for fnetwork in filteredarray:
                fnetworkArray = fnetwork.nodes()
                logger.debug("----fnetworkArray: %s", fnetworkArray)
                tnetworkArray = tnetwork.nodes()
                logger.debug("----tnetworkArray: %s", tnetworkArray)
                minus = list(set(tnetworkArray) - set(fnetworkArray))
                logger.debug("difference between: %s ", minus)
                change = len(minus)
                logger.debug("Change: %d", change)
                if change > 0 and len(list(set(minus) - set(fnetworkArray))):
                    newOne += 1
                else:
                    exists += 1
            if exists == 0:
                logger.debug("pushing tnetwork to list of filtered arrays")
                filteredarray.append(tnetwork)
            exists = 0

        logger.debug("End with %d solutions.", len(filteredarray))
        filterCount = len(filteredarray)
        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(11, 1, "End with filterCount solutions.")

        return filteredarray

    def checkMinimumRequirements(self):
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError(
                "This function needs Graphviz and either PyGraphviz or Pydot. Please install GraphViz from http://www.graphviz.org/")
        if self.args['inputfile'] in self.FalseList:
            sys.exit("Inputfile is a required input value: --inputfile=../testdata/testdata.txt")


    def isSeriation(self, nnetwork):

        if self.args['screen'] not in self.FalseList:
            self.scr.addstr(1, 40, "STEP: Testing to see if this is a valid solution  ....      ")
            self.scr.refresh()

        error=0
        # just start from one side (shouldnt matter which).
        assemblages = nnetwork.nodes()
        startAssemblage=assemblages[0]
        testAssemblage=assemblages[1]

        ## first comaprison is different since its always okay and there are no incorrect values

        oneToColumns = range(len(self.assemblages[testAssemblage]))

        seriationList=[] ## hold all the seriations
        for i in oneToColumns:
            startVal = self.assemblages[startAssemblage][i]
            testVal = self.assemblages[testAssemblage][i]
            if startVal<testVal:
                compareVal="U"
            elif startVal>testVal:
                compareVal="D"
            else:
                compareVal="M"

            seriationList.append(compareVal)

        ## move to next assemblage
        startAssemblage = testAssemblage

        ## now start adding on assemblages and looking for failures.
        ######################################################################################
        for testAssemblage in assemblages[2:]:

            for i in oneToColumns:
                ser = seriationList[i]
                startVal = self.assemblages[startAssemblage][i]
                testVal = self.assemblages[testAssemblage][i]

                if startVal<testVal:
                    ser += "U"
                elif startVal>testVal:
                    ser += "D"
                else:
                    ser += "M"

                seriationList[i]=ser
            startAssemblage=testAssemblage

        ## now look for errors
        error = 0
        for s in seriationList:
            #print "Seriation test: ", s
            test = re.compile('DU|DM*U').search(s)
            if test not in self.FalseList:
                error +=1

        if error > 0:
            return False
        else:
            #print "Good seriation solution found:  ", nnetwork
            return True

    def findValidSeriations(self,graph):
        seriationList=[]
        nodes = graph.nodes()
        #print "Nodes: ", nodes
        pairs = self.all_pairs(nodes)
        for pair in pairs:
            #print "pair1 %s pair2 %s" % (pair[0],pair[1])
            try:
                paths= nx.all_simple_paths(graph,source=pair[0],target=pair[1])
                for p in paths:
                    graphID = uuid.uuid4().urn
                    newGraph = nx.Graph(End1=pair[0],End2=pair[1])
                    newGraph.graph["GraphID"] = graphID
                    newGraph.graph["name"] = graphID
                    newGraph.add_path(p)
                    for n in newGraph.nodes():
                        newGraph.add_node(n,name=n,label=n)
                    test = self.isSeriation(newGraph)
                    if test is True:
                        #print "%s is a seriation!" % newGraph
                        seriationList.append(newGraph)
            except nx.NetworkXNoPath:
                continue
                #print 'No path'

        return seriationList

    def calculateGeographicSolutionPValue(self,graph):

        bootstrap=1000
        solutionDistance=0
        assemblagesInSolution=[]
        edges=0
        for e in graph.edges_iter():
            d = graph.get_edge_data(*e)
            edges +=1
            fromAssemblage = e[0]
            toAssemblage = e[1]
            solutionDistance += math.sqrt(pow((int(self.xAssemblage[fromAssemblage])-int(self.xAssemblage[toAssemblage])),2)
                                +pow((int(self.yAssemblage[fromAssemblage])-int(self.yAssemblage[toAssemblage])),2))
            assemblagesInSolution.append(fromAssemblage)
            assemblagesInSolution.append(toAssemblage)
        assemblageSet=set(assemblagesInSolution)

        rnd.seed() # uses system time to initialize random number generator, or you can pass in a deterministic seed as an argument if you want
        x=[]
        pvalueScore=0.000
        for b in range(0,bootstrap):
            # code to use to generate K pairs
            list1 = self.labels
            list2 = self.labels

            testDistance=0
            for p in range(0,edges-1):
                test = False
                p1 = p2 = ""
                while test is False:
                    p1 = rnd.choice(list1)
                    p2 = rnd.choice(list2)
                    if p1 != p2:
                        test = True
                #print "Pair: ", p1, "-", p2
                testDistance += math.sqrt(pow((int(self.xAssemblage[p1])-int(self.xAssemblage[p2])),2)
                            +pow((int(self.yAssemblage[p1])-int(self.yAssemblage[p2])),2))
                #print "Test Distance: ", testDistance

            if testDistance <= solutionDistance:
                #print "TEST is less than solutionDistance: ",testDistance
                pvalueScore += 1
            x.append(testDistance)
        filename=self.outputDirectory+self.inputFile[0:-4]+"-geographic-distance.png"
        f=plt.figure(filename, figsize=(8, 8))
        plt.rcParams['font.family']='sans-serif'
        #f=plt.figure("Geographic Distance", figsize=(8, 8))
        num_bins = 20
        # the histogram of the data
        n, bins, patches = plt.hist(x, num_bins, facecolor='green', alpha=0.5)

        plt.axvline(solutionDistance, color='r', linestyle='dashed', linewidth=2)
        figure_label = self.inputFile[0:-4]
        plt.xlabel(figure_label)
        plt.ylabel('Count')
        plt.title(r'Histogram of Summed Geographic Distance')
        plt.savefig(filename, dpi=75)

        # Tweak spacing to prevent clipping of ylabel
        plt.subplots_adjust(left=0.15)
        minx =min(x)
        maxx=max(x)
        pvalue = pvalueScore/bootstrap
        x1,x2,y1,y2 = plt.axis()
        text="p-value: "+ str(pvalue)
        plt.text(maxx/3, (y2-y1)*2/3, text, style='italic')

        if pvalue == 0:
            pvalue ="0.000"
        return pvalue, solutionDistance, mean(x), std(x)

    #Prints everything in set b that's not in set a
    def difference(self, a, b):
        """This function compares two strings, a and b
        :param a: A string to be compared
        :param b: A string to be compared
        :returns: a list of the things that are different between the two strings
        """
        return list(set(b).difference(set(a)))

    def filterInclusiveSolutions(self, array):
        """This function filters out subsets in a list of solutions
        :param array: an array of graphs
        :returns: an array of the graphs that are unique and not subsets of each other
        """
        solutionSet = []
        for graph in sorted(array, key=lambda x: x.number_of_nodes(), reverse=True):
            #if len(solutionSet) == 0:
            #    solutionSet.append(graph)
            addSolution = True
            for sol in solutionSet:
                s = sol.nodes()
                g = graph.nodes()
                if set(g).issubset(s):
                    #print "Difference is: ", self.difference(s,g)
                    addSolution=False
            if addSolution is not False:
                solutionSet.append(graph)
        #print "length of final solution set: ", len(solutionSet)
        return solutionSet

    def seriate(self,args):

        self.checkMinimumRequirements()
        #####################################DEBUG OUTPUT#############################################################
        if self.args['debug'] is not None:
            ## Logging
            logger.basicConfig(stream=sys.stderr, level=logger.DEBUG)
            self.args['screen'] = None
        else:
            logger.basicConfig(stream=sys.stderr, level=logger.ERROR)
            self.args['screen'] = True

        logger.debug("Arguments: %s", self.args)

        ##################################################################################################
        if (self.args['screen'] is not None) and (self.args['debug'] is None ):
            ## Set up the screen display (default).
            ## the debug option should not use this since it gets messy
            try:
                # Initialize curses
                self.scr = curses.initscr()
                # Turn off echoing of keys, and enter cbreak mode,
                # where no buffering is performed on keyboard input
                curses.noecho()
                curses.cbreak()
                self.scr.nodelay(1)
                self.scr.addstr(0, 0, "Iterative Seriation Program V.2.0", curses.A_BOLD)
                self.scr.addstr(20, 35, "Hit <q> to quit.")
                self.scr.refresh()
            except:
                # In event of error, restore terminal to sane state.
                self.scr.keypad(0)
                curses.echo()
                curses.nocbreak()
                curses.endwin()
                curses.resetty()
                traceback.print_exc()           # Print the exception
                os.system("reset")

        if self.args['continuity'] in self.FalseList and self.args['frequency'] in self.FalseList and self.args['occurrence'] in self.FalseList:
            sys.exit(
                "You must specify --continuity=1 and/or frequency=1 to set the kind(s) of seriations you would like.")

        ######################################FILE INPUT#############################################################
        filename = self.args['inputfile']
        if filename is "":
            logger.error("You must enter a filename to continue.")
            print "You must enter a filename to continue."
            sys.exit("Quitting due to errors.")

        try:
            logger.debug("Going to try to open and load: %s", filename)
            self.openFile(filename)
        except IOError as e:
            logger.error("Cannot open %s. Error: %s", filename, e.strerror)

            print("Cannot open %s. Error. %s ", filename, e.strerror)
            if self.args['screen'] not in self.FalseList:
                curses.endwin()
                curses.resetty()
            sys.exit("Quitting due to errors.")

        try:
            inputparts = map(str, self.args['inputfile'].split("/"))
            self.inputFile = inputparts[len(inputparts) - 1]
        except:
            sys.exit("There was a problem with parsing the input file. Check it and try again.")

        ############################################################################################################
        if self.args['outputdirectory'] not in self.FalseList:
            self.outputDirectory = self.args['outputdirectory']
        else:
            self.outputDirectory = "../output/"


        ###########################################################################################################
        ### If occurrence seriation, collapse all the identical solutions (shouldnt be a problem for the frequency set but I suppose it could be
        if self.args['occurrence'] not in self.FalseList:
            self.aggregateIdenticalAssemblages()
        if len(self.assemblages)<3:
            print "Problem: you only have 2 assemblages (perhaps after consolidation). There is no point in continuing. "
            sys.exit()

        ############################################################################################################
        logger.debug("Going to open pairwise file it is exists.")
        if self.args['pairwisefile'] not in self.FalseList:
            self.openPairwiseFile(self.args['pairwisefile'])

        ############################################################################################################
        logger.debug("Going to open XY file if it exists.")
        if self.args['xyfile'] not in self.FalseList:
            self.openXYFile(self.args['xyfile'])
        else:
            for ass in self.assemblages:
                self.xAssemblage[ass] = 0.0
                self.yAssemblage[ass] = 0.0
            allp = self.all_pairs(self.assemblages)
            for pr in allp:
                name = pr[0] + "*" + pr[1]
                self.distanceBetweenAssemblages[name] = 0

        ############################################################################################################
        logger.debug("Assume threshold is 1.0 unless its specified in arguments.")
        threshold = 1.0
        if self.args['threshold'] is not None:
            threshold = float(self.args['threshold'])

        logger.debug("Going to create list of valid pairs for comparisons.")
        self.thresholdDetermination(threshold)

        ###########################################################################################################
        logger.debug("Now calculate the bootstrap comparisons based ")
        logger.debug("on specified confidence interval, if in the arguments.")

        if self.args['bootstrapCI'] is not None:
            if self.args['bootstrapSignificance'] not in self.FalseList:
                confidenceInterval = self.args['bootstrapSignificance']
            else:
                confidenceInterval = 0.95
            self.bootstrapCICalculation( 100, float(confidenceInterval))

        ###########################################################################################################
        ### setup the output files. Do this now so that if it fails, its not AFTER all the seriation stuff
        OUTFILE, OUTPAIRSFILE, OUTMSTFILE, OUTMSTDISTANCEFILE = self.setupOutput()



        ###########################################################################################################
        logger.debug("Now pre-calculating all the combinations between pairs of assemblages. ")
        logger.debug("This returns a graph with all pairs and the comparisons as weights.")
        pairGraph = self.preCalculateComparisons()

        #####################################

        logger.debug("Now calculate sum of differences between all pairs")
        self.preCalculateSumOfDifferencesBetweenPairs()

        #####################################

        frequencyArray = []
        continuityArray = []
        maxNodes = 3
        notPartOfSeriationsList = []

        if self.args['frequency'] not in self.FalseList or self.args['occurrence'] not in self.FalseList:
            ###########################################################################################################
            logger.debug("Calculate all the valid triples.")
            triples = self.findAllValidTriples()
            ###########################################################################################################
            stepcount = 0
            currentMaxSeriationSize = 2
            newNetworks = []
            self.solutionCount = len(triples)     ## the current # of solutions that will be checked
            self.solutionsChecked = len(triples)  ## set the total solutions checked value
            currentTotal = len(triples)
            solutions = []
            all_solutions = []
            all_solutions = all_solutions + triples  ## add the triples to the initial solution

            ## create a directory for the processing of the pickle files if it doesnt exist already
            if not os.path.exists(".p"):
                os.makedirs(".p")
            ## pickle the stuff I need for parallel processing
            ch = open('.p/validComparisonsHash.p', 'wb')
            pickle.dump(self.validComparisonsHash,open('.p/validComparisonsHash.p', 'wb'))
            #ch.close()
            pg=open('.p/pairGraph.p','wb')
            pickle.dump(self.pairGraph,open('.p/pairGraph.p','wb'))
            #pg.close
            ass=open('.p/assemblages.p','wb')
            pickle.dump(self.assemblages,open('.p/assemblages.p','wb'))
            #ass.close()
            a=open('.p/args.p','wb')
            pickle.dump(self.args,open('.p/args.p','wb'))
            #a.close()
            fuci=open('.p/typeFrequencyUpperCI.p','wb')
            pickle.dump(self.typeFrequencyUpperCI,open('.p/typeFrequencyUpperCI.p','wb'))
            #fuci.close()
            flci=open('.p/typeFrequencyLowerCI.p','wb')
            pickle.dump(self.typeFrequencyLowerCI,open('.p/typeFrequencyLowerCI.p','wb'))
            #flci.close()
            while currentMaxSeriationSize <= self.maxSeriationSize:
                currentMaxSeriationSize += 1
                ### first time through copy the triples, else get the previous new ones.
                if currentMaxSeriationSize == 3:  ## first time through. Just copy the triples to the working arrays
                    networks = triples
                    solutions = triples # clear the
                else:
                    i = 0
                    #print("Currently have %d solutions at step %d"%( len(newNetworks), currentMaxSeriationSize))
                    if len(newNetworks) > 0:
                        # there were no networks the previous times so nothing to do.
                        logger.debug("These solutions are ---  ")
                        for sol in newNetworks:
                            logger.debug("solution %d: %s", i, nx.shortest_path(sol, sol.graph["End1"], sol.graph["End2"]))
                            i += 1
                        networks = []
                        networks += newNetworks  # copy the array of previous new ones for this round
                        solutions.append(newNetworks) # append the new list to the previous one
                        newNetworks = []         # clear the array of new solutions

                stepcount += 1
                logger.debug("_______________________________________________________________________________________")
                logger.debug("Step number:  %d", currentMaxSeriationSize)
                logger.debug("_______________________________________________________________________________________")

                if self.args['screen'] not in self.FalseList:
                    self.scr.addstr(4, 0, "Step number:                                    ")
                    msg = "Step number:   %d" % currentMaxSeriationSize
                    self.scr.addstr(4, 0, msg)
                    self.scr.addstr(5, 0, "Number of solutions from previous step:         ")
                    msg = "Number of solutions from previous step: %d" % len(networks)
                    self.scr.addstr(5, 0, msg)
                    self.scr.refresh()

                logger.debug("Number of solutions from previous step: %d", len(networks))
                match = 0      ## set the current match to zero for this step (sees if there are any new solutions for step)
                ## look through the set of existing valid networks.
                validNewNetworks = []

                try:
                    cpus = multiprocessing.cpu_count()
                except NotImplementedError:
                    cpus = 2   # arbitrary default

                # Each process will get 'chunksize' nums and a queue to put his out
                # dict into

                out_q = multiprocessing.Queue()
                chunksize = int(math.ceil(len(networks) / float(cpus)))
                procs = []

                for i in range(cpus):
                    p = multiprocessing.Process(
                    target=seriationEvaluation.worker,
                    args=(networks[chunksize * i:chunksize * (i + 1)],out_q))
                    procs.append(p)
                    p.start()

                # Collect all results into a single result dict. We know how many dicts
                # with results to expect.
                resultdict = []
                for i in range(cpus):
                    resultdict.append(out_q.get())

                # Wait for all worker processes to finish
                for p in procs:
                    p.join()

                #validNewNetworks = [x for x in result if not x is False]

                for s in resultdict:
                    if s is not False:
                        #for sol in s:
                        newNetworks += s
                        all_solutions += s
                        self.solutionCount += len(s)
                        logger.debug("Added %d new solutions. Solution count is now:  %d", len(validNewNetworks),
                                 self.solutionCount)
                        if len(s) > maxNodes:
                            maxNodes = len(s)
                        currentTotal = len(newNetworks)

                if self.args['screen'] not in self.FalseList:
                    msg = "Current Max Nodes:  %d " % maxNodes
                    self.scr.addstr(6, 0, msg)
                    msg = "Total number of seriation solutions and sub-solutions: %d" % self.solutionCount
                    self.scr.addstr(7, 0, msg)
                    self.scr.addstr(8, 43, "                                           ")
                    msg = "Number of seriation solutions at this step: %d" % currentTotal
                    self.scr.addstr(8, 0, msg)
                    if os.name != "nt":
                        msg = "Memory used:        " + str(self.mem.memory())
                    self.scr.addstr(9, 0, msg)
                    self.scr.refresh()

                if len(newNetworks) > 0:
                    end_solutions = newNetworks
                else:
                    end_solutions = networks

            logger.debug("Process complete at seriation size %d with %d solutions before filtering.",
                         self.maxSeriationSize, len(end_solutions))

            all_solutions += end_solutions
            ###########################################################################################################
            frequencyArray = self.filterSolutions(end_solutions, all_solutions)
            #if self.args['allsolutions'] not in self.FalseList:
            #frequencyArray = self.filterInclusiveSolutions(all_solutions,end_solutions)
            #else:
            #    frequencyArray = self.filterInclusiveSolutions(end_solutions)
            #filteredarray = all_solutions

            logger.debug("Process complete at seriation size %d with %d solutions after filtering.",
                         self.maxSeriationSize, len(frequencyArray))

            if self.args['verbose'] not in self.FalseList:
                ## determine time elapsed
                #time.sleep(5)
                timeNow = time.time()
                timeElapsed = timeNow - self.start
                print "Time elapsed for frequency seriation processing: %d seconds" % timeElapsed

            #################################################### OUTPUT SECTION ####################################################
            self.output(frequencyArray, OUTFILE, OUTPAIRSFILE, OUTMSTFILE, OUTMSTDISTANCEFILE, maxNodes)

            if self.args['atlas'] not in self.FalseList:
                self.createAtlasOfSolutions(frequencyArray, "frequency")

            sumGraphByWeight = self.sumGraphsByWeight(frequencyArray)
            self.sumGraphOutput(sumGraphByWeight, self.outputDirectory + self.inputFile[0:-4] + "-sumgraph-by-weight")

            sumGraphByCount = self.sumGraphsByCount(frequencyArray)
            self.sumGraphOutput(sumGraphByCount, self.outputDirectory + self.inputFile[0:-4] + "-sumgraph-by-count")

            if self.args['excel'] not in self.FalseList:
                excelFileName,textFileName=self.outputExcel(frequencyArray, self.outputDirectory+self.inputFile[0:-4], "frequency")

            if self.args['frequencyseriation'] not in self.FalseList:
                excelFileName,textFileName=self.outputExcel(frequencyArray, self.outputDirectory+self.inputFile[0:-4], "frequency")
                seriation = frequencySeriationMaker()
                #self.args={'inputfile':textFileName,'pdf':1}
                self.args['inputfile']=textFileName
                self.args['multiple']=1
                seriation.makeGraph(self.args)

            if self.args['occurrenceseriation'] not in self.FalseList:
                excelFileName,textFileName=self.outputExcel(frequencyArray, self.outputDirectory+self.inputFile[0:-4], "occurrence")
                seriation = occurrenceSeriationMaker()
                #self.args={'inputfile':textFileName,'pdf':1}
                self.args['multiple']=1
                self.args['inputfile']=textFileName
                seriation.makeGraph(self.args)
            #################################################### MinMax Graph ############################################
            #print self.args
            minMaxGraphByWeight = self.createMinMaxGraphByWeight(input_graph=sumGraphByWeight, weight='weight')
            if self.args['xyfile'] not in self.FalseList:
                pscore, distance, geodistance, sd_geodistance = self.calculateGeographicSolutionPValue(minMaxGraphByWeight)
                print "Geographic p-value for the frequency seriation minmax solution: ", pscore
                filename=self.outputDirectory + "geography.txt"
                with open(filename, "a") as myfile:
                    text=self.inputFile[0:-4]+"\t"+str(pscore)+"\t"+str(distance)+"\t"+str(geodistance)+"\t" \
                         + str(sd_geodistance)+"\t"+str(self.totalAssemblageSize)+"\n"
                    myfile.write(text)

            minMaxGraphByCount = self.createMinMaxGraphByCount(input_graph=sumGraphByCount, weight='weight')
            #if self.args['graphs'] not in self.FalseList:
            self.graphOutput(minMaxGraphByWeight,
                        self.outputDirectory + self.inputFile[0:-4] + "-minmax-by-weight.png")
            self.graphOutput(minMaxGraphByCount,
                        self.outputDirectory + self.inputFile[0:-4] + "-minmax-by-count.png")

            #################################################### MST SECTION ####################################################
            if self.args['mst'] not in self.FalseList:
                outputFile = self.outputDirectory + self.inputFile[0:-4] + ".vna"
                # Need to have the shapefile flag and the XY file in order to create a valid shapefile.
                if self.args['shapefile'] is not None and self.args['xyfile'] is not None:
                    shapefile = 1
                else:
                    shapefile = None
                mst = MST.MST(outputFile, self.outputDirectory, shapefile)
                mst.createMST()
                #minimumSpanningTree(all_solutions,xAssemblage,yAssemblage,distanceBetweenAssemblages,assemblageSize,outputDirectory,inputFile)
            #################################################### END SECTION ####################################################

            if self.args['verbose'] not in self.FalseList:
                print "Seriation complete."
                print "Maximum size of seriation: %d" % maxNodes
                print "Number of frequency seriation solutions at last step: %d" % len(frequencyArray)
                print "Assemblages not part of final solution:"
                nodeList = sumGraphByWeight.nodes()
                for a in self.assemblages:
                    if a not in nodeList:
                        notPartOfSeriationsList.append(a)
                        print a
                if len(notPartOfSeriationsList) == 0:
                    print "*** All assemblages used in seriations.***"

        if self.args['continuity'] not in self.FalseList:
            # experimental
            continuityArray = self.continunityMaximizationSeriation()
            #self.outputGraphArray(array)
            sGraphByCount = self.sumGraphsByCount(continuityArray)
            sGraphByWeight = self.sumGraphsByWeight(continuityArray)
            self.graphOutput(sGraphByCount, self.inputFile[0:-4] + "-continuity-sumgraph.png")
            self.MST(sGraphByCount, self.outputDirectory + self.inputFile[0:-4] + "-mst-of-min.png")
            minMaxGraphByWeight = self.createMinMaxGraphByWeight(input_graph=sGraphByWeight, weight='weight')
            self.graphOutput(minMaxGraphByWeight, self.outputDirectory +  self.inputFile[0:-4] + "-continuity-minmax-by-weight.png")
            minMaxGraphByCount = self.createMinMaxGraphByCount(input_graph=sGraphByCount, weight='weight')
            self.graphOutput(minMaxGraphByCount, self.outputDirectory + self.inputFile[0:-4] + "-continuity-minmax-by-count.png")
            if self.args['atlas'] not in self.FalseList:
                self.createAtlasOfSolutions(continuityArray, "continuity")

            if self.args['excel'] not in self.FalseList:
                self.outputExcel(continuityArray, self.outputDirectory+self.inputFile[0:-4], "continuity")

            if self.args['frequencyseriation'] not in self.FalseList:
                excelFileName,textFileName=self.outputExcel(continuityArray, self.outputDirectory+self.inputFile[0:-4], "continuity")
                seriation = frequencySeriationMaker()
                argument={'inputfile':textFileName, 'multiple':True}
                seriation.makeGraph(argument)

            if self.args['verbose'] not in self.FalseList:
                ## determine time elapsed
                #time.sleep(5)
                timeNow = time.time()
                timeElapsed = timeNow - self.start
                print "Number of continuity seriation solutions at end: %d " % len(continuityArray)
                print "Time elapsed for continuity seriation processing: %d seconds" % timeElapsed

            validSeriations = self.findValidSeriations(minMaxGraphByWeight)

            self.createAtlasOfSolutions(validSeriations, "Valid_Seriations")
            filteredSet=self.filterInclusiveSolutions(validSeriations)
            #print "Filtered set:", filteredSet
            self.createAtlasOfSolutions(filteredSet, "Unique_Valid_Seriations")
            if self.args['excel'] not in self.FalseList:
                excelFileName,textFileName=self.outputExcel(filteredSet, self.outputDirectory+self.inputFile[0:-4], "valid_continuity")
                seriation = frequencySeriationMaker()
                argument={'inputfile':textFileName, 'multiple':True}
                seriation.makeGraph(argument)

            if self.args['xyfile'] not in self.FalseList:
                pscore ,distance, geodistance, sd_geodistance = self.calculateGeographicSolutionPValue(minMaxGraphByWeight)
                print "Geographic p-value for the continuity seriation minmax solution: ", pscore
                filename=self.outputDirectory + "geography.txt"
                with open(filename, "a") as myfile:
                    text=self.inputFile[0:-4]+"\t"+str(pscore)+"\t"+str(distance)+"\t"+str(geodistance)+"\t" \
                         + str(sd_geodistance)+"\t"+str(self.totalAssemblageSize)+"\n"
                    myfile.write(text)

        ## determine time elapsed
        #time.sleep(5)
        timeNow = time.time()
        timeElapsed = timeNow - self.start
        if self.args['verbose'] not in self.FalseList:
            print "Time elapsed for completion of program: %d seconds" % timeElapsed

        if self.args['graphs'] not in self.FalseList:
            plt.show() # display

        ## say goodbye and clean up the screen stuff #########################
        self.finalGoodbye()

        return frequencyArray, continuityArray, notPartOfSeriationsList

    def addOptions(self):
        self.args = {'debug': None, 'bootstrapCI': None, 'bootstrapSignificance': None,
                'filtered': None, 'largestonly': None, 'individualfileoutput': None, 'xyfile':None,
                'excel': None, 'threshold': None, 'noscreen': None, 'xyfile': None, 'pairwisefile': None, 'mst': None,
                'stats': None, 'screen': None, 'allsolutions': None, 'inputfile': None, 'outputdirectory': None,
                'shapefile': None, 'frequency': None, 'continuity': None, 'graphs': None, 'graphroot': None, ''
                'continuityroot': None, 'verbose':None, 'occurrenceseriation':None,'occurrences':None,'frequency':None,
                'occurrence':None,'frequencyseriation':None, 'pdf':None, 'atlas':None}

    def parse_arguments(self):
        self.addOptions()
        parser = argparse.ArgumentParser(description='Conduct an iterative deterministic seriation analysis')
        parser.add_argument('--debug', '-d', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
        parser.add_argument('--bootstrapCI', '-b', default=None,
                            help="Sets whether you want to use the bootstrap confidence intervals for the comparisons between assemblage type frequencies. Set's to on or off.")
        parser.add_argument('--bootstrapSignificance', '-bs', default=0.95, type=float,
                            help="The significance to which the confidence intervals are calculated. Default is 0.95.")
        parser.add_argument('--filtered','-f', default=1,
                            help="The script will complete by checking to see if smaller valid solutions are included in the larger sets. If not, they are added to the final set. Default is true. ")
        parser.add_argument('--largestonly','-lo', default=None,
                            help="If set, the results will only include the results from the last and largest successful series of solutions. Smaller solutions will be excluded. Default is false.")
        parser.add_argument('--individualfileoutput', default=None,
                            help="If true, a .VNA files will be created for every solution.")
        parser.add_argument('--threshold', default=None,
                            help="Sets the maximum difference between the frequencies of types that will be examine. This has the effect of keeping one from evaluating trivial solutions or solutions in which there is limited warrant for establishing continuity. Default is false.")
        parser.add_argument('--noscreen', default=None,
                            help="If true, there will be no text output (i.e., runs silently). Default is false.")
        parser.add_argument('--xyfile', default=None,
                            help="Enter the name of the XY file that contains the name of the assemblage and the X and Y coordinates for each.")
        parser.add_argument('--pairwisefile', default=None,
                            help="If you have precalculated the bootstrap comparative p-value, enter the name of the file here and it will be used as the basis of the graphical output for showing significance of comparisons. Default is false.")
        parser.add_argument('--mst', default=None,
                            help="If true, will produce a minimum spanning tree diagram from the set of final solutions.")
        parser.add_argument('--stats', default=None,
                            help="(Not implemented). If true, a histogram of the solutions will be shown in terms of the #s of time pairs are included. Default is false.")
        parser.add_argument('--screen', default=True,
                            help="Sets whether the output will be sent all to the screen or not. Default is false. When true, the screen output is all captured through curses.")
        parser.add_argument('--allsolutions', default=None,
                            help="If set, all of the valid solutions are produced even if they are subsets of larger solutions.")
        parser.add_argument('--inputfile',
                            help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
        parser.add_argument('--outputdirectory', default=None,
                            help="If you want the output to go someplace other than the /output directory, specify that here.")
        parser.add_argument('--shapefile', default=None,
                            help="Produces a shapefile as part of the output. You must have specified the --xyfile (coordinates for each point) as well.")
        parser.add_argument('--graphs', default=None,
                            help="If true, the program will display the graphs that are created. If not, the graphs are just saved as .png files.")
        parser.add_argument('--frequency', default=None,
                            help="Conduct a standard frequency seriation analysis. Default is None.")
        parser.add_argument('--continuity', default=None, help="Conduct a continuity seriation analysis. Default is None.")
        parser.add_argument('--graphroot', default=None,
                            help="The root of the graph figures (i.e., name of assemblage you want to treat as one end in the graphs.")
        parser.add_argument('--continuityroot', default=None,
                            help="If you have a outgroup or root of the graph, set that here.")
        parser.add_argument('--atlas', default=None,
                            help="If you want to have a figure that shows all of the results independently, set that here.")
        parser.add_argument('--excel', default=None,
                            help="Will create excel files with the assemblages in seriation order.")
        parser.add_argument('--noheader',default=None,
                            help="If you do not use type names as the first line of the input file, use this option to read the data.")
        parser.add_argument('--frequencyseriation', default=None, help="Generates graphical output for the results in a frequency seriation form.")
        parser.add_argument('--verbose',default=True, help='Provides output for your information')
        parser.add_argument('--occurrence', default=None, help="Treats data as just occurrence information and produces valid occurrence solutions.")
        parser.add_argument('--occurrenceseriation', default=None, help="Generates graphical output for occurrence seriation.")
        try:
            self.args = vars(parser.parse_args())
        except IOError, msg:
            parser.error(str(msg))
            sys.exit()
        return self.args

if __name__ == "__main__":
    seriation = IDSS()
    args = seriation.parse_arguments()
    frequencyResults, continuityResults, exceptionList = seriation.seriate(args)

''''
From the command line:

python ./IDSS.py --inputfile=../testdata/pfg.txt --xyfile=../testdata/pfgXY.txt --largestonly=1 --mst=1 --graphs=1"


As a module:

from IDSS import IDSS

seriation = IDSS()

args={}
args{'inputfile'}="../testdata/testdata-5.txt"
args{'screen'}=1
args{'debug'}=1
args('graphs'}=1

frequencyResults,continuityResults,exceptions = seriation.seriate(args)

'''''


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

