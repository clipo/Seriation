#!/usr/bin/env python
# Copyright (c) 2013.  Carl P. Lipo <clipo@csulb.edu>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.
__author__ = 'carllipo'

import MST
import shapefile
import csv
from datetime import datetime
import argparse
import sys
import logging as logger
import itertools
import math
import random
import curses
from itertools import chain
import numpy as np
import scipy as sp
import scipy.stats
import networkx as nx
import traceback
import memory
import operator
import time
from datetime import datetime
import os
from pylab import *
import matplotlib.pyplot as plt
import re
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class IDSS():

    color=[ "b","r","m","y","k","w",(0.976,0.333,0.518),(0.643,0.416,0.894),
            (0.863,0.66,0.447),(0.824, 0.412, 0.118) ]

    def __init__(self):
        self.inputfile = ""
        self.outputDirectory=""
        self.mem=memory.Memory()
        self.start = time.time()
        self.assemblageSize={}
        self.assemblageFrequencies={}
        self.assemblages={}
        self.countOfAssemblages=0
        self.assemblageValues={}
        self.labels={}
        self.numberOfClasses =0
        self.maxSeriationSize=0
        self.xAssemblage = {}
        self.yAssemblage = {}
        self.xyAssemblages=[]
        self.distanceBetweenAssemblages={}
        self.largestX =0
        self.largestY =0
        self.distanceBetweenAssemblages={}
        self.validComparisonsHash={}
        self.typeFrequencyLowerCI = {}
        self.typeFrequencyUpperCI = {}
        self.typeFrequencyMeanCI = {}
        self.pairwise={}
        self.pairwiseError={}
        logger.debug("Start time:  %s ", self.start)
        self.scr = None

    def saveGraph(self,graph,filename,args):
        nx.write_gml(graph,filename)

    def all_pairs(self,lst):
        return list((itertools.permutations(lst, 2)))

    def all_tuples(self,lst):
        tuples=list(itertools.combinations(lst, 3))
        useable_tuples=[]
        for e in tuples:
            useable_tuples.append(e)
        return useable_tuples

    def openFile(self, filename, args):
        try:
            logger.debug("trying to open: %s ", filename)
            file=open(filename,'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", filename, e)
            sys.exit('file %s does not open: %s') %( filename, e)

        reader = csv.reader(file, delimiter='\t', quotechar='|')
        values=[]
        for row in reader:
            row = map(str, row)
            label=row[0]
            self.labels[ label ] = label
            row.pop(0)
            row = map(float, row)
            self.numberOfClasses = len(row)
            freq=[]
            rowtotal=sum(row)
            for r in row:
                freq.append(float(float(r)/float(rowtotal)))
                values.append(float(r))
            self.assemblages[ label ] = freq
            self.assemblageFrequencies[ label ]  = freq
            self.assemblageValues[ label ] = values
            self.assemblageSize[ label ]= rowtotal
            self.countOfAssemblages +=1
        self.maxSeriationSize=self.countOfAssemblages
        return True

    def preCalculateComparisons(self,args):
        logger.debug("Precalculating the comparisons between all pairs of assemblages...")
        pairs = self.all_pairs(self.assemblages)
        pairGraph = nx.Graph()
        for pair in pairs:
            pairGraph.add_node(pair[0])
            pairGraph.add_node(pair[1])
            columns=range(len(self.assemblages[pair[0]]))
            ass1 = self.assemblages[pair[0]]
            ass2 = self.assemblages[pair[1]]
            comparison=""
            for i in columns:
                val1 = ass1[i]
                val2 = ass2[i]
                logger.debug( "\t\tComparing Assemblage: %s  and    Assemblage: %s  ########",pair[0],pair[1])
                logger.debug( "\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,i,i,i,i,i,i)

                if args['bootstrapCI'] not in (None, ""):
                    upperCI_test = self.typeFrequencyUpperCI[pair[0]][i]
                    lowerCI_test  = self.typeFrequencyLowerCI[pair[0]][i]
                    upperCI_end =  self.typeFrequencyUpperCI[pair[1]][i]
                    lowerCI_end=  self.typeFrequencyLowerCI[pair[1]][i]

                    if upperCI_test < lowerCI_end:
                        comparison +=  "D"
                    elif lowerCI_test > upperCI_end:
                        comparison +=  "U"
                    else:
                        comparison +=  "M"
                else:
                    if val1 < val2 :
                        comparison +=  "D"
                    if val1 > val2:
                        comparison +=  "U"
                    if val1==val2:
                        comparison +=  "M"
                logger.debug( "Type %d: - comparison is: %s ",i, comparison[i])

            logger.debug("Comparison for %s and %s is: %s ",pair[0],pair[1],comparison)
            pairGraph.add_edge(pair[0],pair[1],weight=comparison)
        return pairGraph

    def openPairwiseFile(self,filename ):
        logger.debug("Opening pairwise file %", filename)
        try:
            pw = open(filename,'r' )
        except csv.Error as e:
            sys.exit('pairwise file %s does not open: %s') %( filename, e)

        reader = csv.reader(pw, delimiter='\t', quotechar='|')
        for row in reader:
            pair = row[0]+"#"+row[1]
            self.pairwise[pair]=row[2]
            self.pairwiseError[pair]=row[3]
        return True

    def openXYFile(self,filename ):
        logger.debug("Opening pairwise file %s", filename)
        ## open the xy file
        try:
            xyf= open( filename,'r')
        except csv.Error as e:
            sys.exit('file %s does not open: %s') %( filename, e)

        reader = csv.reader(xyf, delimiter='\t', quotechar='|')

        for row in reader:
            label = row[0]
            self.xyAssemblages.append(label)
            self.yAssemblage[label]=row[1]
            self.xAssemblage[label]=row[2]

        assemblagePairs = self.all_pairs(self.xyAssemblages)
        ## Go through all of the combinations
        for combo in assemblagePairs:
            pairname = combo[0]+"*"+combo[1]
            xdistance = float(self.xAssemblage[combo[0]]) - float(self.xAssemblage[combo[1]])
            xxdistance = xdistance * xdistance
            ydistance = float(self.yAssemblage[combo[0]]) - float(self.yAssemblage[combo[1]])
            yydistance = ydistance * ydistance
            distance = math.sqrt( xxdistance + yydistance)
            self.distanceBetweenAssemblages[ pairname ] = distance
        largestXname = max(self.xAssemblage.iterkeys(), key=(lambda key: self.xAssemblage[key]))
        largestYname= max(self.yAssemblage.iterkeys(), key=(lambda key: self.yAssemblage[key]))
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

    def thresholdDetermination(self,threshold,args):
        assemblageComparison={}
        ##  get all the combinations of 2
        pairs = self.all_pairs(self.assemblages)

        ## Go through all of the combinations
        for combo in pairs:
            logger.debug("comparing combination of %s and %s ", combo[0] , combo[1] )
            pairname1  = combo[0] + "*" + combo[1]
            pairname2  = combo[1] + "*" + combo[0]
            maxDifference = 0
            assemblage1 = self.assemblages[combo[0]]
            assemblage2 = self.assemblages[combo[1]]
            i=-1
            columns= len(self.assemblages[combo[0]])
            logger.debug("Number of columns: %d", columns)
            ## calculate the maximum differences between the pairs of assemblages (across all types)
            for i in (0, columns-1):
                logger.debug("i: %d ",i)
                ass1 = float(assemblage1[i])
                ass2 = float(assemblage2[i])
                diff = abs( ass1 - ass2 )
                logger.debug("assemblage1: %f assemblage2: %f diff: %f",ass1,ass2,diff)
                if diff > maxDifference :
                    maxDifference = diff
            assemblageComparison[ pairname1 ] = maxDifference
            assemblageComparison[ pairname2 ] = maxDifference

        ############## pre calculate the valid pairs of assemblages and stuff into hash ############################
        for assemblage1 in self.assemblages:
            cAssemblages=[]
            for assemblage2 in self.assemblages:
                if not assemblage1 == assemblage2:
                    testpair = assemblage1 + "*" + assemblage2
                    logger.debug("Pairs: %s and %s", assemblage1,assemblage2)
                    logger.debug("Comp value of pairs: %s:  %f and threshold is: %f",testpair, assemblageComparison[testpair],threshold)
                    if assemblageComparison[ testpair ] <= threshold:
                        logger.debug("Appending %s to the list of valid comparisons for %s ", assemblage1, assemblage2)
                        cAssemblages.append( assemblage2 )
            self.validComparisonsHash[ assemblage1]  = cAssemblages
        return True

    def confidence_interval(self,data, confidence=0.05):
        a = 1.0*np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
        return m, m-h, m+h

    ########################################### BOOTSTRAP CI SECTION ####################################
    def bootstrapCICalculation(self,args, bootsize=1000,confidenceInterval=0.05):

        if args['screen']:
            self.scr.addstr(1,40, "STEP: Bootstrap CIs...        ")
            self.scr.refresh()

        ## for each assemblage
        logger.debug("Calculating bootstrap confidence intervals")
        # for each assemblage

        for currentLabel in sorted( self.assemblages.iterkeys()):
            assemblage =  self.assemblages[ currentLabel ]
            types=len(self.assemblages[ currentLabel ])
            currentAssemblageSize = self.assemblageSize[ currentLabel ]

            ## create an array of arrays - one for each type
            arrayOfStats=[]
            for c in range(0,types):
                array=[]
                arrayOfStats.append([])

            ## size of bootstrapping (how many assemblages to create)
            loop = bootsize
            for counter in range(0,bootsize):

                assemsize = currentAssemblageSize
                # clear and set the array
                cumulate=[]
                for d in range(0,types):
                    cumulate.append(0.0)

                index = 0
                count=0
                ## now count through the classes and set up the frequencies
                for typeFrequency in assemblage:
                    index += typeFrequency
                    cumulate[count] = index  ## this is the index of the frequency for this class
                    ## index should be total # of types at end
                    count += 1

                ## set new_assemblage
                new_assemblage=[]
                for c in range(0,types):
                    new_assemblage.append(0.0)

                for sherd in range(0,int(currentAssemblageSize)):
                    rand = random()             ## random number from 0-1
                    classVar = 0
                    typeIndex=types-1
                    found=0
                    total = sum(cumulate)
                    for t in reversed(cumulate):
                        if rand <= t:
                            found=typeIndex
                        typeIndex -=1
                    new_assemblage[found] += 1

                ## count new assemblage frequencies
                counter=0
                new_assemblage_freq = []
                for g in new_assemblage:
                    new_assemblage_freq.append(float(g/float(bootsize)))
                    arrayOfStats[counter].append(float(g/float(bootsize)))
                    counter += 1
                ## this should result in a new assemblage of the same size

            lowerCI=[]
            upperCI=[]
            meanCI=[]
            for freq in arrayOfStats:
                upper=0.0
                lower=0.0
                mean=0.0
                if sum(freq) > 0.0:
                    mean, lower, upper = self.confidence_interval(freq, confidence=float(confidenceInterval))
                else:
                    mean=lower=upper=0
                if math.isnan(lower) is True:
                    lower=0.0
                if math.isnan(upper) is True:
                    upper=0.0
                if math.isnan(mean) is True:
                    mean=0.0
                lowerCI.append( lower)
                upperCI.append( upper )
                meanCI.append( mean )

            self.typeFrequencyLowerCI[ currentLabel ] = lowerCI
            self.typeFrequencyUpperCI[ currentLabel ] = upperCI
            self.typeFrequencyMeanCI[ currentLabel ] = meanCI

        return True

    ########################################### FIND ALL THE VALID TRIPLES  ####################################
    ########################################### #################################### ###########################
    def findAllValidTriples(self,args):
        triples=[]
        error = 0
        numberOfTriplets = 0

        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40, "STEP: Find valid triples....      ")
            self.scr.refresh()
        permutations = self.all_tuples(self.assemblages)

        for permu in permutations:
            if args['screen'] not in (None, ""):
                c = self.scr.getch()
                if c == ord('q'):
                    curses.endwin()
                    curses.resetty()
                    sys.exit("Quitting as requested.\n\r")
            logger.debug("Triple test: %s * %s * %s", permu[0],permu[1],permu[2])

            comparison12 = ""
            comparison23 = ""
            error = 0
            columns=len( self.assemblages[ permu[0] ])
            logger.debug("Columns: %d", columns)
            difscore=0
            difscore2=0
            comparison12=""
            comparison23=""
            for i in range(0,columns):
                ass1 = self.assemblages[ permu[0] ][i]
                ass2 = self.assemblages[ permu[1] ][i]
                ass3 = self.assemblages[ permu[2] ][i]
                logger.debug( "ass1: %f ass2: %f ass3: %f",ass1,ass2,ass3)

                if args['bootstrapCI'] not in (None, ""):
                    low1 = self.typeFrequencyLowerCI[permu[0]][i]
                    low2 = self.typeFrequencyLowerCI[permu[1]][i]
                    low3 = self.typeFrequencyLowerCI[permu[2]][i]
                    high1 = self.typeFrequencyUpperCI[permu[0]][i]
                    high2 = self.typeFrequencyUpperCI[permu[1]][i]
                    high3 = self.typeFrequencyUpperCI[permu[2]][i]
                    mean1 = self.typeFrequencyMeanCI[permu[0]][i]
                    mean2 = self.typeFrequencyMeanCI[permu[1]][i]
                    mean3 = self.typeFrequencyMeanCI[permu[2]][i]

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

                        logger.debug("\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f" % ass1,ass2,ass3)
                        print "\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f \n\r" % ass1,ass2,ass3
                        print "I must quit. Debugging required.\n\r"
                        sys.exit()

                    logger.debug("Comparison12: %s Comparison23: %s", comparison12,comparison23)

            comparison = comparison12 + comparison23
            test = re.compile('DU').search(comparison)
            if test not in (None, ""):
                error +=1

            if error == 0:
                # uses NetworkX
                net = nx.Graph(name=str(numberOfTriplets), GraphID=str(numberOfTriplets), End1=permu[0], End2=permu[2], Middle=permu[1])
                net.add_node(permu[0], name=permu[0], site="end", end=1, connectedTo=permu[1] )
                net.add_node(permu[1], name=permu[1], site="middle", end=0, connectedTo="middle")
                net.add_node(permu[2], name=permu[2], site="end", end=1, connectedTo=permu[1])
                net.add_edge(permu[0], permu[1],weight=comparison12, GraphID=numberOfTriplets,end=1)
                net.add_edge(permu[2], permu[1],weight=comparison23, GraphID=numberOfTriplets,end=1)
                logger.debug("VALID TRIPLE SOLUTION: %s * %s * %s " , permu[0],permu[1], permu[2])
                logger.debug("VALID TRIPLE SOLUTION: %s  <--->   %s", comparison12, comparison23)
                logger.debug("VALID TRIPLE SOLUTION: %s ", net.adjacency_list())
                path = nx.shortest_path(net, source=permu[0], target=permu[2])
                logger.debug("VALID TRIPLE SOLUTION: Ends are: %s and %s",permu[0],permu[2])
                logger.debug("VALID TRIPLE SOLUTION: Shortest Path: %s ", path)
                triples.append( net )
                numberOfTriplets += 1
                logger.debug("Current number of triplets: %d", numberOfTriplets)
        return triples

    def filter_list(self,full_list, excludes):
        s = set(excludes)
        return (x for x in full_list if x not in s)


    def checkForValidAdditionsToNetwork(self,nnetwork, pairGraph, solutionCount, args):

        logger.debug("######################Starting check for solution %s with %s nodes ######################################",nnetwork.graph['GraphID'],len(nnetwork))
        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40, "STEP: Testing for addition to seriation ....      ")
            self.scr.refresh()

        logger.debug("The end of assemblages of network %d are: %s and %s", nnetwork.graph['GraphID'], nnetwork.graph["End1"] , nnetwork.graph["End2"])
        logger.debug("Network:  %s", nnetwork.adjacency_list())
        logger.debug("Seriation %d to evaluate: Shortest Path: %s ", nnetwork.graph['GraphID'], nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
        array_of_new_networks=[]  ## a list of all the valid new networks that we run into
        maxnodes=len(nnetwork.nodes())

        for assEnd in ("End1","End2"):
            if assEnd=="End1":
                otherEnd="End2"
            else:
                otherEnd="End1"

            endAssemblage=nnetwork.graph[assEnd]
            logger.debug(">>>>>> Checking ends of seriation %d:  %s is %s", nnetwork.graph['GraphID'], assEnd,endAssemblage)
            list1 = self.validComparisonsHash[ endAssemblage ]
            list2 = nnetwork.nodes()
            logger.debug("List 1 (valid comparisons): %s", list1)
            logger.debug("List 2 (existing nodes): %s", list2)

            validAssemblages = list(self.filter_list(list1, list2))
            logger.debug("Valid assemblages: %s", validAssemblages)
            logger.debug("The list of valid comparative assemblages for %s is %s",endAssemblage,validAssemblages)

            ######################################################################################
            for testAssemblage in validAssemblages:
                logger.debug(" Now checking %s to see if we can be put next to %s",testAssemblage,endAssemblage)
                if args['screen'] not in (None, ""):
                    msg = "Now checking %s against %s." % (testAssemblage, endAssemblage)
                    self.scr.addstr(3,0, msg)
                    self.scr.refresh()
                    c = self.scr.getch()
                    if c == ord('q'):
                        curses.endwin()
                        curses.resetty()
                        os.system("reset")
                        sys.exit("Quitting as requested.\n\r")

                ## now see if the test assemblages fits on the end.
                logger.debug("Checking assemblage %s to see if it fits on the end of the current solution.", testAssemblage )

                #### FIND INNER EDGE RELATIVE TO THE EXISTING END ASSEMBLAGE ##############
                #neighbors = nnetwork.neighbors(endAssemblage)

                logger.debug("Seriation %d with this %s: %s has this many nodes: %d", nnetwork.graph['GraphID'],assEnd,nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]), len(nnetwork.nodes()))
                logger.debug("End assemblage for this seriation: %s",nnetwork.graph[assEnd])
                logger.debug("Which end: %s", assEnd)
                innerNeighbor = None
                if assEnd=="End1":
                    path = nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                    innerNeighbor=path[1]
                    logger.debug("End1: %s Neighbor: %s", endAssemblage, innerNeighbor)
                elif assEnd=="End2":
                    path = nx.shortest_path(nnetwork, nnetwork.graph["End2"] , nnetwork.graph["End1"])
                    innerNeighbor=path[1]
                    logger.debug("End2: %s Neighbor: %s", endAssemblage, innerNeighbor)
                else: ## sanity check
                    print "\r\n\r\n\r\nSomething is wrong finding the next assemblage over.. Error!\n\r"
                    print "\r\nWe are testing endAssemblage: %s "% endAssemblage
                    print "\r\n with neighbors:", innerNeighbor
                    print "For this network:  ", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                    #print nx.write_adjlist(nnetwork,sys.stdout) # write adjacency list to screen
                    sys.exit("Quitting due to errors.")
                ## Sanity check
                if innerNeighbor is None:
                    print "\r\n\r\n\r\nSomething is wrong finding the next assemblage over.. Error!\n\r"
                    print "\r\nWe are testing endAssemblage: %s "% endAssemblage
                    print "\r\n with neighbor:", innerNeighbor
                    print "For this network:  ", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                    #print nx.write_adjlist(nnetwork,sys.stdout) # write adjacency list to screen
                    sys.exit("Quitting due to errors.")
                logger.debug( "\t\t\tThere should be just 1 neighbor to %s and that is: %s", endAssemblage, innerNeighbor )
                c = pairGraph.get_edge_data(innerNeighbor,endAssemblage )
                comparison=c['weight']
                logger.debug( "\t\t\tCompare current pair with previous comparison: %s", comparison)
                ##########################################################################
                comparisonMap =""
                oneToColumns=range(len(self.assemblages[testAssemblage]))
                logger.debug("Number of columns to check: %d", len(oneToColumns))

                error = 0  ## set the error check to 0
                for i in oneToColumns:
                    logger.debug( "\t\t\tComparing Assemblage: %s  and    Assemblage: %s  ########",testAssemblage,endAssemblage)
                    logger.debug( "\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,i,i,i,i,i,i)
                    c=""
                    p=nx.shortest_path(nnetwork, nnetwork.graph[assEnd] , nnetwork.graph[otherEnd])
                    logger.debug( "Working on path: %s",p)
                    newVal=self.assemblages[testAssemblage][i]
                    logger.debug("Start comparison with %s",testAssemblage)
                    previousAssemblage=testAssemblage
                    for compareAssemblage in p:
                        oldVal=self.assemblages[compareAssemblage][i]
                        logger.debug("Compare %s with %s ", previousAssemblage,compareAssemblage)
                        logger.debug("Old value: %f  vs new value: %f",oldVal,newVal)
                        if args['bootstrapCI'] not in (None, ""):
                            upperCI_test = self.typeFrequencyUpperCI[previousAssemblage][i]
                            lowerCI_test = self.typeFrequencyLowerCI[previousAssemblage][i]
                            upperCI_end = self.typeFrequencyUpperCI[compareAssemblage][i]
                            lowerCI_end = self.typeFrequencyLowerCI[compareAssemblage][i]
                            mean_test = self.typeFrequencyMeanCI[previousAssemblage][i]
                            mean_end = self.typeFrequencyMeanCI[compareAssemblage][i]

                            if upperCI_test < lowerCI_end:
                                c += "D"
                            elif lowerCI_test > upperCI_end:
                                c += "U"
                            else:
                                c += "M"
                        else:
                            logger.debug("Outer value: %f Inner value: %f", oldVal, newVal)
                            if newVal<oldVal:
                                c += "U"
                                c1="U"
                            elif newVal>oldVal:
                                c += "D"
                                c1 = "U"
                            elif newVal == oldVal:
                                c += "M"
                                c1="U"
                            else:
                                logger.debug("Error. Quitting.")
                                sys.exit("got null value in comparison of value for type %d in the comparison of %s", i, compareAssemblage)
                            logger.debug("Comparison %s is now %s", c1,c)
                            newVal=oldVal

                        previousAssemblage=compareAssemblage

                    test = re.compile('DU|DM*U').search(c)
                    if test not in (None, ""):
                        logger.debug("Comparison is %s. Error!",c)
                        error +=1

                logger.debug("Checked out %s. Found %d total errors.", testAssemblage, error)
                if error == 0:
                    logger.debug("Found no errors!  Going to add %s to end of existing network at %s", testAssemblage, endAssemblage)
                    logger.debug( "Original network: %s ",nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
                    logger.debug( "New comparison map is: %s ", comparisonMap)

                    new_network = nnetwork.copy()
                    new_network.graph["GraphID"]= str(solutionCount + 1)
                    new_network.graph["name"]=str(solutionCount + 1)
                    logger.debug( "Here's the new network (before addition): %s", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
                    logger.debug("From %s the ends of the seriation are %d (before): %s and %s",assEnd, nnetwork.graph['GraphID'],nnetwork.graph["End1"],nnetwork.graph["End2"] )
                    path = nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                    logger.debug(" New network shortest path (before): %s ", path)

                    ## mark this vertice as the new "END"
                    new_network.add_node(testAssemblage, name=testAssemblage,end=1,site="end")
                    ## mark the interior vertice as not "END
                    new_network.add_node(endAssemblage, name=endAssemblage,site="middle", end=0)

                    #### This adds the comparison to the new edge that has been added.
                    new_network.add_edge( testAssemblage, endAssemblage, weight=comparisonMap, end=1, site="end", GraphID=solutionCount )
                    logger.debug("Ends of the seriation %d (before): %s and %s ",new_network.graph['GraphID'], new_network.graph["End1"],new_network.graph["End2"] )

                    logger.debug("Reassigning the new end %s from %s to %s", assEnd,new_network.graph[assEnd],testAssemblage )
                    new_network.graph[assEnd]=testAssemblage
                    logger.debug("From %s end of the seriation %s (after): %s and %s",assEnd, new_network.graph['GraphID'],new_network.graph["End1"],new_network.graph["End2"] )
                    logger.debug("Here's the new network %s (with addition): %s", new_network.graph['GraphID'], new_network.adjacency_list())
                    path = nx.shortest_path(new_network, new_network.graph["End1"] , new_network.graph["End2"])
                    logger.debug("New network %d shortest path (after): %s ", new_network.graph['GraphID'], path)

                    ## copy this solution to the new array of networks
                    array_of_new_networks.append(new_network)

                    if len(new_network)> maxnodes:
                        maxnodes = len(new_network)
                logger.debug( "----------------#############-------End of check for %s ---------#############-----------------",testAssemblage)
            logger.debug("--------------------------------------Finished with %s-----------------------------------------------------",assEnd)
        logger.debug("------------------------------- Finished with Both Ends-----------------------------------------------------------------")

        if len(array_of_new_networks)>0:
            return array_of_new_networks,maxnodes,
        else:
            return False,0

    def iso(self,G1, glist):
        """Quick and dirty nonisomorphism checker used to check isomorphisms."""
        for G2 in glist:
            if isomorphic(G1,G2):
                return True
        return False

    def MST(self, sGraph,filename,args):

        plt.rcParams['text.usetex'] = False
        plt.figure(filename,figsize=(8,8))
        M=nx.minimum_spanning_tree(sGraph)

        os.environ["PATH"] += ":/usr/local/bin:"
        pos=nx.graphviz_layout(M)
        #pos=nx.graphviz_layout(M,prog="twopi",root=args['graphroot'])
        edgewidth=[]
        weights = nx.get_edge_attributes(M, 'weight')
        for w in weights:
            edgewidth.append(weights[w])
        maxValue = max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)
        assemblageSizes=[]
        sizes = nx.get_node_attributes(M, 'size')
        for s in sizes:
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(M,pos,alpha=0.3,width=widths)
        sizes = nx.get_node_attributes(M,'size')
        nx.draw_networkx_nodes(M,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
        nx.draw_networkx_edges(M,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
        nx.draw_networkx_labels(M,pos,fontsize=10)
        font = {'fontname'   : 'Helvetica',
            'color'      : 'k',
            'fontweight' : 'bold',
            'fontsize'   : 10}
        plt.axis('off')
        plt.savefig(filename,dpi=75)
        self.saveGraph(sGraph,filename+".gml",args)
        if args['shapefile'] is not None and args['xyfile'] is not None:
            self.createShapefile(M,filename+".shp",args)

    def minimumSpanningTree(self,networks,sumGraph, outputDirectory,inputFile):
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError("This function needs Graphviz and either PyGraphviz or Pydot")

        newfilename=outputDirectory+inputFile[0:-4]+"-mst.png"
        plt.figure(newfilename,figsize=(8,8))

        graphs=[]
        megaGraph = nx.Graph()
        number=0
        graphCount=0
        for net in networks:
            graphCount +=1
            g = nx.Graph()
            number = net.graph['GraphID']
            for nodey in net.nodes(data=True):
                xCoordinate = 0
                yCoordinate = 0
                name = nodey[0]
                xCoordinate = self.xAssemblage[name]
                yCoordinate = self.yAssemblage[name]
                megaGraph.add_node(name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                                   size=self.assemblageSize[name])

            count=0
            for e in net.edges_iter():
                d = net.get_edge_data(*e)
                fromAssemblage = e[0]
                toAssemblage = e[1]
                g.add_node(fromAssemblage, label=fromAssemblage, x=xCoordinate, y=yCoordinate,
                                            name=fromAssemblage, size=self.assemblageSize[name])
                g.add_node(toAssemblage, label=toAssemblage, x=xCoordinate, y=yCoordinate,
                                            name=toAssemblage, size=self.assemblageSize[name])

                weight = d['weight']+1
                distance = self.distanceBetweenAssemblages[fromAssemblage + "*" + toAssemblage]
                #count = megaGraph.get_edge_data(fromAssemblage,toAssemblage,'weight'
                count += 1
                megaGraph.add_path([fromAssemblage, toAssemblage], weight=count,
                                   distance=distance, color=number,
                                   size=(self.assemblageSize[fromAssemblage], self.assemblageSize[toAssemblage]))

                g.add_path([fromAssemblage, toAssemblage],
                                            xy1=(self.xAssemblage[fromAssemblage], self.yAssemblage[fromAssemblage]),
                                            xy2=(self.xAssemblage[toAssemblage], self.yAssemblage[toAssemblage]),
                                            weight=weight,
                                            meanDistance=distance,
                                            size=(self.assemblageSize[fromAssemblage], self.assemblageSize[toAssemblage]))
            graphs.append(g)
        plt.rcParams['text.usetex'] = False
        mst=nx.minimum_spanning_tree(megaGraph,weight='weight')
        os.environ["PATH"] += ":/usr/local/bin:"
        pos=nx.graphviz_layout(mst)
        edgewidth=[]
        weights = nx.get_edge_attributes(mst, 'weight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)

        color = nx.get_edge_attributes(mst, 'color')
        colorList = []
        for c in color:
            colorList.append(color[c])
        colors=[]
        colorMax = max(colorList)
        for c in colorList:
            colors.append(c/colorMax)
        assemblageSizes=[]
        sizes = nx.get_node_attributes(mst, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(mst,pos,alpha=0.3,width=widths, edge_color=colorList)
        sizes = nx.get_node_attributes(mst,'size')
        nx.draw_networkx_nodes(mst,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
        nx.draw_networkx_edges(mst,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
        nx.draw_networkx_labels(mst,pos,fontsize=10)
        font = {'fontname'   : 'Helvetica',
            'color'      : 'k',
            'fontweight' : 'bold',
            'fontsize'   : 10}

        plt.axis('off')
        plt.savefig(newfilename,dpi=75)
        if args['shapefile'] is not None and args['xyfile'] is not None:
            self.createShapefile(mst,outputDirectory+inputFile[0:-4]+"-mst.shp",args)
        self.saveGraph(mst,newfilename+".gml",args)
        atlasFile=outputDirectory+inputFile[0:-4]+"-atlas.png"
        plt.figure(atlasFile,figsize=(8,8))
        UU=nx.Graph()
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist=[] # list of nonisomorphic graphs
        for G in graphs:
            # check against all nonisomorphic graphs so far
            if not self.iso(G, nlist):
                nlist.append(G)

        UU=nx.disjoint_union_all(graphs) # union the nonisomorphic graphs
        pos=nx.graphviz_layout(UU,prog="twopi",root=args['graphroot'])
        # color nodes the same in each connected subgraph
        C=nx.connected_component_subgraphs(UU)
        for g in C:
            c = [random.random()] * nx.number_of_nodes(g) # random color...
            nx.draw(g,
                pos,
                node_size=40,
                node_color=c,
                vmin=0.0,
                vmax=1.0,
                alpha=.2,
                font_size=7,
            )
        plt.savefig(atlasFile,dpi=250)
        self.saveGraph(UU,atlasFile+".gml",args)


    def finalGoodbye(self,maxNodes,frequencyTotal,continuityTotal,args):
        if args['screen'] is not None:
            curses.endwin()
            curses.resetty()
            curses.nl()
            curses.echo()
        ## determine time elapsed
        #time.sleep(5)
        timeNow = time.time()
        timeElapsed = timeNow-self.start
        print "Seriation complete."
        print "Maximum size of seriation: %d" % maxNodes
        print "Number of frequency seriation solutions at last step: %d" % frequencyTotal
        print "Number of continuity seriation solutions at end: %d " % continuityTotal
        print "Time elapsed for calculation: %d seconds" % timeElapsed
        if args['screen'] is not None:
            os.system("reset")

    #################################################### set up all the output files ####################################################
    def setupOutput(self,args):
        outputFile = self.outputDirectory + self.inputFile[0:-4]+".vna"
        OUTMSTFILE=OUTMSTDISTANCEFILE=""
        try:
            OUTFILE = open(outputFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outputFile, e
            sys.exit(msg)

        outpairsFile = self.outputDirectory +self.inputFile[0:-4]+"-pairs.vna"
        try:
            OUTPAIRSFILE = open(outpairsFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outpairsFile, e
            sys.exit(msg)


        outmstFile=  self.outputDirectory + self.inputFile[0:-4] + "-mst.vna"
        outmst2File = self.outputDirectory + self.inputFile[0:-4] + "-mst-distance.vna"

        if args['mst'] not in (None, ""):
            try:
                OUTMSTFILE = open(outmstFile, 'w')
                OUTMSTDISTANCEFILE = open(outmst2File, 'w')
            except csv.Error as e:
                msg = "Can't open file %s to write: %s" % outputFile, e
                sys.exit(msg)

        sumgraphOutputFile = self.outputDirectory + self.inputFile[0:-4]+"-sumgraph.vna"
        try:
            SUMGRAPH = open(sumgraphOutputFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % sumgraphOutputFile, e
            sys.exit(msg)

        return OUTFILE,OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,SUMGRAPH

    #################################################### sort by multiple keys ####################################################
    def multikeysort(self,items, columns):
        from operator import itemgetter
        comparers = [ ((itemgetter(col[1:].strip()), -1) if col.startswith('-') else (itemgetter(col.strip()), 1)) for col in columns]
        def comparer(left, right):
            for fn, mult in comparers:
                result = cmp(fn(left), fn(right))
                if result:
                    return mult * result
            else:
                return 0
        return sorted(items, cmp=comparer)

    def createShapefile(self,graph,shapefilename,args):
        w = shapefile.Writer(shapefile.POLYLINE)  # 3= polylines
        xCoordinates=nx.get_node_attributes(graph,"xCoordinate")
        yCoordinates=nx.get_node_attributes(graph,"yCoordinate")
        for e in graph.edges_iter():
            d = graph.get_edge_data(*e)
            node1 = e[0]
            node2 = e[1]
            x1 = float(xCoordinates[node1])
            y1 = float(yCoordinates[node1])
            x2 = float(xCoordinates[node2])
            y2 = float(yCoordinates[node2])
            w.poly(parts=[[[x1,y1],[x2,y2]]])
        w.save(shapefilename)

    def createAtlasOfSolutions(self,filteredarray,args):
        atlasGraph=nx.disjoint_union_all(filteredarray)
        pos=nx.graphviz_layout(atlasGraph,prog="twopi",root=args['graphroot'])
        atlasFile=self.outputDirectory + self.inputFile[0:-4]+"-new-atlas.png"
        plt.savefig(atlasFile,dpi=250)
        plt.show() # display

    def createAtlas(self,filteredarray,args):
        # remove isolated nodes, only connected graphs are left
        U=nx.Graph() # graph for union of all graphs in atlas
        for G in filteredarray:
            U=nx.disjoint_union(U,G)
        # list of graphs of all connected components
        C=nx.connected_component_subgraphs(U)

        UU=nx.Graph()
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist=[] # list of nonisomorphic graphs
        for G in C:
            # check against all nonisomorphic graphs so far
            if not self.iso(G,nlist):
                nlist.append(G)
                UU=nx.disjoint_union(UU,G) # union the nonisomorphic graphs
        return UU

    def outputGraphArray(self,array,args):
        num=0
        os.environ["PATH"] += ":/usr/local/bin:"
        for g in array:
            num +=1
            pos=nx.graphviz_layout(g,prog="twopi",root=['graphroot'])
            gfile=self.outputDirectory + self.inputFile[0:-4]+"-min-sol-"+str(num)+".png"
            filename=self.outputDirectory + self.inputFile[0:-4]+"-min-sol-"+str(num)+".gml"
            self.saveGraph(g,filename,args)
            edgewidth=[]
            weights = nx.get_edge_attributes(g, 'weight')
            for w in weights:
                edgewidth.append(weights[w])

            maxValue = max(edgewidth)
            widths=[]
            for w in edgewidth:
                widths.append(((maxValue-w)+1)*5)

            assemblageSizes=[]
            sizes = nx.get_node_attributes(g, 'size')
            #print sizes
            for s in sizes:
                assemblageSizes.append(sizes[s])
            nx.draw_networkx_edges(g,pos,alpha=0.3,width=widths)
            sizes = nx.get_node_attributes(g,'size')
            nx.draw_networkx_nodes(g,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
            nx.draw_networkx_edges(g,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
            nx.draw_networkx_labels(g,pos,fontsize=10)
            font = {'fontname'   : 'Helvetica',
                'color'      : 'k',
                'fontweight' : 'bold',
                'fontsize'   : 10}
            plt.axis('off')
            plt.savefig(gfile,dpi=75)
            plt.figure(gfile,figsize=(8,8))

    def sumGraphsByWeight(self,filteredarray,args):
        sumGraph=nx.Graph()
        ## go through all the graphs
        for g in filteredarray:
            ## go through all the edges for each graph
            for node in g.nodes(data=True):
                xCoordinate = 0
                yCoordinate = 0
                name = node[0]
                if args['xyfile'] is not None:
                    xCoordinate = self.xAssemblage[name]
                    yCoordinate = self.yAssemblage[name]
                sumGraph.add_node(name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                                   size=self.assemblageSize[name])

            maxWeight=0
            for e in g.edges_iter():
                d = g.get_edge_data(*e)
                fromAssemblage = e[0]
                toAssemblage = e[1]
                exists = False
                currentWeight=1
                for e in sumGraph.edges():
                    dd = sumGraph.get_edge_data(*e)
                    if fromAssemblage in e and toAssemblage in e:   ## if exists
                        exists = True
                    currentWeight=1
                    if exists is True:
                        currentWeight = int(dd['weight']) + 1

                if currentWeight > maxWeight:
                    maxWeight=currentWeight
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight)

            for e in sumGraph.edges_iter():
                d = sumGraph.get_edge_data(*e)
                currentWeight=int(d['weight'])
                inverseWeight=(maxWeight+1)-currentWeight
                fromAssemblage = e[0]
                toAssemblage = e[1]
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight,inverseweight=inverseWeight )

        return sumGraph

    def sumGraphsByCount(self,filteredarray,args):
        sumGraph=nx.Graph()
        ## go through all the graphs
        for g in filteredarray:
            ## go through all the edges for each graph
            for node in g.nodes(data=True):
                xCoordinate = 0
                yCoordinate = 0
                name = node[0]
                if args['xyfile'] is not None:
                    xCoordinate = self.xAssemblage[name]
                    yCoordinate = self.yAssemblage[name]
                sumGraph.add_node(name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                                   size=self.assemblageSize[name])

            maxWeight=0
            for e in g.edges_iter():
                d = g.get_edge_data(*e)
                fromAssemblage = e[0]
                toAssemblage = e[1]
                exists = False
                currentWeight=1
                for e in sumGraph.edges():
                    dd = sumGraph.get_edge_data(*e)
                    if fromAssemblage in e and toAssemblage in e:   ## if exists
                        exists = True
                    currentWeight=1
                    if exists is True:
                        currentWeight += 1

                if currentWeight > maxWeight:
                    maxWeight=currentWeight
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight)

            for e in sumGraph.edges_iter():
                d = sumGraph.get_edge_data(*e)
                currentWeight=int(d['weight'])
                inverseWeight=(maxWeight+1)-currentWeight
                fromAssemblage = e[0]
                toAssemblage = e[1]
                sumGraph.add_path([fromAssemblage, toAssemblage], weight=currentWeight,inverseweight=inverseWeight )

        return sumGraph

    def calculateSumOfDifferences(self,assemblage1,assemblage2,args):
        diff=0
        for type in range(0,self.numberOfClasses):
                diff += abs(float(self.assemblageFrequencies[assemblage1][type]) - float(self.assemblageFrequencies[assemblage2][type]))
        return diff

    def continunityMaximizationSeriation(self,args):
        graphList=[]
        numGraphs=0
        if args['continuityroot'] is not None:
            numGraphs +=1
            g=nx.Graph(startAssemblage=args['continuityroot'], End1=args['continuityroot'])
            g.add_node(args['continuityroot'],size=self.assemblageSize[args['continuityroot']],
                       xCoordinate=self.xAssemblage[args['continuityroot']],yCoordinate=self.yAssemblage[args['continuityroot']])
            graphList.append(g)   ## create a starting graph for each of assemblage put into an array
        else:
            for ass in self.assemblages:
                numGraphs +=1
                g=nx.Graph(startAssemblage=ass, End1=ass)
                g.add_node(ass,size=self.assemblageSize[ass], xCoordinate=self.xAssemblage[ass],yCoordinate=self.yAssemblage[ass])
                graphList.append(g)   ## create a starting graph for each of assemblage put into an array

        ## special case for the first time through
        for g in graphList:
            minMatch = 10000000
            currentMinimumMatch=""
            nodelist = g.nodes()
            for node in nodelist:
                for b in self.assemblages:
                    if b is not node:
                        diff = self.calculateSumOfDifferences(node,b,args)
                        if diff < minMatch:
                            minMatch = diff
                            currentMinimumMatch=b

            g.add_node(currentMinimumMatch, xCoordinate=self.xAssemblage[currentMinimumMatch],
                       yCoordinate=self.xAssemblage[currentMinimumMatch], size=self.assemblageSize[currentMinimumMatch])
            if minMatch==0:
                minMatch=10000000

            g.add_path([node, currentMinimumMatch], weight=minMatch, inverseweight=(1/minMatch ))
            g.graph['End2']=currentMinimumMatch

        ## Now go through list looking at each one and increasing as long as I can. Add graphs when there are equivalent solutions
        for g in graphList:
            for assemblage in self.assemblages:
                globalMinMatch = 100000
                endMinMatch={"End1":10000,"End2":10000}
                currentMinimumMatch={}
                matchEnd={}
                matchEndAssemblage={}   ## this will contain the end assemblages and the differences
                match=False
                smallestMatchEnd=[]
                assemblagesMatchedToEnd=[]

                ## examine both ends to see which is the smallest summed difference.
                for assEnd in ("End1","End2"):
                    if assEnd=="End1":
                        otherEnd="End2"
                    else:
                        otherEnd="End1"

                    ## only go from one side if you use declare a root
                    if args['continuityroot'] is not None:
                        assEnd=="End2"

                    ## set the current end assemblages
                    endAssemblage=g.graph[assEnd]

                    for a in self.assemblages:
                        if a not in g.nodes():
                            diff = self.calculateSumOfDifferences(endAssemblage,a,args)
                            if diff < globalMinMatch:
                                match=True
                                globalMinMatch = diff
                                currentMinimumMatch[assEnd] = a
                                matchEnd=assEnd
                                matchEndAssemblage[assEnd]=endAssemblage

                ## at this point we should have the minimum distance match for each end.
                ## we then need to compare each end to find which one is the smallest
                ## three possibilities -- end1, end2 and both (i.e., the diff is the same)

                if minMatch['End1'] < minMatch['End2']:
                    smallestMatchEnd.append('End1')
                    assemblagesMatchedToEnd.append(matchEndAssemblage['End1'])

                elif minMatch['End2']< minMatch['End1']:
                    smallestMatchEnd.append('End2')
                    assemblagesMatchedToEnd.append(matchEndAssemblage['End2'])
                else:
                    smallestMatchEnd.append('End1')
                    smallestMatchEnd.append('End2')
                    assemblagesMatchedToEnd.append(matchEndAssemblage['End1'])
                    assemblagesMatchedToEnd.append(matchEndAssemblage['End2'])

                ## find out if there are others that have the same minimum value
                for a in self.assemblages:
                    if a not in g.nodes() and a is not endAssemblage and a not in assemblagesMatchedToEnd:
                        diff = self.calculateSumOfDifferences(a,endAssemblage,args)
                        if diff == globalMinMatch:
                            ## add this as a matched equivalent assemblage. We will then deal with more than one match
                            assemblagesMatchedToEnd.append(a)

                firstOne=True
                for match in assemblagesMatchedToEnd:
                    for endAss in smallestMatchEnd:
                        # for the first time we need to simply add it to the right end but after this we copy...
                        if firstOne == True:
                            firstOne=False
                            g.add_node(match, xCoordinate=self.xAssemblage[match],
                               yCoordinate=self.xAssemblage[match],
                                size=self.assemblageSize[match])
                            if globalMinMatch==0:
                                globalMinMatch=10000000
                            g.add_path([matchEndAssemblage[endAss], match], weight=globalMinMatch, inverseweight=(1/globalMinMatch ))
                        ## if there are more than one we need to copy first before adding node
                        else:
                            new_network = g.copy()
                            new_network.add_node(match, xCoordinate=self.xAssemblage[match],
                                                 yCoordinate=self.xAssemblage[match],size=self.assemblageSize[match])
                            if globalMinMatch==0:
                                globalMinMatch=10000000
                            new_network.add_path([matchEndAssemblage[endAss], match], weight=globalMinMatch, inverseweight=(1/globalMinMatch ))
                            graphList.append(new_network)
                            numGraphs += 1
        return graphList

    ## Output to file and to the screen
    def graphOutput(self,sumGraph,sumgraphfilename, args):
        ## Now make the graphic for set of graphs
        plt.rcParams['text.usetex'] = False
        newfilename=self.outputDirectory+sumgraphfilename
        gmlfilename=self.outputDirectory+sumgraphfilename+".gml"
        self.saveGraph(sumGraph,gmlfilename,args)
        if args['shapefile'] is not None and args['xyfile'] is not None:
            self.createShapefile(sumGraph,newfilename+".shp",args)
        plt.figure(newfilename,figsize=(8,8))
        os.environ["PATH"] += ":/usr/local/bin:"
        pos=nx.graphviz_layout(sumGraph)
        #pos=nx.graphviz_layout(sumGraph,prog="twopi",root=['graphroot'])
        #pos=nx.spring_layout(mst,iterations=500)
        edgewidth=[]
        weights = nx.get_edge_attributes(sumGraph, 'weight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)

        assemblageSizes=[]
        sizes = nx.get_node_attributes(sumGraph, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(sumGraph,pos,alpha=0.3,width=widths)
        sizes = nx.get_node_attributes(sumGraph,'size')
        nx.draw_networkx_nodes(sumGraph,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
        nx.draw_networkx_edges(sumGraph,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
        nx.draw_networkx_labels(sumGraph,pos,fontsize=10)
        font = {'fontname'   : 'Helvetica',
            'color'      : 'k',
            'fontweight' : 'bold',
            'fontsize'   : 10}
        plt.axis('off')
        plt.savefig(newfilename,dpi=75)
        self.saveGraph(sumGraph,newfilename+".gml",args)


    ## Output to file and to the screen
    def sumGraphOutput(self,sumGraph,SUMGRAPH,sumgraphfilename, args):


        nodeList = sumGraph.nodes()
        for a in self.assemblages:
            if a not in nodeList:
                sumGraph.add_node(a,  xCoordinate=self.xAssemblage[a], yCoordinate=self.yAssemblage[a],
                                   size=self.assemblageSize[a])
        SUMGRAPH.write( "*Node data\n")
        SUMGRAPH.write("ID AssemblageSize X Y Easting Northing\n")

        for node in sumGraph.nodes(data=True):
            nodeName =node[0]
            x = 0
            y = 0
            northing = 0
            easting = 0
            if args['xyfile'] is not None:
                x = float(self.xAssemblage[ nodeName ]) / 1000000.0
                y = (float(self.largestY)- float(self.yAssemblage[nodeName]))/100000.0
                easting = self.xAssemblage[nodeName]
                northing = self.yAssemblage[nodeName]
            msg = nodeName + " "+ str(self.assemblageSize[ nodeName ])+" "+ str(x)+" "+str(y)+" "+str(easting)+" "+str(northing)+"\n"
            SUMGRAPH.write(msg)
        SUMGRAPH.write("*Tie data\nFrom To Edge Weight InverseWeight\n")
        edgeCount=0
        for e in sumGraph.edges_iter():
            d = sumGraph.get_edge_data(*e)
            edgeCount += 1
            text = e[0]+ " "+ e[1]+" "+ str(edgeCount)+ " "+ str(d['weight'])+ " "+ str(d['inverseweight'])+ "\n"
            SUMGRAPH.write(text)

        ## Now make the graphic for the sumgraph
        newfilename=self.outputDirectory+self.inputFile[0:-4]+"-sumgraph.png"
        self.saveGraph(sumGraph,newfilename+".gml",args)
        plt.figure(newfilename,figsize=(8,8))
        plt.rcParams['text.usetex'] = False
        os.environ["PATH"] += ":/usr/local/bin:"
        #pos=nx.graphviz_layout(sumGraph,prog="twopi",root=['graphroot'])
        pos=nx.graphviz_layout(sumGraph)

        edgewidth=[]
        weights = nx.get_edge_attributes(sumGraph, 'weight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)

        assemblageSizes=[]
        sizes = nx.get_node_attributes(sumGraph, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(sumGraph,pos,alpha=0.3,width=widths)
        sizes = nx.get_node_attributes(sumGraph,'size')
        nx.draw_networkx_nodes(sumGraph,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
        nx.draw_networkx_edges(sumGraph,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
        nx.draw_networkx_labels(sumGraph,pos,fontsize=10)
        font = {'fontname'   : 'Helvetica',
            'color'      : 'k',
            'fontweight' : 'bold',
            'fontsize'   : 10}
        plt.axis('off')
        plt.savefig(newfilename,dpi=75)
        if args['shapefile'] is not None and args['xyfile'] is not None:
            self.createShapefile(sumGraph,self.outputDirectory+self.inputFile[0:-4]+"-sumgraph.shp",args)

        newfilename=self.outputDirectory+sumgraphfilename
        plt.figure(newfilename,figsize=(8,8))
        mst=nx.minimum_spanning_tree(sumGraph,weight='inverseweight')

        plt.rcParams['text.usetex'] = False

        pos=nx.graphviz_layout(mst)
        #pos=nx.graphviz_layout(mst,prog="twopi",root=['graphroot'])
        edgewidth=[]
        weights = nx.get_edge_attributes(mst, 'weight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)

        assemblageSizes=[]
        sizes = nx.get_node_attributes(mst, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(mst,pos,alpha=0.3,width=widths)
        sizes = nx.get_node_attributes(mst,'size')
        nx.draw_networkx_nodes(mst,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
        nx.draw_networkx_edges(mst,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
        nx.draw_networkx_labels(mst,pos,fontsize=10)
        font = {'fontname'   : 'Helvetica',
            'color'      : 'k',
            'fontweight' : 'bold',
            'fontsize'   : 10}
        plt.axis('off')
        plt.savefig(newfilename,dpi=75)
        self.saveGraph(mst,newfilename+".gml",args)
        if args['shapefile'] is not None and args['xyfile'] is not None:
            self.createShapefile(mst,self.outputDirectory+sumgraphfilename+".shp",args)

        # layout graphs with positions using graphviz neato

    #################################################### OUTPUT SECTION ####################################################
    def output(self,filteredArray,OUTFILE,OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,maxEdges,args):
        if args['screen'] not in (None, ""):
            self.scr.addstr(13,1, "Now printing output file... ")
            self.scr.addstr(1,40,"STEP: Output files...         ")
            self.scr.refresh()

        OUTFILE.write( "*Node data\n")
        OUTFILE.write("ID AssemblageSize X Y Easting Northing\n")
        OUTPAIRSFILE.write("*Node data\n")
        OUTPAIRSFILE.write("ID AssemblageSize X Y Easting Northing\n")
        count = 0
        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40,"STEP: Printing list of nodes....     ")
            self.scr.refresh()
        ## note this assumes the use of UTM coordinates (northing and easting)
        for l in self.assemblages:
            x = 0
            y = 0
            northing = 0
            easting = 0
            if args['xyfile'] not in (None, ""):
                x = float(self.xAssemblage[ l ]) / 1000000.0
                y = (float(self.largestY)- float(self.yAssemblage[l]))/100000.0
                easting = self.xAssemblage[l]
                northing = self.yAssemblage[l]

            msg = l + " "+ str(self.assemblageSize[ l ])+" "+ str(x)+" "+str(y)+" "+str(easting)+" "+str(northing)+"\n"
            OUTFILE.write(msg)
            OUTPAIRSFILE.write(msg)
            if args['mst'] not in (None, ""):
                OUTMSTFILE.write(msg)
                OUTMSTDISTANCEFILE.write(msg)

        OUTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
        OUTPAIRSFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

        if args['mst'] not in (None, ""):
            OUTMSTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
            OUTMSTDISTANCEFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40,"STEP: Printing list of nodes attributes... ")
            self.scr.refresh()
        for l in self.assemblages:
            easting = 0
            northing = 0
            x = 0
            y = 0
            if args['xyfile'] not in (None, ""):
                x = float(self.xAssemblage[l])/1000000
                y = (float(self.largestY)-float(self.yAssemblage[l]))/100000
                easting = self.xAssemblage[l]
                northing = self.yAssemblage[l]
            msg = l +" "+ str(self.assemblageSize[ l])+" "+str(x)+" "+str(y)+" "+str(easting)+" "+str(northing)+"\n"
            OUTFILE.write(msg)
            OUTPAIRSFILE.write(msg)
            if args['mst'] not in (None, ""):
                OUTMSTFILE.write( msg )
                OUTMSTDISTANCEFILE.write(msg)

        ## This prints out counts of the edges as they appear in ALL of the solutions
        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40,"STEP: Going through and counting pairs...     ")
            self.scr.refresh()
        OUTPAIRSFILE.write("*Tie data\nFrom To Edge Count\n")
        if args['mst'] not in (None, ""):
            OUTMSTFILE.write( "*Tie data\nFrom To Edge End Weight ID\n")
            OUTMSTDISTANCEFILE.write("*Tie data\nFrom To Edge End Weight ID\n")

        ## first count up all of the edges by going through the solutions and each edge
        ## put the edge count in a hash of edges
        edgeHash={}

        for network in filteredArray:
            for e in network.edges_iter():
                pairname= e[0]+"*"+e[1]
                edgeHash[ pairname ] = 0

            for e in network.edges_iter():
                pairname= e[0]+"*"+e[1]
                edgeHash[ pairname ] += 1

        ## now go through the edgeHash and print out the edges
        ## do this is sorted order of the counts. For fun.
        if args['screen'] is not None:
            self.scr.addstr(1,40,"STEP: Doing the pair output...                ")
            self.scr.refresh()

        sorted_pairs = sorted(edgeHash.iteritems(), key=operator.itemgetter(1))

        for key,value in sorted_pairs:
            ass1,ass2=key.split("*")
            msg = ass1+"\t"+ass2+"\t" + " 1 "+ str(value) +"\n"
            OUTPAIRSFILE.write(msg)

        OUTFILE.write("*Tie data\nFrom To Edge Weight Network End pValue pError meanSolutionDistance\n")
        if args['screen'] not in (None, ""):
            self.scr.addstr(1,40,"STEP: Eliminating duplicates...     ")
            self.scr.addstr(1,40,"STEP: Printing edges...     ")
            self.scr.refresh()

        uniqueArray = set(filteredArray)
        distanceHash={}
        seriationHash ={}
        ## only print unique ones...
        pairwise={}
        pairwiseError={}
        for network in filteredArray:
            if args['screen'] not in (None, ""):
                self.scr.addstr(14,1, "Now on solution: ")
                self.scr.addstr(14,18,str(network.graph["GraphID"]) )
                #print "now on solution: ", network["GraphID"],"\n"
            if args['largestonly'] not in (None, "") and len(network.edges()) == maxEdges-1:
                edgeCount = len(network.edges())
                groupDistance=0
                meanDistance=0.0
                eCount=0
                if args['xyfile'] not in (None, ""):
                    for e in network.edges_iter():
                      pairname= e[0]+"*"+e[1]
                      groupDistance += self.distanceBetweenAssemblages[ pairname ]
                      eCount += 1
                    meanDistance = groupDistance/eCount         ##use the average distance as the metric
                else:
                    meanDistance = "0"

                ## initialize edges
                for e in network.edges_iter():
                    pairname= e[0]+"#"+e[1]
                    pairwise[ pairname ] = 0
                    pairwiseError[ pairname ] = 0

                for e in network.edges_iter():
                    pVal=0.0
                    pErr=0.0
                    if args['pairwisefile'] is not None:
                        pairname = e[0]+"#"+e[1]
                        pVal = pairwise[ pairname ]
                        pErr = pairwiseError[ pairname ]
                    else:
                        pVal = 0.0
                        pErr = 0.0
                    text = e[0]+ " "+ e[1]+" 1 "+str(edgeCount)+ " "+ str(network.graph["GraphID"])+ " "+\
                           str(pVal)+" "+ str(pErr)+ " "+str(meanDistance)+"\n"
                    OUTFILE.write(text)

                network.graph["meanDistance"]=meanDistance
                distanceHash[ text] = meanDistance

            else:  ## not just the largest, but ALL seriation solutions
                edgeCount = len(network.edges())
                groupDistance=0
                meanDistance=0.0
                eCount=0
                if args['xyfile'] is not None:
                    for e in network.edges_iter():
                      pairname= e[0]+"*"+e[1]
                      groupDistance += self.distanceBetweenAssemblages[ pairname ]
                      eCount += 1
                    meanDistance = groupDistance/eCount         ##use the average distance as the metric
                else:
                    meanDistance = "0"

                ## initialize edges
                for e in network.edges_iter():
                    pairname= e[0]+"#"+e[1]
                    pairwise[ pairname ] = 0
                    pairwiseError[ pairname ] = 0

                for e in network.edges_iter():
                    pVal=0.0
                    pErr=0.0
                    if args['pairwisefile'] is not None:
                        pairname= e[0]+"#"+e[1]
                        pVal = pairwise[ pairname ]
                        pErr = pairwiseError[ pairname ]
                    else:
                        pVal = 0.0
                        pErr = 0.0
                    text = e[0]+ " "+ e[1]+" 1 "+str(edgeCount)+ " "+ str(network.graph["GraphID"])+ " "+\
                           str(pVal)+" "+ str(pErr)+ " "+str(meanDistance)+"\n"
                    OUTFILE.write(text)

                network.graph["meanDistance"]=meanDistance
                distanceHash[ text] = meanDistance

    def filterSolutions(self,end_solutions,all_solutions,args):
        ################################################# FILTERING  ####################################
        # now do some weeding. Basically start with the last network ( largest), and work backwards to smaller and smaller solutions. Ignore any
        # network that is already represented larger since these are trivial (e.g., A->B->C->D already covers
        # A->C->D.) That should then leave all the unique maximum solutions (or so it seems)
        ################################################# FILTERING  ####################################

        filteredarray =[]
        if args['filtered'] not in (None, ""):  ## only get the largest set that includes ALL
            if args['screen'] not in (None, ""):
                self.scr.addstr(1,40,"STEP: Filter to get uniques... ")
            logger.debug("--- Filtering solutions so we only end up with the unique ones.")
            logger.debug("--- Start with %d solutions.", len(end_solutions))
            filteredarray.append(end_solutions[-1])
            newOne=0
            for tnetwork in reversed(all_solutions):
                exists=0
                for fnetwork in filteredarray:
                    fnetworkArray= fnetwork.nodes()
                    logger.debug("----fnetworkArray: %s", fnetworkArray)
                    tnetworkArray = tnetwork.nodes()
                    logger.debug("----tnetworkArray: %s", tnetworkArray)
                    minus = list(set(tnetworkArray) - set(fnetworkArray))
                    logger.debug("difference between: %s ", minus)
                    change= len(minus)
                    logger.debug("Change: %d", change)
                    if change > 0 and len(list(set(minus)-set(fnetworkArray))):
                        newOne += 1
                    else:
                        exists += 1
                if exists==0:
                    logger.debug("pushing tnetwork to list of filtered arrays")
                    filteredarray.append(tnetwork)
                exists=0

            logger.debug("End with %d solutions.", len(filteredarray))
            filterCount= len(filteredarray)
            if args['screen'] not in (None, ""):
                self.scr.addstr(11,1,"End with filterCount solutions.")
        elif args['allsolutions'] not in (None, ""):
            filteredarray = all_solutions  ## all possible networks
        else:
            filteredarray = end_solutions ## just the largest ones (from the last round)
        return filteredarray

    def checkMinimumRequirements(self,args):
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError("This function needs Graphviz and either PyGraphviz or Pydot. Please install GraphViz from http://www.graphviz.org/")
        if args['inputfile'] in (None,""):
            sys.exit("Inputfile is a required input value: --inputfile=../testdata/testdata.txt")

    def addOptions(self,oldargs):

        args={'debug':None ,'bootstrapCI':None,'bootstrapSignificance':None,
              'filtered':None,'largestonly':None,'individualfileoutput':None,
              'excel':None,'threshold':None,'noscreen':None,'xyfile':None,'pairwisefile':None,'mst':None,
              'stats':None,'screen':None,'allsolutions':None,'inputfile':None,'outputdirectory':None,
              'shapefile':None,'frequency':None,'continuity':None,'graphs':None,'graphroot':None,'continuityroot':None}
        for a in oldargs:
            args[a]=oldargs[a]
        return args

    def seriate(self, args):
        args=self.addOptions(args)
        self.checkMinimumRequirements(args)
        #####################################DEBUG OUTPUT#############################################################
        if args['debug'] is not None:
            ## Logging
            logger.basicConfig(stream=sys.stderr, level=logger.DEBUG)
            args['screen']= None
        else:
            logger.basicConfig(stream=sys.stderr, level=logger.ERROR)
            args['screen'] = True

        logger.debug("Arguments: %s", args)

        ##################################################################################################
        if (args['screen'] is not None) and (args['debug'] is None ):
            ## Set up the screen display (default).
            ## the debug option should not use this since it gets messy
            try:
                # Initialize curses
                self.scr=curses.initscr()
                # Turn off echoing of keys, and enter cbreak mode,
                # where no buffering is performed on keyboard input
                curses.noecho()
                curses.cbreak()
                self.scr.nodelay(1)
                self.scr.addstr(0,0,"Iterative Seriation Program V.2.0", curses.A_BOLD)
                self.scr.addstr(20,35,"Hit <q> to quit.")
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

        if args['continuity'] in (None,False,0) and args['frequency'] in (None,False,0):
            sys.exit("You must specify --continuity=1 and/or frequency=1 to set the kind(s) of seriations you would like.")

        ######################################FILE INPUT#############################################################
        filename=args['inputfile']
        if filename is "":
            logger.error("You must enter a filename to continue.")
            print "You must enter a filename to continue."
            sys.exit("Quitting due to errors.")

        try:
            logger.debug("Going to try to open and load: %s", filename)
            self.openFile(filename,args)
        except IOError as e:
            logger.error("Cannot open %s. Error: %s", filename, e.strerror)

            print("Cannot open %s. Error. %s ", filename, e.strerror)
            if args['screen'] not in (None, ""):
                curses.endwin()
                curses.resetty()
            sys.exit("Quitting due to errors.")


        try:
            inputparts =map(str,args['inputfile'].split("/"))
            self.inputFile= inputparts[len(inputparts)-1]
        except:
            sys.exit("There was a problem with parsing the input file. Check it and try again.")

        ############################################################################################################
        if args['outputdirectory']  not in (None, ""):
            self.outputDirectory = args['outputdirectory']
        else:
            self.outputDirectory = "../output/"
        ############################################################################################################
        logger.debug("Going to open pairwise file it is exists.")
        if args['pairwisefile']  not in (None, ""):
            self.openPairwiseFile(args['pairwisefile'])

        ############################################################################################################
        logger.debug("Going to open XY file if it exists.")

        if args['xyfile'] not in (None, ""):
            self.openXYFile(args['xyfile'])
        else:
            for ass in self.assemblages:
                self.xAssemblage[ass]=0.0
                self.yAssemblage[ass]=0.0
            allp=self.all_pairs(self.assemblages)
            for pr in allp:
                name = pr[0]+"*"+pr[1]
                self.distanceBetweenAssemblages[name]=0

        ############################################################################################################
        logger.debug("Assume threshold is 1.0 unless its specified in arguments.")
        threshold=1.0
        if args['threshold'] is not None :
            threshold=float(args['threshold'])

        logger.debug("Going to create list of valid pairs for comparisons.")
        self.thresholdDetermination(threshold,args)

        ###########################################################################################################
        logger.debug("Now calculate the bootstrap comparisons based ")
        logger.debug("on specified confidence interval, if in the arguments.")

        if args['bootstrapCI'] is not None:
            if args['bootstrapSignificance'] not in (None, ""):
                confidenceInterval= args['bootstrapSignificance']
            else:
                confidenceInterval=0.95
            self.bootstrapCICalculation(args, 100, float(confidenceInterval))

        ###########################################################################################################
        ### setup the output files. Do this now so that if it fails, its not AFTER all the seriation stuff
        OUTFILE,OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,SUMGRAPH=self.setupOutput(args)

        ###########################################################################################################
        logger.debug("Now pre-calculating all the combinations between pairs of assemblages. ")
        logger.debug("This returns a graph with all pairs and the comparisons as weights.")
        pairGraph = self.preCalculateComparisons(args)

        frequencyArray=[]
        continuityArray=[]
        maxNodes=3
        notPartOfSeriationsList=[]

        if args['frequency'] not in (None,False,0):
            ###########################################################################################################
            logger.debug("Calculate all the valid triples.")
            triples = self.findAllValidTriples(args)
            ###########################################################################################################
            stepcount = 0
            currentMaxSeriationSize = 2
            newNetworks=[]
            solutionCount=len(triples)

            currentTotal = len(triples)
            solutions=[]
            all_solutions=[]
            all_solutions= all_solutions + triples  ## add the triples to the intial solution

            while currentMaxSeriationSize < self.maxSeriationSize:
                currentMaxSeriationSize += 1
                ### first time through copy the triples, else get the previous new ones.
                if currentMaxSeriationSize==3:  ## first time through. Just copy the triples to the working arrays
                    networks = triples
                    solutions = triples # clear the
                else:
                    i = 0
                    logger.debug("Currently have %d solutions at step %d", len(newNetworks),currentMaxSeriationSize)
                    if len(newNetworks)==0:
                        # there were no networks the previous times so nothing to do.
                        break
                    logger.debug("These solutions are ---  ")
                    for sol in newNetworks:
                        logger.debug("solution %d: %s", i, nx.shortest_path(sol, sol.graph["End1"] , sol.graph["End2"]))
                        i += 1
                    networks=[]
                    networks += newNetworks  # copy the array of previous new ones for this round
                    solutions.append(newNetworks ) # append the new list to the previous one
                    #print "Number of new Networks:", len(newNetworks)
                    newNetworks=[]         # clear the array of new solutions

                stepcount += 1
                logger.debug("_______________________________________________________________________________________")
                logger.debug("Step number:  %d", currentMaxSeriationSize)
                logger.debug("_______________________________________________________________________________________")

                if args['screen'] not in (None, ""):
                    self.scr.addstr(4,0,"Step number:                                    ")
                    msg = "Step number:   %d" % currentMaxSeriationSize
                    self.scr.addstr(4,0,msg)
                    self.scr.addstr(5,0,"Number of solutions from previous step:         ")
                    msg= "Number of solutions from previous step: %d" % len(networks)
                    self.scr.addstr(5,0,msg)
                    self.scr.refresh()

                logger.debug("Number of solutions from previous step: %d", len(networks))
                match = 0      ## set the current match to zero for this step (sees if there are any new solutions for step)
                ## look through the set of existing valid networks.
                validNewNetworks=[]

                ##pool.map(seriationCheck, networks)
                for nnetwork in networks:
                    logger.debug("-----------------------------------------------------------------------------------")
                    logger.debug("Network: %s", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
                    logger.debug("-----------------------------------------------------------------------------------")
                    ## find the ends
                    ## given the ends, find the valid set of assemblages that can be potentially added
                    ## this list is all assemblages meet the threshold requirements
                    validNewNetworks,currentMaxNodes = self.checkForValidAdditionsToNetwork(nnetwork, pairGraph,solutionCount,args)
                    if validNewNetworks is not False:
                        newNetworks += validNewNetworks
                        all_solutions += validNewNetworks
                        solutionCount += len(validNewNetworks)
                        logger.debug("Added %d new solutions. Solution count is now:  %d", len(validNewNetworks),solutionCount)
                        if currentMaxNodes > maxNodes:
                            maxNodes = currentMaxNodes
                        currentTotal = len(newNetworks)

                if args['screen'] not in (None, ""):
                    msg = "Current Max Nodes:  %d " % maxNodes
                    self.scr.addstr(6, 0, msg)
                    msg = "Total number of seriation solutions and sub-solutions: %d" % solutionCount
                    self.scr.addstr(7, 0, msg)
                    self.scr.addstr(8, 43, "                                           ")
                    msg = "Number of seriation solutions at this step: %d" % currentTotal
                    self.scr.addstr(8, 0, msg)
                    msg = "Memory used:        " + str(self.mem.memory())
                    self.scr.addstr(9, 0, msg)
                    self.scr.refresh()

                if len(newNetworks)>0:
                    end_solutions = newNetworks
                else:
                    end_solutions = networks

            logger.debug("Process complete at seriation size %d with %d solutions before filtering.",self.maxSeriationSize,len(end_solutions))

            ###########################################################################################################
            frequencyArray = self.filterSolutions(end_solutions,all_solutions,args)

            #filteredarray = all_solutions

            logger.debug("Process complete at seriation size %d with %d solutions after filtering.",self.maxSeriationSize,len(frequencyArray))

            #################################################### OUTPUT SECTION ####################################################
            self.output(frequencyArray,OUTFILE,OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,maxNodes,args)

            sumGraph=self.sumGraphsByWeight(frequencyArray,args)
            self.sumGraphOutput(sumGraph,SUMGRAPH,self.inputFile[0:-4]+"-mst-sumgraph.png",args)
            #self.createAtlasOfSolutions(frequencyArray,args)

            #################################################### MST SECTION ####################################################
            if args['mst'] not in (None,False,0):
                outputFile = self.outputDirectory + self.inputFile[0:-4]+".vna"
                # Need to have the shapefile flag and the XY file in order to create a valid shapefile.
                if args['shapefile'] is not None and args['xyfile'] is not None:
                    shapefile = 1
                else:
                    shapefile= None
                mst = MST.MST(outputFile,self.outputDirectory,shapefile)
                mst.createMST()
                #minimumSpanningTree(all_solutions,xAssemblage,yAssemblage,distanceBetweenAssemblages,assemblageSize,outputDirectory,inputFile)
            #################################################### MST SECTION ####################################################

            print "Assemblages not part of final solution:"
            nodeList=sumGraph.nodes()
            for a in self.assemblages:
                if a not in nodeList:
                    notPartOfSeriationsList.append(a)
                    print a

        if args['continuity'] not in (None,False,0):
            # experimental
            continuityArray=self.continunityMaximizationSeriation(args)
            #self.outputGraphArray(array,args)
            sGraph=self.sumGraphsByCount(continuityArray,args)
            self.graphOutput(sGraph,self.inputFile[0:-4]+"-minimum-sumgraph.png", args)
            self.MST(sGraph,self.inputFile[0:-4]+"-mst-of-min.png",args)

        if args['graphs'] not in (None,False,0):
            plt.show() # display

        ## say goodbye and clean up the screen stuff #########################
        self.finalGoodbye(maxNodes,len(frequencyArray),len(continuityArray), args)

        return frequencyArray, continuityArray, notPartOfSeriationsList

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Conduct an iterative deterministic seriation analysis')
    parser.add_argument('--debug', default=None, help='Sets the DEBUG flag for massive amounts of annoated output.')
    parser.add_argument('--bootstrapCI', default=None, help="Sets whether you want to use the bootstrap confidence intervals for the comparisons between assemblage type frequencies. Set's to on or off.")
    parser.add_argument('--bootstrapSignificance', default=0.95, type=float, help="The significance to which the confidence intervals are calculated. Default is 0.95.")
    parser.add_argument('--filtered',default=1,help="The script will complete by checking to see if smaller valid solutions are included in the larger sets. If not, they are added to the final set. Default is true. ")
    parser.add_argument('--largestonly',default=None, help="If set, the results will only include the results from the last and largest successful series of solutions. Smaller solutions will be excluded. Default is false.")
    parser.add_argument('--individualfileoutput',default=None,help="If true, a .VNA files will be created for every solution.")
    parser.add_argument('--excel',default=None, help="Not implemented.")
    parser.add_argument('--threshold',default=None,help="Sets the maximum difference between the frequencies of types that will be examine. This has the effect of keeping one from evaluating trivial solutions or solutions in which there is limited warrant for establishing continuity. Default is false.")
    parser.add_argument('--noscreen',default=None, help="If true, there will be no text output (i.e., runs silently). Default is false.")
    parser.add_argument('--xyfile',default=None,help="Enter the name of the XY file that contains the name of the assemblage and the X and Y coordinates for each.")
    parser.add_argument('--pairwisefile',default=None, help="If you have precalculated the bootstrap comparative p-value, enter the name of the file here and it will be used as the basis of the graphical output for showing significance of comparisons. Default is false.")
    parser.add_argument('--mst', default=None, help="If true, will produce a minimum spanning tree diagram from the set of final solutions.")
    parser.add_argument('--stats', default=None, help="(Not implemented). If true, a histogram of the solutions will be shown in terms of the #s of time pairs are included. Default is false.")
    parser.add_argument('--screen', default=True, help="Sets whether the output will be sent all to the screen or not. Default is false. When true, the screen output is all captured through curses." )
    parser.add_argument('--allsolutions', default=None,help="If set, all of the valid solutions are produced even if they are subsets of larger solutions.")
    parser.add_argument('--inputfile',help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
    parser.add_argument('--outputdirectory',default=None, help="If you want the output to go someplace other than the /output directory, specify that here.")
    parser.add_argument('--shapefile',default=None,help="Produces a shapefile as part of the output. You must have specified the XYfile as well.")
    parser.add_argument('--graphs',default=None,help="If true, the program will display the graphs that are created. If not, the graphs are just saved as .png files.")
    parser.add_argument('--frequency',default=None,help="Conduct a standard frequency seriation analysis. Default is None.")
    parser.add_argument('--continuity',default=None,help="Conduct a continuity seriation analysis. Default is None.")
    parser.add_argument('--graphroot',default=None,help="The root of the graph figures (i.e., name of assemblage you want to treat as one end in the graphs.")
    parser.add_argument('--continuityroot',default=None,help="If you have a outgroup or root of the graph, set that here.")
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    seriation = IDSS()

    frequencyResults,continuityResults,exceptionList=seriation.seriate(args)

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