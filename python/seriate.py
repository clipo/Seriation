__author__ = 'carllipo'

import csv
from datetime import datetime
import pprint
import argparse
import sys
import logging
import itertools
import math
import random
import curses
import numpy as np
import scipy as sp
from scipy import stats
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
from copy import copy, deepcopy

# start prettyprint (python Dumper)
pp = pprint.PrettyPrinter(indent=4)

### yields a<->b and b<->a
#def all_pairs(lst):
#    for p in itertools.permutations(lst):
#        i = iter(p)
#        yield zip(i,i)

def all_pairs(lst):
    return list((itertools.permutations(lst, 2)))

def all_tuples(lst):
    tuples=list(itertools.combinations(lst, 3))
    useable_tuples=[]
    for e in tuples:
        useable_tuples.append(e)
    return useable_tuples

def openFile(filename):
    assemblageValues={}
    assemblageSize={}
    assemblageFrequencies={}
    assemblages={}
    countOfAssemblages=0
    labels={}
    if screenFlag>0:
        msg1 = "Filename: %s " % filename
        scr.addstr(1,0,msg1)
        scr.refresh()

    ## Read in the data
    # the input of the classes -- note these are absolute counts, not frequencies
    # might want to change that...\
    if screenFlag:
        scr.addstr(1,40,"STEP: Read in data...")
        scr.refresh()
    try:
        logging.debug("trying to open: %s ", filename)
        file=open(filename,'r')
    except csv.Error as e:
        logging.error("Cannot open %s. Error: %s", filename, e)
        sys.exit('file %s does not open: %s') %( filename, e)

    reader = csv.reader(file, delimiter='\t', quotechar='|')

    values=[]
    for row in reader:
        row = map(str, row)
        label=row[0]
        labels[ label ] = label
        row.pop(0)
        #row=[map(int,x.split()) for x in f if x.strip()]
        row = map(float, row)
        numberOfClasses = len(row)
        #rowtotal = sum(row)
        freq=[]
        rowtotal=sum(row)
        for r in row:
            freq.append(float(float(r)/float(rowtotal)))
            values.append(float(r))
        assemblages[ label ] = freq
        assemblageFrequencies[ label ]  = freq
        assemblageValues[ label ] = values
        assemblageSize[ label ]= rowtotal
        countOfAssemblages +=1

    return len(assemblages), assemblages, assemblageFrequencies,assemblageValues,assemblageSize,numberOfClasses


def preCalculateComparisons(assemblages,bootstrapCI,typeFrequencyUpperCI,typeFrequencyLowerCI):
    logging.debug("Precalculating the comparisons between all pairs of assemblages...")
    pairs = all_pairs(assemblages)
    pairGraph = nx.Graph()
    for pair in pairs:
        pairGraph.add_node(pair[0])
        pairGraph.add_node(pair[1])
        columns=range(len(assemblages[pair[0]]))
        ass1 = assemblages[pair[0]]
        ass2 = assemblages[pair[1]]
        comparison=""
        for i in columns:
            val1 = ass1[i]
            val2 = ass2[i]
            logging.debug( "\t\tComparing Assemblage: %s  and    Assemblage: %s  ########",pair[0],pair[1])
            logging.debug( "\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,i,i,i,i,i,i)

            if bootstrapCI > 0:
                upperCI_test = typeFrequencyUpperCI[pair[0]][i]
                lowerCI_test  = typeFrequencyLowerCI[pair[0]][i]
                upperCI_end =  typeFrequencyUpperCI[pair[1]][i]
                lowerCI_end=  typeFrequencyLowerCI[pair[1]][i]

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
            logging.debug( "Type %d: - comparison is: %s ",i, comparison[i])

        logging.debug("Comparison for %s and %s is: %s ",pair[0],pair[1],comparison)
        pairGraph.add_edge(pair[0],pair[1],weight=comparison)
    return pairGraph

def openPairwiseFile( filename ):
    logging.debug("Opening pairwise file %", filename)
    try:
        pw = open(filename,'r' )
    except csv.Error as e:
        sys.exit('pairwise file %s does not open: %s') %( filename, e)

    reader = csv.reader(pw, delimiter='\t', quotechar='|')
    for row in reader:
        pair = row[0]+"#"+row[1]
        pairwise[pair]=row[2]
        pairwiseError[pair]=row[3]


def openXYFile( filename ):
    logging.debug("Opening pairwise file %", filename)
    ## open the xy file
    try:
        xyf= open( filename,'r')
    except csv.Error as e:
        sys.exit('file %s does not open: %s') %( filename, e)

    xAssemblage={}
    yAssemblage={}
    reader = csv.reader(xyf, delimiter='\t', quotechar='|')

    for row in reader:
        label = row[0]
        xyAssemblages.append(label)
        yAssemblage[ label ]= row[1]
        xAssemblage[ label ]= row[2]

    assemblagePairs = all_pairs(xyAssemblages)
    ## Go through all of the combinations
    for combo in assemblagePairs:
        pairname = combo[0]+"*"+combo[1]
        distance = math.sqrt( (xAssemblage[ combo[0] ] -xAssemblage[ combo[1]])^2 + (yAssemblage[combo[0]] -yAssemblage[combo[1]])^2)
        distanceBetweenAssemblages[ pairname ] = distance

    largestX = max(xAssemblage.iterkeys(), key=(lambda key: xAssemblage[key]))
    largestY= max(yAssemblage.iterkeys(), key=(lambda key: yAssemblage[key]))
    return largestX,largestY,distanceBetweenAssemblages,xAssemblage,yAssemblage


#############################################  THRESHOLD DETERMINATION ####################################
## first compare each assemblage and see if the threshold is exceeded.
## By this I mean the maximum difference between frequencies of any Type %ds greater than what is specified.
## When threshold = 0, all combinations are used. Later the "ends" of solutions are not evaluated if the
## difference between the last assemblage and any other free assemblage is > the threshold.
## This arbitrary setting is to keep from arbitrary solutions being stuck on that come from the "ends"
## of solutions.
##
## Precalculate all of the max differences between types in assembalge pairs.

def thresholdDetermination(threshold, assemblages):
    assemblageComparison={}
    validComparisonsArray={}
    ##  get all the combinations of 2
    pairs = all_pairs(assemblages)

    ## Go through all of the combinations
    for combo in pairs:
        logging.debug("comparing combination of %s and %s ", combo[0] , combo[1] )
        pairname1  = combo[0] + "*" + combo[1]
        pairname2  = combo[1] + "*" + combo[0]
        maxDifference = 0
        assemblage1 = assemblages[combo[0]]
        assemblage2 = assemblages[combo[1]]
        i=-1
        columns= len(assemblages[combo[0]])
        logging.debug("Number of columns: %d", columns)
        ## calculate the maximum differences between the pairs of assemblages (across all types)
        for i in (0, columns-1):
            logging.debug("i: %d ",i)
            ass1 = float(assemblage1[i])
            ass2 = float(assemblage2[i])
            diff = abs( ass1 - ass2 )
            logging.debug("assemblage1: %f assemblage2: %f diff: %f",ass1,ass2,diff)
            if diff > maxDifference :
                maxDifference = diff
        assemblageComparison[ pairname1 ] = maxDifference
        assemblageComparison[ pairname2 ] = maxDifference

    ############## pre calculate the valid pairs of assemblages and stuff into hash ############################
    for assemblage1 in assemblages:
        cAssemblages=[]
        for assemblage2 in assemblages:
            if not assemblage1 == assemblage2:
                testpair = assemblage1 + "*" + assemblage2
                logging.debug("Pairs: %s and %s", assemblage1,assemblage2)
                logging.debug("Comp value of pairs: %s:  %f and threshold is: %f",testpair, assemblageComparison[testpair],threshold)
                if assemblageComparison[ testpair ] <= threshold:
                    logging.debug("Appending %s to the list of valid comparisons for %s ", assemblage1, assemblage2)
                    cAssemblages.append( assemblage2 )

        validComparisonsArray[ assemblage1]  = cAssemblages
    return validComparisonsArray

def confidence_interval(data, confidence=0.95):
  ci2 = ( 1.0 - float(confidence))*.5
  d = 1.0*np.array(data)
  low_idx = int(ci2* d.size)
  high_idx = int((1-ci2)*d.size)
  d.sort()
  return d.mean(), d[low_idx], d[high_idx]

########################################### BOOTSTRAP CI SECTION ####################################
def bootstrapCICalculation(assemblages, assemblageSize, bootsize=100, confidenceInterval=0.95):

    typeFrequencyLowerCI = {}
    typeFrequencyUpperCI = {}
    typeFrequencyMeanCI = {}

    # now do ALL the pairwise assemblage comparisons
    # go to sleep and come back later.

    if screenFlag:
        scr.addstr(1,40, "STEP: Bootstrap CIs...        ")
        scr.refresh()

    ## for each assemblage
    logging.debug("Calculating bootstrap confidence intervals")
    # for each assemblage

    for currentLabel in sorted( assemblages.iterkeys()):
        assemblage =  assemblages[ currentLabel ]
        types=len(assemblages[ currentLabel ])
        currentAssemblageSize = assemblageSize[ currentLabel ]

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
            mean, lower, upper = confidence_interval(freq, confidence=float(confidenceInterval))
            if math.isnan(lower) is True:
                lower=0.0
            if math.isnan(upper) is True:
                upper=0.0
            if math.isnan(mean) is True:
                mean=0.0
            lowerCI.append( lower)
            upperCI.append( upper )
            meanCI.append( mean )

        typeFrequencyLowerCI[ currentLabel ] = lowerCI
        typeFrequencyUpperCI[ currentLabel ] = upperCI
        typeFrequencyMeanCI[ currentLabel ] = meanCI

    return typeFrequencyLowerCI, typeFrequencyUpperCI, typeFrequencyMeanCI

########################################### FIND ALL THE VALID TRIPLES  ####################################
########################################### #################################### ###########################
def findAllValidTriples(assemblages,pairGraph,validAssemblagesForComparisons,
                        bootstrapCI,typeFrequencyLowerCI, typeFrequencyUpperCI,typeFrequencyMeanCI):
    triples=[]
    error = 0
    numberOfTriplets = 0

    if screenFlag > 0:
        scr.addstr(1,40, "STEP: Find valid triples....      ")
        scr.refresh()
    permutations =all_tuples(assemblages)

    for permu in permutations:
        if screenFlag >0:
            c = scr.getch()
            if c == ord('q'):
                curses.endwin()
                curses.resetty()
                sys.exit("Quitting as requested.\n\r")
        logging.debug("Triple test: %s * %s * %s", permu[0],permu[1],permu[2])

        comparison12 = ""
        comparison23 = ""
        error = 0
        columns=len( assemblages[ permu[0] ])
        logging.debug("Columns: %d", columns)
        difscore=0
        difscore2=0
        comparison12=""
        comparison23=""
        for i in range(0,columns):
            ass1 = assemblages[ permu[0] ][i]
            ass2 = assemblages[ permu[1] ][i]
            ass3 = assemblages[ permu[2] ][i]
            logging.debug( "ass1: %f ass2: %f ass3: %f",ass1,ass2,ass3)

            if bootstrapCI >0 :
                low1 = typeFrequencyLowerCI[permu[0]][i]
                low2 = typeFrequencyLowerCI[permu[1]][i]
                low3 = typeFrequencyLowerCI[permu[2]][i]
                high1 = typeFrequencyUpperCI[permu[0]][i]
                high2 = typeFrequencyUpperCI[permu[1]][i]
                high3 = typeFrequencyUpperCI[permu[2]][i]
                mean1 = typeFrequencyMeanCI[permu[0]][i]
                mean2 = typeFrequencyMeanCI[permu[1]][i]
                mean3 = typeFrequencyMeanCI[permu[2]][i]

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

                    logging.debug("\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f" % ass1,ass2,ass3)
                    print "\n\rNo match to our possibility of combinations. ass1: %f ass2: %f  ass3: %f \n\r" , ass1,ass2,ass3
                    print "I must quit. Debugging required.\n\r"
                    sys.exit()

                logging.debug("Comparison12: %s Comparison23: %s", comparison12,comparison23)

        comparison = comparison12 + comparison23
        test = re.compile('DU').search(comparison)
        if test is not None:
            error +=1

        if error == 0:
            # uses NetworkX
            net = nx.Graph(name=numberOfTriplets, GraphID=numberOfTriplets, End1=permu[0], End2=permu[2], Middle=permu[1])
            net.add_node(permu[0], name=permu[0], site="end", end=1, connectedTo=permu[1] )
            net.add_node(permu[1], name=permu[1], site="middle", end=0, connectedTo="middle")
            net.add_node(permu[2], name=permu[2], site="end", end=1, connectedTo=permu[1])
            net.add_edge(permu[0], permu[1],weight=comparison12, GraphID=numberOfTriplets,end=1)
            net.add_edge(permu[2], permu[1],weight=comparison23, GraphID=numberOfTriplets,end=1)
            logging.debug("VALID TRIPLE SOLUTION: %s * %s * %s " , permu[0],permu[1], permu[2])
            logging.debug("VALID TRIPLE SOLUTION: %s  <--->   %s", comparison12, comparison23)
            logging.debug("VALID TRIPLE SOLUTION: %s ", net.adjacency_list())
            path = nx.shortest_path(net, source=permu[0], target=permu[2])
            logging.debug("VALID TRIPLE SOLUTION: Ends are: %s and %s",permu[0],permu[2])
            logging.debug("VALID TRIPLE SOLUTION: Shortest Path: %s ", path)
            triples.append( net )
            numberOfTriplets += 1
            logging.debug("Current number of triplets: %d", numberOfTriplets)

    return triples

def filter_list(full_list, excludes):
    s = set(excludes)
    return (x for x in full_list if x not in s)


def checkForValidAdditionsToNetwork(nnetwork, pairGraph, validAssemblagesForComparisons, assemblages,
                                    typeFrequencyLowerCI, typeFrequencyUpperCI, bootstrapCI,
                                    typeFrequencyMeanCI,
                                    solutionCount):
    logging.debug(" ######################Starting check for solution %s with %s nodes ######################################",nnetwork.graph['GraphID'],len(nnetwork))
    if screenFlag > 0:
        scr.addstr(1,40, "STEP: Testing for addition to seriation ....      ")
        scr.refresh()

    logging.debug("The end of assemblages of network %d are: %s and %s", nnetwork.graph['GraphID'], nnetwork.graph["End1"] , nnetwork.graph["End2"])
    logging.debug("Network:  %s", nnetwork.adjacency_list())
    logging.debug("Seriation %d to evaluate: Shortest Path: %s ", nnetwork.graph['GraphID'], nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
    array_of_new_nodes=[]
    maxnodes=len(nnetwork.nodes())

    for assEnd in ("End1","End2"):
        if assEnd=="End1":
            otherEnd="End2"
        else:
            otherEnd="End1"

        endAssemblage=nnetwork.graph[assEnd]
        logging.debug(">>>>>> Checking ends of seriation %d:  %s is %s", nnetwork.graph['GraphID'], assEnd,endAssemblage)
        list1 = validAssemblagesForComparisons[ endAssemblage ]
        list2 = nnetwork.nodes()
        logging.debug("List 1 (valid comparisons): %s", list1)
        logging.debug("List 2 (existing nodes): %s", list2)

        validAssemblages = list(filter_list(list1, list2))
        logging.debug("Valid assemblages: %s", validAssemblages)
        logging.debug("The list of valid comparative assemblages for %s is %s",endAssemblage,validAssemblages)

        ######################################################################################
        for testAssemblage in validAssemblages:
            logging.debug(" Now checking %s to see if we can be put next to %s",testAssemblage,endAssemblage)
            if screenFlag >0:
                msg = "Now checking %s against %s." % (testAssemblage, endAssemblage)
                scr.addstr(3,0, msg)
                scr.refresh()
                c = scr.getch()
                if c == ord('q'):
                    curses.endwin()
                    curses.resetty()
                    sys.exit("Quitting as requested.\n\r")

            ## now see if the test assemblages fits on the end.
            logging.debug("Checking assemblage %s to see if it fits on the end of the current solution.", testAssemblage )

            #### FIND INNER EDGE RELATIVE TO THE EXISTING END ASSEMBLAGE ##############
            #neighbors = nnetwork.neighbors(endAssemblage)

            logging.debug("Seriation %d with this %s: %s has this many nodes: %d", nnetwork.graph['GraphID'],assEnd,nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]), len(nnetwork.nodes()))
            logging.debug("End assemblage for this seriation: %s",nnetwork.graph[assEnd])
            logging.debug("Which end: %s", assEnd)
            innerNeighbor = None
            if assEnd=="End1":
                path = nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                innerNeighbor=path[1]
                logging.debug("End1: %s Neighbor: %s", endAssemblage, innerNeighbor)
            elif assEnd=="End2":
                path = nx.shortest_path(nnetwork, nnetwork.graph["End2"] , nnetwork.graph["End1"])
                innerNeighbor=path[1]
                logging.debug("End2: %s Neighbor: %s", endAssemblage, innerNeighbor)
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
            logging.debug( "\t\t\tThere should be just 1 neighbor to %s and that is: %s", endAssemblage, innerNeighbor )
            c = pairGraph.get_edge_data(innerNeighbor,endAssemblage )
            comparison=c['weight']
            logging.debug( "\t\t\tCompare current pair with previous comparison: %s", comparison)
            ##########################################################################
            comparisonMap =""
            oneToColumns=range(len(assemblages[testAssemblage]))
            logging.debug("Number of columns to check: %d", len(oneToColumns))

            error = 0  ## set the error check to 0
            for i in oneToColumns:
                logging.debug( "\t\t\tComparing Assemblage: %s  and    Assemblage: %s  ########",testAssemblage,endAssemblage)
                logging.debug( "\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,i,i,i,i,i,i)
                c=""
                p=nx.shortest_path(nnetwork, nnetwork.graph[assEnd] , nnetwork.graph[otherEnd])
                logging.debug( "Working on path: %s",p)
                newVal=assemblages[testAssemblage][i]
                logging.debug("Start comparison with %s",testAssemblage)
                previousAssemblage=testAssemblage
                for compareAssemblage in p:
                    oldVal=assemblages[compareAssemblage][i]
                    logging.debug("Compare %s with %s ", previousAssemblage,compareAssemblage)
                    logging.debug("Old value: %f  vs new value: %f",oldVal,newVal)
                    if bootstrapCI > 0:
                        upperCI_test = typeFrequencyUpperCI[previousAssemblage][i]
                        lowerCI_test = typeFrequencyLowerCI[previousAssemblage][i]
                        upperCI_end = typeFrequencyUpperCI[compareAssemblage][i]
                        lowerCI_end = typeFrequencyLowerCI[compareAssemblage][i]
                        mean_test = typeFrequencyMeanCI[previousAssemblage][i]
                        mean_end = typeFrequencyMeanCI[compareAssemblage][i]

                        if upperCI_test < lowerCI_end:
                            c += "D"
                        elif lowerCI_test > upperCI_end:
                            c += "U"
                        else:
                            c += "M"
                    else:
                        logging.debug("Outer value: %f Inner value: %f", oldVal, newVal)
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
                            logging.debug("Error. Quitting.")
                            sys.exit("got null value in comparison of value for type %d in the comparison of %s", i, compareAssemblage)
                        logging.debug("Comparison %s is now %s", c1,c)
                        newVal=oldVal

                    previousAssemblage=compareAssemblage

                test = re.compile('DU|DM*U').search(c)
                if test is not None:
                    logging.debug("Comparison is %s. Error!",c)
                    error +=1

            logging.debug("Checked out %s. Found %d total errors.", testAssemblage, error)
            if error == 0:
                logging.debug("Found no errors!  Going to add %s to end of existing network at %s", testAssemblage, endAssemblage)
                logging.debug( "Original network: %s ",nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
                logging.debug( "New comparison map is: %s ", comparisonMap)

                new_network = nnetwork.copy()
                new_network.graph["GraphID"]= solutionCount + 1
                logging.debug( "Here's the new network (before addition): %s", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
                logging.debug("From %s the ends of the seriation are %d (before): %s and %s",assEnd, nnetwork.graph['GraphID'],nnetwork.graph["End1"],nnetwork.graph["End2"] )
                path = nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"])
                logging.debug(" New network shortest path (before): %s ", path)

                ## mark this vertice as the new "END"
                new_network.add_node(testAssemblage, name=testAssemblage,end=1,site="end")
                ## mark the interior vertice as not "END
                new_network.add_node(endAssemblage, name=endAssemblage,site="middle", end=0)

                #### This adds the comparison to the new edge that has been added.
                new_network.add_edge( testAssemblage, endAssemblage, weight=comparisonMap, end=1, site="end", GraphID=solutionCount )
                logging.debug("Ends of the seriation %d (before): %s and %s ",new_network.graph['GraphID'], new_network.graph["End1"],new_network.graph["End2"] )

                logging.debug("Reassigning the new end %s from %s to %s", assEnd,new_network.graph[assEnd],testAssemblage )
                new_network.graph[assEnd]=testAssemblage
                logging.debug("From %s end of the seriation %s (after): %s and %s",assEnd, new_network.graph['GraphID'],new_network.graph["End1"],new_network.graph["End2"] )
                logging.debug("Here's the new network %s (with addition): %s", new_network.graph['GraphID'], new_network.adjacency_list())
                path = nx.shortest_path(new_network, new_network.graph["End1"] , new_network.graph["End2"])
                logging.debug("New network %d shortest path (after): %s ", new_network.graph['GraphID'], path)
                ## copy this solution to the new array of networks
                array_of_new_nodes.append(new_network)
                if len(new_network)> maxnodes:
                    maxnodes = len(new_network)
            logging.debug( "----------------#############-------End of check for %s ---------#############-----------------",testAssemblage)
        logging.debug("--------------------------------------Finished with %s-----------------------------------------------------",assEnd)
    logging.debug("------------------------------- Finished with Both Ends-----------------------------------------------------------------")

    if len(array_of_new_nodes)>0:
        return new_network,maxnodes,
    else:
        return False,0

def minimumSpanningTree(networks,xAssemblage,yAssemblage,distanceBetweenAssemblages,assemblageSize,outputDirectory,inputFile):
    try:
        from networkx import graphviz_layout
    except ImportError:
        raise ImportError("This function needs Graphviz and either PyGraphviz or Pydot")

    graphs=[]
    megaGraph = nx.Graph(name="MST")
    number=0
    graphCount=0
    for net in networks:
        graphCount += 1
        number = net.graph['GraphID']
        for nodey in net.nodes(data=True):
            xCoordinate = 0
            yCoordinate = 0
            name = nodey[0]
            xCoordinate = xAssemblage[name]
            yCoordinate = yAssemblage[name]
            megaGraph.add_node(name, name=name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                               size=assemblageSize[name])
            #graphs[graphCount].add_node(fromAssemblage, label=fromAssemblage, x=xCoordinate, y=yCoordinate,
            #                            name=fromAssemblage, size=assemblageSize[name])
            #graphs[graphCount].add_node(toAssemblage, label=toAssemblage, x=xCoordinate, y=yCoordinate,
            #                            name=toAssemblage, size=assemblageSize[name])

        count=0
        for e in net.edges_iter():   ### no need to go in order -- just look at all the other edges to see if there is an X
            d = net.get_edge_data(*e)
            fromAssemblage = e[0]
            toAssemblage = e[1]
            weight = d['weight']
            distance = distanceBetweenAssemblages[fromAssemblage + "*" + toAssemblage]
            #count = megaGraph.get_edge_data(fromAssemblage,toAssemblage,'weight'
            count += 1
            megaGraph.add_path([fromAssemblage, toAssemblage], weight=count,
                               distance=distance, color=number,
                               size=(assemblageSize[fromAssemblage], assemblageSize[toAssemblage]))

            #graphs[graphCount].add_path([fromAssemblage], [toAssemblage],
            #                            xy1=(xAssemblage[fromAssemblage], yAssemblage[fromAssemblage]),
            #                            xy2=(xAssemblage[toAssemblage], yAssemblage[toAssemblage]),
            #                            weight=weight,
            #                            meanDistance=distance,
            #                            size=(assemblageSize[fromAssemblage], assemblageSize[toAssemblage]))

    plt.rcParams['text.usetex'] = False
    plt.figure(0,figsize=(8,8))
    mst=nx.minimum_spanning_tree(megaGraph,weight='weight')
    os.environ["PATH"]=os.environ["PATH"]+":/usr/local/bin:"
    #print os.environ["PATH"]
    pos=nx.graphviz_layout(mst)
    #pos=nx.spring_layout(mst,iterations=500)
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
        'fontsize'   : 14}
    edgelist = list(mst) # make a list of the edges
    #print edgelist
    #nx.draw(mst)
    #plt.savefig("path.png")
    plt.axis('off')
    newfilename=outputDirectory+inputFile[:4]+"-mst.png"
    plt.savefig(newfilename,dpi=75)
    plt.figure(1,figsize=(30,20))
    # layout graphs with positions using graphviz neato

    return
    UU=nx.Graph()
    # do quick isomorphic-like check, not a true isomorphism checker
    nlist=[] # list of nonisomorphic graphs
    for G in graphs:
        # check against all nonisomorphic graphs so far
        if not nx.iso(G, nlist):
            nlist.append(G)
    UU=nx.union_all(graphs) # union the nonisomorphic graphs
    #UU=nx.disjoint_union_all(nlist) # union the nonisomorphic graphs
    pos=nx.spring_layout(UU,iterations=50)

    #pos=nx.graphviz_layout(UU)
    #pos=nx.graphviz_layout(UU,prog="twopi",root=0)
    ##labels=nx.draw_networkx_labels(UU,pos)
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
    atlasFile=outputDirectory+inputFile[:4]+"-atlas.png"
    plt.savefig(atlasFile,dpi=250)
    plt.show() # display

def finalGoodbye(start,maxNodes,currentTotal):
    if screenFlag >0:
        curses.endwin()
        curses.resetty()
        curses.nl()
        curses.echo()
    ## determine time elapsed
    #time.sleep(5)
    timeNow = time.time()
    timeElapsed = timeNow-start
    print "Seriation complete.\r"
    print "Maximum size of seriation: %d\r" % maxNodes
    print "Number of solutions at last step: %d\r" % currentTotal
    print "Time elapsed for calculation: %d seconds\r" % timeElapsed


#################################################### set up all the output files ####################################################
def setupOutput(filename, pairwiseFlag,mstFlag, outputDirectory,inputFile):
    outputFile = outputDirectory + inputFile[0:-4]+".vna"
    try:
        OUTFILE = open(outputFile, 'w')
    except csv.Error as e:
        msg = "Can't open file %s to write: %s" % outputFile, e
        sys.exit(msg)

    outpairsFile = outputDirectory +inputFile[0:-4]+"-pairs.vna"
    if pairwiseFlag is not None:
        try:
            OUTPAIRSFILE = open(outpairsFile, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outputFile, e
            sys.exit(msg)

    outmstFile=  outputDirectory + inputFile[0:-4] + "-mst.vna"
    outmst2File = outputDirectory + inputFile[0:-4] + "-mst-distance.vna"
    
    if mstFlag is not None:
        try:
            OUTMSTFILE = open(outmstFile, 'w')
            OUTMSTDISTANCEFILE = open(outmst2File, 'w')
        except csv.Error as e:
            msg = "Can't open file %s to write: %s" % outputFile, e
            sys.exit(msg)
            
    return OUTFILE,OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE

#################################################### sort by multiple keys ####################################################
def multikeysort(items, columns):
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

#################################################### OUTPUT SECTION ####################################################
def output(assemblages,assemblageSize,distanceBetweenAssemblages,xAssemblage,yAssemblage,largestX,largestY,filteredArray,
             OUTFILE, OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,mstFlag,largestonlyFlag,maxEdges,xyfileFlag,pairwiseFileFlag):
    if screenFlag:
        scr.addstr(13,1, "Now printing output file... ")
        scr.addstr(1,40,"STEP: Output files...         ")
        scr.refresh()
    OUTFILE.write( "*Node data\n")
    OUTFILE.write("ID AssemblageSize X Y Easting Northing\n")
    OUTPAIRSFILE.write("*Node data\n")
    OUTPAIRSFILE.write("ID AssemblageSize X Y Easting Northing\n")
    count = 0
    if screenFlag:
        scr.addstr(1,40,"STEP: Printing list of nodes....     ")
        scr.refresh()
    ## note this assumes the use of UTM coordinates (northing and easting)
    for l in assemblages:
        x = 0
        y = 0
        northing = 0
        easting = 0
        if xyfileFlag:
            x = xAssemblage[l]/1000000
            y = (largestY-yAssemblage[ l ])/100000
            easting = xAssemblage[ l ]
            northing = yAssemblage[ l ]

        msg = l + " "+ str(assemblageSize[ l ])+" "+ str(x)+" "+str(y)+" "+str(easting)+" "+str(northing)+"\n"
        OUTFILE.write(msg)
        OUTPAIRSFILE.write(msg)
        if mstFlag is not None:
            OUTMSTFILE.write(msg)
            OUTMSTDISTANCEFILE.write(msg)

    OUTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
    OUTPAIRSFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

    if mstFlag is not None:
        OUTMSTFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")
        OUTMSTDISTANCEFILE.write("*Node properties\nID AssemblageSize X Y Easting Northing\n")

    if screenFlag:
        scr.addstr(1,40,"STEP: Printing list of nodes attributes... ")
        scr.refresh()
    for l in assemblages:
        easting = 0
        northing = 0
        x = 0
        y = 0
        if xyfileFlag:
            x = xAssemblage[l]/1000000
            y = (largestY-yAssemblage[l])/100000
            easting = xAssemblage[l]
            northing = yAssemblage[l]
        msg = l +" "+ str(assemblageSize[ l])+" "+str(x)+" "+str(y)+" "+str(easting)+" "+str(northing)+"\n"
        OUTFILE.write(msg)
        OUTPAIRSFILE.write(msg)
        if mstFlag >0:
            OUTMSTFILE.write( msg )
            OUTMSTDISTANCEFILE.write(msg)

    ## This prints out counts of the edges as they appear in ALL of the solutions
    if screenFlag:
        scr.addstr(1,40,"STEP: Going through and counting pairs...     ")
        scr.refresh()
    OUTPAIRSFILE.write("*Tie data\nFrom To Edge Count\n")
    if mstFlag is not None:
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
    if screenFlag:
        scr.addstr(1,40,"STEP: Doing the pair output...                ")
        scr.refresh()

    sorted_pairs = sorted(edgeHash.iteritems(), key=operator.itemgetter(1))

    for key,value in sorted_pairs:
        ass1,ass2=key.split("*")
        msg = ass1+"\t"+ass2+"\t" + " 1 "+ str(value) +"\n"
        OUTPAIRSFILE.write(msg)

    OUTFILE.write("*Tie data\nFrom To Edge Weight Network End pValue pError meanSolutionDistance\n")
    if screenFlag:
        scr.addstr(1,40,"STEP: Eliminating duplicates...     ")
        scr.addstr(1,40,"STEP: Printing edges...     ")
        scr.refresh()

    uniqueArray = set(filteredArray)
    distanceHash={}
    seriationHash ={}
    ## only print unique ones...
    pairwise={}
    pairwiseError={}
    for network in uniqueArray:
        if screenFlag:
            scr.addstr(14,1, "Now on solution: ")
            scr.addstr(14,18,str(network.graph["GraphID"]) )
            #print "now on solution: ", network["GraphID"],"\n"

        if largestonlyFlag>0:
            if len(network.edges()) == maxEdges:
                groupDistance=0
                edgeCount = len(network.edges())
                meanDistance=0.0
                eCount=0
                if xyfileFlag is not None:
                    for e in network.edges_iter():
                        pairname= e[0]+"*"+e[1]
                        groupDistance += distanceBetweenAssemblages[ pairname ]
                        eCount += 1

                    meanDistance = groupDistance/eCount      ## use the average for the group for now
                    #print "\n\rMean distance for this group is: ", meanDistance, "\n\r"
                    network["meanDistance"]= meanDistance
                else:
                    meanDistance="0"
                    network["meanDistance"]= "0"

                for e in network.edges_iter():
                    pVal=0.0
                    pErr=0.0
                    d = e.get
                    if pairwiseFileFlag is not None:
                        pairname= e[0]+"#"+e[1]
                        pVal = pairwise[ pairname ]
                        pErr = pairwiseError[pairname]
                    text = e[0]+" "+e[1]+" 1 "+str(edgeCount)+ " "+network["GraphID"]+ " "\
                            +e[0]+ " End "+str(pVal)+" "+ str(pErr) + " " +str(meanDistance)+"\n"
                    OUTFILE.write(text)
                network['meanDistance'] = meanDistance
                distanceHash[ network["GraphID"] ]= meanDistance
                #seriationHash[ network["GraphID"] ]['meanDistance']= meanDistance
                #seriationHash[ network["GraphID"] ]['ID']=network["GraphID"]
                #seriationHash[ network["GraphID"] ]['size']=edgeCount
        else:  ## not just the largest, but ALL seriation solutions
            edgeCount = len(network.edges())
            groupDistance=0
            meanDistance=0.0
            eCount=0
            if xyfileFlag > 0:
                for e in network.edges_iter():
                  pairname= e[0]+"*"+e[1]
                  groupDistance += distanceBetweenAssemblages[ pairname ]
                  eCount += 1
                meanDistance = groupDistance/eCount         ##use the average distance as the metric
            else:
                meanDistance = "0"

            ## initialize ashes
            for e in network.edges_iter():
                pairname= e[0]+"#"+e[1]
                pairwise[ pairname ] = 0
                pairwiseError[ pairname ] = 0

            for e in network.edges_iter():
                pVal=0.0
                pErr=0.0
                #print "e0: ", e[0],"\n"
                #print "e1: ", e[1],"\n"
                if pairwiseFileFlag >0 :
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
            #seriationHash[network.graph["GraphID"] ]['meanDistance']= meanDistance
            #seriationHash[network.graph["GraphID"] ]['ID']=network["GraphID"]
            #seriationHash[network.graph["GraphID"] ]['size']=edgeCount

def main():
    mem=memory.Memory()
    parser = argparse.ArgumentParser(description='Conduct seriation analysis')
    parser.add_argument('--debug')
    parser.add_argument('--bootstrapCI')
    parser.add_argument('--bootstrapSignificance', type=float)
    parser.add_argument('--filtered')
    parser.add_argument('--largestonly')
    parser.add_argument('--individualfileoutput')
    parser.add_argument('--excel')
    parser.add_argument('--threshold')
    parser.add_argument('--noscreen')
    parser.add_argument('--xyfile')
    parser.add_argument('--pairwisefile')
    parser.add_argument('--mst')
    parser.add_argument('--stats')
    parser.add_argument('--nosum')
    parser.add_argument('--screen')
    parser.add_argument('--allsolutions')
    parser.add_argument('--memusage')
    parser.add_argument('--inputfile')
    parser.add_argument('--outputdirectory')
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    ##################################################################################################
    global scr
    global screenFlag
    screenFlag=0
    pairwiseFlag=0
    mstFlag=0
    if args['screen'] is not None:
        screenFlag=1
        ## Set up the screen display (default).
        ## the debug option should not use this since it gets messy
        try:
            # Initialize curses
            scr=curses.initscr()
            # Turn off echoing of keys, and enter cbreak mode,
            # where no buffering is performed on keyboard input
            curses.noecho()
            curses.cbreak()
            scr.nodelay(1)
            scr.addstr(0,0,"Iterative Seriation Program V.2.0", curses.A_BOLD)
            scr.addstr(20,35,"Hit <q> to quit.")
            scr.refresh()
        except:
            # In event of error, restore terminal to sane state.
            scr.keypad(0)
            curses.echo()
            curses.nocbreak()
            curses.endwin()
            curses.resetty()
            traceback.print_exc()           # Print the exception

    ##################################################################################################
    if args['debug'] is not None:
        ## Logging
        logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    else:
        logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

    # start the clock to track how long this run takes
    global start
    start = time.time()
    logging.debug("Start time:  %s ", start)
    logging.debug("Arguments: %s", args)
    bootstrapCI=0
    filename=args['inputfile']

    if filename is "":
        logging.error("You must enter a filename to continue.")
        print "You must enter a filename to continue."
        sys.exit("Quitting due to errors.")
    try:
        logging.debug("Going to try to open and load: %s", filename)
        maxSeriationSize, assemblages, assemblageFrequencies,assemblageValues,assemblageSize,numberOfClasses = openFile(filename)
    except IOError as e:
        logging.error("Cannot open %s. Error: %s", filename, e.strerror)

        print("Cannot open %s. Error. %s ", filename, e.strerror)
        if screenFlag is not None:
            curses.endwin()
            curses.resetty()
        sys.exit("Quitting due to errors.")
    inputFile=""
    try:
        inputparts =map(str,args['inputfile'].split("/"))
        inputFile= inputparts[len(inputparts)-1]
    except:
        sys.exit("There was a problem with parsing the input file. Check it and try again.")

    ############################################################################################################
    outputDirectory = ""
    if args['outputdirectory'] is not None:
        outputDirectory = args['outputdirectory']
    else:
        outputDirectory = "../output/"

    ############################################################################################################
    if args['largestonly'] is not None:
        largestonlyFlag=1

    if args['xyfile'] is not None:
        xyfileFlag=1

    ############################################################################################################
    logging.debug("Going to open pairwise file it is exists.")
    if args['pairwisefile'] is not None:
        openPairwiseFile(args['pairwisefile'])
        pairwiseFlag = 1

    ############################################################################################################
    logging.debug("Going to open XY file if it exists.")
    largestX=0
    largestY=0
    distanceBetweenAssemblages={}
    xAssemblage={}
    yAssemblage={}
    largestonlyFlag=0
    xyfileFlag=0
    if args['xyfile'] is not None:
        largestX,largestY,distanceBetweenAssemblages,xAssemblage,yAssemblage=openXYFile(args['xyfile'])
    else:
        for ass in assemblages:
            xAssemblage[ass]=0.0
            yAssemblage[ass]=0.0

        allp=all_pairs(assemblages)
        for pr in allp:
            name = pr[0]+"*"+pr[1]
            distanceBetweenAssemblages[name]=0

    ############################################################################################################
    logging.debug("Assume threshold is 1.0 unless its specified in arguments.")
    threshold=1.0
    if args['threshold']>0 :
        threshold=float(args['threshold'])

    logging.debug("Going to create list of valid pairs for comparisons.")
    validAssemblagesForComparisons={}
    validAssemblagesForComparisons = thresholdDetermination(threshold, assemblages)
    typeFrequencyLowerCI={}
    typeFrequencyUpperCI={}
    typeFrequencyMeanCI ={}
    ###########################################################################################################
    logging.debug("Now calculate the bootstrap comparisons based ")
    logging.debug("on specified confidence interval, if in the arguments.")
    if args['bootstrapCI'] is not None:
        bootstrapCI = 1
        confidenceInterval=0.95
        if args['bootstrapSignificance'] is not None:
            confidenceInterval= args['bootstrapSignificance']
        else:
            confidenceInterval=0.95
        typeFrequencyLowerCI, typeFrequencyUpperCI, typeFrequencyMeanCI= bootstrapCICalculation(assemblages,
                                                                            assemblageSize, 100,
                                                                            float(confidenceInterval))
    if args['mst'] is not None:
        mstFlag=1

    ###########################################################################################################
    ### setup the output files. Do this now so that if it fails, its not AFTER all the seriation stuff
    OUTFILE, OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE=setupOutput(filename,pairwiseFlag,mstFlag,outputDirectory,inputFile)

    ###########################################################################################################
    logging.debug("Now precalculating all the combinations between pairs of assemblages. ")
    logging.debug("This returns a graph with all pairs and the comparisons as weights.")
    pairGraph = preCalculateComparisons(assemblages, bootstrapCI, typeFrequencyUpperCI, typeFrequencyLowerCI)

    ###########################################################################################################
    logging.debug("Calculate all the valid triples.")
    triples = []
    triples = findAllValidTriples(assemblages, pairGraph, validAssemblagesForComparisons, bootstrapCI,
                                  typeFrequencyLowerCI, typeFrequencyUpperCI,typeFrequencyMeanCI)

    ###########################################################################################################
    stepcount = 0
    currentMaxSeriationSize = 2
    newNetworks=[]
    solutionCount=len(triples)
    maxNodes=3
    currentTotal = len(triples)
    solutions=[]
    networks=[]

    while currentMaxSeriationSize < maxSeriationSize:
        currentMaxSeriationSize += 1
        ### first time through copy the triples, else get the previous new ones.
        #print "current step: ", currentMaxSeriationSize
        if currentMaxSeriationSize==3:  ## first time through. Just copy the triples to the working arrays
            networks = triples
            solutions = triples # clear the
        else:
            i = 0
            logging.debug("Currently have %d solutions at step %d", len(newNetworks),currentMaxSeriationSize)

            if len(newNetworks)==0:
                # there were no networks the previous times so nothing to do.
                break

            logging.debug("These solutions are ---  ")
            for sol in newNetworks:
                logging.debug("solution %d: %s", i, nx.shortest_path(sol, sol.graph["End1"] , sol.graph["End2"]))
                i += 1
            networks=[]
            networks += newNetworks  # copy the array of previous new ones for this round
            solutions.append(newNetworks ) # append the new list to the previous one
            #print "Number of new Networks:", len(newNetworks)
            newNetworks=[]         # clear the array of new solutions


        #print "Number of networks:", len(networks)
        #print "Number of solutions:", len(solutions)

        stepcount += 1
        logging.debug("_______________________________________________________________________________________")
        logging.debug("Step number:  %d", currentMaxSeriationSize)
        logging.debug("_______________________________________________________________________________________")

        if screenFlag>0:
            scr.addstr(4,0,"Step number:                                    ")
            msg = "Step number:   %d" % currentMaxSeriationSize
            scr.addstr(4,0,msg)
            scr.addstr(5,0,"Number of solutions from previous step:         ")
            msg= "Number of solutions from previous step: %d" % len(networks)
            scr.addstr(5,0,msg)
            scr.refresh()

        logging.debug("Number of solutions from previous step: %d", len(networks))
        match = 0      ## set the current match to zero for this step (sees if there are any new solutions for step)
        ## look through the set of existing valid networks.
        validNewNetworks=[]
        for nnetwork in networks:
            logging.debug("-----------------------------------------------------------------------------------")
            logging.debug("Network: %s", nx.shortest_path(nnetwork, nnetwork.graph["End1"] , nnetwork.graph["End2"]))
            logging.debug("-----------------------------------------------------------------------------------")
            ## find the ends
            ## given the ends, find the valid set of assemblages that can be potentially added
            ## this list is all assemblages meet the threshold requirements
            validNewNetworks,currentMaxNodes = checkForValidAdditionsToNetwork(nnetwork, pairGraph, validAssemblagesForComparisons,
                                                              assemblages, typeFrequencyLowerCI, typeFrequencyUpperCI,
                                                              bootstrapCI, typeFrequencyMeanCI,solutionCount)
            if  validNewNetworks is not False:
                newNetworks.append(validNewNetworks)
                solutionCount += len(validNewNetworks)
                logging.debug("Added %d new solutions. Solution count is now:  %d", len(validNewNetworks),solutionCount)
                if currentMaxNodes > maxNodes:
                    maxNodes = currentMaxNodes
                currentTotal = len(newNetworks)

        if screenFlag > 0:
            msg = "Current Max Nodes:  %d " % maxNodes
            scr.addstr(6, 0, msg)
            msg = "Total number of seriation solutions and sub-solutions: %d" % solutionCount
            scr.addstr(7, 0, msg)
            scr.addstr(8, 43, "                                           ")
            msg = "Number of seriation solutions at this step: %d" % currentTotal
            scr.addstr(8, 0, msg)
            msg = "Memory used:        " + str(mem.memory())
            scr.addstr(9, 0, msg)
            scr.refresh()

    if len(newNetworks)>0:
        end_solutions = newNetworks
    else:
        end_solutions = networks

    logging.debug("Process complete at seriation size %d with %d solutions.",maxSeriationSize,len(end_solutions))

    ###########################################################################################################
    #if len(networks):
    #    print "\n\r\n\r\n\r\n\r\n\rNo solutions Found!!\n\r"
    #    finalGoodbye(start,maxNodes,currentTotal)

    ###########################################################################################################
    if args['mst'] is not None:
        minimumSpanningTree(end_solutions,xAssemblage,yAssemblage,distanceBetweenAssemblages,assemblageSize,outputDirectory,inputFile)

    if args['largestonly'] is not None:
        largestonlyFlag = 1

    ################################################# FILTERING  ####################################
    # now do some weeding. Basically start with the first network that is the largest, and work backwards. Ignore any
    # network that is already represented in the smaller ones since these are trivial (e.g., A->B->C->D already covers
    # A->C->D.) That should then leave all the unique maximum solutions (or so it seems)
    ################################################# FILTERING  ####################################
    ## first need to sort the networks by size
    filteredarray = []
    if args['filtered'] is not None:  ## only get the largest set that includes ALL
        if screenFlag:
            scr.addstr(1,40,"STEP: Filter to get uniques... ")
        logging.debug("---Filtering solutions so we only end up with the unique ones.")
        logging.debug("---Start with % solutions.", len(solutions))
        for i in range(0,len(solutions),-1):
            exists=0
            for tnetwork in filteredarray:
                fnetworkArray = solutions[i].nodes()
                tnetworkArray = tnetwork.nodes()
                minus = fnetworkArray - tnetworkArray
                if len(minus)== len(fnetworkArray)-len(tnetworkArray):
                    exists += 1
            if exists==0:
                ##print "pushing fnetwork to list\n\r"
                filteredarray.append(solutions[i])

        logging.debug("End with %d solutions.", len(filteredarray))
        filterCount= len(filteredarray)
        scr.addstr(11,1,"End with filterCount solutions.")
    elif args['allsolutions'] is not None:
        filteredarray = solutions  ## all possible networks
    else:
        filteredarray = end_solutions ## just the largest ones (from the last round)


    #################################################### OUTPUT SECTION ####################################################
    output(assemblages,assemblageSize,distanceBetweenAssemblages,xAssemblage,yAssemblage,largestX,largestY,filteredarray,
             OUTFILE, OUTPAIRSFILE,OUTMSTFILE,OUTMSTDISTANCEFILE,mstFlag,largestonlyFlag,maxNodes,xyfileFlag,pairwiseFlag)

    ## say goodbye
    finalGoodbye(start,maxNodes,currentTotal)

if __name__ == "__main__":
    main()


