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
import networkx as nx
import traceback
import matplotlib.pyplot as pltc

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
        label=row[0]
        labels[ label ] = label
        row.pop(0)
        numberOfClasses = len(row)
        row = map(float, row)
        rowtotal = sum(row)
        freq=[]
        for r in row:
            freq.append(float(float(r)/float(rowtotal)))
            values.append(float(r))
        assemblages[ label ] = freq
        assemblageFrequencies[ label ]  = freq
        assemblageValues[ label ] = values
        assemblageSize[ label ]= rowtotal
        countOfAssemblages +=1

    return countOfAssemblages, assemblages, assemblageFrequencies,assemblageValues,assemblageSize,numberOfClasses

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
                     ##  COMBINATIONS of VALUES
                       #   dif	comparison	result	comparisonMap
                       #   1	U	      okay	U
                       #   1	M	      okay	U
                       #   1	X	      bad	--
                       #   1	D	      bad	--
                       #   0	U	      okay	U
                       #   0	M	      okay	M
                       #   0	D	      okay	D
                       #   0	X	      okay	X
                       #   -1	U	      bad	--
                       #   -1	M	      okay	M
                       #   -1	D	      okay	D
                       #   -1	X	      okay	D

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
                dif1 = val1 - val2
                if dif1 < 0:
                    comparison +=  "D"
                if dif1 > 0:
                    comparison +=  "M"
                if dif1 == 0:
                    comparison +=  "U"

            logging.debug( "Type %d: - comparison is: %s  a score of: %d",i, comparison[i], dif1)

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
    return largestX,largestY,distanceBetweenAssemblages


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
        pairname  = combo[0] + "*" + combo[1]

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

        assemblageComparison[ pairname ] = maxDifference

    ############## pre calculate the valid pairs of assemblages and stuff into hash ############################
    for assemblage1 in assemblages:
        cAssemblages=[]
        for assemblage2 in assemblages:
            if not assemblage1 == assemblage2:
                testpair = assemblage1 + "*" + assemblage2
                logging.debug("Pairs: %s and %s", assemblage1,assemblage2)
                logging.debug("Comp value:  %f and threshold is: %f",assemblageComparison[testpair],threshold)
                if assemblageComparison[ testpair ] <= threshold:
                    logging.debug("Appending %s to the list of valid comparisons for %s ", assemblage1, assemblage2)
                    cAssemblages.append( assemblage2 )

        validComparisonsArray[ assemblage1]  = cAssemblages
    return validComparisonsArray

def confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

########################################### BOOTSTRAP CI SECTION ####################################
def bootstrapCICalculation(assemblages, assemblageSize, bootsize=1000, confidenceInterval=0.95):
    random.seed(start)
    perror ={}
    pvalue={}
    results = 0
    ptr1=0
    classes = 0
    typeFrequencyLowerCI = {}
    typeFrequencyUpperCI = {}

    # now do ALL the pairwise assemblage comparisons
    # go to sleep and come back later.

    if screenFlag:
        scr.addstr(1,40, "STEP: Bootstrap CIs...        ")
        scr.refresh()
    countup=0
    ## for each assemblage
    logging.debug("Calculating bootstrap confidence intervals")
    # for each assemblage

    for currentLabel in sorted( assemblages.iterkeys()):
        #label = assemblages[countup]
        a =  assemblages[ currentLabel ]
        columns=len(assemblages[ currentLabel ])
        currentAssemblageSize = assemblageSize[ currentLabel ]

        ## create an array of arrays - one for each type
        arrayOfStats = []
        for c in assemblages:
            array=[]
            arrayOfStats.append([])

        ## size of bootstrapping (how many assemblages to create)
        loop = bootsize

        while loop:
            assemsize = currentAssemblageSize
            cumulate=[]
            classes = columns
            index = 0
            total = 0.0
            count   = 0

            ## now count through the classes and set up the frequencies
            for count in range (0,classes):
                cumulate[index] = a[0][count]  ## this is the index of the frequency for this class
                total += a[0][count]            ## should ultimate equal 100
                index += 1                      ## index should be total # of types at end

            new_assemblage=[]
            while assemsize>0:
                rand = random()              ## random number from 0-1
                classVar = 0
                while (classVar < index) and (rand > cumulate[classVar] ):
                    rand -= cumulate[classVar]
                new_assemblage.append(classVar)
                assemsize -= 1

            ## this should result in a new assemblage of the same size
            ahat=[]
            aholder={}
            bholder={}
            aholder = {}

            ## initialize arrauy
            indexN = 0
            for indexN in range(0,classes):
                ahat[indexN] = 0

            ## count up the classes
            for assem in new_assemblage:
                aholder[assem] +=1

            classCount=0
            for stat in arrayOfStats.iterkeys():
                results = aholder[classCount]
                arrayOfStats.append(results / currentAssemblageSize)
                classCount += 1
            loop -= 1

        lowerCI=[]
        upperCI=[]
        for stat in arrayOfStats:
            upper=0
            lower=0
            mean=0
            mean, upper, lower = confidence_interval(stat, confidence=confidenceInterval)
            lowerCI.append( lower)
            upperCI.append( upper )

        typeFrequencyLowerCI[ currentLabel ] = lowerCI
        typeFrequencyUpperCI[ currentLabel ] = upperCI
        results = 0
        countup += 1

    return typeFrequencyLowerCI, typeFrequencyUpperCI


########################################### FIND ALL THE VALID TRIPLES  ####################################
########################################### #################################### ###########################
def findAllValidTriples(assemblages,pairGraph,validAssemblagesForComparisons,bootstrapCI,typeFrequencyLowerCI, typeFrequencyUpperCI):
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
            # first compare assemblages 1 and 2
            if bootstrapCI:
                upperCI_1 = typeFrequencyUpperCI[ assemblages[ permu[0] ] ][i]
                lowerCI_1 = typeFrequencyUpperCI[ assemblages[ permu[0] ] ][i]
                upperCI_2 = typeFrequencyUpperCI[ assemblages[ permu[1] ] ][i]
                lowerCI_2 = typeFrequencyUpperCI[ assemblages[ permu[1] ] ][i]
                dif1 = ass1 - ass2
                if upperCI_1 < lowerCI_2:
                    difscore = -1
                elif lowerCI_1 > upperCI_2:
                    difscore = 1
                else:
                    difscore = 0
            else:    #if the bootstrapCI is not being used
                dif1 = ass1 - ass2
                if dif1 < 0:
                    difscore = -1
                if dif1 > 0:
                    difscore = 1
                if dif1 == 0:
                    difscore = 0
            logging.debug("Difscore between ass1 and ass2:  %d", difscore)
            # now compare assemblages 2 and 3
            if bootstrapCI:   # boostrap confidence intervals
                upperCI_2 = typeFrequencyUpperCI[ permu[1] ][i]
                lowerCI_2 = typeFrequencyUpperCI[ permu[1] ][i]
                upperCI_3 = typeFrequencyUpperCI[ permu[2] ][i]
                lowerCI_3 = typeFrequencyUpperCI[ permu[2] ][i]
                if upperCI_2 < lowerCI_3:
                    difscore = -1
                elif lowerCI_2 > upperCI_3:
                    difscore = 1
                else:
                    difscore = 0
            else:          ## if the bootstrapCI is not being used
                dif2 = ass3 - ass2
                if dif2 < 0:
                    difscore2 = -1
                if dif2 > 0:
                    difscore2 = 1
                if dif2 == 0:
                    difscore2 = 0
            logging.debug("Difscore2 (between ass2 and ass3:  %d", difscore2)
            if difscore == 1 and difscore2 == 1:     # F1 > F2 < F3  criteria not met
                error += 1
                continue
            elif difscore == 1  and difscore2 == -1:  #  F1 > F2 > F3 OK
                comparison12 += "U"
                comparison23 += "D"
            elif   difscore == -1  and  difscore2 == -1: #   F1 < F2 >F3 OK
                comparison12 += "X"
                comparison23 += "X"
            elif difscore == -1  and difscore2 == 1:   # F1 < F2 < F3
                comparison12 += "D"
                comparison23 += "U"
            elif difscore == 0  and  difscore2 == 1  :
                comparison12 += "M"
                comparison23 += "U"
            elif difscore == 1 and  difscore2 == 0  :
                comparison12 += "U"
                comparison23 += "M"
            elif   difscore == 0  and  difscore2 == -1  :
                comparison12 += "M"
                comparison23 += "D"
            elif   difscore == -1  and  difscore2 == 0  :
                comparison12 += "D"
                comparison23 += "M"
            elif   difscore == 0  and  difscore2 == 0  :
                comparison12 += "M"
                comparison23 += "M"
            else:
                error += 1
                print "\n\rNo match to our possibility of combinations. Difscore 1: %d Difscore 2: %d \n\r" % difscore,difscore2
                print "I must quit. Debugging required.\n\r"
                sys.exit()

            logging.debug("Comparison12: %s Comparison23: %s", comparison12,comparison23)

        if error == 0:
            # uses NetworkX
            net = nx.Graph(name=numberOfTriplets, GraphID=numberOfTriplets,End1=permu[0],End2=permu[2])

            net.add_node(permu[0], name=permu[0], site="end", end=1, connectedTo=permu[1] )
            net.add_node(permu[1], name=permu[1], site="middle", end=0, connectedTo="middle")
            net.add_node(permu[2], name=permu[2], site="end", end=1, connectedTo=permu[1])

            net.add_path([permu[0], permu[1]],weight=comparison12, GraphID=numberOfTriplets,end=1)
            net.add_path([permu[2], permu[1]],weight=comparison23, GraphID=numberOfTriplets,end=1)

            logging.debug("VALID SOLUTION: %s * %s * %s " , permu[0],permu[1], permu[2])
            logging.debug("VALID SOLUTION: %s  <--->   %s", comparison12, comparison23)
            logging.debug("VALID SOLUTION: %s ", net.adjacency_list())
            triples.append( net )
            numberOfTriplets += 1
            logging.debug("Current number of triplets: %d", numberOfTriplets)
        error = 0

    return triples

def filter_list(full_list, excludes):
    s = set(excludes)
    return (x for x in full_list if x not in s)

def checkForValidAdditionsToNetwork(nnetwork,pairGraph,validAssemblagesForComparisons,assemblages,typeFrequencyLowerCI, typeFrequencyUpperCI, bootstrapCI,solutionCount):

    logging.debug(" ######################Starting Seriation ######################################")
    if screenFlag > 0:
        scr.addstr(1,40, "STEP: Adding to seriations ....      ")
        scr.refresh()

    whichEnd = 0
    logging.debug("The end of assemblages of network are: %s and %s", nnetwork.graph["End1"] , nnetwork.graph["End2"])
    logging.debug("Network:  %s", nnetwork.adjacency_list())
    new_network=nx.Graph()
    for endAssemblage in (nnetwork.graph["End1"],nnetwork.graph['End2']):
        whichEnd += 1 ## either a 1 or a 2
        list1 = validAssemblagesForComparisons[ endAssemblage ]
        list2 = nnetwork.nodes()
        validAssemblages = list(filter_list(list1, list2))
        for testAssemblage in validAssemblages:
            if screenFlag >0:
                msg = "Checking %s against %s." % (testAssemblage, endAssemblage)
                scr.addstr(3,0, msg)
                scr.refresh()

                c = scr.getch()
                if c == ord('q'):
                    curses.endwin()
                    sys.exit("Quitting as requested.\n\r")
            # We dont want to include the assemblage twice or anything already in network
            if testAssemblage in nnetwork.nodes():
                continue

            ## now see if the test assemblages fits on the end.
            logging.debug("\t\tChecking assemblage %s to see if it fits on the end of the current solution.", testAssemblage )

            #### FIND INNER EDGE RELATIVE TO THE EXISTING END ASSEMBLAGE ##############
            neighbors = nnetwork.neighbors(endAssemblage)
            if len(neighbors) > 1 or len(neighbors)==0:
                print "\r\n\r\n\r\nThere are too many or too few neighbors (should only be 1!). Error!\n\r"
                print "\r\nWe are testing endAssemblage: %s "% endAssemblage
                print "\r\n with neighbors:", neighbors
                print "For this network:  ", nnetwork.adjacency_list()
                #print nx.write_adjlist(nnetwork,sys.stdout) # write adjacency list to screen
                sys.exit("Quitting due to errors.")

            logging.debug( "\t\t\t The number of neighbors at endAssemblage is %d (should be just one).", len(neighbors))
            logging.debug( "\t\t\tThere should be just 1 neighbor to endAssemblage and that is: %s", neighbors[0])
            c = pairGraph.get_edge_data(neighbors[0],endAssemblage )
            comparison=c['weight']
            n = pairGraph.get_edge_data(endAssemblage, testAssemblage )
            newComparison=n['weight']

            ##########################################################################
            comparisonMap =""
            oneToColumns=range(len(assemblages[testAssemblage]))
            logging.debug("One to Columns: %s", oneToColumns)
            newassemblage=assemblages[testAssemblage]
            oldassemblage=assemblages[endAssemblage]
            error = 0  ## set the error check to 0
            for i in oneToColumns:
                logging.debug( "\t\tComparing Assemblage: %s  and    Assemblage: %s  ########",testAssemblage,endAssemblage)
                logging.debug( "\t\t\t\tType %d- Type %d - Type %d - Type %d - Type %d - Type %d - Type %d  ########", i,i,i,i,i,i,i)
                logging.debug( "\t\t\t\tType %d:  testAssemblage 1: %d  endAssemblage 2: %d ", i, newassemblage[i],oldassemblage[i])
                         ##  COMBINATIONS of VALUES
                           #   dif      comparison      result  comparisonMap
                           #   1        U             okay      U
                           #   1        M             okay      U
                           #   1        X             bad       --
                           #   1        D             bad       --
                           #   0        U             okay      U
                           #   0        M             okay      M
                           #   0        D             okay      D
                           #   0        X             okay      X
                           #   -1       U             bad       --
                           #   -1       M             okay      M
                           #   -1       D             okay      D
                           #   -1       X             okay      D

                if bootstrapCI > 0:
                    upperCI_test = typeFrequencyUpperCI[testAssemblage][i]
                    lowerCI_test  = typeFrequencyLowerCI[testAssemblage][i]
                    upperCI_end =  typeFrequencyUpperCI[endAssemblage][i]
                    lowerCI_end=  typeFrequencyLowerCI[endAssemblage][i]

                    if upperCI_test < lowerCI_end:
                        difscore = -1
                    elif lowerCI_test > upperCI_end:
                        difscore = 1
                    else:
                        difscore = 0

                else:
                    dif1 = newassemblage[i] - oldassemblage[i]
                    if dif1 < 0:
                        difscore = -1
                    if dif1 > 0:
                        difscore = 1
                    if dif1 == 0:
                        difscore = 0

                logging.debug( "\t\t\t\t#### Type %d: - comparison is: %s  a score of: %d",i, comparison[i], difscore)
                #################################################################################       #### 1 U  #############
                if difscore == 1  and comparison[i] is "U":
                    comparisonMap += "U"
                    logging.debug(  "\t\t\t\t#### Type %d: Got a difscore of 1 and a comparison of a U. This works.",i)
                    logging.debug( " \t\t\t\tAdding %s to vertices %s", testAssemblage,endAssemblage)

                #################################################################################           ### 1 M   #############
                elif difscore == 1 and comparison[i] is "M":
                    # this is okay - its a match and the new value is greater. New value shoudl be U
                    # need to find what was happening on the previous comparison to know whether this needs
                    # to be up or down.
                    logging.debug( "\t\t\t\t#### Type %d: Got a difscore of 1 and a comparison of a M. This could be okay.", i)
                    logging.debug( "\t\t\t\tType %d:   Matching case A (1, M)", i)
                    logging.debug( "\t\t\t\tThis will only work if there no Xs anywhere previously OR if the opposite end doesnt ALSO go up!")
                    logging.debug( "\t\t\t\t\t %s", nnetwork.adjacency_list())
                    xerror =0
                    ccount=0
                    for e in nnetwork.edges_iter(): ### no need to go in order -- jsut look at all the other edges to see if there is an X
                        ccount+=1
                        d = nnetwork.get_edge_data(*e)
                        logging.debug( "\t\t\t\t\ Now on %s => weight: %s",e,d['weight'])
                        newComparison = d['weight']
                        if newComparison[i] is None:
                            print "Comparison is empty. Error! Stopping.\n\r\n\r"
                            sys.exit("Quitting due to errors.")

                        logging.debug( "\t\t\t\t#### Type %d: Here is what we get for comparison # %d ",i, ccount)  ## note that i is the current type
                        logging.debug( " \t\t\t\t\t inwardEdge - outwardEdge: %s ->  %s", comparison[i],newComparison[i])
                        if newComparison[i] is "X" or newComparison[i] is "U":
                            xerror += 1  ### BLARGH a previous X or an UP ! This will not be tolerated!
                            logging.debug( "\t\t\t\t\t Since I got %s my potential new value is still X.",newComparison[i])
                            logging.debug( "\t\t\t\t\t Now going to get the continue pair of assemblages to examine in the chain")
                            ccount +=1
                    if xerror > 0:
                        error +=1
                        break
                    else:
                        comparisonMap += "U"
                        logging.debug( "\t\t\t\t ####Type %d: For this type, OK to add %s to vertices %s ",i,testAssemblage,endAssemblage)
                        logging.debug( "\t\t\t\t\t No previous X values anywhere. ")
                        logging.debug( "\t\t\t\t Type %d: Adding an U to the comparisons.", i)
                        logging.debug( "\t\t\t\t\t Comparison map is now comparisonMap")
                #################################################################################    ## 1 D   #############
                elif difscore == 1 and comparison[i] is  "D" :
                    ## error the new value is greater but should be less. Error!
                    error += 1
                    continue
                    logging.debug( "\t\t\t\t####Type %d: Value 1:  %d value 2: %d", i,newassemblage[i],oldassemblage[i])
                    logging.debug( "\t\t\t\tType %d: Comparison is: %s a score of: %d  ", i, comparison[i], difscore)
                    logging.debug( "\t\t\t\tType %d: Rejecting %s from %s", testAssemblage,endAssemblage)
                    logging.debug( "\t\t\t\t\t because value is 1 and comparison is D.")
                #################################################################################     # -1 U #############
                elif difscore == -1  and comparison[i] is  "U":
                    ## new score is less and the Comparison is up .. Error!
                    ## note -- need to check to see if there is a previous change in direction because
                    ## its possible that the this is a mode change.
                    ## do this by logging all modes in the original triplet -- keep track of modes per type
                    ## first check to see if there is already and X in this column somewhere else.
                    xerror= 0
                    logging.debug( "\t\t\t\t####Type %d:  Case B (-1, U). Potentially can add %s and vert %s",i,testAssemblage,endAssemblage)
                    logging.debug( "\t\t\t\tType %d:  But need to check the rest of the chain for X's (can't already be an X).",i)
                    ccount=0
                    for e in nnetwork.edges_iter():   ### no need to go in order -- just look at all the other edges to see if there is an X
                        ccount +=1
                        d =nnetwork.get_edge_data(*e)
                        newComparison = d['weight']
                        if newComparison[i] is None:
                            print "Comparison is empty. Error! Stopping.\n\r\n\r"
                            sys.exit("Quitting due to errors.")
                        logging.debug( "\t\t\t\t####Type %d: Here is what we get for comparison # %s ", i, ccount) ## note that i is the current type
                        logging.debug( " \t\t\t\t\t inwardEdge - outwardEdge: %s -> %s ", comparison[i],newComparison[i])
                        if newComparison[i] is  "X":
                            xerror += 1  ### BLARGH a previous X! This will not be tolerated!

                    if xerror > 0:
                        error += 1
                        logging.debug( "\t\t\t\tType %d: Rejecting %s from %s) because there was an X ", i, testAssemblage, endAssemblage)
                        logging.debug( "\t\t\t\t\t  This would make it multimodal - so error.")
                        break
                    else:
                        comparisonMap += "X"   ## this is an X unless there is an error....
                        logging.debug( "\t\t\t\t#### Type %d:Definitely OK to add %s to vertices %s because score", i,testAssemblage,endAssemblage)
                        logging.debug( "\t\t\t\t\tis -1 and the comparison is U but no other Xs in the previous linkages.")
                        logging.debug( "\t\t\t\tType %d: Adding an X to the comparisons for type %d. ",i,i)
                        logging.debug( "\t\t\t\t\tComparison map is now %s", comparisonMap)
                         #end if if check error (xerror)
                #################################################################################  ## -1   D #############
                elif difscore == -1 and comparison[i] is  "D":
                    ## new score is less and the comparison is down .. Yes!
                    comparisonMap += "D"
                    logging.debug( "\t\t\t\t#### Type %d: Adding a D to the comparisons for type %d. Comparison map is now: %s", i, i, comparisonMap)

                #################################################################################  ## ## -1 M #############
                elif difscore == -1 and  comparison[i] is  "M":
                    # new score is less but comparison is Match. Okay
                    #vHere    = endAssemblage
                    xerror   = 0  ## count for errors
                    logging.debug( "\t\t\t\t#### For type %d we have a matching Case C (-1, M)", i)
                    logging.debug( "\t\t\t\t\tWe can potentially add %s and vert %s but need to check further", testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\t\tbecause score is -1 and the comparison is M.")
                    logging.debug(" \t\t\t\t\tCould be X or U or M or D")
                    change        =   "M"  ## the new comparison variable to use
                    stopFlag      =   0   ## use this flag to determine if one needs to keep checking through the pairs.
                    ## now get the continue set of comparisons
                    ccount=0
                    potential_change=""
                    for e in nnetwork.edges_iter():     ### no need to go in order -- jsut look at all the other edges to see if there is an X
                        ccount +=1
                        d = nnetwork.get_edge_data(*e)
                        compArray = d['weight']
                        if compArray[i] is None:
                            print "Comparison is empty. Error! Stopping.\n\r\n\r"
                            sys.exit("Quitting due to errors.")
                        logging.debug( "\t\t\t\tType %d: Here is what we get for comparison # %d: ", i, ccount)  ## note that i is the current type
                        logging.debug( " \t\t\t\t\t inwardEdge:%s - outwardEdge: %s ", comparison, compArray)
                        if compArray[i] is  "U":
                            potential_change += "X"
                        ############################################
                        elif compArray[i] is  "X" or compArray[i] is  "D":
                            potential_change += "D"
                        ############################################
                        elif compArray[i] is  "M":
                            potential_change += "M"
                             ## in this case we have to keep going
                        else:
                            print "ERROR: missing value -- comparison is compArray[i]. Must have a value.\n\r\n\r" % compArray[i]
                            sys.exit("Quitting due to errors.")
                        logging.debug( "\t\t\t\t\t Since I got %s my potential new value is change.", compArray[i])
                        logging.debug( "\t\t\t\t\t Now going to continue to check pair of assemblages in the chain to look for problems.")

                    ## now decide what the change should be. Here are cases:
                    ## X exists, then it must be an X, Otherwise, D
                    if "X" in potential_change:
                        change = "X"
                    else:
                        change = "D"

                    ## in this case I dont think there are any errors possible. types can always go down from any other value
                    comparisonMap +=  change      ## use the value from above.
                    logging.debug( "\t\t\t\t#### Type %d: OK to add %s to vertices %s because ", i, testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\t score is -1 and the comparison is D. ComparisonMap is now %s ", comparisonMap)
                    if comparisonMap=="":
                        print "\n\rERROR: comparisonMap can't be empty. Bug here. \n\r\n\r"
                        sys.exit("Quitting due to errors.")

                #################################################################################     ## 0  U  #############
                elif difscore == 0  and comparison[i] is  "U":
                    # new score is match but comparison is Match. Okay
                    comparisonMap += "U"
                    logging.debug( "\t\t\t\#### tType %d:  Ok to add  %s to vertices %s because its a match.",i, testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\tType %d: ComparisonMap is now: %s ", i, comparisonMap)

                #################################################################################  ## 0 D  #############
                elif difscore == 0 and comparison[i] is  "D":
                      # new score is match but comparison is Match. Okay
                    comparisonMap += "D"
                    logging.debug( "\t\t\t\t#### Type %d:  Ok to add  %s to vertices %s because its a match.", i,testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\tType %d: ComparisonMap is now: %s ", i, comparisonMap)
                #################################################################################     ## 0 M #############
                elif  difscore == 0 and comparison[i] is  "M":
                    # new score is match but comparison is Match. Okay
                    comparisonMap += "M"
                    logging.debug( "\t\t\t\t##### Type %d:  Ok to add  %s to vertices %s because its a match.",i,  testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\tType %d: ComparisonMap is now: %s", i, comparisonMap)
                #################################################################################  ## -1 X #############
                elif difscore == -1  and comparison[i] is  "X":
                    # newscore is down but comparison is X. This means that there was already a peak
                    ## this is okay since it is down from a mode peak
                    logging.debug( "\t\t\t\t#### Type %d:  Ok to add  %s to vertices %s because ", i, testAssemblage, endAssemblage)
                    logging.debug( " \t\t\t\tscore is -1 and the comparison is D. ComparisonMap is now %s ", comparisonMap)
                    comparisonMap += "D"
                #################################################################################    ## 1  X #############
                elif difscore == 1 and comparison[i] is  "X":
                    ## new score is up but comparison is X.. no cant work because past peak
                    error += 1
                    break
                    logging.debug( "\t\t\t\t#### Type %d: Rejecting %s from %s]. We can't go up ", testAssemblage, endAssemblage)
                    logging.debug( " \t\t\t\t after a peak. so error. Error now error")
                ################################################################################# ## 0  X #############
                elif difscore == 0 and comparison[i] is  "X":
                   # newscore is down but comparison is X. This means that there was already a peak
                    ## this is okay since it is down from a mode peak
                    comparisonMap += "X"
                    logging.debug( "\t\t\t\t#### Type %d:  Ok to add  %s to vertices %s because ", i, testAssemblage, endAssemblage)
                    logging.debug( "\t\t\t\t is 0 and the comparison is X. ComparisonMap is now %s ", comparisonMap)
                else:
                    print "\t\t\t\tERROR!!!! Not found match combination! MUST FIX! Some combination is not being caught correctly... Exiting."
                    print "\t\t\t\tHere is the score of the differences in  for Type: %d %s:" % (i,difscore)
                    print "\t\t\t\tHere is the comparison value: %s " % comparison[i]
                    sys.exit("Quitting due to errors.")

                logging.debug( "\t\t\t\t#### Type %d:  Errors so far error")

            logging.debug("Checked out %s. Found %d errors.", testAssemblage, error)
            if error == 0:
                logging.debug("Found no errors!  Going to add %s to end of existing network at %s", testAssemblage, endAssemblage)
                logging.debug( "Original network: %s ", nnetwork.adjacency_list())
                logging.debug( "New comparison map is: %s ", comparisonMap)
                #first make a copy
                new_network = nnetwork.to_directed()
                new_network.graph["GraphID"]= solutionCount+1
                ## now add to this *new one*
                new_network.add_node(testAssemblage, name=testAssemblage,end=1,site="end")
                ## mark this vertice as the new "END"
                ## mark the interior vertice as not "END"
                new_network.add_node(endAssemblage, name=endAssemblage,site="middle", end=0)
                #### This adds the comparison to the new edge that has been added.
                new_network.add_path( [testAssemblage, endAssemblage], weight=comparisonMap, end=1, site="end", GraphID=solutionCount )
                if whichEnd==1:
                    new_network.graph["End1"]=testAssemblage
                    whichEnd += 1
                else:
                    new_network.graph["End2"]=testAssemblage
                    whichEnd = 1

                logging.debug( "Here's the new network (with addition): %s", new_network.adjacency_list())
                ## copy this solution to the new array of networks

                logging.debug( "----------------#############-------END CHECK LOOP---------#############-----------------" )

    if len(new_network.nodes())>0:
        return new_network
    else:
        return False

class CursesWindow(object):
    def __enter__(self):
        curses.initscr()

    def __exit__(self):
        curses.endwin()

def main():
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
    parser.add_argument('--pairwiseFile')
    parser.add_argument('--mst')
    parser.add_argument('--stats')
    parser.add_argument('--nosum')
    parser.add_argument('--screen')
    parser.add_argument('--allSolutions')
    parser.add_argument('--memusage')
    parser.add_argument('--inputfile')
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    ##################################################################################################

    global scr
    global screenFlag
    screenFlag=0
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
            traceback.print_exc()           # Print the exception

    ##################################################################################################
    if args['debug'] is not None:
        ## Logging
        logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    else:
        logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

    # start the clock to track how long this run takes
    global start
    start = datetime.now().time()
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
    except IOError, msg:
        logging.error("Cannot open %s. Error: %s", filename, msg)
        print("Cannot open %s. Error. %s " % filename, msg)
        sys.exit("Quitting due to errors.")

    ###############################################################################
    logging.debug("Going to open pairwise file it is exists.")
    if args['pairwiseFile'] is not None:
        openPairwiseFile(args['pairwiseFile'])

    ###############################################################################
    logging.debug("Going to open XY file if it exists.")
    largestX=0
    largestY=0
    distanceBetweenAssemblages={}
    if args['xyfile'] is not None:
        largestX,largestY,distanceBetweenAssemblages=openXYFile(args['xyfile'])

    ###############################################################################
    logging.debug("Assume threshold is 1.0 unless its specified in arguments.")
    threshold=1.0
    if args['threshold'] is not None:
        threshold=args[threshold]
    logging.debug("Going to create list of valid pairs for comparisons.")
    validAssemblagesForComparisons={}
    validAssemblagesForComparisons = thresholdDetermination(threshold, assemblages)
    typeFrequencyLowerCI={}
    typeFrequencyUpperCI={}

    ###############################################################################
    logging.debug("Now calculate the bootstrap comparisons based on specified confidence interval, if in the arguments.")
    if args['bootstrapCI'] is not None:
        bootstrapCI=1
        typeFrequencyLowerCI, typeFrequencyUpperCI = bootstrapCICalculation(assemblages, assemblageSize,1000,args['bootstrapSignificance'])

    ###############################################################################
    logging.debug("Now precalculating all the combinations between pairs of assemblages. This returns a graph with all pairs and the comparisons as weights.")
    pairGraph = preCalculateComparisons(assemblages,bootstrapCI,typeFrequencyUpperCI,typeFrequencyLowerCI)

    ###############################################################################
    logging.debug("Calculate all the valid triples.")
    triples=[]
    triples = findAllValidTriples(assemblages,pairGraph,validAssemblagesForComparisons,bootstrapCI,typeFrequencyLowerCI, typeFrequencyUpperCI)

    ###############################################################################
    stepcount = 0
    currentMaxSeriationSize = 3
    newNetworks=[]
    solutionCount=len(triples)
    maxEdges=2
    currentTotal = len(triples)
    solutions=[]
    while currentMaxSeriationSize <= maxSeriationSize:
        currentMaxSeriationSize += 1
        ### first time through copy the triples...
        if currentMaxSeriationSize==4:
            networks = triples
            solutions = triples
        else:
            networks = newNetworks

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
        for nnetwork in networks:
            logging.debug( "-----------------------------------------------------------------------------------")
            logging.debug( "Network: %s", nnetwork.adjacency_list())
            logging.debug( "-----------------------------------------------------------------------------------")
            ## find the ends
            ## given the ends, find the valid set of assemblages that can be potentially added
            ## this list is all assemblages meet the threshold requirements
            validNewNetwork = checkForValidAdditionsToNetwork(nnetwork,pairGraph,validAssemblagesForComparisons, assemblages, typeFrequencyLowerCI, typeFrequencyUpperCI,bootstrapCI,solutionCount)
            if not validNewNetwork is False:
                newNetworks.append(validNewNetwork)
                solutions.append(validNewNetwork)
                solutionCount +=1
                logging.debug("Solution count is now:  %d", solutionCount)
                if len(validNewNetwork.edges()) > maxEdges:
                    maxEdges = len(validNewNetwork.edges())
                currentTotal = len(newNetworks)

            if screenFlag>0:
                msg="Current Max Edges:  %d"% maxEdges
                scr.addstr(6,0,msg)
                msg="Sum of all solutions up to this step: %d" % solutionCount
                scr.addstr(7,0,msg)
                scr.addstr(8,43,"                                           ")
                msg="Current number of seriation linkages at this step: %d"% currentTotal
                scr.addstr(8,0,msg)
                scr.refresh()

if __name__ == "__main__":

    main()
    print "\n\rDone!\n\r"
    scr.nodelay(0)
    scr.getch()
    curses.endwin()