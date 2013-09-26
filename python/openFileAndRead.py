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

import matplotlib.pyplot as pltc

# screen
scr = curses.initscr()

## Logging
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

# start the clock to track how long this run takes
start = datetime.now().time()
print start

# start prettyprint (python Dumper)
pp = pprint.PrettyPrinter(indent=4)

## the numerous global arrays and values
assemblageNumber={}
assemblageValues={}
assemblageSize={}
assemblageFrequencies={}
assemblages={}
pairwise={}
pairwiseError={}
xAssemblage={}
yAssemblage={}
xyAssemblages=[]
distanceBetweenAssemblages={}
largestX = 0.0
largestY = 0.0
labels=[]
assemblageComparison={}
columns=0
validComparisonsArray=[]
validComparisonAssemblages=()
screenFlag = 0
numrows = len(assemblages)
triples =[]
triplettype =[]
tripletNames = []
tripletArray = []
typeFrequencyLowerCI={}
typeFrequencyUpperCI={}
solutions=[]
networks=[]

def all_pairs(lst):
    for p in itertools.permutations(lst):
        i = iter(p)
        yield zip(i,i)

def all_tuples(lst):
    tuples=itertools.combinations(lst, 3)
    useable_tuples=[]
    for e in tuples:
        useable_tuples.append(e)
    return useable_tuples

def openFile(filename):
    if screenFlag:
        msg1 = "Filename: %s " % filename
        scr.addstr(1,1,msg1)

    ## Read in the data
    # the input of the classes -- note these are absolute counts, not frequencies
    # might want to change that...\
    if screenFlag:
        scr.addstr(1,40,"STEP: Read in data...")

    try:
        print "trying to open ", filename
        file=open(filename,'r')
    except csv.Error as e:
        sys.exit('file %s does not open: %s') %( filename, e)
    reader = csv.reader(file, delimiter='\t', quotechar='|')
    count=0
    values=[]
    for row in reader:
        label=row[0]
        labels.append(labels)
        row.pop(0)
        columns = len(row)
        row = map(float, row)
        rowtotal = sum(row)
        freq=[]
        rowtotals=[]
        for r in row:
            freq.append(float(float(r)/float(rowtotal)))
            values.append(float(r))
        assemblages[ label ] = freq
        assemblageFrequencies[ label ]  = freq
        assemblageValues[ label ] = values
        assemblageSize[ label ]= rowtotal
        count +=1
    return count

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


#############################################  THRESHOLD DETERMINATION ####################################
## first compare each assemblage and see if the threshold is exceeded.
## By this I mean the maximum difference between frequencies of any type is greater than what is specified.
## When threshold = 0, all combinations are used. Later the "ends" of solutions are not evaluated if the
## difference between the last assemblage and any other free assemblage is > the threshold.
## This arbitrary setting is to keep from arbitrary solutions being stuck on that come from the "ends"
## of solutions.
##
## Precalculate all of the max differences between types in assembalge pairs.

def threshholdDetermination(threshhold):
    ## We use the Math::Combinatorics to get all the combinations of 2
    pairs = all_pairs(assemblageNumber)
    ## Go through all of the combinations
    for combo in pairs:
        logging.debug("comparing combination of %s and %s ", labels[ combo[0] ], labels[ combo[1] ] )
        pairname  = labels( combo[0] ) + " * " + labels(combo[1] )

        maxDifference = 0
        assemblage1 = assemblageFrequencies[combo[0]]
        assemblage2 = assemblageFrequencies[combo[1]]
        i=0
        while i < columns:
            ass1 = float(assemblage1[i])
            ass2 = float(assemblage2[i])
            diff = abs( ass1 - ass2 )
            if diff > maxDifference :
                maxDifference = diff

        assemblageComparison[ pairname ] = maxDifference
    cAssemblages=()

    ############## pre calcuate the valid pairs of assemblages and stuff into hash ############################
    for label in labels:
        cAssemblages=()
        for comparativeLabel in labels:
            test = label+"*"+comparativeLabel
            if (assemblageComparison[ test ]  <= threshhold) and (comparativeLabel != label):
                cAssemblages.append(comparativeLabel)
        validComparisonsArray[ label ] = cAssemblages


def confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

########################################### BOOTSTRAP CI SECTION ####################################
def bootstrapCI(bootsize=1000, confidenceInterval=0.95):
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
    countup=0
    ## for each assemblage
    logging.debug("Calculating bootstrap confidence intervals")
    # for each assemblage

    for currentLabel in sorted( assemblageFrequencies.iterkeys()):
        label = labels[countup]
        a =  assemblageFrequencies[ currentLabel ]
        currentAssemblageSize = assemblageSize[ currentLabel ]

        ## create an array of arrays - one for each type
        arrayOfStats = []
        for c in assemblageFrequencies:
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
                rand  = random()              ## random number from 0-1
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
                classCount +=1
            loop -=1

        lowerCI=[]
        upperCI=[]
        for stat in arrayOfStats:
            upper=0
            lower=0
            mean=0
            mean, upper, lower = confidence_interval(stat, confidence=confidenceInterval)
            lowerCI.append( lower)
            upperCI.append( upper )

        typeFrequencyLowerCI[ label ] = lowerCI
        typeFrequencyUpperCI[ label ] = upperCI
        results = 0
        countup += 1


########################################### FIND ALL THE VALID TRIPLES  ####################################
########################################### #################################### ###########################
def findAllValidTriples(bootstrapCI):

    directionstate=""
    error=0
    numberOfTriplets = 0

    if screenFlag>0:
        scr.addstr(1,40,"STEP: Find valid triples....      ")

    permutations =all_tuples(assemblageNumber)

    for permu in permutations:
    tripletname = labels[ permu[0] ] + " * "+ labels[ permu[1] ] + " * "+ labels[ permu[2] ]
    comparison12 = ""
    comparison23 = ""
    error= 0

    for i in range(0,columns):

        ass1 = assemblages[ permu[0] ][i]
        ass2 = assemblages[ permu[1] ][i]
        ass3 = assemblages[ permu[2] ][i]

        ## first compare assemblages 1 and 2
        if bootstrapCI:
            upperCI_1 = typeFrequencyUpperCI[ labels[ permu[0] ] ][i]
            lowerCI_1 = typeFrequencyUpperCI[ labels[ permu[0] ] ][i]
            upperCI_2 = typeFrequencyUpperCI[ labels[ permu[1] ] ][i]
            lowerCI_2 = typeFrequencyUpperCI[ labels[ permu[1] ] ][i]
            dif1 = ass1 - ass2
            if upperCI_1 < lowerCI_2:
                difscore = -1
            elif lowerCI_1 > upperCI_2:
                difscore = 1
            else:
                difscore = 0

        else:   ### if the bootstrapCI is not being used
            dif1 = ass1 - ass2
            if dif1 < 0:
                difscore = -1
            if dif1 > 0:
                difscore = 1
            if dif1 == 0:
                difscore = 0

        ## now compare assemblages 2 and 3
        if bootstrapCI:  ## boostrap confidence intervals
            upperCI_2 = typeFrequencyUpperCI[ labels[ permu[1] ] ][i]
            lowerCI_2 = typeFrequencyUpperCI[ labels[ permu[1] ] ][i]
            upperCI_3 = typeFrequencyUpperCI[ labels[ permu[2] ] ][i]
            lowerCI_3 = typeFrequencyUpperCI[ labels[ permu[2] ] ][i]
            if upperCI_2 < lowerCI_3:
                difscore = -1
            elif lowerCI_2 > upperCI_3:
                difscore = 1
            else:
                difscore = 0

        else:         ### if the bootstrapCI is not being used
            dif2 = ass3 - ass2
            if dif2 < 0:
                difscore2 = -1
            if dif2 > 0:
                difscore2 = 1
            if dif2 == 0:
                difscore2 = 0

        if ( ( difscore == 1 ) and ( difscore2 == 1 ) ):     ## F1 > F2 < F3 ## criteria not met
            error += 1
            next
        elif ( ( difscore == 1 ) and ( difscore2 == -1 ) ) :   ## F1 > F2 > F3 OK
            comparison12 += "U"
            comparison23 += "D"
        elif ( ( difscore == -1 ) and ( difscore2 == -1 ) ) :  #  F1 < F2 >F3 OK
            comparison12 += "X"
            comparison23 += "X"
        elif ( ( difscore == -1 ) and ( difscore2 == 1 ) ):   # F1 < F2 < F3
            comparison12 += "D"
            comparison23 += "U"
        elif ( ( difscore == 0 ) and ( difscore2 == 1 ) ):
            comparison12 += "M"
            comparison23 += "U"
        elif ( ( difscore2 == 0 ) and ( difscore == 1 ) ):
            comparison12 += "U"
            comparison23 += "M"
        elif ( ( difscore == 0 ) and ( difscore2 == -1 ) ):
            comparison12 += "M"
            comparison23 += "D"
        elif ( ( difscore == -1 ) and ( difscore2 == 0 ) ):
            comparison12 += "D"
            comparison23 += "M"
        elif ( ( difscore == 0 ) and ( difscore2 == 0 ) ):
            comparison12 += "M"
            comparison23 += "M"
        else:
            error += 1
            print "\n\rNo match to our possibility of combinations. Difscore 1: difscore Difscore 2: difscore2 \n\r"
            print "I must quit. Debugging required.\n\r"
            exit()

    if error == 0:

        # uses NetworkX
        net = nx.Graph(GraphID=numberOfTriplets,End1=labels[permu[0]],End2=labels[permu[2]])

        net.add_node( name=labels[ permu[0] ],site="end",end=1 )
        net.add_node( name=labels[ permu[1] ],site="middle",end=0)
        net.add_node( name=labels[ permu[2] ], site="end",end=1)
        net.add_edge( labels[ permu[1] ],labels[ permu[0] ],weight=comparison12, GraphID=numberOfTriplets,end=1)
        net.add_edge( labels[ permu[1] ],labels[ permu[2] ],weight=comparison23, GraphID=numberOfTriplets,end=1)
        logging.debug("VALID SOLUTION: %s * %s * %s " , (labels[ permu[0] ],labels[ permu[1] ], labels[ permu[2] ]))
        logging.debug("VALID SOLUTION: \t  %st  ---   comparison23\n",  (comparison12, comparison23))
        triples.append(net)
        numberOfTriplets+=1
    error = 0

def checkForValidAdditionsToNetwork(nnetwork):

    whichEnd = 0

    for endAssemblage in nnetwork.nodes(data=True):
        whichEnd += 1 ## either a 1 or a 2
        ##print "The end of assemblages of nnetwork are: ". nnetwork->get_graph_attribute("End_1") . "  and ". nnetwork->get_graph_attribute("End_2"), "\n"

        if (endAssemblage[])

        foreach testAssemblage ( @ { validComparisonsArray{ endAssemblage } }) {

            DEBUG  and print "\t\tChecking assemblage: ", testAssemblage, " to see if it fits on the end of the current solution.\n"
            DEBUG  and print "\t\tFirst check to see if it is included already. If it has, move on.\n"
            my @vertices = nnetwork->vertices

            if ( (! grep { _ eq testAssemblage} @vertices) ) {   ## if the assemblage is NOT in the list of existing vertices.

                # get the exterior vertices (should be 2)
                DEBUG  and print "\t\tFind the ends of the network. Do this by getting all the vertices \r\n"
                DEBUG  and print " \t\tand looking for the ones with only 1 connection. There should be just 2 here.\r\n"
                ## loop through all of the edges to see if they can be stuck on the ends of the networks.

                comparisonMap

                logging.debug( "\t\t", endAssemblage, " is on the edge since it only has one vertice.\r\n"

                #################### THRESHOLDING STUFF ####################################
                #first determine if the pairs are within the threshold value (0 = all assemblages)
                #pairname = endAssemblage. " * " . testAssemblage
                #diff     = assemblageComparison{pairname}
                #logging.debug( "\t\t\tFor pairname the max frequency difference is diff.\n"
                #if (diff eq undef) {
                #       print "\n\rError: pairname: pairname not found in hash lookup!\n\r"
                #       print "Probably a problem with the names of the assemblages. Check for weird characters. Exiting.\n\r"
                #       exit()
                #}
                ## go through process only if threshold is 0 or difference value is below threshold
                ## this should mean that network will not grow unless the above conditions are met.
                error = 0
                ##if (  (threshold>0 ) and (diff > threshold ) )  {
                ##   logging.debug( "\t\t\tThreshold = threshold and Diff = diff. Since diff < threshold, continue.\n"
                ##   error++ # this should ensure future failure....
                ##    next
                ##}
                ############################################################################

                my @newassemblage = @{ assemblageFrequencies{ testAssemblage } }
                my @oldassemblage = @{ assemblageFrequencies{ endAssemblage } }

                #### FIND INNER EDGE RELATIVE TO THE EXISTING END ASSEMBLAGE ##############
                my @neighbors = nnetwork->neighbors(endAssemblage)
                if (scalar(@neighbors) > 1 || scalar(@neighbors)==0){
                    print "\r\n\r\n\r\nThere are too many or two few neighbors (should only be 1!). Error!\n\r"
                    print "\r\nWe are testing endAssemblage and got ", Dumper(\@neighbors)
                    print Dumper(nnetwork)
                    print nnetwork
                    exit()
                }
                logging.debug( "\t\t\t The number of neighbors at endAssemblage is ", scalar(@neighbors), " (should be just one).\n"
                g = nnetwork->get_edge_weight( neighbors[0], endAssemblage )
                logging.debug( "\t\t\tThere should be just 1 neighbor to endAssemblage and that is: neighbors[0]\r\n"
                logging.debug( "\t\t\t\t it has a relation of g\r\n"
                outerEdge= endAssemblage
                innerEdge= neighbors[0]
                ##########################################################################


                my @comparison = split //, g
                cols = scalar(@newassemblage)
                comparisonMap =""

                # go through the columns
                if (!error) {
                       for ( i = 0  i < cols  i++ ) {
                            my ( difscore, difscore2 )
                            val1 = newassemblage[i]
                            val2 = oldassemblage[i]
                            logging.debug( "\t\tComparing Assemblage: testAssemblage    and    Assemblage: endAssemblage            ########\n "
                            logging.debug( "\t\t\t\tType i - Type i - Type i - Type i - Type i - Type i - Type i  ########  \n"
                            logging.debug( "\t\t\t\tType i:  testAssemblage 1: ", newassemblage[i], " endAssemblage 2: ", oldassemblage[i], "\n"
                             ## ALL COMBINAATIONS
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
                            if (bootstrapCI ) {
                                upperCI_test = typeFrequencyUpperCI{ testAssemblage  }->[i]
                                lowerCI_test = typeFrequencyUpperCI{ testAssemblage  }->[i]
                                upperCI_end = typeFrequencyUpperCI{ endAssemblage }->[i]
                                lowerCI_end = typeFrequencyUpperCI{ endAssemblage }->[i]
                                if ( upperCI_test < lowerCI_end )  {
                                    difscore = -1
                                } elsif (lowerCI_test > upperCI_end ) {
                                    difscore = 1
                                } else {
                                    difscore = 0
                                }
                            } else {
                                dif1 = newassemblage[i] - oldassemblage[i]
                                if ( dif1 < 0 )  { difscore = -1 }
                                if ( dif1 > 0 )  { difscore = 1 }
                                if ( dif1 == 0 ) { difscore = 0 }
                            }
                            logging.debug( "\t\t\t\tType i: - comparison is:  ", comparison[i], " a score of: ", difscore, "\n"
                             if (   ( difscore == 1 ) && ( comparison[i] =~ "U" ) ) {                                                 #### 1 U
                                comparisonMap .= "U"
                                logging.debug(  "\t\t\t\tType i: Got a difscore of 1 and a comparison of a U. This works. \n"
                                logging.debug( " \t\t\t\tAdding testAssemblage to vertices endAssemblage\n"

                            } elsif (( difscore == 1 ) && ( comparison[i] =~ "M" ) ) {                                                ### 1 M
                                # this is okay - its a match and the new value is greater. New value shoudl be U
                                # need to find what was happening on the previous comparison to know whether this needs
                                # to be up or down.
                                logging.debug( "\t\t\t\tType i: Got a difscore of 1 and a comparison of a M. This could be okay.\n"
                                xerror        =     0
                                logging.debug( "\t\t\t\tType i:   Matching case A (1, M)  \n"
                                logging.debug( "\t\t\t\t\ This will only work if there no Xs anywhere previously OR if the opposite end doesnt ALSO go up!\n"
                                stopFlag      =     0   ## use this flag to determine if one needs to keep checking through the pairs.
                                my (outwardEdge, inwardEdge)
                                old_inner= outerEdge
                                logging.debug( "\t\t\t\t\t ", nnetwork
                                my @currentEdges  = nnetwork->edges05(innerEdge)
                                ccount
                                foreach checkEdge (@currentEdges) {     ### no need to go in order -- jsut look at all the other edges to see if there is an X
                                   logging.debug( "\t\t\t\t\ Now on @checkEdge[0] @checkEdge[1] \n"
                                   inwardEdge= @checkEdge[0]
                                   outwardEdge = @checkEdge[1]
                                   comparison = nnetwork->get_edge_weight( outwardEdge, inwardEdge) || nnetwork->get_edge_weight( inwardEdge, outwardEdge )
                                   my @compArray = split //, comparison
                                   if (!compArray[i]) {
                                      print "Comparison is empty. Error! Stopping.\n\r\n\r"
                                      exit()
                                   }
                                   logging.debug( "\t\t\t\tType i: Here is what we get for comparison # ccount \n"  ## note that i is the current type
                                   logging.debug( " \t\t\t\t\t inwardEdge - outwardEdge: ", comparison, "->", compArray[i], "\n"
                                   if (compArray[i] =~ "X" || compArray[i] =~ "U" )  {
                                         xerror = 1  ### BLARGH a previous X or an UP ! This will not be tolerated!
                                         stopFlag = 1 ## We can stop
                                   }
                                   logging.debug( "\t\t\t\t\t Since I got compArray[i] my potential new value is still X.\n"
                                   logging.debug( "\t\t\t\t\t Now going to get the next pair of assembalges to examine in the chain\n"
                                   ccount++
                                }
                                if ( xerror ) {
                                    error++
                                    last
                                } else {
                                    comparisonMap .= "U"
                                    logging.debug( "\t\t\t\t Type i: For this type, OK to add testAssemblage to vertices endAssemblage \n"
                                    logging.debug( "\t\t\t\t\t No previous X values anywhere. \n"
                                    logging.debug( "\t\t\t\t Type i: Adding an U to the comparisons for type i. \n"
                                    logging.debug( "\t\t\t\t\t Comparison map is now comparisonMap\n"
                                }
                            }
                            elsif
                              ## error the new value is greater but shoudl be less. Error!
                              ( ( difscore == 1 ) &  ( comparison[i] =~ "D" ) ) {                                                  ## 1 D
                                #print "mismatch!\n"
                                error++
                                next
                                logging.debug( "\t\t\t\tType i: Value 1:  ", newassemblage[i], " value 2: ", oldassemblage[i], "\n"
                                logging.debug( "\t\t\t\tType i: Comparison is:  ", comparison[i], " a score of: ", difscore, "\n"
                                logging.debug( "\t\t\t\tType i: Rejecting testAssemblage from endAssemblage \n"
                                logging.debug( "\t\t\t\t\t because value is 1 and comparison is D.\n"
                            }
                            elsif
                              ## new score is less and the Comparison is up .. Error!
                              ## note -- need to check to see if there is a previous change in direction because
                              ## its possible that the this is a mode change.
                              ## do this by logging all modes in the original triplet -- keep track of modes
                              ## per type
                              (    ( difscore == -1 ) && ( comparison[i] =~ "U" ) )                                               # -1 U
                            {
                                ## first check to see if there is already and X in this column somewhere else.
                                xerror   = 0
                                logging.debug( "\t\t\t\tType i:  Case B (-1, U). Potentially can add testAssemblage and vert endAssemblage\n "
                                logging.debug( "\t\t\t\tType i:  But need to check the rest of the chain for X's (can't already be an X). \n"
                                stopFlag      =   0   ## use this flag to determine if one needs to keep checking through the pairs.
                                my (outwardEdge, inwardEdge)
                                my @currentEdges  = nnetwork->edges05(innerEdge)
                                ccount
                                foreach checkEdge (@currentEdges) {     ### no need to go in order -- just look at all the other edges to see if there is an X
                                   logging.debug( "\t\t\t\t\ Now on @checkEdge[0] @checkEdge[1] \n"
                                   inwardEdge= @checkEdge[0]
                                   outwardEdge = @checkEdge[1]
                                   comparison = nnetwork->get_edge_weight( outwardEdge, inwardEdge) || nnetwork->get_edge_weight( inwardEdge, outwardEdge )
                                   my @compArray = split //, comparison
                                   if (!compArray[i]) {
                                      print "Comparison is empty. Error! Stopping.\n\r\n\r"
                                      exit()
                                   }
                                   logging.debug( "\t\t\t\tType i: Here is what we get for comparison # ccount \n"  ## note that i is the current type
                                   logging.debug( " \t\t\t\t\t inwardEdge - outwardEdge: ", comparison, "->", compArray[i], "\n"
                                   if (compArray[i] =~ "X")  {
                                         xerror += 1  ### BLARGH a previous X! This will not be tolerated!
                                         stopFlag = 1 ## We can stop
                                   }
                                }
                                if (xerror > 0) {
                                    error += 1
                                    last
                                    logging.debug( "\t\t\t\tType i: Rejecting testAssemblage from endAssemblage) because there was an X \n"
                                    logging.debug( "\t\t\t\t\t  This would make it multimodal - so error.\n"
                                } else {
                                    comparisonMap .= "X"   ## this is an X unless there is an error....
                                    logging.debug( "\t\t\t\tType i:Definitely OK to add testAssemblage to vertices endAssemblage because score \n"
                                    logging.debug( "\t\t\t\t\tis -1 and the comparison is U but no other Xs in the previous linkages.\n"
                                    logging.debug( "\t\t\t\tType i: Adding an X to the comparisons for type i. \n"
                                    logging.debug( "\t\t\t\t\tComparison map is now comparisonMap\n"
                                }    #end if if check error (xerror)
                            }
                            elsif ## new score is less and the comparison is down .. Yes!
                              (    ( difscore == -1 ) && ( comparison[i] =~ "D" ) )                                          ## -1   D
                            {
                                comparisonMap .= "D"
                                logging.debug( "\t\t\t\tType i: Adding a D to the comparisons for type i. Comparison map is now comparisonMap\n"

                            } elsif  (    ( difscore == -1 )  || ( comparison[i] =~ "M" ) ) {                               ## -1 M

                                # new score is less but comparison is Match. Okay
                                #vHere    = endAssemblage
                                xerror   = 0  ## count for errors
                                logging.debug( "\t\t\t\t #### For type i we have a matching Case C (-1, M) \n"
                                logging.debug( "\t\t\t\t\t We can potentially add testAssemblage and vert endAssemblage but need to check further\n"
                                logging.debug( "\t\t\t\t\t because score is -1 and the comparison is M.\n"
                                logging.debug(" \t\t\t\t\t Could be X or U or M or D\n"
                                change        =   "M"  ## the new comparison variable to use
                                stopFlag      =   0   ## use this flag to determine if one needs to keep checking through the pairs.
                                ## now get the next set of comparisons
                                logging.debug( "\t\t\t\t ", nnetwork
                                my @currentEdges  = nnetwork->edges05(innerEdge)
                                ccount
                                my (outwardEdge, inwardEdge)
                                foreach checkEdge (@currentEdges) {     ### no need to go in order -- jsut look at all the other edges to see if there is an X
                                   logging.debug( "\t\t\t\t\ Now on @checkEdge[0] @checkEdge[1] \n"
                                   inwardEdge= @checkEdge[0]
                                   outwardEdge = @checkEdge[1]
                                   comparison = nnetwork->get_edge_weight( outwardEdge, inwardEdge) || nnetwork->get_edge_weight( inwardEdge, outwardEdge )
                                   my @compArray = split //, comparison
                                   if (!compArray[i]) {
                                      print "Comparison is empty. Error! Stopping.\n\r\n\r"
                                      exit()
                                   }
                                   logging.debug( "\t\t\t\tType i: Here is what we get for comparison # ccount \n"  ## note that i is the current type
                                   logging.debug( " \t\t\t\t\t inwardEdge - outwardEdge: ", comparison, "->", compArray[i], "\n"
                                   if ( compArray[i] =~ "U" ) {
                                         change = "X"
                                         stopFlag=1 ## we can stop
                                   } elsif ((compArray[i] =~ "X") or (compArray[i] =~ "D"))  {
                                         change = "D"
                                         stopFlag = 1 ## We can stop
                                   } elsif (compArray[i] =~ "M") {
                                         change = "M"
                                         ## in this case we have to keep going
                                   } else {
                                      print "ERROR: missing value -- comparison is compArray[i]. Must have a value.\n\r\n\r"
                                      exit()
                                   }
                                   logging.debug( "\t\t\t\t\t Since I got compArray[i] my potential new value is change.\n"
                                   logging.debug( "\t\t\t\t\t Now going to get the next pair of assembalges to examine in the chain\n"

                                   ccount++
                                }
                                ## in this case I dont think there are any errors possible. types can always go down from any other value
                                comparisonMap = comparisonMap . change      ## use the value from above.
                                logging.debug( "\t\t\t\tType i: OK to add testAssemblage to vertices endAssemblage because \n"
                                logging.debug( "\t\t\t\t score is -1 and the comparison is D. ComparisonMap is now comparisonMap\n"
                                if (!comparisonMap ) {
                                   print "\n\rERROR: comparisonMap can't be empty. Bug here. \n\r\n\r"
                                   exit()
                                }
                            } elsif  # new score is match but comparison is Match. Okay
                              ( ( difscore == 0 ) && (comparison[i] =~ "U") ) {                                              ## 0  U
                                comparisonMap .= "U"
                                logging.debug( "\t\t\t\tType i:  Ok to add  testAssemblage to vertices endAssemblage because its a match. \n"
                                logging.debug( "\t\t\t\tType i: ComparisonMap is now comparisonMap\n"
                            } elsif  # new score is match but comparison is Match. Okay
                              ( ( difscore == 0 ) && (comparison[i] =~ "D") ) {                                              ## 0 D
                                comparisonMap .= "D"
                                logging.debug( "\t\t\t\tType i: Ok to add testAssemblage to vertices endAssemblage because its a match. \n"
                                logging.debug( "\t\t\t\tType i: ComparisonMap is now comparisonMap\n"
                            } elsif  # new score is match but comparison is Match. Okay
                              ( ( difscore == 0 ) && (comparison[i] =~ "M") ) {                                              ## 0 M
                                comparisonMap .= "M"
                                logging.debug( "\t\t\t\tType i:  Ok to add  testAssemblage to vertices endAssemblage because its a match. \n"
                                logging.debug( "\t\t\t\tType i: ComparisonMap is now comparisonMap\n"
                            } elsif # newscore is down but comparison is X. This means that there was already a peak
                              (    ( difscore == -1 ) && ( comparison[i] =~ "X" ) )                                          ## -1 X
                            {
                                ## this is okay since it is down from a mode peak
                                logging.debug( "\t\t\t\tType i:  Ok to add  testAssemblage to vertices endAssemblage because \n"
                                logging.debug( " \t\t\t\tscore is -1 and the comparison is D. ComparisonMap is now comparisonMap\n"
                                comparisonMap .= "D"
                            } elsif (( difscore == 1 ) && ( comparison[i] =~ "X" ) )                                         ## 1  X
                            {
                                ## new score is up but comparison is X.. no cant work because past peak
                                error++
                                last
                                logging.debug( "\t\t\t\tType i: Rejecting testAssemblage from endAssemblage]. We can't go up \n"
                                logging.debug( " \t\t\t\tafter a peak. so error. Error now error\n"
                            } elsif # newscore is down but comparison is X. This means that there was already a peak
                              (    ( difscore == 0 ) && ( comparison[i] =~ "X" ) )                                           ## 0  X
                            {
                                ## this is okay since it is down from a mode peak
                                comparisonMap .= "X"
                                logging.debug( "\t\t\t\tType i:  Ok to add  testAssemblage to vertices endAssemblage because score \n"
                                logging.debug( "\t\t\t\t is 0 and the comparison is X. ComparisonMap is now comparisonMap\n"
                            } else {
                                print "\n\r\t\t\t\tERROR!!!! Not found match combination! MUST FIX! Some combination\n\r "
                                print "\t\t\t\t is not being caught correctly... Exiting.\n\r"
                                print "\t\t\t\tHere is the score of the differences in  for Type i: difscore\n\r"
                                print "\t\t\t\tHere is the comparison value: ", comparison[i], "\n\r"
                                error++
                                exit()
                            }
                            logging.debug( "\t\t\t\tType i:  Error so far error\n\r\n\r"
                        }
                    }
                    if ( error==0 )  {
                        logging.debug( "--------------------------------------------------\n\r\n\r"
                        logging.debug( "Original network: ", nnetwork, "\n\r"
                        logging.debug( "New comparison map is: comparisonMap\n\r"

                        ## no errors so add vertice added to the new network
                        #oldedgenum = nnetwork->edges
                        my @vertices = nnetwork->vertices
                        if ( ! grep { _ eq testAssemblage} @vertices) {
                            #first make a copy
                            #print Dumper(nnetwork)
                            new_network = nnetwork->deep_copy_graph
                            #print Dumper(new_network)
                            solutionCount++    ## increment the # of solutions
                            new_network->set_graph_attribute("GraphID", solutionCount)
                            ## now add to this *new one*
                            new_network->add_vertex(testAssemblage)
                            ## mark this vertice as the new "END"
                            new_network->set_vertex_attribute(testAssemblage,"End", "1" )
                            ## mark the interior vertice as not and "END"
                            new_network->set_vertex_attribute(endAssemblage,"End", "0" )
                            #### This adds the comparison to the new edge that has been added.
                            new_network->add_weighted_edge( endAssemblage, testAssemblage,  comparisonMap )
                            ## Mark this edge and a new end edge
                            new_network->set_edge_attribute(endAssemblage , testAssemblage, "End", "1")
                            new_network->set_edge_attribute(endAssemblage , testAssemblage, "GraphID", solutionCount)
                            ## mark the previous edge as no longer the end
                            my @n = new_network->neighbours(endAssemblage)  ## assume 0 is the only neighbor (should be only one!)
                            new_network->set_edge_attribute(n[0], endAssemblage , "End", "0")

                            #print "Which end: ", whichEnd, "\n"
                            if (whichEnd==1) {
                                new_network->set_graph_attribute("End_1", testAssemblage)
                            } else {
                                new_network->set_graph_attribute("End_2", testAssemblage)
                            }
                             ## increment the ends (we start with 0, then 1)

                            logging.debug( "New network (with addition): ", new_network, "\n\r"
                            ## copy this solution to the new array of networks
                            #print Dumper(new_network)
                            #print new_network, "\n"
                            push @newnets, new_network   ## contains solutions for just this step - add to list
                            #push @networks, new_network
                            if (nosum==0 ) {
                                push @solutions, new_network ## this is a list of all the solutions (shallow a)
                            }

                            currentTotal =  scalar(@newnets)
                            if ((new_network->unique_edges) > maxEdges) {
                                maxEdges = new_network->unique_edges
                                scr.addstr(6,1,"Current Max Edges: maxEdges   ")
                            }
                            scr.addstr(7,1,"Sum of all solutions up to this step: solutionCount")
                            scr.addstr(8,43,"                                           ")
                            scr.addstr(8,1,"Current number of seriation linkages at this step: currentTotal")
                            if memusage:
                                text = "Memory used for array @ newnets: ". total_size(\@newnets)
                                scr.addstr(9,1,text)
                            }
                            logging.debug( "-------------------------------------------------\n\r"
                        }

                    }
            }    # end of iterate through the existing network link
        }    # end of if assemblage not already in netowrk check
    }


def main():
    parser = argparse.ArgumentParser(description='Conduct seriation analysis')
    parser.add_argument('--debug')
    parser.add_argument('--bootstrapCI')
    parser.add_argument('--bootstrapSignificance', type=float)
    parser.add_argument('--filtered')
    parser.add_argument('--largestonly')
    parser.add_argument('--individualfileoutput')
    parser.add_argument('--excel')
    parser.add_argument('--threshhold')
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
    args = vars(parser.parse_args())

    filename=args['inputfile']
    if filename is "":
        print "You must enter a filename to continue."
        exit()
    #file = filename[0,-4]
    print "Trying to open: ", filename
    maxSeriations = openFile(filename)

    print pp.pprint(assemblageFrequencies)
    print pp.pprint(assemblageSize)
    print pp.pprint(assemblages)

    if args['screen']:
        ## Set up the screen display (default).
        screenFlag = 1
        ## the debug option should not use this since it gets messy
        scr.refresh()  # clear the screen

    if args['pairwiseFile']:
        openPairwiseFile(args['pairwiseFile'])

    if args['xyfile']:
        openXYFile(args['xyfile'])

    if args['threshhold']:
        threshholdDetermination(args['threshhold'])

    if args['boostrapCI']:
        bootstrapCI(1000,args['bootstrapSignificance'])

    findAllValidTriples(args['bootstrapCI'])
    stepcount = 0
    #currentMaxSeriations = 3
    currentMaxSeriationSize = 4
    while currentMaxSeriationSize <= maxSeriations:
        ### first time through copy the triples...
        if currentMaxSeriationSize==4:
            solutions = triples
            networks = triples
        stepcount += 1
        logging.debug("_______________________________________________________________________________________")
        logging.debug("Step number:  currentMaxSeriationSize")
        logging.debug("_______________________________________________________________________________________")
        if screenFlag>0:
            scr.addstr(4,1,"Step number:                     ")
            scr.addstr(4,1,"Step number:  currentMaxSeriationSize ")
            scr.addstr(5,1,"Number of solutions from previous step:         ")
            scr.addstr(5,1,"Number of solutions from previous step: netnum")
        netnum=len(networks)
        logging.debug("Number of solutions from previous step: netnum")
        match = 0      ## set the current match to zero for this step (sees if there are any new solutions for step)
        ## look through the set of existing valid networks.
        for nnetwork in networks:
            logging.debug( "-----------------------------------------------------------------------------------")
            logging.debug( "Network: %s", nnetwork)
            logging.debug( "-----------------------------------------------------------------------------------")
            ## find the ends
            ## given the ends, find the valid set of assemblages that can be potentially added
            ## this list is all assemblages meet the threshold requirements

            validNewNetwork = checkForValidAdditionsToNetwork(nnetwork)


if __name__ == "__main__":
    main()