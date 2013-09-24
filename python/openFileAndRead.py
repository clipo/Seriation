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
    error        = 0

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
        net = nx.Graph(GraphID=numberOfTriplets,End1=labels[permu[0]],End2=labels[permu[2]])

        net.add_node( name=labels[ permu[0] ],site="end",end=1 )
        net.add_node( name=labels[ permu[1] ],site="middle",end=0)
        net.add_node( name=labels[ permu[2] ], site="end",end=1)


        net.add_edge( labels[ permu[1] ],labels[ permu[0] ],weight=comparison12, GraphID=numberOfTriplets,end=1)
        net.add_edge( labels[ permu[1] ],labels[ permu[2] ],weight=comparison23, GraphID=numberOfTriplets,end=1)



        logger.debug("VALID SOLUTION: %s * %s * %s " , (labels[ permu[0] ],labels[ permu[1] ], labels[ permu[2] ]))
        logger.debug("VALID SOLUTION: \t  %st  ---   comparison23\n",  (comparison12, comparison23))
        triples.append(net)
        numberOfTriplets+=1
    error = 0



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
    openFile(filename)

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

    findAllValidTriples()

if __name__ == "__main__":
    main()