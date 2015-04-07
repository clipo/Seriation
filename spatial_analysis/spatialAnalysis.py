__author__ = 'carllipo'

import csv
import logging as logger
import itertools
from itertools import chain
from pylab import *
import matplotlib.pyplot as plt
import argparse
import numpy as np


class spatialAnalysis():
    color = ["b", "r", "m", "y", "k", "w", (0.976, 0.333, 0.518), (0.643, 0.416, 0.894),
             (0.863, 0.66, 0.447), (0.824, 0.412, 0.118)]
    def __init__(self):

        self.inputFile = ""
        self.outputDirectory = ""

        self.assemblageSize = {}
        self.assemblageFrequencies = {}
        self.assemblages = {}
        self.countOfAssemblages = 0
        self.assemblageValues = {}
        self.labels = []
        self.numberOfClasses = 0
        self.nodeSizeFactor =0
        self.xAssemblage = {}
        self.yAssemblage = {}
        self.xyAssemblages = []
        self.distanceBetweenAssemblages = {}
        self.largestX = 0
        self.largestY = 0
        self.distanceBetweenAssemblages = {}
        self.args={}
        self.totalAssemblageSize=0
        self.FalseList=[None,0,False,"None","0","False"]
        self.sumOfDifferencesBetweenPairs={}
        self.validComparisonsHash={}
        self.vmr_for_blocks = {}


    def openFile(self):

        try:
            logger.debug("trying to open: %s ", self.args['inputfile'])
            file = open(self.args['inputfile'], 'rU')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s",  self.args['inputfile'], e)
            sys.exit('file %s does not open: %s') % (  self.args['inputfile'], e)

        reader = csv.reader(file, delimiter='\t')
        values = []
        rowcount=0
        for row in reader:
            row = map(str, row)
            if rowcount==0:
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

    def preCalculateDistancesBetweenPairs(self):
        logger.debug("calculate differences between pairs")
        list_of_distances = []
        pairs = self.all_pairs(self.assemblages)
        for pair in pairs:
            diff = self.calculateRootSumDiffSquares(pair[0], pair[1])
            list_of_distances.append(diff)
            key1 = pair[0] + "*" + pair[1]
            key2 = pair[1] + "*" + pair[0]
            self.sumOfDifferencesBetweenPairs[key1] = diff
            self.sumOfDifferencesBetweenPairs[key2] = diff
            if diff == 0:
                logger.info("Potential problem:  %s and %s are identical in their frequencies. ", pair[0],pair[1])
        ## graph results
        figure_label = self.inputFile[0:-4]
        plt.xlabel("Square Root of Sum of Differences Squared")
        plt.ylabel('Count')
        plt.title(figure_label)
        hist, bins = np.histogram(list_of_distances, bins=10)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()

    def VMR(self, assemblages):

        varlist = []

        pairs = self.all_pairs(assemblages)
        for pair in pairs:
            diff = self.calculateRootSumDiffSquares(pair[0], pair[1])
            varlist.append( diff)

        variance = var(varlist)
        mean = average(varlist)
        return variance/mean

    def calculateSumOfDifferences(self, assemblage1, assemblage2):
        diff = 0
        for type in range(0, self.numberOfClasses):
            diff += abs(float(self.assemblageFrequencies[assemblage1][type]) - float(
                self.assemblageFrequencies[assemblage2][type]))
        return diff

    def calculateRootSumDiffSquares(self, assemblage1, assemblage2):
        diff = 0
        for type in range(0, self.numberOfClasses):
            diff += pow((float(self.assemblageFrequencies[assemblage1][type]) - float(
                self.assemblageFrequencies[assemblage2][type])),2)
        return pow(diff,0.5)


    def openXYFile(self):
        logger.debug("Opening XY file %s", self.args['xyfile'])
        ## open the xy file
        try:
            xyf = open(self.args['xyfile'], 'r')
        except csv.Error as e:
            sys.exit('file %s does not open: %s') % ( self.args['xyfile'], e)

        reader = csv.reader(xyf, delimiter='\t', quotechar='|')
        next(reader, None)
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

    def all_pairs(self, lst):
        return list((itertools.permutations(lst, 2)))


    def addOptions(self):
        self.args = {'debug': None, 'xyfile':None,
                 'xyfile': None, 'inputfile': None, 'outputdirectory': None }


    def parse_arguments(self):
        parser = argparse.ArgumentParser(description='Conduct an iterative deterministic seriation analysis')
        parser.add_argument('--debug', '-d', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
        parser.add_argument('--xyfile', default=None,
                            help="Enter the name of the XY file that contains the name of the assemblage and the X and Y coordinates for each.")
        parser.add_argument('--inputfile',
                            help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
        parser.add_argument('--graphs', default=0,
                            help="If true, the program will display the graphs that are created. If not, the graphs are just saved as .png files.")
        parser.add_argument('--outputdirectory', default=None,
                            help="If you want the output to go someplace other than the /output directory, specify that here.")
        self.args= parser.parse_args()
        try:
            self.args = vars(parser.parse_args())
        except IOError, msg:
            parser.error(str(msg))
            sys.exit()
        return self.args


    def analyze(self):

        self.openFile()
        self.openXYFile()

        self.preCalculateDistancesBetweenPairs()

        xmax=ymax=0
        xmin=ymin=100000000
        vmr_list=[]
        for assemblage in self.labels:
            if self.xAssemblage[ assemblage ] > xmax:
                xmax = self.xAssemblage[ assemblage ]

            if self.yAssemblage[ assemblage ] > ymax:
                ymax = self.yAssemblage[ assemblage ]

            if self.xAssemblage[ assemblage ] < xmin:
                xmin = self.xAssemblage[ assemblage ]

            if self.yAssemblage[ assemblage ] < ymin:
                ymin = self.yAssemblage[ assemblage ]
            #print "x-min: ", xmin
            #print "y-min: ", ymin

        #print "x-min: ", xmin, "x-max: ", xmax
        #print "y-min: ", ymin, "y-max: ", ymax

        for blocks in range(1,10):  ## this is the # of blocks to do


            quadrat_members=[]
            x_size = int((xmax - xmin)/blocks)
            for xscale in range(1,blocks):
                x_size = int((xmax - xmin)/xscale)
                print "x_size is: ", x_size
                xstart = xmin
                xend= xstart + x_size

                for yscale in range(1,blocks):

                    y_size = int((ymax - ymin)/yscale)
                    #print "y_size is: ", y_size
                    ystart = ymin
                    yend = ystart + y_size



                    #print "now on x: ", xstart, " to: ", xend, " y: ", ystart, " to: ", yend

                    for a in self.labels:
                        if self.xAssemblage[a] > xstart and self.xAssemblage[a] < xend and \
                                        self.yAssemblage[a] > ystart and self.yAssemblage[a] < yend:
                            #print "adding: ", a
                            quadrat_members.append(a)

                    ystart += y_size
                    yend += y_size



                xstart += x_size
                xend += x_size

            vmr = self.VMR(quadrat_members)
            self.vmr_for_blocks[ x_size ] = vmr

        return

    def plot(self):

        colors = list("rgbcmyk")
        x_vals=[]
        y_vals=[]
        for d in self.vmr_for_blocks:
            x = d
            y = self.vmr_for_blocks[d ]
            plt.scatter(x,y,color='b')

        plt.show()

if __name__ == "__main__":

    spatial_analysis = spatialAnalysis()
    spatial_analysis.parse_arguments()

    results = spatial_analysis.analyze()
    spatial_analysis.plot()
