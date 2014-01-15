__author__ = 'carllipo'

import pysvg
import logging as logger
import itertools
import math
import numpy as np
import xlsxwriter
import sys
import csv
import argparse
from pysvg.structure import *
from pysvg.core import *
from pysvg.text import *
from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *


class frequencySeriationMaker():
    color = ["b", "r", "m", "y", "k", "w", (0.976, 0.333, 0.518), (0.643, 0.416, 0.894),
             (0.863, 0.66, 0.447), (0.824, 0.412, 0.118)]

    def __init__(self):
        self.inputfile = ""
        self.outputDirectory = ""
        self.assemblageSize = {}
        self.assemblageFrequencies = {}
        self.assemblages = {}
        self.countOfAssemblages = 0
        self.assemblageValues = {}
        self.labels = {}
        self.numberOfClasses = 0
        self.maxSeriationSize = 0
        self.validComparisonsHash = {}
        self.typeFrequencyLowerCI = {}
        self.typeFrequencyUpperCI = {}
        self.typeFrequencyMeanCI = {}
        self.typeNames = []

    def openFile(self, filename, args):
        try:
            logger.debug("trying to open: %s ", filename)
            file = open(filename, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", filename, e)
            sys.exit('file %s does not open: %s') % ( filename, e)

        reader = csv.reader(file, delimiter='\t', quotechar='|')
        values = []
        count =0
        index=0
        for row in reader:
            if count==0:
                row = map(str,row)
                row.pop(0)
                for r in row:
                    self.typeNames.append(r)
                count +=1
            else:
                if len(row) > 1:
                    index += 1
                    row = map(str, row)
                    label = row[0]
                    self.labels[label+"_"+str(index)] = label
                    row.pop(0)
                    row = map(float, row)
                    self.numberOfClasses = len(row)
                    freq = []
                    rowtotal = sum(row)
                    for r in row:
                        freq.append(float(float(r) / float(rowtotal)))
                        values.append(float(r))
                    self.assemblages[label+"_"+str(index)] = freq
                    self.assemblageFrequencies[label+"_"+str(index)] = freq
                    self.assemblageValues[label+"_"+str(index)] = values
                    self.assemblageSize[label+"_"+str(index)] = rowtotal
                    self.countOfAssemblages += 1
        self.maxSeriationSize = self.countOfAssemblages
        return True

    def makeGraph(self, args):
        self.openFile(args['inputfile'],args)
        outputFile= args['inputfile'][0:-4]+".svg"
        oh = ShapeBuilder()
        s = svg("seriation")

        myStyle = StyleBuilder()
        myStyle.setFontFamily(fontfamily="Verdana")
        myStyle.setFontSize('4')
        maxAssemblageLength=0
        # find longest name in assemblages
        for a in self.assemblages:
            length=len(a)
            if length>maxAssemblageLength:
                maxAssemblageLength=length

        colsize = 100
        colindex = maxAssemblageLength*20
        count=0
        xpos={}
        for t in self.typeNames:
            count+=1
            t1 = text(t, colindex, 100)
            t1.set_style(myStyle.getStyle())
            s.addElement(t1)
            colindex += colsize
            xpos[count]=colindex
        rowIndex=30
        rowPosition=150
        rowHeight = 10
        maxWidth = self.numberOfClasses * 20
        col_spacing = maxWidth/self.numberOfClasses
        sb =ShapeBuilder()

        for assemblage in self.labels:
            t1 = text(assemblage, 0, rowPosition)
            t1.set_style(myStyle.getStyle())
            s.addElement(t1)
            rowPosition += rowIndex
            count = 0
            xposition = xpos[1]
            for typeFreq in self.assemblageFrequencies[assemblage]:
                count += 1
                x = xpos[count]
                width = int(typeFreq*200)
                xposition += rowIndex
                #print "width: ", width
                bar = sb.createRect(x-60-(width*0.5),  rowPosition-50, width,20, strokewidth=1)
                s.addElement(bar)



        s.save(outputFile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='create seriation graph')
    parser.add_argument('--debug', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
    parser.add_argument('--bootstrapCI', default=None,
                        help="Sets whether you want to use the bootstrap confidence intervals for the comparisons between assemblage type frequencies. Set's to on or off.")
    parser.add_argument('--bootstrapSignificance', default=0.95, type=float,
                        help="The significance to which the confidence intervals are calculated. Default is 0.95.")
    parser.add_argument('--inputfile',
                        help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
    parser.add_argument('--outputdirectory', default=None,
                        help="If you want the output to go someplace other than the /output directory, specify that here.")
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    seriation = frequencySeriationMaker()

    seriation.makeGraph(args)

''''
From the command line:

python ./frequencySeriationMaker.py --inputfile=../testdata/pfg.txt"


As a module:

from frequencySeriationMaker import frequencySeriationMaker

seriation = frequencySeriationMaker()

args={}
args{'inputfile'}="../testdata/testdata-5.txt"
args{'debug'}=1


seriation.makeGraph(args)

'''''