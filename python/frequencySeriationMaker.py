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
        self.labels = []
        self.numberOfClasses = 0
        self.maxSeriationSize = 0
        self.validComparisonsHash = {}
        self.typeFrequencyLowerCI = {}
        self.typeFrequencyUpperCI = {}
        self.typeFrequencyMeanCI = {}
        self.typeNames = []
        self.canvas = pysvg.svg()
        self.rowIndex=30  # pixels between rows
        self.fontSize=4
        self.fontFamily="Verdana"
        self.typePositions=[]  ## center of each type column
        self.myStyle = StyleBuilder()
        self.myStyle.setFontFamily(fontfamily=self.fontFamily)
        self.myStyle.setFontSize(self.fontSize)
        self.maxAssemblageLabelLength=0

    def processSeriationData(self, args):
        try:
            logger.debug("trying to open: %s ", self.openFile)
            file = open(self.openFile, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", self.openFile, e)
            sys.exit('file %s does not open: %s') % ( self.openFile, e)


        reader = csv.reader(file, delimiter='\t', quotechar='|')
        rowcount = 0
        for row in reader:
            if rowcount>0:      ## ignore the first row here. we just need the assemblage info
                row = map(str, row)
                sernum=row[0]
                label = row[1]
                self.labels.append(label)
                labelLength=len(label)
                if labeLength>self.maxAssemblageLabelLength:
                    self.maxAssemblageLabelLength=labelLength
            rowcount += 1
        file.close()

        try:
            logger.debug("trying to open: %s ", self.openFile)
            file = open(self.openFile, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", self.openFile, e)
            sys.exit('file %s does not open: %s') % ( self.openFile, e)
        values = []
        count =0
        index=0
        for row in reader:
            ## first row is the header row with the type names.
            if count==0:
                row = map(str,row)
                row.pop(0)
                for r in row:
                    self.typeNames.append(r)
                count +=1
                self.outputHeaderRow(typeNames,args)
            else:
                if len(row) > 1:
                    index += 1
                    row = map(str, row)
                    sernum=row[0]
                    label = row[1]
                    self.labels[label] = label
                    row.pop(0)
                    row.pop(0)
                    row = map(float, row)
                    self.numberOfClasses = len(row)
                    self.countOfAssemblages += 1
                    self.outputAssemblageRow(label,row,rowPosition,args)
                    rowPosition += rowIndex
                else:
                    rowPosition += rowIndex

        self.maxSeriationSize = self.countOfAssemblages
        return True

    def outputHeaderRow(self, typeNames,args):

        maxTypeNameLength=0
        totalLength=0
        # find longest name in assemblages
        for a in self.typeNames:
            length=len(a)
            totalLength+=length
            if length>maxTypeNameLength:
                maxTypeNameLength=length

        colsize = 100
        colindex = maxTypeNameLength*20
        count=0
        xpos={}
        for t in self.typeNames:
            count+=1
            t1 = text(t, colindex, 100)
            t1.set_style(self.myStyle.getStyle())
            self.canvas.addElement(t1)
            colindex += colsize
            self.typePositions[count]=colindex

        rowPosition=150
        rowHeight = 10
        maxWidth = self.numberOfClasses * 20
        col_spacing = maxWidth/self.numberOfClasses

    def outputAssemblageRow(self, assemblageName, row,rowPosition, args):
        rowPosition = 150
        freq = []
        values=[]
        rowtotal = sum(row)
        for r in row:
            freq.append(float(float(r) / float(rowtotal)))
            values.append(float(r))

        t1 = text(assemblageName, 0, rowPosition)
        t1.set_style(self.myStyle.getStyle())
        self.canvas.addElement(t1)
        rowPosition += self.rowIndex
        count = 0
        xposition = self.typePositions[1]
        for typeFreq in freq:
            count += 1
            x = self.typePositions[count]
            width = int(typeFreq*200)
            xposition += self.rowIndex
            #print "width: ", width
            bar = self.sb.createRect(x-60-(width*0.5),  rowPosition-50, width,20, strokewidth=1)
            s.addElement(bar)



        s.save(self.outputFile)


    def setupOutput(self, args):
        self.openFile(args['inputfile'],args)
        self.outputFile= args['inputfile'][0:-4]+".svg"
        self.sb = ShapeBuilder()
        self.canvas = svg("Frequency Seriation")


    def makeGraph(self, args):
        self.setupOutput(args)
        self.processSeriationData(args)


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