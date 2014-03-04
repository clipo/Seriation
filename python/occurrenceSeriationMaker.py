__author__ = 'carllipo'

import pysvg
from pysvg import *
import logging as logger
import itertools
import math
import numpy as np
import xlsxwriter
import sys
import csv
import argparse
import numpy as np
import scipy as sp
import scipy.stats
import svgwrite
from svgwrite import cm, mm
import random
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

class occurrenceSeriationMaker():
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
        self.rowIndex=10  # pixels between rows
        self.fontSize=4
        self.fontFamily="Verdana"
        self.typePositions=[]  ## center of each type column
        self.maxAssemblageLabelLength=0
        self.dwg=None
        self.rowPosition=200
        self.rowHeight=15
        self.columnSize=25
        self.barScale=1000
        self.seriationNumber=0
        self.args={}
        self.outputRowAssemblage={}
        self.occurrenceValues={}
        self.simplifiedList = dict()
        self.occurrenceSeriationList={}
        self.newList={}
        self.columnStart=0

    def processSeriationData(self):
        try:
            logger.debug("trying to open: %s ", self.openFile)
            file = open(self.openFile, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", self.openFile, e)
            sys.exit('file %s does not open: %s') % ( self.openFile, e)

        reader = csv.reader(file, delimiter='\t', quotechar='|')
        rowcount = 0
        for row in reader:
            if rowcount==0:
                row = map(str,row)
                row.pop(0)
                self.numberOfClasses=len(row)
            else:      ## ignore the first row here. we just need the assemblage info
                row = map(str, row)
                if len(row)>0:
                    if self.args['multiple'] in (None, 'False', '0'):
                        sernum=1
                        label=row[0]
                    else:
                        sernum=row[0]
                        label = row[1]
                    self.labels.append(label)
                    labelLength=len(label)
                    if labelLength>self.maxAssemblageLabelLength:
                        self.maxAssemblageLabelLength=labelLength
            rowcount += 1
        file.close()
        self.rowPosition += self.maxAssemblageLabelLength
        ## initialize positions
        for n in range(0,self.numberOfClasses+1):
            self.typePositions.append(n)

        try:
            logger.debug("trying to open: %s ", self.openFile)
            file = open(self.openFile, 'r')
        except csv.Error as e:
            logger.error("Cannot open %s. Error: %s", self.openFile, e)
            sys.exit('file %s does not open: %s') % ( self.openFile, e)
        values = []
        count =0
        index=2
        reader = csv.reader(file, delimiter='\t', quotechar='|')
        self.rowPosition = 200
        startBlock=0
        startcorner=self.rowPosition
        for row in reader:
            ## first row is the header row with the type names.
            if count==0:
                row = map(str,row)
                if self.args['multiple'] in (None, 'False', "0"):
                    row.pop(0)
                else:
                    row.pop(0)
                    row.pop(0)
                for r in row:
                    self.typeNames.append(r)
                count +=1
                self.headerRow=(self.typeNames)
                self.headerRowPosition=self.rowPosition
            else:
                if len(row) > 1:
                    index += 1
                    row = map(str, row)
                    if self.args['multiple'] in (None, 'False', '0'):
                        self.seriationNumber=1
                        label=row[0]
                        self.labels.append(label)
                        row.pop(0)
                    else:
                        self.seriationNumber=float(row[0])
                        label = row[1]
                        self.labels.append(label)
                        row.pop(0)
                        row.pop(0)
                    #print row
                    row = map(float, row)
                    newrow=[]
                    rowtext=""
                    for r in row:
                        if r>0:
                            r = 1.0
                            rowtext += "1"
                        else:
                            r = 0.0
                            rowtext += "0"
                        newrow.append(r)
                    self.numberOfClasses = len(newrow)
                    self.countOfAssemblages += 1
                    self.outputRowAssemblage[label]=newrow
                    self.occurrenceValues[label]=rowtext
                else:
                    endcorner=self.rowPosition+self.rowHeight
                    self.rowPosition += self.rowIndex*2
                    if self.seriationNumber % 2 <> 0:
                        self.createBlock(startcorner,endcorner)
                    startcorner=self.rowPosition

        if self.args['multiple'] in (None, 'False', "0"):
            ## eliminate duplicates (if single)
            #print "length before: ",len(self.outputRowAssemblage)
            self.aggregateIdenticalAssemblages()
            for e in self.outputRowAssemblage:
                if len(e) > self.maxAssemblageLabelLength:
                    self.maxAssemblageLabelLength=len(e)
            #print "length after: ",len(self.outputRowAssemblage)

        self.outputHeaderRow(self.typeNames)
        self.rowPosition += self.rowIndex
        startcorner=self.rowPosition
        for r in self.outputRowAssemblage:
            self.outputAssemblageRow(r,self.outputRowAssemblage[r])
            self.rowPosition += self.rowIndex

        self.maxSeriationSize = len(self.labels)

        #if args['pdf'] not in (None, False, 0):
        #    drawing = svg2rlg(self.outputFile)
        #    renderPDF.drawToFile(drawing, self.outputFile[0:-4]+".pdf")

        return True

    def calculateSumOfDifferences(self, assemblage1, assemblage2):
        diff = 0
        for type in range(0, self.numberOfClasses):
            diff += pow((float(self.outputRowAssemblage[assemblage1][type]) - float(
                self.outputRowAssemblage[assemblage2][type])),2)
        return pow(diff,0.5)

    def createBlock(self,startcorner,endcorner):
        shapes = self.dwg.add(self.dwg.g(id='background', fill='grey', opacity=0.3))
        shapes.add(self.dwg.rect(insert=(20,  startcorner), size=(self.numberOfClasses*self.columnSize+self.maxAssemblageLabelLength*8+280,
                                                                  endcorner-startcorner),
                        fill='grey',opacity=0.3, stroke='none', stroke_width=1))

    def aggregateIdenticalAssemblages(self):
        tempList={}
        for e in self.occurrenceValues:
            combined = ''.join(self.occurrenceValues[e])
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

        self.occurrenceSeriationList={}
        self.outputRowAssemblage={}
        self.labels=[]
        for t in tempList:
            newval = str("/".join(tempList[t]))
            self.occurrenceSeriationList[newval]=t
            #print "NewKey:", newval, "Value:", t
            self.outputRowAssemblage[newval]=t
            self.labels.append(newval)

        return True

    def outputHeaderRow(self, typeNames):

        maxTypeNameLength=0
        totalLength=0
        # find longest name in assemblages
        for a in self.typeNames:
            length=len(a)
            totalLength+=length
            if length>maxTypeNameLength:
                maxTypeNameLength=length

        colindex = self.maxAssemblageLabelLength*10+self.columnSize+self.columnStart
        count=0
        xpos={}
        for t in self.typeNames:
            count+=1
            #self.dwg.add(self.dwg.textArea(text=t,insert=(colindex, self.rowPosition)))
            rotationtext = "rotate(-90,"+str(colindex)+","+str(self.rowPosition)+")"
            self.dwg.add(self.dwg.text(t, insert=(colindex, self.rowPosition),transform=str(rotationtext)  ))
            colindex += self.columnSize
            self.typePositions[count]=colindex

        #self.dwg.add(self.dwg.text("Assemblage Size", insert=(colindex+50,self.rowPosition)))
        #self.assemblageSize=colindex+50

    def outputAssemblageRow(self, assemblageName, row):
        freq=[]
        values=[]
        #rowtotal = sum(row)
        for r in row:
            ## turn into occurrence
            if int(r)>0:
                r=1.0
            else:
                r=0.0
            freq.append(float(r))
            values.append(float(r))
        self.rowPosition += self.rowIndex
        assemblageNameXPosition=self.maxAssemblageLabelLength*10-len(assemblageName)*7.5+self.columnStart
        self.dwg.add(self.dwg.text(assemblageName, insert=(assemblageNameXPosition, self.rowPosition+self.rowHeight)))
        count = 0
        xposition = self.typePositions[1]
        lowerCI,upperCI,meanCI=self.bootstrapCICalculation(freq, sum(values))
        for typeFreq in freq:
            count += 1
            x = self.typePositions[count]
            if typeFreq>0:
                width=10
            else:
                width=0
            xposition += self.rowIndex
            #print "width: ", width
            shapes = self.dwg.add(self.dwg.g(id='freqbar', fill='white'))
            leftx = x-90-(width*0.5)+60
            shapes.add(self.dwg.rect(insert=(leftx,  self.rowPosition), size=(width,self.rowHeight),
                        fill='black', stroke='black', stroke_width=1))
            #self.errorBars(typeFreq,width,x,leftx,lowerCI[count-1],upperCI[count-1],meanCI[count-1])
        #self.dwg.add(self.dwg.text(int(sum(values)), insert=(self.assemblageSize+25,self.rowPosition+5)))
        self.dwg.save()



    def errorBars(self,freq,width,x,original_left,lowerCI,upperCI,meanCI):

        newWidth=meanCI+(meanCI-lowerCI)+(upperCI-meanCI)
        newWidthSize=int(newWidth*self.barScale)
        leftx = x-90-(newWidthSize*0.5)
        #print "freq: ",freq, "lowerCI:", lowerCI, "upperCI", upperCI
        errorBarHeight=self.rowHeight/3
        shapes = self.dwg.add(self.dwg.g(id='errorbar', fill='white'))
        ## left side
        shapes.add(self.dwg.rect(insert=(leftx,self.rowPosition+(0.3*self.rowHeight)), size=((original_left-leftx),errorBarHeight),
                        fill="black",stroke="black", stroke_width=0.5))
        ## right side
        rightx = original_left+width
        rightCIend= rightx + newWidthSize*0.5
        shapes.add(self.dwg.rect(insert=(rightx,self.rowPosition+(0.3*self.rowHeight)), size=((original_left-leftx),errorBarHeight),
                        fill="black",stroke="black", stroke_width=0.5))

        #shapes.add(self.dwg.rect(insert=(leftx,self.rowPosition+(0.5*self.rowHeight)), size=(newWidthSize,errorBarHeight),
                        #fill="black",stroke="black", stroke_width=0.5))

        ## add meanCI bar
        meanWidthSize=int(meanCI*self.barScale)
        leftx=x-90-(meanWidthSize*0.5)
        #shapes.add(self.dwg.rect(insert=(leftx,self.rowPosition+(0.4*self.rowHeight)),size=(meanWidthSize,errorBarHeight*3),
        #                     fill="blue",stroke="blue", stroke_width=0.5))


   ########################################### BOOTSTRAP CI SECTION ####################################
    def bootstrapCICalculation(self, freqs, currentAssemblageSize, bootsize=100, confidenceInterval=0.999):
        types = len(freqs)
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
            for typeFrequency in freqs:
                index += typeFrequency
                cumulate[count] = index  ## this is the index of the frequency for this class
                ## index should be total # of types at end
                count += 1

            ## set new_assemblage
            new_assemblage = []
            for c in range(0, types):
                new_assemblage.append(0.0)

            for sherd in range(0, int(currentAssemblageSize)):
                rand = random.random()             ## random number from 0-1
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
            boot_assem_size=sum(new_assemblage)
            for g in new_assemblage:
                new_assemblage_freq.append(float(g / float(boot_assem_size)))
                arrayOfStats[counter].append(float(g / float(boot_assem_size)))
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

        #self.typeFrequencyLowerCI[currentLabel] = lowerCI
        #self.typeFrequencyUpperCI[currentLabel] = upperCI
        #self.typeFrequencyMeanCI[currentLabel] = meanCI

        return lowerCI,upperCI,meanCI

    def confidence_interval(self, data, confidence=0.05):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * sp.stats.t._ppf((1 + confidence) / 2., n - 1)
        return m, m - h, m + h

    def setupOutput(self):
        self.openFile=self.args['inputfile']
        self.outputFile= self.args['inputfile'][0:-4]+".svg"
        self.dwg = svgwrite.Drawing(self.outputFile, profile='tiny')

    def makeGraph(self,args):
        self.addOptions(args)
        self.setupOutput()
        self.processSeriationData()

    def addOptions(self, oldargs):
        self.args = {'debug': None, 'bootstrapCI': None, 'bootstrapSignificance': None,
                'inputfile': None, 'outputdirectory': None,
                'multiple': None}
        for a in oldargs:
            self.args[a] = oldargs[a]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='create seriation graph')
    parser.add_argument('--pdf',default=None,help='[UNIMPLEMENTED] Create PDF files in addition to SVG')
    parser.add_argument('--debug', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
    parser.add_argument('--bootstrapCI', default=None,
                        help="Sets whether you want to use the bootstrap confidence intervals for the comparisons between assemblage type frequencies. Set's to on or off.")
    parser.add_argument('--bootstrapSignificance', default=0.95, type=float,
                        help="The significance to which the confidence intervals are calculated. Default is 0.95.")
    parser.add_argument('--inputfile',
                        help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
    parser.add_argument('--outputdirectory', default=None,
                        help="If you want the output to go someplace other than the /output directory, specify that here.")
    parser.add_argument('--multiple', default=True, help="If you are creating just one seriation set (with the first column being the assemblage name), set this to FALSE. Otherwise, set the seriation number as the first column and set to TRUE. Default is TRUE. ")
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()

    seriation = occurrenceSeriationMaker()
    seriation.makeGraph(args)

''''
From the command line:

python ./occurrenceSeriationMaker.py --inputfile=../testdata/pfg.txt"

if you want to use this based on the original data that doesn't have a first column with the seriation number (as is true for the IDSS output), you
need to specify that its a single input.

python ./occurrenceSeriationMaker.py --multiple=0 --inputfile=../testdata/testdata-two-branches.txt

As a module:

from occurrenceySeriationMaker import occurrenceSeriationMaker

seriation = occurrenceSeriationMaker()

args={'inputfile':'../testdata/pfg-cpl-seriations.txt','debug':1 }

seriation.makeGraph(args)




'''''