__author__ = 'clipo'

import logging as logger
import itertools
import math
import sys
import csv
import argparse
import numpy as np
import scipy as scipy
import scipy.stats
import random as rnd
import scikits.bootstrap as bootstrap
import numpy as np
#from sklearn.gaussian_process import GaussianProcess
from matplotlib import pyplot as plt
from scipy import stats
import operator

class sampleSizeEvaluation():

    def __init__(self):
        self.cell_text =[]
        self.data={}
        self.upperCI={}
        self.lowerCI={}
        self.mean={}
        self.sampleSizeArray=[]
        self.lowerCIs=[]
        self.upperCIs=[]
        self.means=[]
        self.numberOfClasses=0
        self.typeNames=[]
        self.labels=[]
        self.assemblageTypeCounts={}
        self.args={}
        self.fignum=0
        self.assemblageSampleSize=0
        self.slope=self.intercept= self.r_value=self.p_value=self.std_err=0.0

    def confidence_interval(self, data, alpha=0.05):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * scipy.stats.t._ppf((1 + alpha) / 2., n - 1)
        return m, m - h, m + h

       ######################################
       ## types = array of frequencies
       ## assemblageLabel = name of assemblage
       ########################################### BOOTSTRAP CI SECTION ####################################
    def bootstrapRichness(self, assemblageLabel, types, bootsize=3, alpha=0.05):


        ## original sample size
        currentAssemblageSize = sum(types)
        self.assemblageSampleSize=currentAssemblageSize
        #print "assemblage size: ", currentAssemblageSize
        numberOfTypes = len(types)
        typeFrequencies=[]
        t=0.0
        ## now count through the classes and set up the frequencies
        for type in types:
            t = float(float(type)/float(self.assemblageSampleSize))
            #print "type: ", type, "t: ", t
            typeFrequencies.append(t)
        #print "type freqs: ", typeFrequencies

        ## iterate through bootstrap (in 50s through 2x of the original sample)
        maxBootstrapSampleSize = 2 * self.assemblageSampleSize
        step=maxBootstrapSampleSize/20
        for currentSampleSize in range(10,maxBootstrapSampleSize,step):
            richness=[]
            ## size of bootstrapping (how many assemblages to create)
            for counter in range(0, bootsize):
                # clear and set the array
                cumulate = []
                for d in range(0, numberOfTypes):
                    cumulate.append(0.0)

                index = 0
                count = 0
                ## now count through the classes and set up the frequencies
                for typeFrequency in typeFrequencies:
                    index += typeFrequency
                    cumulate[count] = index  ## this is the index of the frequency for this class
                    ## index should be total # of types at end
                    count += 1

                ## set new_assemblage
                new_assemblage = []
                for c in range(0, self.numberOfClasses):
                    new_assemblage.append(0.0)

                for sherd in range(0, int(currentSampleSize)):
                    rand = rnd.random()             ## random number from 0-1
                    classVar = 0
                    found = 0
                    typeIndex = self.numberOfClasses - 1
                    total = sum(cumulate)
                    for t in reversed(cumulate):
                        if rand <= t:
                            found = typeIndex
                        typeIndex -= 1
                    new_assemblage[found] += 1
                #print "New assemblage: ", new_assemblage
                r = sum(x > 0 for x in new_assemblage)
                #print "Richness: ",r
                richness.append(r)

                # compute 95% confidence intervals around the mean
            #CIs = bootstrap.ci(data=richness, statfunction=scipy.mean, alpha=alpha)
            #mean = scipy.stats.mean(richness, axis=0)

            meanValue,lowerCIValue,upperCIValue=self.confidence_interval(richness,alpha)

            self.sampleSizeArray.append(currentSampleSize)
            self.means.append(meanValue)
            self.upperCIs.append(upperCIValue)
            self.lowerCIs.append(lowerCIValue)
            self.data[currentSampleSize]=richness
            self.lowerCI[currentSampleSize]=lowerCIValue
            self.upperCI[currentSampleSize]=upperCIValue
            self.mean[currentSampleSize]=meanValue
            if currentSampleSize-step+1>self.assemblageSampleSize<currentSampleSize+step-1:
                arrayPosition=len(self.sampleSizeArray)-1

        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = stats.linregress(self.sampleSizeArray[arrayPosition-5:arrayPosition],self.means[arrayPosition-5:arrayPosition])


        #print "SLOPE: ", self.slope, "std_err:", self.std_err, " p-value: ", self.p_value
        return True

    def plotData(self,assemblageName,row,col):
        font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 14,
        }
        self.fignum += 1
        plt.subplot(row, col, self.fignum)
        lowerDiffs=map(operator.sub, self.means,self.lowerCIs)
        #lowerDiffs=[self.means-self.lowerCIs for self.means, self.lowerCIs in zip(self.means, self.lowerCIs)]
        upperDiffs=map(operator.sub,self.upperCIs,self.means)
        #upperDiffs=[self.upperCIs-self.means for self.upperCIs, self.means in zip(self.upperCIs, self.means)]
        CIIntervals=[lowerDiffs,upperDiffs]
        #print CIIntervals
        plt.errorbar(self.sampleSizeArray,self.means,yerr=CIIntervals, fmt='o',label=assemblageName)

        #yerr=[self.lowerCIs,self.upperCIs]
        plt.plot(self.sampleSizeArray,self.means, 'k--')
        text="%s\t%f\t%f\t%f" % (assemblageName, self.slope, self.std_err, self.p_value)
        print text
        self.cell_text.append([assemblageName,self.slope,self.std_err,self.p_value])
        #plt.text(self.assemblageSampleSize*.5, self.numberOfClasses-3,slopetext, fontdict=font)
        plt.axvline(self.assemblageSampleSize, color='r', linestyle='dashed', linewidth=2)
        plt.title(assemblageName, fontdict=font)
        plt.xlabel('Sample Size (N)', fontdict=font)
        plt.ylabel('Bootstrap Richness', fontdict=font)

    def openFile(self,filename):
        try:
            file = open(filename, 'r')
        except csv.Error as e:
            sys.exit('file %s does not open: %s') % ( filename, e)

        reader = csv.reader(file, delimiter='\t', quotechar='|')
        rowcount = 0
        for row in reader:
            if rowcount==0:
                row = map(str,row)
                row.pop(0)
                self.numberOfClasses=len(row)
                self.typeNames=row
                rowcount += 1
            else:
                row = map(str, row)
                if len(row)>0:
                    assemblage = row[0]
                    self.labels.append(assemblage)
                    row.pop(0)
                    row = map(int, row)
                    self.assemblageTypeCounts[assemblage]=row

        file.close()

    def addOptions(self, oldargs):
        self.args = {'debug': None, 'alpha': None,
                'inputfile': None, 'outputdirectory': None }
        for a in oldargs:
            self.args[a] = oldargs[a]

    def conductAnalysis(self,args):
        self.addOptions(args)
        self.checkMinimumRequirements()
        self.openFile(self.args['inputfile'])
        #fig = plt.figure()
        print "Assemblage\tSlope\t Std_Dev\tp-value"
        self.cell_text.append(["Assemblage","Slope","Std.Dev.","p-value"])
        rows=int(len(self.labels)/3)+1
        cols=3
        for ass in self.labels:
            self.sampleSizeArray =[]
            self.means=[]
            self.upperCIs=[]
            self.lowerCIs=[]
            self.data={}
            #self.lowerCI=[]
            #self.upperCI=[]
            #self.mean=[]
            self.bootstrapRichness( ass, self.assemblageTypeCounts[ass], bootsize=1000, alpha=self.args['alpha'])

            self.plotData(ass,rows,cols)
        #the_table = plt.table(cellText=self.cell_text')
        plt.show()

    def checkMinimumRequirements(self):
        if self.args['inputfile'] in (None, ""):
            sys.exit("Inputfile is a required input value: --inputfile=../testdata/testdata.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='bootstrap analysis')
    parser.add_argument('--debug', default=None, help='Sets the DEBUG flag for massive amounts of annotated output.')
    parser.add_argument('--alpha', default=0.05, type=float,
                        help="The alpha significance to which the confidence intervals are calculated. Default is 0.05.")
    parser.add_argument('--inputfile',
                        help="<REQUIRED> Enter the name of the data file with the assemblage data to process.")
    parser.add_argument('--outputdirectory', default=None,
                        help="If you want the output to go someplace other than the /output directory, specify that here.")
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()


    sampleSizeEval = sampleSizeEvaluation()
    sampleSizeEval.conductAnalysis(args)

''''
From the command line:

python ./frequencySeriationMaker.py --inputfile=../testdata/pfg.txt"

if you want to use this based on the original data that doesn't have a first column with the seriation number (as is true for the IDSS output), you
need to specify that its a single input.

python ./frequencySeriationMaker.py --multiple=0 --inputfile=../testdata/testdata-two-branches.txt

As a module:

from frequencySeriationMaker import frequencySeriationMaker

seriation = frequencySeriationMaker()

args={'inputfile':'../testdata/pfg-cpl-seriations.txt','debug':1 }

seriation.makeGraph(args)




'''''