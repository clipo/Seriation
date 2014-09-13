__author__ = 'clipo'

from frequencySeriationMaker import frequencySeriationMaker

seriation = frequencySeriationMaker()

args={'inputfile':'../testdata/pfg-cpl-seriations.txt','debug':1 }


seriation.makeGraph(args)
