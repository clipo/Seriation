__author__ = 'carllipo'
import shapefile
import csv
import getopt
import sys
import networkx as nx
import pprint
import Dumper
from pylab import *
import matplotlib.pyplot as plt
import time

pp = pprint.PrettyPrinter(indent=4)
dumper = Dumper.Dumper(max_depth=5)

def usage():
    print "seriation2GIS -file <filename.csv>\n";

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hf:d", ["help", "file="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    inputFile = None
    verbose = False
    for opt, arg in opts:
        if opt == "-v":
            verbose = True
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--file"):
            inputFile = arg
        else:
            assert False, "unhandled option"

    #Set up blank lists for data
    x,y,id_no,date,target=[],[],[],[],[]
    nodeX={}
    nodeY={}
    nodes=[]
    graphs=[]

    #read data from csv file and store in lists
    block = ""
    count = 0
    old_network = 0
    G = None
    graphCount = 0
    edgeCount = 0
    graphHash = {}
    row=()
    file=open(inputFile)

    while 1:
        line=file.readline()
        #print line
        if not line:
            break
        if "*Node data" in line:
            block = "nodes"
            count = 0
        elif "*Node properties" in line:
            block = "properties"
            count= 0
        elif "*Tie data" in line:
            block = "ties"
            count = 0
        if line in ['\n', '\r\n']:
            break
        row =line.split()
        if row is None:
            break
        if count > 1 and block == "nodes":
            nodename = row[0]
            nodes.append(row[0])
            nodeX[nodename] = row[2]
            nodeY[nodename] = row[3]
        if count > 1 and block == "ties":
            node1 = row[0]
            node2 = row[1]
            weight = row[3]
            network = row[4]
            print "now on network: ", network, " edge: ", edgeCount
            pvalue = row[5]
            pError = row[6]
            meanDistance = row[7]
            if network > old_network:
                old_network = network
                if G is not None:
                    graphs.append(G)
                graphHash[graphCount] = G
                G = nx.Graph()
                graphCount += 1
                edgeCount = 0
            G.add_edge(node1, node2, weight = weight)
            edgeCount += 1
        count += 1


    print count, " graphs "
    c=0
    pp.pprint(graphs)
    while c < count:

        print graphs[ c ].number_of_edges()
        c += 1

if __name__ == "__main__":
    main(sys.argv[1:])

