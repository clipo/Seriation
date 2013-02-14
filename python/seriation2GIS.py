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
import shapefile

pp = pprint.PrettyPrinter(indent=4)
dumper = Dumper.Dumper(max_depth=5)

def usage():
    print "seriation2Network -file <filename.txt> -bootstrap <filename.txt>\n";

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hb:f:d", ["help", "bootstrap=","file="])
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
        elif opt in ("-b", "--bootstrap"):
            bootstrapFile=arg
        else:
            assert False, "unhandled option"

    try:
        file = open(bootstrapFile)
    except IOError:
        print "can't open ", bootstrapFile, ". I will now exit."
        sys.exit()

    bootstrapValues = []
    bootstrapError = []
    while 1:
        line = file.readline()
        if line in ['\n', '\r\n']:
            break
        row = line.split()
        bootstrapValues[row[0]][row[1]] = row[2]
        bootstrapValues[row[1]][row[0]] = row[2]
        bootstrapError[row[1]][row[0]] = row[3]
        bootstrapError[row[0]][row[1]] = row[3]

    #Set up blank lists for data
    x,y,id_no,date,target = [], [], [], [], []
    nodeX = {}
    nodeY = {}
    nodes = []
    graphs = []
    nodeSize = {}

    #read data from csv file and store in lists
    block = ""
    count = 0
    old_network = 0
    G = None
    graphCount = -1
    edgeCount = 0
    graphHash = {}
    row = ()
    ## Read in all of the data from the .vna file Reconstruct the graphs.
    try:
        file = open(inputFile)
    except IOError:
        print "can't open ", inputFile, ". I will now exit."
        sys.exit()

    megaGraph = nx.Graph()
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
            nodeX[nodename] = float(row[2])
            nodeY[nodename] = float(row[3])
            nodeSize[nodename] = float(row[1])
        if count > 1 and block == "ties":
            node1 = row[0]
            node2 = row[1]
            node1x,node1y = nodeX[node1],nodeY[node1]
            node2x,node2y = nodeX[node2],nodeY[node2]
            node1Size=nodeSize[node1]
            node2Size=nodeSize[node2]

            weight = float(row[3])
            network = int(row[4])
            #print "now on network: ", network, " edge: ", edgeCount, " oldnetwork: ", old_network
            pvalue = float(row[5])
            pError = float(row[6])
            meanDistance = float(row[7])
            if network > old_network:
                #print "adding network..."
                old_network = network
                graphs.append(nx.Graph())
                graphCount += 1
                edgeCount = 0
            graphs[graphCount].add_node(node1, x = node1x, y = node1y, name=node1, size=node1Size )
            graphs[graphCount].add_node(node2, x = node2x, y = node2y, name=node2, size=node2Size )
            graphs[graphCount].add_edge(node1,node2,
                                    xy1 = (node1x,node1y),
                                    xy2 = (node2x,node2y),
                                    weight = weight,
                                    meanDistance = meanDistance,
                                    pvalue = pvalue,pError=pError,color='black')
            megaGraph.add_edge(node1,node2, xy1=(node1x,node1y), xy2=(node2x,node2y),
                               weight = bootstrapValues[node1][node2],
                               pvalueError = bootstrapError[node1][node2],
                               meanDistance = meanDistance,
                               pvalue = pvalue,
                               pError = pError,
                               color ='black')
            edgeCount += 1
        count += 1

    w = shapefile.Writer(shapefile.POLYLINE)  # 3= polylines

    print count, " graphs "
    c=0
    #pp.pprint(graphs)
    for g in graphs:
        edges = g.edges()
        for e in edges:
            node1=e[0]
            node2=e[1]

            print g[node1]
            print "--------"
            #w.poly(parts=[[n1x,n1y],[n2x,n2y]], shapeType=shapefile.POLYLINE)
        c += 1

    mst=nx.minimum_spanning_edges(megaGraph,data=True)
    edgelist = list(mst) # make a list of the edges

    for edge in edgelist:
        print edge[0], " - ", edge[1]


if __name__ == "__main__":
    main(sys.argv[1:])

