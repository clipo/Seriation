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
    print "seriation2Network --file <filename.txt> --bootstrap <filename.txt> --gps\n";

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hb:f:dg", ["help", "bootstrap=","file="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    filename = None
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
            gpsfilename = inputFile[0:-4]
        elif opt in ("-b", "--bootstrap"):
            bootstrapFile = arg
        elif opt in ("-g", "--gps"):
            gpsFlag = arg
        else:
            assert False, "unhandled option"

    if bootstrapFile is not None:
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
    nodeEasting = {}
    nodeNorthing = {}
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
            nodeEasting[nodename] = float(row[4])
            nodeNorthing[nodename] = float(row[5])
            nodeSize[nodename] = float(row[1])
        if count > 1 and block == "ties":
            node1 = row[0]
            node2 = row[1]
            node1x,node1y = nodeEasting[node1], nodeNorthing[node1]
            node2x,node2y = nodeEasting[node2], nodeNorthing[node2]
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
            graphs[graphCount].add_edge(node1,node2, xy1=(node1x, node1y), xy2=(node2x, node2y),
                                                weight=weight,
                                                meanDistance=meanDistance,
                                                pvalue=pvalue,pError=pError,color='black')
            if bootstrapFile is not None:
                megaGraph.add_edge(node1,node2, xy1=(node1x,node1y), xy2=(node2x,node2y),
                               weight = bootstrapValues[node1][node2],
                               pvalueError = bootstrapError[node1][node2],
                               meanDistance = meanDistance,
                               pvalue = pvalue,
                               pError = pError,
                               color ='black')
            edgeCount += 1
        count += 1

    if gpsFlag is not None:
        w = shapefile.Writer(shapefile.POLYLINE)  # 3= polylines
        #print count, " graphs "
        c=0
        #pp.pprint(graphs)
        for g in graphs:
            edges = g.edges()
            for e in edges:
                node1 = e[0]
                node2 = e[1]
                print g[node1][node2]
                x1 = g[node1][node2]['xy1'][0]
                y1 = g[node1][node2]['xy1'][1]
                x2 = g[node2][node1]['xy2'][0]
                y2 = g[node2][node1]['xy2'][1]
                #print x1, "-", y1
                #print x2, "-", y2
                w.poly(parts=[[[x1,y1],[x2,y2]]])
            c += 1
        w.save(gpsfilename)

    if bootstrapFile is not None:
        outputFile = inputFile[0,-4]+"-bootstrap.vna"
        f = open(outputFile, 'w')
        f.write("*node data\n")
        f.write("ID X Y")
        vertices = megaGraph.vertices()
        for v in vertices:
            output = v[0]+" "+nodeX[v[0]]+" "+nodeY[v[0]]+"\n"
            f.write(output)
        mst=nx.minimum_spanning_edges(megaGraph,data=True)
        edgelist = list(mst) # make a list of the edges
        f.write("*tie data\n*from to weight error\n")
        for edge in edgelist:
            output = edge[0]+" "+edge[1]+" "+edge[0][1]['weight']+" "+edge[0][1]['pvalueError']+"\n"
            print edge[0], " - ", edge[1], " ", edge[0][1]['weight']," ", edge[0][1]['pvalueError'],"\n"
            f.write(output)

if __name__ == "__main__":
    main(sys.argv[1:])
    print "done!"
