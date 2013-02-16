__author__ = 'carllipo'

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
    print "minimumSpanningTree --file <filename.vna> --shapefile\n";

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hb:f:ds", ["help", "bootstrap=","file="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    filename = None
    inputFile = None
    verbose = False
    shapefileFlag = None
    for opt, arg in opts:
        if opt == "-v":
            verbose = True
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--file"):
            inputFile = arg
            shapefilename = inputFile[0:-4]
        elif opt in ("-s", "--shapefile"):
            shapefileFlag = arg
        else:
            assert False, "unhandled option"

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
    size=[]
    pairwise={}
    pairwiseError={}


    ## Read in all of the data from the .vna file Reconstruct the graphs.
    try:
        file = open(inputFile)
    except IOError:
        print "can't open ", inputFile, ". I will now exit."
        sys.exit()
    except TypeError:
        print "No file specified (use --file=<filename>). I will now exit."
        sys.exit()

    outputFile = inputFile[0:-4]+"-bootstrap.vna"
    f = open(outputFile, 'w')
    f.write("*node data\n")
    f.write("ID Size X Y Easting Northing\n")

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
            output = nodename+" " + str(nodeSize[nodename]) + str(nodeX[nodename])+ " " + str(nodeY[nodename])+ " " + str(nodeEasting[nodename]) + " " + str(nodeNorthing[nodename]) + "\n"
            f.write(output)

        if count > 1 and block == "ties":
            node1 = row[0]
            node2 = row[1]
            node1x,node1y = nodeEasting[node1], nodeNorthing[node1]
            node2x,node2y = nodeEasting[node2], nodeNorthing[node2]
            node1Size=nodeSize[node1]
            node2Size=nodeSize[node2]
            weight = 1-float(row[5])
            weightError = float(row[6])
            pair = node1 + "#" + node2
            pair2 = node2 + "#" + node1
            pairwise[ pair ]=weight
            pairwise[ pair2 ]=weight
            #print weight
            network = int(row[4])
            #print "now on network: ", network, " edge: ", edgeCount, " oldnetwork: ", old_network
            pvalue = float(row[5])
            pError = float(row[6])
            pairwiseError[ pair ]=pError
            pairwiseError[ pair2 ] = pError
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
                                        pvalue=weight,pError=pError,color='black',
                                        size=(nodeSize[node1],nodeSize[node2]))
            megaGraph.add_node(node1, x = node1x, y = node1y, name=node1, size=node1Size )
            megaGraph.add_node(node2, x = node2x, y = node2y, name=node2, size=node2Size )
            megaGraph.add_edge(node1,node2, xy1=(node1x,node1y), xy2=(node2x,node2y),
                                   weight = weight,
                                   pvalueError = pError,
                                   meanDistance = meanDistance,
                                   pvalue = pvalue,
                                   pError = pError,
                                   color ='black',
                                   size=(node1Size,node2Size),
                                    )
            edgeCount += 1
        count += 1

    if shapefileFlag is not None:
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
        w.save(shapefilename)

    mst=nx.minimum_spanning_tree(megaGraph,weight='weight')

    pos=nx.spring_layout(mst,iterations=200)

    # edge width is proportional number of games played
    edgewidth=[]
    weights = nx.get_edge_attributes(mst, 'weight')
    #print weights
    maxval = max(weights)

    for w in weights:
        #print weights[w]
        edgewidth.append(weights[w]*10)
    maxValue = max(edgewidth)
    #print "MAX: ", maxValue
    widths=[]
    for w in edgewidth:
        widths.append(((maxValue-w)+1)*5)

    assemblageSizes=[]
    sizes = nx.get_node_attributes(mst, 'size')
    #print sizes
    for s in sizes:
        #print sizes[s]
        assemblageSizes.append(sizes[s])
    nx.draw_networkx_edges(mst,pos,alpha=0.3,width=widths, edge_color='m')
    sizes = nx.get_node_attributes(mst,'size')
    nx.draw_networkx_nodes(mst,pos,node_size=assemblageSizes,node_color='w',alpha=0.4)
    nx.draw_networkx_edges(mst,pos,alpha=0.4,node_size=0,width=1,edge_color='k')
    nx.draw_networkx_labels(mst,pos,fontsize=10)
    font = {'fontname'   : 'Helvetica',
        'color'      : 'k',
        'fontweight' : 'bold',
        'fontsize'   : 14}
    edgelist = list(mst) # make a list of the edges
    #print edgelist
    #nx.draw(mst)
    #plt.savefig("path.png")
    plt.axis('off')
    pngfile=shapefilename+"-mst.png"
    plt.savefig(pngfile,dpi=75)
    print(pngfile)


    f.write("*tie data\n*from to weight error distance\n")
    for (u,v,d) in mst.edges_iter(data=True):
        output = u +" "+ v + " "+str(d['weight'])+" "+str(d['pError'])+" "+str(d['meanDistance'])+"\n"
        #print output
        f.write(output)
    plt.show() # display


if __name__ == "__main__":
    main(sys.argv[1:])
    print "done!"
