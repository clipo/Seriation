#!/usr/bin/env python
# Copyright (c) 2013.  Carl P. Lipo <clipo@csulb.edu>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.
__author__ = 'carllipo'

# a standalone script for producing a minimum spanning tree solution based on the results of the seriation analysis.
# this requires GraphViz to be installed.

import shapefile
import csv
import sys
import networkx as nx
from pylab import *
import numpy
import matplotlib.pyplot as plt
import argparse

import time
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic
import random


''''
    usage:
        mst = MST()
        mst.createMST("../output/pfg.vna","output/",1)
'''''
class MST():
    def __init__(self,filename,outputdirectory,shapefile):
        self.filename=filename
        self.outputdirectory=outputdirectory
        self.shapefile=1
        self.shapefilename=""
        if shapefile is not None:
            if self.outputdirectory is not None:
                self.shapefilename = self.outputdirectory+self.filename[0:-4]+"-mst.shp"
            else:
                self.shapefilename = self.filename[0:-4]+"-mst.shp"
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")

    def usage(self):
        print "MST.createMST(filename,outputdirectory,shapefileFlag)"

    def iso(self,G1, glist):
        """Quick and dirty nonisomorphism checker used to check isomorphisms."""
        for G2 in glist:
            if isomorphic(G1,G2):
                return True
        return False

    def openFileHandles(self):
        ## Read in all of the data from the .vna file Reconstruct the graphs.
        try:
            file = open(self.filename)
        except IOError:
            print "can't open %s. I will now exit." % self.filename
            sys.exit()
        except TypeError:
            print "No file specified (use --file=<filename>). I will now exit."
            sys.exit()

        if self.outputdirectory is not None:
            outputFile = self.outputdirectory+self.filename[0:-4]+"-mst.vna"
        else:
            outputFile = self.filename[0:-4]+"-mst.vna"

        OUTPUTFILE = open(outputFile, 'w')
        return file,OUTPUTFILE

    def createMST(self):
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

        file,f=self.openFileHandles()

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
                weight = float(row[5])
                weightError = float(row[6])
                pair = node1 + "#" + node2
                pair2 = node2 + "#" + node1
                pairwise[ pair ]=weight
                pairwise[ pair2 ]=weight
                #print weight
                network = int(row[4])
                print "now on network: ", network, " edge: ", edgeCount, " oldnetwork: ", old_network
                pvalue = float(row[5])
                pError = float(row[6])
                pairwiseError[ pair ]=pError
                pairwiseError[ pair2 ] = pError
                meanDistance = float(row[7])
                if network > old_network:
                    #print "adding network..."
                    old_network = network
                    graphs.append(nx.Graph(ID=graphCount))
                    graphCount += 1
                    edgeCount = 0
                node1name = node1+"_"+str(network)
                node2name = node2+"_"+str(network)
                graphs[graphCount].add_node(node1name, label= node1, x = node1x, y = node1y,  size=node1Size )
                graphs[graphCount].add_node(node2name, label= node2, x = node2x, y = node2y,  size=node2Size )
                graphs[graphCount].add_edge(node1name, node2name, xy1=(node1x, node1y), xy2=(node2x, node2y),
                                            weight=weight,
                                            meanDistance=meanDistance,
                                            pvalue=weight,pError=pError,
                                            size=(nodeSize[node1],nodeSize[node2]))
                megaGraph.add_node(node1, x = node1x, y = node1y,  size=node1Size )
                megaGraph.add_node(node2, x = node2x, y = node2y, size=node2Size )
                megaGraph.add_edge(node1,node2, xy1=(node1x,node1y), xy2=(node2x,node2y),
                                       weight = weight,
                                       pvalueError = pError,
                                       meanDistance = meanDistance,
                                       pvalue = pvalue,
                                       pError = pError,
                                       color =network,
                                       size=(node1Size,node2Size),
                                        )
                edgeCount += 1
            count += 1

        if self.shapefile is not None:
            w = shapefile.Writer(shapefile.POLYLINE)  # 3= polylines
            c=0
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
                    w.poly(parts=[[[x1,y1],[x2,y2]]])
                c += 1
            w.save(self.shapefilename)

            time.sleep(15)

        plt.rcParams['text.usetex'] = False
        plt.figure(0,figsize=(8,8))
        mst=nx.minimum_spanning_tree(megaGraph,weight='weight')
        pos=nx.graphviz_layout(mst,prog="neato")
        #pos=nx.spring_layout(mst,iterations=500)

        # edge width is proportional number of games played
        edgewidth=[]
        weights = nx.get_edge_attributes(mst, 'weight')
        for w in weights:
            edgewidth.append(weights[w]*10)

        maxValue = np.max(edgewidth)
        widths=[]
        for w in edgewidth:
            widths.append(((maxValue-w)+1)*5)

        color = nx.get_edge_attributes(mst, 'color')
        colorList = []
        for c in color:
            colorList.append(color[c])
        colors=[]
        colorMax = max(colorList)
        for c in colorList:
            colors.append(c/colorMax)
        assemblageSizes=[]
        sizes = nx.get_node_attributes(mst, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])

        nx.draw_networkx_edges(mst,pos,alpha=0.3,width=widths, edge_color=colorList)
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
        pngfile=self.filename+"-mst.png"
        plt.savefig(pngfile,dpi=75)
        print(pngfile)


        f.write("*tie data\n*from to weight error distance\n")
        for (u,v,d) in mst.edges_iter(data=True):
            output = u +" "+ v + " "+str(d['weight'])+" "+str(d['pError'])+" "+str(d['meanDistance'])+"\n"
            #print output
            f.write(output)

        plt.figure(1,figsize=(30,20))
        # layout graphs with positions using graphviz neato

        UU=nx.Graph()
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist=[] # list of nonisomorphic graphs
        for G in graphs:
            # check against all nonisomorphic graphs so far
            if not self.iso(G, nlist):
                nlist.append(G)

        UU=nx.union_all(graphs) # union the nonisomorphic graphs
        #UU=nx.disjoint_union_all(nlist) # union the nonisomorphic graphs
        #pos=nx.spring_layout(UU,iterations=50)

        ##pos=nx.graphviz_layout(UU,prog="neato")
        pos=nx.graphviz_layout(UU,prog="twopi",root=0)
        ##labels=nx.draw_networkx_labels(UU,pos)
        # color nodes the same in each connected subgraph
        C=nx.connected_component_subgraphs(UU)
        for g in C:
            c = [random.random()] * nx.number_of_nodes(g) # random color...
            nx.draw(g,
                pos,
                node_size=40,
                node_color=c,
                vmin=0.0,
                vmax=1.0,
                alpha=.2,
                font_size=7,
            )
        atlasFile=self.outputdirectory+self.filename[0:-4]+"-atlas.png"
        plt.savefig(atlasFile,dpi=250)
        plt.show() # display
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NST')
    parser.add_argument('--test')
    parser.add_argument('--inputfile',default="../output/pfg.vna",help="The input .vna file.")
    parser.add_argument('--outputdirectory',default="../output/",help="If you want to specify an output directory do it here. The default is ../output/")
    parser.add_argument('--shapefile',default=None,help="flag for creating a shapefile.")
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()
    if args['test'] is not None:
        mst = MST("../output/pfg.vna","",None)
        out = mst.createMST()
    else:
        mst = MST(args['inputfile'],args['outputdirectory'],args['shapefile'])
        out = mst.createMST()

