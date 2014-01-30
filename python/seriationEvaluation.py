__author__ = 'clipo'


import csv
import re
import sys
import logging as logger
import random
from datetime import datetime
import pickle
import multiprocessing
import numpy as np
import scipy as sp
import scipy.stats
import networkx as nx
from random import randint
from pylab import *

def filter_list(full_list, excludes):
    s = set(excludes)
    return (x for x in full_list if x not in s)

def worker(networks, out_q):
    """ The worker function, invoked in a process. 'nums' is a
        list of numbers to factor. The results are placed in
        a dictionary that's pushed to a queue.
    """
    outdict = {}
    for n in networks:
        outdict[n] = checkForValidAdditions(n)

    out_q.put(outdict)

def checkForValidAdditions(nnetwork):
    validComparisonsHash = {}
    filter_list = {}
    pairGraph = {}
    assemblages={}
    args={}
    args={}
    typeFrequencyUpperCI={}
    typeFrequencyLowerCI={}

    ## pickle the stuff I need for parallel processing
    validComparisonsHash=pickle.load(open('validComparisonsHash.p','rb'))
    pairGraph=pickle.load(open('validComparisonsHash.p','rb'))
    assemblages=pickle.load(open('assemblages.p','rb'))
    args=pickle.load(open('args.p','rb'))
    typeFrequencyUpperCI=pickle.load(open('typeFrequencyUpperCI.p','rb'))
    typeFrequencyLowerCI=pickle.load(open('typeFrequencyLowerCI.p','rb'))

    array_of_new_networks = []  ## a list of all the valid new networks that we run into
    maxnodes = len(nnetwork.nodes())

    for assEnd in ("End1", "End2"):
        if assEnd == "End1":
            otherEnd = "End2"
        else:
            otherEnd = "End1"

        endAssemblage = nnetwork.graph[assEnd]

        list1 = [validComparisonsHash[endAssemblage]]
        list2 = [nnetwork.nodes()]
        #validAssemblages = list(filter_list(list1, list2))
        validAssemblages=validComparisonsHash[endAssemblage]
        #print "valid assemblages: ", validAssemblages

        ######################################################################################
        for testAssemblage in validAssemblages:
            if assEnd == "End1":
                path = nx.shortest_path(nnetwork, nnetwork.graph["End1"], nnetwork.graph["End2"])
                innerNeighbor = path[1]
            elif assEnd == "End2":
                path = nx.shortest_path(nnetwork, nnetwork.graph["End2"], nnetwork.graph["End1"])
                innerNeighbor = path[1]
            else: ## sanity check
                sys.exit("Quitting due to errors.")
            ## Sanity check
            if innerNeighbor is None:
                sys.exit("Quitting due to errors.")
            c = pairGraph.get_edge_data(innerNeighbor, endAssemblage)
            comparison = c['weight']
            comparisonMap = ""
            oneToColumns = range(len(assemblages[testAssemblage]))

            error = 0  ## set the error check to 0
            for i in oneToColumns:
                c = ""
                p = nx.shortest_path(nnetwork, nnetwork.graph[assEnd], nnetwork.graph[otherEnd])
                newVal = assemblages[testAssemblage][i]
                previousAssemblage = testAssemblage
                for compareAssemblage in p:
                    oldVal = assemblages[compareAssemblage][i]
                    if args['bootstrapCI'] not in (None, ""):
                        upperCI_test = typeFrequencyUpperCI[previousAssemblage][i]
                        lowerCI_test = typeFrequencyLowerCI[previousAssemblage][i]
                        upperCI_end = typeFrequencyUpperCI[compareAssemblage][i]
                        lowerCI_end = typeFrequencyLowerCI[compareAssemblage][i]

                        if upperCI_test < lowerCI_end:
                            c += "D"
                        elif lowerCI_test > upperCI_end:
                            c += "U"
                        else:
                            c += "M"
                    else:
                        #print("Outer value: %f Inner value: %f"%(oldVal, newVal))
                        if newVal < oldVal:
                            c += "U"
                            c1 = "U"
                        elif newVal > oldVal:
                            c += "D"
                            c1 = "U"
                        elif newVal == oldVal:
                            c += "M"
                            c1 = "U"
                        else:
                            logger.debug("Error. Quitting.")
                            sys.exit("got null value in comparison of value for type %d in the comparison of %s", i,
                                     compareAssemblage)
                        logger.debug("Comparison %s is now %s", c1, c)
                        newVal = oldVal


                    previousAssemblage = compareAssemblage

                test = re.compile('DU|DM*U').search(c)
                if test not in (None, ""):
                    error += 1

            if error == 0:
                name = str(randint(1,1000000000000000))
                new_network = nnetwork.copy()
                new_network.graph["GraphID"] = name
                new_network.graph["name"] = name
                path = nx.shortest_path(nnetwork, nnetwork.graph["End1"], nnetwork.graph["End2"])

                new_network.add_node(testAssemblage, name=testAssemblage, end=1, site="end")
                new_network.add_node(endAssemblage, name=endAssemblage, site="middle", end=0)
                new_network.add_edge(testAssemblage, endAssemblage, weight=comparisonMap, end=1, site="end",
                                     GraphID=name)

                new_network.graph[assEnd] = testAssemblage
                path = nx.shortest_path(new_network, new_network.graph["End1"], new_network.graph["End2"])

                ## copy this solution to the new array of networks
                array_of_new_networks.append(new_network)

                if len(new_network) > maxnodes:
                    maxnodes = len(new_network)

    return array_of_new_networks



