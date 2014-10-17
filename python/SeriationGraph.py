#!/usr/bin/env python
# Copyright (c) 2013.  Carl P. Lipo <clipo@csulb.edu>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.
__author__ = 'carllipo'

# An graph/network abstraction for seriation graphs. Using this will allow us to potentially swap out graph libraries
# as there are a number of options -- networkX, graph-tool, igraph with vastly different performances


#import networkx as nx
import numpy
import uuid
from graph_tool.all import *


class SeriationGraph():

    def __index__(self):
        self.ID = uuid.uuid4().urn
        self.graph = gt.Graph(directed=False)
        self.nodeHash = {}
        self.edgeHash = {}
        nameProp = self.graph.new_graph_property("string")
        end1Prop = self.graph.new_graph_property("string")
        end2Prop = self.graph.new_graph_property("string")
        graphIDProp = self.graph.new_graph_property("string")
        middleProp = self.graph.new_graph_property("string")
        self.graph.graph_properties["Name"] = nameProp
        self.graph.graph_properties["End1"] = end1Prop
        self.graph.graph_properties["End2"] = end2Prop
        self.graph.graph_properties["Middle"] = middleProp
        self.graph.graph_properties["GraphID"] = graphIDProp
        self.graph.graph_properties["GraphID"] = str(self.ID)

    def get_number_of_nodes(self):
        return len(list(self.nodeHash))

    def get_graph_ID(self):
        return self.ID

    def create_graph_properties(self,propertyTypeHash):
        for p in propertyTypeHash:
            prop = self.graph.new_graph_properties[p]=propertyTypeHash[p]

    def set_graph_properties(self,propertyHash):
        for p in propertyHash:
            prop = self.graph.graph_properties[p]=propertyHash[p]

    def get_graph_property(self,property):
        return self.graph.graph_properties[property]

    def add_node(self,ID,name):
        n = self.graph.add_node(ID)
        n.new_vertex_properties["ID"]="string"
        n.new_vertex_properties["Name"]="string"
        n.vertex_properties["ID"]=ID
        n.vertex_properties["Name"]=name
        self.nodeHash[ID]=n

        typeHash=[]
        typeHash["end"]="string"
        typeHash["site"]="string"
        typeHash["connectedTo"]="string"
        self.set_graph_properties(typeHash)

    def get_node_property(self,v,property):
        return v.vertex_properties[property]

    def create_node_properties(self,v,propertyTypeHash):
        for p in propertyTypeHash:
            v.new_vertex_properties[p]=propertyTypeHash[p]

    def set_node_properties(self,v,propertyHash):
        for p in propertyHash:
            v.vertex_properties[p]=propertyHash[p]

    def get_node_by_ID(self,ID):
        return self.nodeHash[ID]

    def get_node_name(self,v):
        return v.vertex_properties["Name"]

    def get_nodes(self):
        return self.nodeHash

    def get_node_name_list(self):
        nodes = []
        for n in self.nodeHash:
            nodes.append(n)
        return n

    def add_edge(self,v1,v2):
        e = self.graph.add_edge(v1,v2)
        v1Name=self.get_node_name(v1)
        v2Name=self.get_node_name(v2)
        ID = str(v1Name)+"*"+str(v2Name)
        self.edgeHash[ID]=e
        typeHash={"weight":"float","GraphID":"string","end":"string","Name":"string"}
        self.create_edge_properties(e,typeHash)
        self.edge_properties["Name"]=ID

    def get_edge_by_ID(self,ID):
        return self.edgeHash[ID]

    def get_edge_name(self,e):
        return e.edge_properties["Name"]

    def create_edge_properties(self,e, propertyTypeHash):
        for p in propertyTypeHash:
            prop = e.new_edge_properties[p]=propertyTypeHash[p]

    def set_edge_properties(self,e,propertyHash):
        for p in propertyHash:
            e.edge_properties[p]=propertyHash[p]
    def get_edge(self,v1,v2):
        return self.edgeHash(str(v1)+"*"+str(v2))

    def get_edge_property(self,e,property):
        return e.edge_properties[property]

    def find_path(self, v1,v2):
        n1 = self.nodeHash[v1]
        n2 = self.nodeHash[v2]
        vlist,elist = gt.shortest_path(self.graph,n1,n2)
        return vlist, elist

    def get_minimum_spanning_tree(self):
        spanningWeight = self.graph.new_edge_property("double")
        for e in self.graph.edges():
            spanningWeight[e]=e.edge_property("Weight")

        return gt.min_spanning_tree(self.graph,weights=spanningWeight)

