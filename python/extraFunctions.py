__author__ = 'clipo'


    def outputGraphArray(self, array):
        num = 0
        os.environ["PATH"] += ":/usr/local/bin:"
        for g in array:
            num += 1
            pos = nx.graphviz_layout(g, prog="twopi", root=['graphroot'])
            gfile = self.outputDirectory + self.inputFile[0:-4] + "-min-sol-" + str(num) + ".png"
            filename = self.outputDirectory + self.inputFile[0:-4] + "-min-sol-" + str(num) + ".gml"
            self.saveGraph(g, filename)
            edgewidth = []
            weights = nx.get_edge_attributes(g, 'weight')
            for w in weights:
                edgewidth.append(weights[w])

            maxValue = max(edgewidth)
            widths = []
            for w in edgewidth:
                widths.append(((maxValue - w) + 1) * 5)

            assemblageSizes = []
            sizes = nx.get_node_attributes(g, 'size')
            #print sizes
            for s in sizes:
                assemblageSizes.append((sizes[s]/self.totalAssemblageSize)*10)
            nx.draw_networkx_edges(g, pos, alpha=0.3, width=widths)
            sizes = nx.get_node_attributes(g, 'size')
            nx.draw_networkx_nodes(g, pos, node_size=assemblageSizes, node_color='w', alpha=0.4)
            nx.draw_networkx_edges(g, pos, alpha=0.4, node_size=0, width=1, edge_color='k')
            nx.draw_networkx_labels(g, pos, fontsize=10)
            font = {'fontname': 'Helvetica',
                    'color': 'k',
                    'fontweight': 'bold',
                    'fontsize': 10}
            plt.axis('off')
            plt.savefig(gfile, dpi=75)
            plt.figure(gfile, figsize=(8, 8))



    def minimumSpanningTree(self, networks, sumGraph, outputDirectory, inputFile):
        try:
            from networkx import graphviz_layout
        except ImportError:
            raise ImportError("This function needs Graphviz and either PyGraphviz or Pydot")

        newfilename = outputDirectory + inputFile[0:-4] + "-mst.png"
        plt.figure(newfilename, figsize=(8, 8))

        graphs = []
        megaGraph = nx.Graph(is_directed=False)
        graphCount = 0
        for net in networks:
            graphCount += 1
            g = nx.Graph(is_directed=False)
            for nodey in net.nodes(data=True):
                xCoordinate = 0
                yCoordinate = 0
                name = nodey[0]
                xCoordinate = self.xAssemblage[name]
                yCoordinate = self.yAssemblage[name]
                megaGraph.add_node(name, name=name, xCoordinate=xCoordinate, yCoordinate=yCoordinate,
                                   size=(self.assemblageSize[name]/self.totalAssemblageSize)*10)

            count = 0
            for e in net.edges_iter():
                d = net.get_edge_data(*e)
                fromAssemblage = e[0]
                toAssemblage = e[1]
                g.add_node(fromAssemblage, label=fromAssemblage, x=xCoordinate, y=yCoordinate,
                           name=fromAssemblage, size=self.assemblageSize[fromAssemblage]/(self.totalAssemblageSize)*10)
                g.add_node(toAssemblage, label=toAssemblage, x=xCoordinate, y=yCoordinate,
                           name=toAssemblage, size=self.assemblageSize[toAssemblage]/(self.totalAssemblageSize)*10)

                weight = d['weight'] + 1
                distance = self.distanceBetweenAssemblages[fromAssemblage + "*" + toAssemblage]
                #count = megaGraph.get_edge_data(fromAssemblage,toAssemblage,'weight'
                count += 1
                megaGraph.add_path([fromAssemblage, toAssemblage], weight=count,
                                   distance=distance,
                                   size=(self.assemblageSize[fromAssemblage], self.assemblageSize[toAssemblage]))

                g.add_path([fromAssemblage, toAssemblage],
                           xy1=(self.xAssemblage[fromAssemblage], self.yAssemblage[fromAssemblage]),
                           xy2=(self.xAssemblage[toAssemblage], self.yAssemblage[toAssemblage]),
                           weight=weight,
                           meanDistance=distance,
                           size=(self.assemblageSize[fromAssemblage], self.assemblageSize[toAssemblage]))
            graphs.append(g)
        plt.rcParams['text.usetex'] = False
        mst = nx.minimum_spanning_tree(megaGraph, weight='weight')
        os.environ["PATH"] += ":/usr/local/bin:"
        pos = nx.graphviz_layout(mst)
        edgewidth = []
        weights = nx.get_edge_attributes(mst, 'inverseweight')
        for w in weights:
            edgewidth.append(weights[w])

        maxValue = max(edgewidth)
        widths = []
        for w in edgewidth:
            widths.append(((maxValue - w) + 1) * 5)

        color = nx.get_edge_attributes(mst, 'color')
        colorList = []
        for c in color:
            colorList.append(color[c])
        colors = []
        colorMax = max(colorList)
        for c in colorList:
            colors.append(c / colorMax)
        assemblageSizes = []
        sizes = nx.get_node_attributes(mst, 'size')
        #print sizes
        for s in sizes:
            #print sizes[s]
            assemblageSizes.append(sizes[s])
        nx.draw_networkx_edges(mst, pos, alpha=0.3, width=widths, edge_color=colorList)
        sizes = nx.get_node_attributes(mst, 'size')
        nx.draw_networkx_nodes(mst, pos, node_size=assemblageSizes, node_color='w', alpha=0.4)
        nx.draw_networkx_edges(mst, pos, alpha=0.4, node_size=0, width=1, edge_color='k')
        nx.draw_networkx_labels(mst, pos, fontsize=10)
        font = {'fontname': 'Helvetica',
                'color': 'k',
                'fontweight': 'bold',
                'fontsize': 10}

        plt.axis('off')
        plt.savefig(newfilename, dpi=75)
        if self.args['shapefile'] is not None and self.args['xyfile'] is not None:
            self.createShapefile(mst, outputDirectory + inputFile[0:-4] + "-mst.shp")
        self.saveGraph(mst, newfilename + ".gml")
        atlasFile = outputDirectory + inputFile[0:-4] + "-atlas.png"
        plt.figure(atlasFile, figsize=(8, 8))
        UU = nx.Graph(is_directed=False)
        # do quick isomorphic-like check, not a true isomorphism checker
        nlist = [] # list of nonisomorphic graphs
        for G in graphs:
            # check against all nonisomorphic graphs so far
            if not self.iso(G, nlist):
                nlist.append(G)

        UU = nx.disjoint_union_all(graphs) # union the nonisomorphic graphs
        pos = nx.graphviz_layout(UU, prog="twopi", root=self.args['graphroot'])
        # color nodes the same in each connected subgraph
        C = nx.connected_component_subgraphs(UU)
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
        plt.savefig(atlasFile, dpi=250)
        self.saveGraph(UU, atlasFile + ".gml")
