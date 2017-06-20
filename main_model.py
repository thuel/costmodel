#!/usr/bin/env python

"""Main cost model script
"""

from pltgrph import *
from netwrkgrph import *

graph=Graph()
graph.add_edge(2,5,1,1)
graph.add_edge(2,5,7,5)
graph.add_edge(2,5,3,6)
graph.add_edge(3,3,1,1)
graph.add_edge(3,6,7,5)
graph.add_edge(0,0,1,1,[(0.25,0.3),(0.5,0.8),(0.75,0.8)])
graph.add_edge(0,2,1,1)
graph.add_edge(1,1,6,1)
graph.add_edge(3,3,4,2)
calc_neighbours(graph.vertices, graph.edges)

print_graph(graph)
