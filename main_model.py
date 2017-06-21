#!/usr/bin/env python

"""Main cost model script
"""

from pltgrph import *
from netwrkgrph import *

graph=Graph()
graph.add_edge(2,5,1,1)
graph.add_edge(2,5,7,5)
graph.add_edge(2,5,3,6)
graph.add_edge(3,3,1,1) #try to add list with sinus wave
"""sinus wave:
startx, step, endx as givens normed to 0
xgses = [ i*step for i in range(int((endx - startx)/step) + 1) ]
yies = []
for x in xgses:
    yies.append(x - sin(x*step*360)*(step/sin(step*360)))

"""
graph.add_edge(3,6,7,5)
graph.add_edge(0,0,1,1,[(0.2,0.375),(0.375,0.625),(0.7,0.875)])
graph.add_edge(0,2,1,1)
graph.add_edge(1,1,6,1)
graph.add_edge(3,3,4,2)
calc_neighbours(graph.vertices, graph.edges)

print_graph(graph)
