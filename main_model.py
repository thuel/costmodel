#!/usr/bin/env python

"""Main cost model script
"""

#import essential modules and packages
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *
from pltgrph import *
from netwrkgrph import *
from math import radians, sin

def sinus_line_points(startx, starty, endx, step=0.25):
    """sinus line points calculation:
    """
    xgses = [ i*step for i in range(int((endx - startx)/step) + 1) ]
    yies = []
    for x in xgses:
        yies.append(x - sin(radians(x*360))*(step/sin(radians(step*360))))
    return lists_to_tuples([x + startx for x in xgses ], [y + starty for y in yies])

def lists_to_tuples(list1, list2):
    """Function to combine the elements of 2 lists into a list of tuples
    """
    if len(list1) == len(list2):
        result = []
        for i in range(len(list1)):
            result.append((list1[i], list2[i]))
        return result

if __name__=="__main__":
    graph=Graph()
    graph.add_edge(2,5,1,1, sinus_line_points(1,1,2,0.2))
    graph.add_edge(2,5,7,5)
    graph.add_edge(2,5,3,6)
    graph.add_edge(3,3,1,1, sinus_line_points(1,1,3,0.1)) 
    graph.add_edge(3,6,7,5)
    graph.add_edge(0,0,1,1,[(0.2,0.375),(0.375,0.625),(0.7,0.875)])
    graph.add_edge(0,0,1,1,[(0.2,0.025),(0.6,0.45),(0.8,0.6)])
    graph.add_edge(0,2,1,1)
    graph.add_edge(1,1,6,1)
    graph.add_edge(3,3,4,2)
    graph.add_edge(7,5,6,1)
    graph.add_edge(7,5,10,8)

    calc_neighbours(graph.vertices, graph.edges)
    
    for vertice in graph.vertices:
        #print graph.get_connecting_edges(graph.vertices[vertice])
        print("neighbours: ", graph.vertices[vertice].neighbours)

    for edge in graph.edges:
        e = graph.edges[edge]
        print(e, e.length())

    print(graph.min_corner_xy())
    print(graph.max_corner_xy())

    print("dimensions: ", graph.dimensions())
    print([edge.id() for edge in graph.edges.values()])

    d2 = dijkstra(graph, '011')
    d1 = dijkstra(graph, '011', '011')
    
    print('distances d2: %s' % d2[0])

    result = shortest_path(graph,'0108', '011')
    print(result)

    print(all_paths(dijkstra(graph, '011')))

    print_graph(graph)
