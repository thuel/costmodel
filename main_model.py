#!/usr/bin/env python

"""Main cost model script
"""

from pltgrph import *
from netwrkgrph import *
from math import radians, sin

def sinus_line_points(startx, starty, endx, endy, step=0.25):
    """sinus line points calculation:
    """
    len_x = endx - startx
    npoints = int(len_x/step) + 1 # number of points to calculate
    xgses = [ i/len_x*step for i in range(npoints) ] # start and end to be removed
    xs = [ ]
    ys = [ ]
    for x in xgses[1:npoints-1]:
        ys.append( abs(endy-starty) * ( (x - sin( radians( x*360 ) ) * ( step / sin( radians( step*360 ) ) ) ) ) )
        xs.append( len_x * x )
    return lists_to_tuples([x + startx for x in xs ], [y + starty for y in ys])

def lists_to_tuples(list1, list2):
    if len(list1) == len(list2):
        result = []
        for i in range(len(list1)):
            result.append((list1[i], list2[i]))
        return result

if __name__=="__main__":
    graph=Graph()
    graph.add_edge(1,1,2,5, sinus_line_points(1,1,2,5,0.05))
    graph.add_edge(2,5,7,5)
    graph.add_edge(2,5,3,6)
    graph.add_edge(3,3,1,1, sinus_line_points(1,1,3,3,0.1)) 
    graph.add_edge(3,6,7,5, sinus_line_points(3,6,7,5,0.1))
    graph.add_edge(0,0,1,1,[(0.2,0.375),(0.375,0.625),(0.7,0.875)])
    graph.add_edge(0,2,1,1)
    graph.add_edge(1,1,6,1)
    graph.add_edge(3,3,4,2)
    calc_neighbours(graph.vertices, graph.edges)

    print sinus_line_points(1,1,3,5,0.1)

    print_graph(graph)
