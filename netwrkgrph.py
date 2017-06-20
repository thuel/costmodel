#!/usr/bin/env python

"""Network graph library
"""

#import essential modules, libraries and methods/functions
from math import sqrt, asin, degrees

"""Define classes of this library
"""

class Vertice(object):
    """Vertice to be used in a network graph. Attributes: x and y coordinates, neigbour vertices,
        type.
    """
    def __init__(self, x=0.0, y=0.0, neighbours={}):
        self.x = x
        self.y = y
        self.neighbours = neighbours

    def __str__(self):
        return "Vertice with coordinates: %f / %f" % (self.x, self.y)

    def id(self):
        return str(0) + str(self.x) + str(self.y)

class Edge(object):
    """Edge to be used in a network graph. Attributes: start, end as Vertice objects.
    """
    def __init__(self, start=Vertice(0.0,0.0), end=Vertice(1.0,1.0), intermediates=list()):
        self.start = start
        self.end = end
        self.sort_coordinates()
        self.intermediates = intermediates
        
    def __str__(self):
        return "Edge with vertices at %f / %f and %f / %f." % (self.start.x, self.start.y, self.end.x, self.end.y)

    def length(self):
        return distance(self.start, self.end, self.intermediates)

    def id(self):
        return str(9) + str(self.start.x) + str(self.start.y) + str(self.end.x) + str(self.end.y)

    def sort_coordinates(self):
        if self.start.x > self.end.x:
            self.start, self.end = (self.end, self.start)

    def angle(self):
        return degrees(asin((self.end.y - self.start.y)/self.length()))

class Graph(object):
    """Graph to be populated by Vertice and Edge objects.
    """
    def __init__(self, vertices={}, edges={}):
        self.vertices = vertices
        self.edges = edges

    def __str__(self):
        return "Graph with vertices %s and edges %s." % (self.vertices, self.edges)

    def add_edge(self, x1, y1, x2, y2, intermediates=list()):
        start = self.vertices.get(str(0) + str(x1) + str(y1), Vertice(x1,y1))
        if not self.check_vertice(start):
            self.add_vertice(start)
        end = self.vertices.get(str(0) + str(x2) + str(y2), Vertice(x2,y2))
        if not self.check_vertice(end):
            self.add_vertice(end)
        new_edge = Edge(start, end, intermediates)
        if self.check_edge(new_edge):
            print "%s already exists. Doing noting" % new_edge
        else:
            self.edges[new_edge.id()] = new_edge

    def add_vertice(self, vertice):
        self.vertices[vertice.id()] = vertice

    def check_vertice(self, vertice):
        if vertice.id() in self.vertices:
            return True
        else:
            return False

    def check_edge(self, edge):
        if edge.id() in self.edges:
            return True
        else:
            return False

    def min_corner_xy(self):
        """Returns the minimum point of the rectangle spanwning the graph
        as tuple of coordinates.
        """
        min_x = float("inf")
        min_y = float("inf")
        d = self.vertices
        for i in d:
            vertice = d[i]
            min_x = min(vertice.x, min_x)
            min_y = min(vertice.y, min_y)
        return (min_x, min_y)

    def max_corner_xy(self):
        """Returns the minimum point of the rectangle spanwning the graph
        as tuple of coordinates.
        """
        max_x = -float("inf")
        max_y = -float("inf")
        d = self.vertices
        for i in d:
            vertice = d[i]
            max_x = max(vertice.x, max_x)
            max_y = max(vertice.y, max_y)
        return (max_x, max_y)

    def dimensions(self):
        """Returns the width and the height of the graph as tuple.
        """
        min_x, min_y = self.min_corner_xy()
        max_x, max_y = self.max_corner_xy()
        return (max_x - min_x, max_y - min_y)
        

def get_connecting_edges(vertice, edges):
    """The edges - given as a dictionary of edge objects - connecting a given vertice
       object to a graph are returned as list.
    """
    connecting_edges = []
    for edge in edges:
        if edges[edge].start is vertice or edges[edge].end is vertice:
            connecting_edges.append(edges[edge])
    return connecting_edges

def get_neighbours(vertice, edges):
    """Get a dictionary of vertice objects of the neighbouring vertices for a given
       vertice object and the edges of the graph. TODO: change order of calling
       connecting edges and the neighbour calculation.
    """
    n = {}
    edges = get_connecting_edges(vertice, edges)
    for edge in edges:
        if edge.start is vertice:
            n[edge.end.id()] = edge.end
        elif edge.end is vertice:
            n[edge.start.id()] = edge.start
    return n

def calc_neighbours(vertices, edges):
    """For a dictionary of vertice objctes and a directory of edge objects, calculate
       the neighbouring vertices of every vertice object.
    """
    for vertice in vertices:
        vertices[vertice].neighbours = get_neighbours(vertices[vertice], edges)

def distance(p1, p2, i=list(), dist=0):
    if len(i) == 0:
        dist += sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)
        return dist
    else:
        x2, y2 = i[0]
        dist += sqrt((p1.x - x2) ** 2 + (p1.y - y2) ** 2)
        return distance(Vertice(x2,y2), p2, i[1:], dist)
        
if __name__ == "__main__":
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
    for vertice in graph.vertices:
        #print graph.get_connecting_edges(graph.vertices[vertice])
        print "neighbours: ", graph.vertices[vertice].neighbours

    for edge in graph.edges:
        e = graph.edges[edge]
        print e, e.length()

    print graph.min_corner_xy()
    print graph.max_corner_xy()

    print "dimensions: ", graph.dimensions()
    
