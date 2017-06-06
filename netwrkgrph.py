#!/usr/bin/env python

"""Network graph library
"""

#import essential modules, libraries and methods/functions
from math import sqrt, asin, degrees

"""Define classes of this library
"""

#network graph point/vertice class
class Vertice(object):
    """Vertice to be used in a network graph. Attributes: x and y coordinatesm, neigbours.
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
    def __init__(self, start=Vertice(0.0,0.0), end=Vertice(1.0,1.0)):
        self.start = start
        self.end = end
        self.sort_coordinates()
        print "angle: %f, start: %s, end: %s" % (self.angle(), self.start, self.end)

    def __str__(self):
        return "Edge with vertices at %f / %f and %f / %f." % (self.start.x, self.start.y, self.end.x, self.end.y)

    def length(self):
        return sqrt((self.start.x - self.end.x) ** 2 + (self.start.y - self.end.y) ** 2)

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

    def add_edge(self, x1, y1, x2, y2):
        start = self.vertices.get(str(0) + str(x1) + str(y1), Vertice(x1,y1))
        if not self.check_vertice(start):
            self.add_vertice(start)
        end = self.vertices.get(str(0) + str(x2) + str(y2), Vertice(x2,y2))
        if not self.check_vertice(end):
            self.add_vertice(end)
        new_edge = Edge(start, end)
        if self.check_edge(new_edge):
            print "%s already exists. Doing noting" % new_edge
        else:
            self.edges[new_edge.id()] = new_edge

    def add_vertice(self, vertice):
        self.vertices[vertice.id()] = vertice

    def check_vertice(self, vertice):
        if vertice.id() in self.vertices:
            print "Vertice exists."
            return True
        else:
            print "Vertice doesn't exist."
            return False

    def check_edge(self, edge):
        if edge.id() in self.edges:
            return True
        else:
            return False

    def get_connecting_edges(self, vertice):
        connecting_edges = []
        for edge in self.edges:
            if self.edges[edge].start is vertice or self.edges[edge].end is vertice:
                connecting_edges.append(self.edges[edge])
        return connecting_edges

    def get_neighbours(self, vertice):
        n = {}
        edges = self.get_connecting_edges(vertice)
        for edge in edges:
            if edge.start is vertice:
                n[edge.end.id()] = edge.end
            elif edge.end is vertice:
                n[edge.start.id()] = edge.start
        return n

    def calc_neighbours(self):
        for vertice in self.vertices:
            self.vertices[vertice].neighbours = self.get_neighbours(self.vertices[vertice])

if __name__ == "__main__":
    graph=Graph()
    graph.add_edge(2,5,1,1)
    graph.add_edge(3,3,1,1)
    graph.add_edge(3,6,7,5)
    graph.add_edge(0,0,1,1)
    graph.add_edge(0,2,1,1)
    graph.add_edge(1,1,6,1)
    graph.add_edge(3,3,4,2)
    graph.calc_neighbours()
    print graph
    for vertice in graph.vertices:
        #print graph.get_connecting_edges(graph.vertices[vertice])
        print "neighbours: ", graph.vertices[vertice].neighbours
