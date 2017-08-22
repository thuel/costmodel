#!/usr/bin/env python

"""Network graph library
"""

#import essential modules, libraries and methods/functions
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *
from math import sqrt, asin, degrees


"""Define classes of this library
"""

class Vertice(object):
    """Vertice to be used in a network graph. Attributes: x and y coordinates, neighbour vertices,
    and assigned graph object.
    """
    def __init__(self, x=0.0, y=0.0, neighbours={}, graph=None):
        self.x = x
        self.y = y
        self.neighbours = neighbours
        self.graph = graph

    def __str__(self):
        return "Vertice with coordinates: %f / %f" % (self.x, self.y)

    def id(self):
        return str(0) + str(self.x) + str(self.y)

    def get_neighbours(self):
        assert self.graph is not None, 'Can not get neighbours: vertice not part of graph!'
        calc_neighbours({self.id(): self}, self.graph.edges)

class Edge(object):
    """Edge to be used in a network graph. Attributes: start, end as Vertice objects, 
    intermediate points as a list of coordinate tuples and the assigned graph object.
    """
    def __init__(self, start=Vertice(0.0,0.0), end=Vertice(1.0,1.0), intermediates=list(), graph=None):
        self.start = start
        self.end = end
        self.sort_coordinates()
        self.intermediates = intermediates
        self.graph = graph
        
    def __str__(self):
        return "Edge with vertices at %f / %f and %f / %f." % (self.start.x, self.start.y, self.end.x, self.end.y)
        
    def length(self):
        return distance(self.start, self.end, self.intermediates)

    def id(self):
        """Create a unique id for an edge object. Set together from coordinates of
        start and end eventually joined by the intermediate points' coordinates.
        """
        interm = ""
        for t in self.intermediates:
            interm = interm + str(t[0]) + str(t[1])
        return str(9) + str(self.start.x) + str(self.start.y) + str(self.end.x) + str(self.end.y) + interm

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
        """Add an edge to the graph. Needs x/y coordinates of start and end vertice.
        Optional: a list of intermediate coordinate tuples.
        """
        start = self.vertices.get(str(0) + str(x1) + str(y1), Vertice(x1,y1))
        self.add_vertice(start)
        """Add the starting vertice to the graph
        """
        end = self.vertices.get(str(0) + str(x2) + str(y2), Vertice(x2,y2))
        self.add_vertice(end)
        """Add the ending vertice to the graph
        """
        new_edge = Edge(start, end, intermediates)
        if self.check_edge(new_edge):
            print("%s already exists. Doing noting" % new_edge)
        else:
            self.edges[new_edge.id()] = new_edge
            new_edge.graph = self
        """Add edge to graph if it doesn't exist.
        """

    def add_vertice(self, vertice):
        """Add vertice to the graph if it doesn't exist already."""
        if not self.check_vertice(vertice):
            self.vertices[vertice.id()] = vertice
            vertice.graph = self

    def check_vertice(self, vertice):
        """Check if a given vertice is already in the graph.
        """
        if vertice.id() in self.vertices:
            return True
        else:
            return False

    def check_edge(self, edge):
        """Check if a given edge is already in the graph.
        """
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
        
"""Functions used by the netwrkgrph classes, but which may eventually be used independent
of those classes are declared hereafter:
"""

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
       vertice object and the edges of the graph.
    """
    n = {}
    edges = get_connecting_edges(vertice, edges)
    for edge in edges:
        if edge.start is vertice:
            n[edge.end.id()] = edge.length()
        elif edge.end is vertice:
            n[edge.start.id()] = edge.length()
    return n

def calc_neighbours(vertices, edges):
    """For a dictionary of vertice objctes and a directory of edge objects, calculate
       the neighbouring vertices of every vertice object.
    """
    for vertice in vertices:
        vertices[vertice].neighbours = get_neighbours(vertices[vertice], edges)

def distance(p1, p2, i=list(), dist=0):
    """Calculate the distance between two points p1 and p2. optionally a list of
    intermediate points between p1 and p2 and a starting distance can be given to
    the function. Intermediate points allow to represent curved lines.
    The optional starting distance is needed in the self referencing part of the
    function with intermediates.
    """
    if len(i) == 0:
        dist += sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2)
        return dist
    else:
        x2, y2 = i[0]
        dist += sqrt((p1.x - x2) ** 2 + (p1.y - y2) ** 2)
        return distance(Vertice(x2,y2), p2, i[1:], dist)

def dijkstra(graph,src,dest,visited=[],distances={},predecessors={}):
    """ calculates a shortest path tree routed in src
    """
    assert src in graph.vertices, 'The root of the shortest path tree cannot be found'
    assert dest in graph.vertices, 'The root of the shortest path tree cannot be found'
    """ check if src and dest are part of the network graph.
    """
    if src == dest:
        """ Ending condition.
        """
        if len(distances) == 0:
            distances[dest] = 0
        path = []
        pred = dest
        while pred != None:
            path.append(pred)
            pred=predecessors.get(pred,None)
        print('shortest path %s is %d long.' % (str(path), distances[dest]))
        """ Build the shortes path and display it.
        """
        visited = []
        distances = {}
        predecessors = {}
        """ Reset the key variables to not influence following calculations
        """
    else:
        if not visited:
            distances[src] = 0
        """ Initialise the distance in the initial run.
        """
        for neighbour in graph.vertices[src].neighbours:
            if neighbour not in visited:
                new_distance = distances[src] + graph.vertices[src].neighbours[neighbour]
                if new_distance < distances.get(neighbour, float('inf')):
                    distances[neighbour] = new_distance
                    predecessors[neighbour] = src
                    print('predecessors: %s' % predecessors)
        """ Visit the neighbours and calculate distances.
        """
        visited.append(src)
        """ Mark the src as visited.
        """
        unvisited = {}
        for vertice in graph.vertices:
            if vertice not in visited:
                unvisited[vertice] = distances.get(vertice, float('inf'))
        """ Prepare for recursion of the dijkstra function.
        """
        new_src = min(unvisited, key=unvisited.get)
        return dijkstra(graph, new_src, dest, visited, distances, predecessors)
        """ Recurse through the dijkstra function with the new "source"
        """
        
    
if __name__ == "__main__":
    graph=Graph()
    graph.add_edge(2,5,1,1)
    graph.add_edge(2,5,7,5)
    graph.add_edge(2,5,3,6)
    graph.add_edge(3,3,1,1)
    graph.add_edge(3,6,7,5)
    graph.add_edge(0,0,1,1,[(0.25,0.3),(0.5,0.8),(0.75,0.8)])
    graph.add_edge(0,0,1,1,[(0.25,0.2),(0.5,0.2),(0.75,0.2)])
    graph.add_edge(0,2,1,1)
    graph.add_edge(1,1,6,1)
    graph.add_edge(3,3,4,2)
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

    dijkstra(graph, '000', '075')
    dijkstra(graph, '011', '011')
    
