#!/usr/bin/env python

"""Network graph library
"""

#import essential modules, libraries and methods/functions
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *
from math import sqrt, asin, degrees
from priodict import priorityDictionary


"""Define classes of this library
"""

class Vertex(object):
    """Vertex to be used in a network graph. Attributes: x and y coordinates, neighbour vertices,
    adjoining edges and assigned graph object.
    """
    def __init__(self, x=0.0, y=0.0, graph=None):
        if graph is None:
            graph = None
        self.x = x
        self.y = y
        self.graph = graph
        self.neighbours = {}
        self.edges = {}

    def __str__(self):
        return "Vertex with coordinates: %f / %f" % (self.x, self.y)

    def id(self):
        return str(0) + str(self.x) + str(self.y)

    def get_neighbours(self):
        assert self.graph is not None, 'Can not get neighbours: vertex not part of graph!'
        calc_neighbours({self.id(): self}, self.graph.edges)

class Edge(object):
    """Edge to be used in a network graph. Attributes: start, end as Vertex objects, 
    intermediate points as a list of coordinate tuples and the assigned graph object.
    """
    def __init__(self, start=Vertex(0.0,0.0), end=Vertex(1.0,1.0), intermediates=list(), graph=None):
        self.start = start
        self.end = end
        self.sort_coordinates()
        self.intermediates = intermediates
        self.graph = graph

        length = self.length()
        if length < start.neighbours.get(end.id(),float('inf')):
            start.neighbours[end.id()] = length
        start.edges[self.id()] = length
        if length < end.neighbours.get(start.id(),float('inf')):
            end.neighbours[start.id()] = length
        end.edges[self.id()] = length
        
        
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
    """Graph to be populated by Vertex and Edge objects.
    """
    def __init__(self, vertices={}, edges={}):
        self.vertices = vertices
        self.edges = edges

    def __str__(self):
        return "Graph with vertices %s and edges %s." % (self.vertices, self.edges)

    def add_edge(self, x1, y1, x2, y2, intermediates=list()):
        """Add an edge to the graph. Needs x/y coordinates of start and end vertex.
        Optional: a list of intermediate coordinate tuples.
        """
        start = self.vertices.get(str(0) + str(x1) + str(y1), Vertex(x1,y1))
        self.add_vertex(start)
        """Add the starting vertex to the graph
        """
        end = self.vertices.get(str(0) + str(x2) + str(y2), Vertex(x2,y2))
        self.add_vertex(end)
        """Add the ending vertex to the graph
        """
        new_edge = Edge(start, end, intermediates)
        if self.check_edge(new_edge):
            print("%s already exists. Doing noting" % new_edge)
        else:
            self.edges[new_edge.id()] = new_edge
            new_edge.graph = self
        """Add edge to graph if it doesn't exist.
        """

    def add_vertex(self, vertex):
        """Add vertex to the graph if it doesn't exist already."""
        if not self.check_vertex(vertex):
            self.vertices[vertex.id()] = vertex
            vertex.graph = self

    def check_vertex(self, vertex):
        """Check if a given vertex is already in the graph.
        """
        if vertex.id() in self.vertices:
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
            vertex = d[i]
            min_x = min(vertex.x, min_x)
            min_y = min(vertex.y, min_y)
        return (min_x, min_y)

    def max_corner_xy(self):
        """Returns the minimum point of the rectangle spanwning the graph
        as tuple of coordinates.
        """
        max_x = -float("inf")
        max_y = -float("inf")
        d = self.vertices
        for i in d:
            vertex = d[i]
            max_x = max(vertex.x, max_x)
            max_y = max(vertex.y, max_y)
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

def get_connecting_edges(vertex, edges):
    """The edges - given as a dictionary of edge objects - connecting a given vertex
       object to a graph are returned as list.
    """
    connecting_edges = []
    for edge in edges:
        if edges[edge].start is vertex or edges[edge].end is vertex:
            connecting_edges.append(edges[edge])
    return connecting_edges

def get_neighbours(vertex, edges):
    """Get a dictionary of vertex objects of the neighbouring vertices for a given
       vertex object and the edges of the graph.
    """
    n = {}
    edges = get_connecting_edges(vertex, edges)
    for edge in edges:
        if edge.start is vertex:
            n[edge.end.id()] = edge.length()
        elif edge.end is vertex:
            n[edge.start.id()] = edge.length()
    return n

def calc_neighbours(vertices, edges):
    """For a dictionary of vertex objctes and a directory of edge objects, calculate
       the neighbouring vertices of every vertex object.
    """
    for vertex in vertices:
        vertices[vertex].neighbours = get_neighbours(vertices[vertex], edges)

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
        return distance(Vertex(x2,y2), p2, i[1:], dist)

def dijkstra(graph, start, end=None):
    """ Find shortest paths from the  start vertex to all vertices nearer
    than or equal to the end.

    The input graph "graph" is assumed to be Graph object. A vertex is a   
    Vertex object. For any vertex v, graph.vertices[v].neighbours is itself 
    a dictionary, indexed by the neighbors of v. For any edge v->w, 
    graph.vertices[v].neighbours[w] is the length of the edge.

    The output is a pair (D,P) where D[v] is the distance from start to
    v and P[v] is the predecessor of v along the shortest path from s to
    v.

    Dijkstra's algorithm is only guaranteed to work correctly when all
    edge lengths are positive. This code does not verify this property
    for all edges (only the edges examined until the end vertex is
    reached), but will correctly compute shortest paths even for some
    graphs with negative edges, and will raise an exception if it
    discovers that a negative edge has caused it to make a mistake.
    """

    D = {}  # dictionary of final distances
    P = {}  # dictionary of predecessors
    Q = priorityDictionary()  # estimated distances of non-final vertices
    Q[start] = 0

    for vertex in Q:
        D[vertex] = Q[vertex]
        if vertex == end:
            break

        for neighbour in graph.vertices[vertex].neighbours:
            length = D[vertex] + graph.vertices[vertex].neighbours[neighbour]
            if neighbour in D:
                if length < D[neighbour]:
                    raise ValueError("Dijkstra: found better path to already-final vertex")
            elif neighbour not in Q or length < Q[neighbour]:
                Q[neighbour] = length
                P[neighbour] = vertex

    return (D, P, start)


def shortest_path(graph, start, end):
    """ Find a single shortest path from the given start vertex to the given
    end vertex. The input has the same conventions as Dijkstra(). The
    output is a list of the vertices in order along the shortest path.
    """

    D, P = dijkstra(graph, start, end)[:2]
    path = []
    while True:
        path.append(end)
        if end == start:
            break
        end = P[end]
    path.reverse()
    return path
        
def all_paths(dijkstra):
    """ Returns a dict of all the paths that were calculated by dijkstra().
    """
    D, P, start = dijkstra
    paths = {}
    for vertex in D:
        p_index = vertex
        path = []
        while True:
            path.append(vertex)
            if vertex == start:
                break
            vertex = P[vertex]
        path.reverse()
        paths[p_index] = {'path': path, 'distance': D[p_index]}
    return paths


    
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
    graph.add_edge(7,5,6,1)
    graph.add_edge(7,5,10,8)

#    calc_neighbours(graph.vertices, graph.edges)
    
    for vertex in graph.vertices:
        #print graph.get_connecting_edges(graph.vertices[vertex])
        print("neighbours: ", graph.vertices[vertex].neighbours)

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
