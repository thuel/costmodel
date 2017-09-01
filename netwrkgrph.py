#!/usr/bin/env python

"""Network graph library
"""

#import essential modules, libraries and methods/functions
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *
from math import sqrt, sin, cos, asin, atan2, acos, degrees, radians
from priodict import priorityDictionary
import copy


"""Define classes of this library
"""

class Vertex(object):
    """Vertex to be used in a network graph. Attributes: x and y coordinates, neighbour vertices,
    adjoining edges and assigned graph object.
    """
    def __init__(self, x=0.0, y=0.0, kind='helper', graph=None):
        if graph is None:
            graph = None
        self.x = x
        self.y = y
        self.graph = graph
        self.neighbours = {}
        self.edges = {}
        self.sub_edges = {}
        self.kind = kind

    def __str__(self):
        return "Vertex of type %s with coordinates: %f / %f" % (self.kind, self.x, self.y)

    def id(self):
        return str(0) + str(self.x) + str(self.y)

    def update_neighbours(self):
        assert self.graph is not None, 'Can not get neighbours: vertex not part of graph!'
        update_neighbours({self.id(): self})

class Edge(object):
    """Edge to be used in a network graph. Attributes: start, end as Vertex objects, 
    intermediate points as a list of Vertex objects and the assigned graph object.
    """
    def __init__(self, start=Vertex(0.0,0.0), end=Vertex(1.0,1.0), intermediates=None, graph=None):
        self.start = start
        self.end = end
        self.sort_coordinates()
        self.vertices = [start.id(), end.id()]
        if intermediates is None:
            intermediates = list()
        self.intermediates = intermediates
        if graph is None:
            graph = None
        self.graph = graph
        self.sub_edges = self.get_sub_edges()

        length = self.length()
        if length < start.neighbours.get(end.id(),float('inf')) and self.graph == start.graph:
            start.neighbours[end.id()] = length
        if self.id() not in start.edges:
            start.edges[self.id()] = self
        if length < end.neighbours.get(start.id(),float('inf')) and self.graph == end.graph:
            end.neighbours[start.id()] = length
        if self.id() not in end.edges:
            end.edges[self.id()] = self
        self.start.kind = 'vertex'
        self.end.kind = 'vertex'
        
        
    def __str__(self):
        return "Edge with vertices at %f / %f and %f / %f." % (self.start.x, self.start.y, self.end.x, self.end.y)
        
    def __delete__(self):
        del self.start.edges[self.id()]
        del self.end.edges[self.id()]
        del self.start.neighbours[self.end.id()]
        del self.end.neighbours[self.start.id()]
        
    def length(self):
        return distance(self.start, self.end, self.intermediates)

    def id(self):
        """Create a unique id for an edge object. Set together from coordinates of
        start and end eventually joined by the intermediate points' coordinates.
        """
        interm = ""
        for t in self.intermediates:
            interm = interm + str(t.x) + str(t.y)
        return str(9) + str(self.start.x) + str(self.start.y) + str(self.end.x) + str(self.end.y) + interm

    def sort_coordinates(self):
        if self.start.x > self.end.x:
            self.start, self.end = (self.end, self.start)

    def angle(self):
        return angle(self.start, self.end)

    def get_sub_edges(self):
        edges = []
        vertices = [self.start] + self.intermediates + [self.end]
        for i in range(len(vertices) - 1):
            edges.append(SubEdge(vertices[i], vertices[i+1], self.id(), self.graph))
        return edges

class SubEdge(Edge):
    """ Sub edge to be used in a network graph to represent edges between intermediate points of
    an Edge object. Attributes: start, end as Vertex objects, the parent Edge object and the
    assigned graph object.
    """
    def __init__(self, start, end, parent, graph=None):
        self.start = start
        self.end = end
        self.sort_coordinates()
        self.vertices = [start.id(), end.id()]
        if graph is None:
            graph = None
        self.graph = graph
        self.parent = parent

        self.start.sub_edges[self.id()] = self
        self.end.sub_edges[self.id()] = self
                
    def __str__(self):
        return "Sub edge with vertices at %f / %f and %f / %f." % (self.start.x, self.start.y, self.end.x, self.end.y)

    def id(self):
        """Create a unique id for an edge object. Set together from coordinates of
        start and end eventually joined by the intermediate points' coordinates.
        """
        return str('sb9') + str(self.start.x) + str(self.start.y) + str(self.end.x) + str(self.end.y)
    
    def delete(self):
        del self.start.sub_edges[self.id()]
        del self.end.sub_edges[self.id()]

class Graph(object):
    """Graph to be populated by Vertex and Edge objects.
    """
    def __init__(self, name, vertices=None, edges=None):
        if vertices is None:
            vertices = dict()
        if edges is None:
            edges = dict()
        self.name = name
        self.vertices = vertices
        self.edges = edges

    def __str__(self):
        return "Graph with vertices %s and edges %s." % (self.vertices, self.edges)

    def add_edge(self, x1, y1, x2, y2, intermediates=None):
        """Add an edge to the graph. Needs x/y coordinates of start and end vertex.
        Optional: a list of intermediate coordinate tuples.
        """
        if intermediates is None:
            intermediates = list()
        start = self.vertices.get(str(0) + str(x1) + str(y1), Vertex(x1,y1,'vertex',self))
        self.add_vertex(start)
        """Add the starting vertex to the graph
        """
        end = self.vertices.get(str(0) + str(x2) + str(y2), Vertex(x2,y2,'vertex',self))
        self.add_vertex(end)
        """Add the ending vertex to the graph
        """
        intermediates = [ Vertex(i[0], i[1], 'helper', self) for i in intermediates ]
        for vertex in intermediates:
            self.add_vertex(vertex)
        """ Generate a list of Vertex objects for the intermediate points and add them to
        the network Graph object.
        """
        new_edge = Edge(start, end, intermediates, self)
        if self.check_edge(new_edge):
            print("%s already exists. Doing noting" % new_edge)
        else:
            self.edges[new_edge.id()] = new_edge
        """Add edge to graph if it doesn't exist.
        """

    def check_edge(self, edge):
        """Check if a given edge is already in the graph.
        """
        if edge.id() in self.edges:
            return True
        else:
            return False

    def delete_edge(self, edge):
        """ Delete the given Edge from the edges dict of the corresponding Graph 
        object and the start and end Vertex objects.
        """            
        del edge.start.edges[edge.id()]
        del edge.start.neighbours[edge.end.id()]
        del edge.end.edges[edge.id()]
        del edge.end.neighbours[edge.start.id()]
        for e in edge.sub_edges:
            e.delete()
        edge.sub_edges = []
        del self.edges[edge.id()]

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

def get_neighbours(vertex):
    """Get a dictionary of vertex objects of the neighbouring vertices for a given
       vertex object and the edges of the graph.
    """
    n = {}
    edges = vertex.edges
    graph_vertices = vertex.graph.vertices
    for e_id, edge in edges.items():
        start, end, length = edge.start, edge.end, edge.length()
        if start is vertex and length < n.get(start.id(), float('inf')) and end.id() in graph_vertices:
            n[end.id()] = length
        elif end is vertex and length < n.get(end.id(), float('inf')) and start.id() in graph_vertices:
            n[start.id()] = length
    return n

def update_neighbours(vertices):
    """For a dictionary of vertex objctes and a directory of edge objects, calculate
       the neighbouring vertices of every vertex object.
    """
    for vertex in vertices:
        vertices[vertex].neighbours = get_neighbours(vertices[vertex])

def distance(v1, v2, i=None, dist=None):
    """ Calculate the distance between two Vertex objects v1 and v2 and optionally a list
    of intermediate points between v1 and v2 and a starting distance can be given to
    the function. Intermediate points allow to represent curved lines.
    The optional starting distance is needed in the self referencing part of the
    function with intermediates.
    """
    if i is None:
        i = list()
    if dist is None:
        dist = 0
    """ Reset optional arguments to be able to call the function more than once in
    the same session.
    """
    if len(i) == 0:
        dist += sqrt((v1.x - v2.x) ** 2 + (v1.y - v2.y) ** 2)
        return dist
    else:
        x2, y2 = (i[0].x, i[0].y)
        dist += sqrt((v1.x - x2) ** 2 + (v1.y - y2) ** 2)
        return distance(Vertex(x2,y2), v2, i[1:], dist)

def angle(v1, v2):
    """ Return the angle of the line connecting to Vertex objects.
    """
    if distance(v1, v2) != 0:
        return atan2(v2.y - v1.y, v2.x - v1.x)
    else:
        return 0
    """ Ensure that there is a distance between the two Vertex objects, otherwise
    return 0.
    """
    #return degrees(asin((v2.y - v1.y)/distance(v1, v2)))

def near_points(vertex, vertices, radius):
    """ For a given id of a Vertex object "vertex_id" in a dictionary of Vertex
    objects "vertices", return a dictionary of vertices which are within a given
    distance "radius" from the Vertex object with id vertex_id.
    """
    near_points = {}
    for key, value in vertices.items():
        """ key: vertex id, value: Vertex object
        """
        dist = distance(vertex, value)
        if dist < radius and dist > 0:
            near_points[key] = value
    return near_points

def near_edges(vertex, edges, radius):
    """ Return a dictionary of the nearest edges for a given Vertex object,
    a dictionary of Edge objects and the search radius.
    """
    near_edges = {}
    graph = vertex.graph
    vertices = get_vertices_from_edges(edges, graph)
    """ The dictionary with the Vertex objects to be used next has to be
    initiated from the edges given to the function.
    """
    if vertex.id() not in vertices:
        vertices[vertex.id()] = vertex
    """ Add the vertex itself to the dictionary of vertices, otherwise
    near_points() won't be able to return a valid value.
    """
    nearest_vertices = near_points(vertex, vertices, radius)
    """ Get the relevant vertices from all the vertices. Eventually the following
    calculations are faster.
    """
    for v_id, v in nearest_vertices.items():
        """ Iterate through every Vertex object in nearest_vertices.
        """
        for e_id, edge in v.sub_edges.items():
            if e_id not in near_edges and v_id in edge.vertices:
                near_edges[e_id] = edge
    """ Iterate over the relevant vertices and get their connecting edges. If those
    edges are not an item of near_edges yet and part of the original edges dict given
    to the function, add them to the dictionary.
    """
    return near_edges

def nearest_point_on_edges(vertex, radius=3, edges=None):
    """ Returns the nearest point to a Vertex object on an existing Edge object as Vertex.
    """
    graph = vertex.graph
    if edges is None:
        edges = graph.edges
    else:
        edges = edges
    relevant_edges = near_edges(vertex, edges, radius)
    """ From the edges in the graph, get a new dictionary of edges within the relevant
    distance and divide edges with intermediates into new Edge objects, referencing their
    parent Edge object.
    """
    intersections = {} # dict with edge id as key and Vertex as value
    distances = {} # dict with edge id as key and distance to vertex as value
    for e_id, edge in relevant_edges.items():
        coords = perpendicular_on_edge(vertex, edge)
        """ For every edge get the point on the edge which is perpendicular to the Vertex object
        the nearest point is searched for.
        """
        if coords is not None:
            new_vertex = Vertex(coords[0], coords[1],'vertex')
            new_vertex = graph.vertices.get(new_vertex.id(), new_vertex)
            """ If the nearest point corresponds to an existing Vertex object in in the
            vertices dictionary of the Graph object, point the new_vertex variable to
            this object.
            """
            intersections[e_id] = new_vertex
            distances[e_id] = distance(vertex, new_vertex)
        else:
            dist_s = distance(vertex, edge.start)
            dist_e = distance(vertex, edge.end)
            if dist_s <= dist_e:
                intersections[e_id] = edge.start
            else:
                intersections[e_id] = edge.end
            distances[e_id] = min([dist_s, dist_e])
        """ If there is a perpendicular point on the edge, create a new Vertex object with the
        coordinates of this point. Add this object to the intersecting vertex dictionary. Add the
        distance from the Vertex "vertex" to the new Vertex object to the dictionary distances.
        Else, return the nearest existing start or end point of an edge within the search radius.
        """
    if len(distances) == 0:
        return None
        """ If no Edge object is within the search radius return None.
        """
    key = min(distances, key=distances.get)
    near_vertex = intersections[key]
    near_vertex.sub_edges[key] = relevant_edges[key]
    near_vertex.kind = 'vertex'
    near_vertex.nearest_edge = relevant_edges[key].parent
    graph.add_vertex(near_vertex)
    return near_vertex
    """ Get the key of the new Vertex object with the minimal distance to the Vertex object
    "vertex" and return this object from the dictionary of new Vertex objects.
    """

def nearest_point_to_point(vertex, vertices):
    """ Returns the nearest Vertex object to a specific Vertex object "vertex" out of dictionary
    of Vertex objects "vertices".
    """
    distances = {}
    for v_id, v in vertices.items():
        distances[v_id] = distance(v, vertex)
    key = min(distances, key=distances.get)
    return vertices[key]

def perpendicular_on_edge(vertex, edge):
    """ Calculate the coordinates on an edge, which are perpendicular to the Vertex
    object given to the function.
    """
    beta = edge.angle()
    b = distance(edge.start, vertex) * cos(angle(edge.start, vertex) - beta)
    xi = cos(beta) * b + edge.start.x
    yi = sin(beta) * b + edge.start.y
    if is_in_range(xi, edge.start.x, edge.end.x) and is_in_range(yi, edge.start.y, edge.end.y):
        return (xi, yi)

def is_in_range(coord, start, end):
    """ Function to check if some coordinate is in the range specified.
    """
    return coord >= start and coord <= end

def split_edge_at_point(edge, vertex):
    """ Replace an edge object with two new ones. The first starting at the old start Vertex object
    and ending at the Vertex object. The second one starting at the Vertex object and ending at
    old edge's end Vertex.
    """
    graph = edge.graph
    if vertex.id() in edge.vertices or vertex in edge.intermediates:
        intermediates = split_intermediates(edge, vertex)
        start = edge.start
        end = edge.end
        graph.delete_edge(edge)
        graph.add_edge(start.x, start.y, vertex.x, vertex.y, [(i.x,i.y) for i in intermediates[0]])
        graph.add_edge(vertex.x, vertex.y, end.x, end.y, [(i.x,i.y) for i in intermediates[1]])
        start.update_neighbours()
        vertex.update_neighbours()
        end.update_neighbours()
            
def split_intermediates(edge, vertex):
    """ Return the list of intermediate points of an Edge object divided at Vertex object "vertex"
    as a tuple.
    """
    if len(edge.intermediates) != 0:
        result_1 = []
        result_2 = []
        split_edge = vertex.sub_edges.values()[0]
        if split_edge.end.x == edge.end.x and split_edge.end.y == edge.end.y:
            limit = len(edge.intermediates) + 1
        else:
            limit = edge.intermediates.index((split_edge.end.x, split_edge.end.y))
        return (edge.intermediates[:limit], edge.intermediates[limit:])
    else:
        return ([],[])
    
def min_distance():
    pass

def get_vertices_from_edges(edges, graph):
    """ Returns a dictionary of unique vertices belonging to the edges given in form
    of a dictionary with edge ids as keys and Edge objects as values.
    """
    vertices = {}
    for edge in edges.values():
        for sub_edge in edge.sub_edges:
            for v_id in sub_edge.vertices:
                vertices.setdefault(v_id, graph.vertices[v_id])
    return vertices

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
        """print()
        print('vertex: ', vertex)
        print('neighbours of %s: %s' % (graph.vertices[vertex], graph.vertices[vertex].neighbours))
        print()"""
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
    graph=Graph('graph1')
    graph.add_edge(2,5,1,1, [(1.3,4)])
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

    e1_bs = copy.deepcopy(graph.edges)

    print("dimensions: ", graph.dimensions())
    #print([edge.id() for edge in graph.edges.values()])

    d2 = dijkstra(graph, '011')
    d1 = dijkstra(graph, '011', '011')

    print('distances d2: %s' % d2[0])
    print()

    result = shortest_path(graph,'0108', '011')
    print(result)
    print()

    print('All paths from dijkstra calculation')
    print(all_paths(dijkstra(graph, '011')))
    print()
    
    for edge in graph.edges.values():
        for sub_edge in edge.sub_edges:
            print(sub_edge.id(), "\t", sub_edge.vertices)
    for v_id in graph.vertices:
        print(v_id)

    newVertex = Vertex(0.39,0.8472)
    graph.add_vertex(newVertex)
    #print(near_points(newVertex, graph.vertices, 3))
    #print(near_edges(newVertex, inter_edges(graph.edges), 3))
 
    np = nearest_point_on_edges(newVertex, 3, graph.edges)
    graph.add_vertex(np)
    print('near_point:', np.x, np.y)
    print()
    
    ne = graph.edges[np.nearest_edge]
    split_edge_at_point(ne, np)

    print('All paths from dijkstra calculation')
    print(all_paths(dijkstra(graph, '011')))
    print()

    """ Test run with second graph """
    graph2=Graph('graph2')
    graph2.add_edge(2,5,1,1, [(1.301,4)])
    graph2.add_edge(2,5,7,5)
    graph2.add_edge(2,5,3,6)
    graph2.add_edge(3,3,1,1)
    graph2.add_edge(3,6,7,5)
    graph2.add_edge(0,0,1,1,[(0.25,0.3),(0.5,0.8),(0.75,0.8)])    
    graph2.add_edge(0,0,1,1,[(0.25,0.2),(0.5,0.2),(0.75,0.2)])
    graph2.add_edge(0,2,1,1)    
    graph2.add_edge(1,1,6,1)
    graph2.add_edge(3,3,4,2)
    graph2.add_edge(7,5,6,1)
    graph2.add_edge(7,5,10,8)

    e2_bs = copy.deepcopy(graph2.edges)

    d2 = dijkstra(graph2, '011')
    d1 = dijkstra(graph2, '011', '011')

    print('All paths from dijkstra calculation')
    print(all_paths(dijkstra(graph2, '011')))
    print()

    newVertex2 = Vertex(0.39,0.8472)
    graph2.add_vertex(newVertex2)
    #print(near_points(newVertex, graph2.vertices, 3))
    #print(near_edges(newVertex, inter_edges(graph2.edges), 3))

    np = nearest_point_on_edges(newVertex2, 3, graph2.edges)
    graph2.add_vertex(np)
    print('near_point:', np.x, np.y)
    print()

    ne = graph2.edges[np.sub_edges.values()[0].parent]
    split_edge_at_point(ne, np)

    print('All paths from dijkstra calculation')
    print(all_paths(dijkstra(graph2, '011')))
    print()
    """
    inte1 = inter_edges(e1_bs)
    inte2 = inter_edges(e2_bs)

    print('edges only in inte1:')
    for e_id in inte1:
        if e_id not in inte2:
            print(e_id)
    print()
    print('edges only in inte2:')
    for e_id in inte2:
        if e_id not in inte1:
            print(e_id)
    print()
    v1 = Vertex(0.39,0.8472)
    v10 = nearest_point_on_edges(v1, 3, e1_bs)
    print(v10.x, v10.y)
    v1 = graph.vertices['00.390.8472']
    v2 = graph2.vertices['00.390.8472']

    intersections = {}
    distances = {}
    for e_id, edge in inte1.items():
        print(edge.start.x, "\t", edge.start.y, "\t", edge.end.x, "\t", edge.end.y, "\t", perpendicular_on_edge(v1, edge))
        coords = perpendicular_on_edge(v1, edge)
        if coords is not None:
            new_vertex = Vertex(coords[0], coords[1])
            new_vertex.edges[e_id] = edge
            intersections[e_id] = new_vertex
            distances[e_id] = distance(v1, new_vertex)
        else:
            dist_s = distance(v1, edge.start)
            dist_e = distance(v1, edge.end)
            if dist_s <= dist_e:
                intersections[e_id] = edge.start
            else:
                intersections[e_id] = edge.end
            distances[e_id] = min([dist_s, dist_e])
    print(intersections, distances)
    key = min(distances, key=distances.get)
    print(key)
    nearp = intersections[key]
    print('nearest point: ', '\tx:', nearp.x, '\ty:', nearp.y)
    print()

    rele1 = near_edges(v1, inte1, 3)
    intersections = {}
    distances = {}
    for e_id, edge in rele1.items():
        print(edge.start.x, "\t", edge.start.y, "\t", edge.end.x, "\t", edge.end.y, "\t", perpendicular_on_edge(v1, edge))
        coords = perpendicular_on_edge(v1, edge)
        if coords is not None:
            new_vertex = Vertex(coords[0], coords[1])
            new_vertex.edges[e_id] = edge
            intersections[e_id] = new_vertex
            distances[e_id] = distance(v1, new_vertex)
        else:
            dist_s = distance(v1, edge.start)
            dist_e = distance(v1, edge.end)
            if dist_s <= dist_e:
                intersections[e_id] = edge.start
            else:
                intersections[e_id] = edge.end
            distances[e_id] = min([dist_s, dist_e])
    print(intersections, distances)
    key = min(distances, key=distances.get)
    print(key)
    nearp = intersections[key]
    print('nearest point: ', '\tx:', nearp.x, '\ty:', nearp.y)
    print()

    for edge in inte2.values():
        print(edge.start.x, "\t", edge.start.y, "\t", edge.end.x, "\t", edge.end.y, "\t", perpendicular_on_edge(v2, edge))
    print()

    rele2 = near_edges(v2, inte2, 3)
    for edge in rele2.values():
        print(edge.start.x, "\t", edge.start.y, "\t", edge.end.x, "\t", edge.end.y, "\t", perpendicular_on_edge(v2, edge))
    print()"""
