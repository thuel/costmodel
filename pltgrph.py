#!/usr/bin/env python

"""Print Network graph library
"""

#import essential modules, libraries and methods/functions
from Tkinter import *
from netwrkgrph import *

zoom = 1
mirror = 150

def print_graph(graph):
    global zoom
    width, height = graph.dimensions()
    zoom_factor((width, height))
    window = init_canvas(width * zoom, height * zoom)
    lines = [graph.edges[i] for i in graph.edges]
    add_lines_from_list(window, lines)
    points = [graph.vertices[i] for i in graph.vertices]
    add_points_from_list(window, points)
    mainloop()

def zoom_factor(dimensions, w_max=500, h_max=500):
    """Adapt dimensions of a graph to fit proportionally to the
    the given maximum width and height. Eventually this could be set in
    a variable.
    """
    global zoom
    width, height = dimensions
    if width > w_max and height > h_max:
        zoom = min(width / w_max, height / h_max)
    else:
        zoom = min(w_max / width, h_max / height)
    

def init_canvas(w=300, h=150):
    """Initialize a canvas with width w and height h.
    """
    global mirror
    mirror = h + 2
    master = Tk()
    window = Canvas(master, width=w, height=h)
    window.pack()
    return window

def add_line(can, edge):
    """Add a line with starting coordinates ax/ay and ending coordinates
    zx/zy to canvas can.
    """
    global zoom, mirror
    can.create_line(edge.start.x * zoom, mirror - edge.start.y * zoom, edge.end.x * zoom, mirror - edge.end.y * zoom)

def add_lines_from_list(can, lst):
    for line in lst:
        add_line(can, line)

def add_point(can, vertice):
    """Add a line with starting coordinates ax/ay and ending coordinates
    zx/zy to canvas can.
    """
    global zoom, mirror
    t = 2 # line thickness
    can.create_oval(vertice.x * zoom - t, mirror - (vertice.y * zoom - t), (vertice.x * zoom + t), mirror - (vertice.y * zoom + t), fill="black")

def add_points_from_list(can, lst):
    for point in lst:
        add_point(can, point)

if __name__=="__main__":
    window = init_canvas(300, 150)
    lines = [Edge(Vertice(20,50),Vertice(0,0)),Edge(Vertice(20,50),Vertice(10,10)),Edge(Vertice(20,50),Vertice(70,50))]
    add_lines_from_list(window, lines)

    mainloop()
