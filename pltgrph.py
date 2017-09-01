#!/usr/bin/env python

"""Print Network graph library
"""

#import essential modules, libraries and methods/functions
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *
from Tkinter import *

zoom = 1
mirror = 150
border = 10

def print_graph(graph):
    global zoom
    width, height = graph.dimensions()
    zoom_factor((width, height))
    window = init_canvas(width * zoom + border * 2, height * zoom + border * 2)
    lines = [graph.edges[i] for i in graph.edges]
    add_lines_from_list(window, lines)
    points = [i for i in graph.vertices.values() if i.kind != 'helper']
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
    global mirror, border
    mirror = h - border * 2 + 2
    master = Tk()
    window = Canvas(master, width=w, height=h)
    window.pack()
    return window

def add_line(can, edge):
    """Add a line representing an Edge() object to canvas can.
    """
    global zoom, mirror, border
    axay = [(edge.start.x, edge.start.y)]
    zxzy = [(edge.end.x, edge.end.y)]
    interlst = [(v.x,v.y) for v in edge.intermediates]
    points = axay + interlst + zxzy
    for i in range(len(points)-1):
        can.create_line(points[i][0] * zoom + border, mirror - points[i][1] * zoom + border, points[i+1][0] * zoom + border, mirror - points[i+1][1] * zoom + border) 

def add_lines_from_list(can, lst):
    for line in lst:
        add_line(can, line)

def add_point(can, vertex):
    """Add a line with starting coordinates ax/ay and ending coordinates
    zx/zy to canvas can.
    """
    global zoom, mirror, border
    t = 2 # line thickness
    can.create_oval(vertex.x * zoom - t + border, mirror - (vertex.y * zoom - t) + border, vertex.x * zoom + t + border, mirror - (vertex.y * zoom + t) + border, fill="black")

def add_points_from_list(can, lst):
    for point in lst:
        add_point(can, point)

if __name__=="__main__":
    window = init_canvas(300, 150)

    mainloop()
