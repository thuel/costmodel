#!/usr/bin/env python

def near_edge(point, edge):
    slop1 = (y1-y2)/(x1-x2)
    operator = 90
    if slop1 > 0:
    	operator = -90
    slop2 = slop1 + operator
