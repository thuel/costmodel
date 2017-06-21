#!/usr/bin/env python
from math import radians, sin

def sinus_line_points(startx, starty, endx, endy, step=0.25):
    """sinus line points calculation:
    """
    npoints = int((endx - startx)/step) + 1 # number of points to calculate
    xgses = [ i*step for i in range(npoints) ] # start and end to be removed
    # ystep = (endy-starty) / npoints
    # yies = [ i*ystep for i in ]
    for x in xgses:
        yies.append((endy-starty) * (x - sin(radians(x*360))*(step/sin(radians(step*360)))))
    return lists_to_tuples([x + startx for x in xgses ], [y + starty for y in yies])

def lists_to_tuples(list1, list2):
    if len(list1) == len(list2):
        result = []
        for i in range(len(list1)):
            result.append((list1[i], list2[i]))
        return result

if __name__=="__main__":