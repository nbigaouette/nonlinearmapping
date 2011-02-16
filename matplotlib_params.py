#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 3.0
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['lines.markeredgewidth'] = 2.0
mpl.rcParams['figure.figsize'] = (12,4)

mpl.rcParams['font.size'] = 24.0

#mpl.rcParams['figure.subplot.left']   = 0.07
mpl.rcParams['figure.subplot.left']   = 0.1
mpl.rcParams['figure.subplot.right']  = 0.985
mpl.rcParams['figure.subplot.top']    = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.2

# Fit fonts with LaTeX
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "cm"

mpl.rcParams['text.usetex'] = True


import pylab, sys
###################################################################
def on_key(event):
    #print 'you pressed', event.key, event.xdata, event.ydata
    if   (event.key == 'q'):
        sys.exit(0)
    #
    elif (event.key == 'w'):
        pylab.close(pylab.gcf())
    #
    elif (event.key == 'd'):
        print "Please click two points to get the distance."
        points = pylab.ginput(2)
        print "points = ", points
        print "distance (x) =", abs(points[0][0] - points[1][0])
        print "distance (y) =", abs(points[0][1] - points[1][1])
        print "distance =", pylab.sqrt( (points[0][0] - points[1][0])**2 + (points[0][1] - points[1][1])**2 )
    #
    elif (event.key == 'a'):
        print "Please click a point to get the position."
        #from matplotlib.widgets import Cursor
        #cursor = Cursor(pylab.gca(), useblit=True, color='red', linewidth=2 )
        print pylab.ginput(1)
#
###################################################################

def savefigure(fig, _filename):
    import os
    extensions = ['png', 'svg', 'pdf', 'eps']
    path = "figures"
    if (not os.path.exists(path)):
        os.makedirs(path)
    #
    for ext in extensions:
        filename = os.path.join(path, _filename + "." + ext)
        print "Saving figure to " + filename
        fig.savefig(filename)
    #
#

