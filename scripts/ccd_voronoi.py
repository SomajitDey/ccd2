#!/usr/bin/env python3
# Creates Periodic Voronoi tesellations from a given set of points
# This uses the Python module Freud, which in turn uses Voro++. Freud is multithreaded.

import sys
import numpy as np
import freud

pts_in = 'points.in' # Input file containing XY coordinates of the points
cells_out = 'cells.out' # Output file containing Voronoi cell vertices
hexop_out = 'psi6.out' # Output file containing psi6 (hexatic order parameter) for each voronoi cell

points = np.loadtxt(pts_in) # Loading the input points (2D)
npoints = np.size(points, 0) # Number of points
z = np.zeros((npoints, 1)) # Prepping null z coordinates for feeding to Freud
points = np.hstack((points, z)) # Turn input 2D points into 3D with z=0

boxlen = float(sys.argv[1]) # Access box length from first command-line argument

# Periodic Voronoi construction using Freud
box = freud.box.Box.square(boxlen)
voro = freud.locality.Voronoi()
voro.compute((box, points))
cells = voro.polytopes
cells = np.delete(cells, -1, 2) # Delete the extraneous z coordinates

# Output Voronoi cell vertices. Cells are separated by a blank line.
f = open(cells_out, 'w')
for cell in cells:
    np.savetxt(f, cell)
    f.write("\n") # Blank line to separate cells
f.close()

# Compute hexatic order
nlist = voro.nlist
hex_order = freud.order.Hexatic(k=6)
hex_order.compute(system=(box, points), neighbors=nlist)
np.savetxt(hexop_out, hex_order.particle_order)
