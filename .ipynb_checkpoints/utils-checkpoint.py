'''
    Utilities for use in flopy: these are generally implementations of
    ModelMuse methods for handling objects
'''
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull
import numpy as np
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

def surface_interpolation(path, L, delL):
    '''
     Synopsis: This function performs triangular interpolation using object vertices as
     interpolation points. The outer points of the grid are set using nearest point
     on the convex hull

     Inputs:
        path: path to vertex data
        L: size of grid (change to array later)
        delL: discretiztion length (change to array later)

     Outputs:
        top: 2D array of the top value of each cell in the grid

     Author: Connor Cleary, University of Canterbury
    '''
    # import matplotlib.pyplot as plt

    N = int(L/delL)
    # read csv data
    vertices = np.genfromtxt(path, skip_header=1, delimiter=',')

    # create array of points to interpolate
    x = np.linspace(delL/2, L-delL/2, N)
    y = np.linspace(-delL/2, -L+delL/2, N)
    X, Y = np.meshgrid(x, y)

    # perform interpolation
    interp = LinearNDInterpolator(vertices[:,:2], vertices[:,2], fill_value = np.inf)
    top = interp(X, Y)

    # create convex hull: this is faster than creating a hull from the triangulation
    hull = ConvexHull(vertices[:,:2])
    simplices = hull.simplices

    # fill in values outside the convex hull
    for j in range(N):
        for i in range(N):
            if top[j][i] == np.inf:
                mind = np.inf
                # find facet which it is closest to (as well as being between simplices)
                x = (i*delL) + (delL/2)
                y =  (-j*delL) - (delL/2)
                p3 = np.asarray([x, y])

                # for each facet check whether it is the closest
                for (ifac, facet) in enumerate(simplices):
                    p1 = np.asarray(vertices[facet[0]][:2])
                    p2 = np.asarray(vertices[facet[1]][:2])
                    # finding normal of point to facet
                    dx = p2[0] - p1[0]
                    dy = p2[1] - p1[1]
                    d2 = dx*dx + dy*dy
                    nx = ((p3[0]-p1[0])*dx + (p3[1]-p1[1])*dy) / d2
                    # this chooses either end of the facet if point is not inside segment
                    nx = min(1, max(0, nx))
                    # projection
                    p = (dx*nx + p1[0], dy*nx + p1[1])
                    # distance from point to projection
                    d = np.sqrt((p3[0]-p[0])**2 + (p3[1]-p[1])**2)

                    if d < mind:
                        minp = p
                        mind = d
                        imin = ifac
                        minnx = nx

                top[j][i] = (1-minnx)*vertices[simplices[imin][0]][2] + minnx*vertices[simplices[imin][1]][2]

    # plt.pcolormesh(X, Y, top, shading='auto', cmap = 'Greens_r')
    # plt.legend()
    # plt.colorbar()
    # plt.title("Heathcote model surface elevation")
    # plt.axis("equal")
    # plt.show()

    return top