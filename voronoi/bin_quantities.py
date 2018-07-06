#!/usr/bin/env python
# -----------------------------------------------------------------------------
# VORONOI.BIN_QUANTITIES
# - converted from IDL code by Michele Cappellari
#   (bin2d_compute_useful_bin_quantities)
# -----------------------------------------------------------------------------

import numpy as np
from .weighted_centroid import weighted_centroid


def bin_quantities(x, y, signal, noise, xnode, ynode, scale):
    
    """
    Recomputes (weighted) voronoi tessellation of the pixels grid to make
    sure that the clas number corresponds to the proper Voronoi generator.
    This is done to take into account possible zero-size Voronoi bins
    in output from the previous CVT (or WVT).
    
    INPUTS
      x      : x-coordinates of pixels
      y      : y-coordinates of pixels
      signal : signal in pixels
      noise  : noise in pixels
      xnode  : x-coordinates of bins
      ynode  : y-coordinates of bins
      scale  : bin scale
    """
    
    clas = np.zeros(x.size, dtype="int")   # will contain bin num of given pixels
    for j in range(x.size):
        clas[j] = (((x[j]-xnode)/scale)**2 + ((y[j]-ynode)/scale)**2).argmin()
    
    # At the end of the computation evaluate the bin luminosity-weighted
    # centroids (xbar,ybar) and the corresponding final S/N of each bin.
    
    area, lim = np.histogram(clas, bins=int(clas.ptp())+1,
        range=(clas.min()-0.5, clas.max()+0.5))
    cent = (lim[:-1]+lim[1:])/2.
    
    xb = np.zeros(xnode.size)
    yb = np.zeros(xnode.size)
    sn = np.zeros(xnode.size)
    for j in range(xnode.size):
        pix = np.where(clas == cent[j])[0]
        xb[j], yb[j] = weighted_centroid(x[pix], y[pix], signal[pix])
        sn[j] = signal[pix].sum()/np.sqrt((noise[pix]**2).sum())
    
    return clas, xb, yb, sn, area
