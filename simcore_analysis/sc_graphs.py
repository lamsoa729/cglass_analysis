#!/usr/bin/env python

"""@package docstring
File: sc_graphs.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def graph_avg_xlink_distr(h5_data, fig, ax):
    """!TODO: Docstring for graph_avg_xlink_distr.

    Parameters
    ----------
    h5_data : TODO

    Returns
    -------
    TODO

    """
    # h5_data = h5py.File('test_data.h5','r+')
    def flatten_dset(l): return [item for sublist in l for item in sublist]
    length = h5_data['filament_data'].attrs['lengths'][0]
    fil_bins = np.linspace(-.5 * length, .5 * length, 120)

    #fil_bins += (fil_bins[1]-fil_bins[0])*.5
    print(fil_bins)
    # Combine all time data to get an average density
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']
    print(dbl_xlink_dset.shape[0])
    #fil1_lambdas = [item for item in dbl_xlink_dset[i,0].tolist() for i in range(dbl_xlink_dset.shape[0])]
    #fil1_lambdas = [item for lst in dbl_xlink_dset[:,0] for item in lst]
    fil0_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 0]))
    fil1_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 1]))
    # print(fil1_lambdas)
    dbl_2D_distr, xedges, yedges = np.histogram2d(
        fil0_lambdas, fil1_lambdas, fil_bins)
    print(dbl_2D_distr)
    ax.set_aspect('equal')
    ax.pcolormesh(fil_bins, fil_bins, dbl_2D_distr.T)
