#!/usr/bin/env python

"""@package docstring
File: sc_analyze_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
import numpy as np


def flatten_dset(dset):
    """!Flattens a time series of positions of crosslinker heads on filament

    @param dset: TODO
    @return: TODO

    """
    return [item for sublist in dset for item in sublist]


def analyze_xlink_moments(h5_data):
    anal_grp = h5_data['analysis']
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']

    dbl_num_arr = [xld_t.size for xld_t in dbl_xlink_dset[:, 0]]
    xl_zeroth_mom_dset = anal_grp.create_dataset(
        'xl_zeroth_moment', data=dbl_num_arr)

    mu10_arr = [np.sum(xld_t) for xld_t in dbl_xlink_dset[:, 0]]
    mu01_arr = [np.sum(xld_t) for xld_t in dbl_xlink_dset[:, 1]]
    xl_first_mom_arr = np.vstack((mu10_arr, mu01_arr))
    xl_first_mom_dset = anal_grp.create_dataset(
        'xl_first_moments', data=xl_first_mom_arr.T)

    mu20_arr = [np.sum(np.power(xld_t, 2)) for xld_t in dbl_xlink_dset[:, 0]]
    mu02_arr = [np.sum(np.power(xld_t, 2)) for xld_t in dbl_xlink_dset[:, 1]]
    mu11_arr = [np.dot(xld_0, xld_1) for xld_0, xld_1 in
                zip(dbl_xlink_dset[:, 0], dbl_xlink_dset[:, 1])]
    xl_second_mom_arr = np.vstack((mu11_arr, mu20_arr, mu02_arr))
    xl_second_mom_dset = anal_grp.create_dataset(
        'xl_second_moments', data=xl_second_mom_arr.T)


def analyze_avg_xlink_distr(h5_data):
    """!TODO: Docstring for analyze_average_xlink_distr.

    @param h5_data: TODO
    @return: TODO

    """
    length = h5_data['filament_data'].attrs['lengths'][0]
    fil_bins = np.linspace(-.5 * length, .5 * length, 120)

    # Combine all time data to get an average density
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']
    fil0_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 0]))
    fil1_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 1]))
    # print(fil1_lambdas)
    dbl_2D_distr, xedges, yedges = np.histogram2d(
        fil0_lambdas, fil1_lambdas, fil_bins)

    n_steps = h5_data.attrs['n_steps']
    n_spec = h5_data['xl_data'].attrs['n_spec']

    dbl_2D_distr *= float(n_spec / n_steps)

    xl_avg_distr_dset = h5_data['analysis'].create_dataset(
        'average_doubly_bound_distr', data=dbl_2D_distr)
    xl_avg_distr_dset.attrs['xedges'] = xedges
    xl_avg_distr_dset.attrs['yedges'] = yedges
