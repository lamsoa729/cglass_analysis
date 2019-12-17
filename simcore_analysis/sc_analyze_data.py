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


def xl_zrl_force(r_i, r_j, u_i, u_j, s_i, s_j, ks):
    """!TODO: Docstring for xl_zrl_force.

    @param r_i: TODO
    @param r_j: TODO
    @param u_i: TODO
    @param u_j: TODO
    @param s_i: TODO
    @param s_j: TODO
    @param ks: TODO
    @return: TODO

    """
    return -ks * (r_j + (u_j * s_j) - r_i - (u_i * s_i))


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


def analyze_xlink_force(h5_data):
    """!Analyze the force on filament_j by crosslinkers attached to filament_i

    @param h5_data: TODO
    @return: TODO

    """
    xl_dbl_dset = h5_data['xl_data/doubly_bound']
    ks = h5_data['xl_data'].attrs['k_spring']
    fil_pos_dset = h5_data['filament_data/filament_position']
    fil_orient_dset = h5_data['filament_data/filament_orientation']
    u_i = fil_orient_dset[:, :, 0]
    u_j = fil_orient_dset[:, :, 1]
    r_i = fil_pos_dset[:, :, 0]
    r_j = fil_pos_dset[:, :, 1]
    xl_si = xl_dbl_dset[:, 0]
    xl_sj = xl_dbl_dset[:, 1]
    force_arr = np.zeros((len(u_i), 3))
    for i in range(len(u_i)):
        for xl in range(len(xl_si[i])):
            force_arr[i, :] += xl_zrl_force(r_i[i], r_j[i],
                                            u_i[i], u_j[i],
                                            xl_si[i][xl], xl_sj[i][xl],
                                            ks)
    xl_force_dset = h5_data['analysis'].create_dataset('xl_forces',
                                                       data=force_arr)