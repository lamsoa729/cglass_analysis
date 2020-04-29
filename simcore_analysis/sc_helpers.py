#!/usr/bin/env python

"""@package docstring
File: sc_helpers.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import numpy as np

SQRT_PI = np.sqrt(np.pi)
sec = .0358  # sec
um = .025  # um
pN = .1644  # pN
nm = 25.  # nm
nM = 106382.978723404  # nM for molecular simulations
nM_fp = 1147.604948473  # nM for PDE and ME simulations


def make_pde_dict_from_sc_h5(h5_data):
    p_dict = {}
    xl_grp = h5_data['xl_data']
    fil_grp = h5_data['filament_data']
    p_dict['L'] = fil_grp.attrs['length'] * nm
    p_dict['ks'] = xl_grp.attrs['k_spring'] * pN / nm
    p_dict['fs'] = xl_grp.attrs['f_stall'] * pN
    p_dict['ko'] = 2. * xl_grp.attrs['k_off_d'] / sec
    p_dict['co'] = xl_grp.attrs['concentration'] * nM / nM_fp
    # (xl_grp.attrs['k_on_s'] * xl_grp.attrs['k_on_d'] *
    #                xl_grp.attrs['bind_site_density']**2 *
    #                xl_grp.attrs['concentration'] /
    #                (xl_grp.attrs['k_off_s'] *
    # xl_grp.attrs['k_off_d'])) * (nM / (nM_fp*nm*nm) #TODO might want to
    # check this
    p_dict['vo'] = xl_grp.attrs['velocity_d'] * nm / sec
    p_dict['beta'] = 0.243309002
    return p_dict


def find_start_time(arr, reps=1):
    """! A function to find when simulations reaches a steady state with
    respect to array, arr.

    @param arr: Array to find steady state in
    @param reps: repetitions of recursion
    @return: st Start time, the index of time array when the simulation
    first reaches a the steady state average

    """
    # Test to make sure correct parameters types were given to function
    if not isinstance(arr, np.ndarray):
        raise TypeError(" Array arr must be numpy.ndarray type ")
    if reps > 0:
        start_time = find_start_time(arr - arr.mean(), reps - 1)
    else:
        # Get array of sign values, ie. sign with respect to mean
        sign_arr = np.sign(arr)
        # Create array of differences from one index to the next
        diff_arr = np.diff(sign_arr)
        # Find the non-zero differences and record the indices
        index_arr = np.where(diff_arr)[0]  # always produces a tuple
        if index_arr.size == 0:  # System was in steady state all along
            start_time = 0
        else:
            start_time = index_arr[0]
    return start_time


##########################################
