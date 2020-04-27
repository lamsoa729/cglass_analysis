#!/usr/bin/env python

"""@package docstring
File: sc_helpders.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import numpy as np


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
