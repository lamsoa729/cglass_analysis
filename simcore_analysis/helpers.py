#!/usr/bin/env python

"""@package docstring
File: helpers.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import re
import numpy as np


def contiguous_regions(condition):
    """By Joe Kington, Finds contiguous True regions of the boolean array
    "condition". Returns a 2D array where the first column is the start
    index of the region and the second column is the end index.
    http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array"""

    # Find the indices of changes in "condition"
    diff = np.diff(condition.astype(int))
    idx, = diff.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll  shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # if the end of the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # if the end of condtion is True, append the length of the array
        idx = np.r_[idx, condition.size]

    idx.shape = (-1, 2)
    return idx


def get_ot_force(h5_data, species_name='filament-trap'):
    """ Get force on bead from optical based on position"""
    k_spring = h5_data['optical_trap_data/' +
                       species_name].attrs['trap_spring']
    trap_pos_arr = h5_data['optical_trap_data/' +
                           species_name + '/trap_position'][:, :, 0]
    bead_pos_arr = h5_data['optical_trap_data/' +
                           species_name + '/bead_position'][:, :, 0]
    time = h5_data['optical_trap_data/' + species_name + '/time'][...]
    return time, k_spring * (trap_pos_arr - bead_pos_arr)


def get_suff(path):
    """ Get the suffix from a parameter directory path. """
    scinote = re.compile(r'[+-]?\d+(?:\.?\d*(?:[eE][+-]?\d+)?)?')
    patterns = re.findall(scinote, str(path))
    return float(patterns[-1])  # only return the last one for now


def collect_contiguous_intervals(arr, delta):
    """ Collect and return different contiguous regions of an array. """
    arr_deriv = np.gradient(arr, delta)

    # Get a matrix of indices. First column is starting point of positive
    # interval, second column is ending point of positive interval
    contig_idx_arr = contiguous_regions(arr_deriv > 0)
    start_pos = contig_idx_arr[:, 0]
    end_pos = contig_idx_arr[:, 1]
    pos_lengths = (end_pos - start_pos) * 1.
    neg_lengths = (start_pos[1:] - end_pos[:-1]) * 1.

    return arr_deriv, contig_idx_arr, pos_lengths, neg_lengths

# End #########################################
