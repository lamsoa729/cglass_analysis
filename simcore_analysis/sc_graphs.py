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
from matplotlib.lines import Line2D
from matplotlib.patches import (Circle, RegularPolygon, FancyArrowPatch,
                                ArrowStyle)

from .sc_analyze_seed import flatten_dset


# def flatten_dset(l): return [item for sublist in l for item in sublist]

def convert_size_units(d, ax, reference='y'):
    """
    convert a linewidth in data units to linewidth in points.

    parameters
    ----------
    linewidth: float
        linewidth in data units of the respective reference-axis
    axis: matplotlib axis
        the axis which is used to extract the relevant transformation
        data (data limits and size must not change afterwards)
    reference: string
        the axis that is taken as a reference for the data width.
        possible values: 'x' and 'y'. defaults to 'y'.

    returns
    -------
    d: float
        linewidth in points
    """
    fig = ax.get_figure()
    if reference == 'x':
        length = fig.bbox_inches.width * ax.get_position().width
        value_range = np.diff(ax.get_xlim())[0]
    elif reference == 'y':
        length = fig.bbox_inches.height * ax.get_position().height
        value_range = np.diff(ax.get_ylim())[0]
    # convert length to points
    length *= 72
    # scale linewidth to value range
    return d * (length / value_range)


class LineDataUnits(Line2D):
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop("linewidth", 1)
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72. / self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data)) - trans((0, 0))) * ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)


def xlink_end_pos(r_vec, u_vec, s):
    return (r_vec + (u_vec * s))


def get_max_min_ends(r_i, r_j, u_i, u_j, L_i, L_j):
    """!TODO: Docstring for get_max_min_ends.

    @param r_i: TODO
    @param r_j: TODO
    @param u_i: TODO
    @param u_j: TODO
    @param L_i: TODO
    @param L_j: TODO
    @return: TODO

    """
    return [np.amax(0.5 * L_i * u_i + r_i), np.amin(0.5 * L_i * u_i + r_i),
            np.amax(-.5 * L_i * u_i + r_i), np.amin(-.5 * L_i * u_i + r_i),
            np.amax(0.5 * L_j * u_j + r_j), np.amin(0.5 * L_j * u_j + r_j),
            np.amax(-.5 * L_j * u_j + r_j), np.amin(-.5 * L_j * u_j + r_j)]


def draw_rod(ax, r_vec, u_vec, L, rod_diam, color='tab:green'):
    line = LineDataUnits((r_vec[1] - .5 * L * u_vec[1],
                          r_vec[1] + .5 * L * u_vec[1]),
                         (r_vec[2] - .5 * L * u_vec[2],
                          r_vec[2] + .5 * L * u_vec[2]),
                         linewidth=rod_diam, solid_capstyle='round',
                         color=color, clip_on=False, )

    tip = Circle((r_vec[1] + .5 * L * u_vec[1], r_vec[2] + .5 * L * u_vec[2]),
                 .5 * rod_diam, color='r', zorder=3)
    ax.add_patch(tip)
    ax.add_line(line)


def draw_xlink(ax, e_i, e_j, lw=10, color='k', alpha=.5):
    line = LineDataUnits((e_i[1], e_j[1]), (e_i[2], e_j[2]),
                         linewidth=lw,  # solid_capstyle='round',
                         color=color, clip_on=False, alpha=alpha)
    ax.add_line(line)


def graph_2d_rod_diagram(ax, sd_data, n=-1):
    """!TODO: Docstring for graph_2d_rod_diagram.

    @param ax: TODO
    @param anal: TODO
    @param n: TODO
    @return: TODO

    """
    # params = anal._params
    L_i, L_j = sd_data.h5_data['filament_data'].attrs['lengths']
    # L_j = params["L2"]
    lw = 1.
    r_i_arr = sd_data.r_i_arr
    r_j_arr = sd_data.r_j_arr
    u_i_arr = sd_data.u_i_arr
    u_j_arr = sd_data.u_j_arr

    draw_rod(ax, r_i_arr[n], u_i_arr[n], L_i, lw, color='tab:green')
    draw_rod(ax, r_j_arr[n], u_j_arr[n], L_j, lw, color='tab:purple')

    # Get all extreme positions of tips in the first dimension to maintain
    # consistent graphing size
    x_ends = get_max_min_ends(
        r_i_arr[:, 1], r_j_arr[:, 1], u_i_arr[:, 1], u_j_arr[:, 1], L_i, L_j)
    # Get all extreme positions of tips in the second dimension to maintain
    # consistent graphing size
    y_ends = get_max_min_ends(
        r_i_arr[:, 2], r_j_arr[:, 2], u_i_arr[:, 2], u_j_arr[:, 2], L_i, L_j)

    max_x = max(x_ends + y_ends)
    max_x = max_x * 1.25 if max_x > 0 else .75 * max_x
    min_x = min(x_ends + y_ends)
    min_x = min_x * 1.25 if min_x < 0 else .75 * min_x

    # Make a square box always
    max_y = max_x
    min_y = min_x

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_xlabel(r'x (nm)')
    ax.set_ylabel(r'y (nm)')
    labels = ["fil$_i$", "fil$_j$", "Plus-end"]
    # if anal.OT1_pos is not None or anal.OT2_pos is not None:
    #     labels += ["Optical trap", "Bead"]
    ax.legend(labels, loc="upper right")


def sc_graph_all_data_2d(n, fig, axarr, sc_data):
    print("frame =", n)
    # Clean up if lines
    if not sc_data.init_flag:
        for ax in axarr.flatten():
            ax.clear()
        for artist in fig.gca().lines + fig.gca().collections:
            artist.remove()
            del artist

    graph_2d_rod_diagram(axarr[0], sc_data, n)

    if sc_data.init_flag:
        axarr[0].set_aspect(1.0)
        axarr[1].set_aspect(1.0)
        sc_data.init_flag = False

    return fig.gca().lines + fig.gca().collections


def graph_avg_xlink_distr(h5_data, fig, ax):
    """!TODO: Docstring for graph_avg_xlink_distr.

    Parameters
    ----------
    h5_data : TODO
    fig : TODO
    ax : TODO

    Returns
    -------
    TODO

    """
    length = h5_data['filament_data'].attrs['lengths'][0]
    fil_bins = np.linspace(-.5 * length, .5 * length, 120)

    print(fil_bins)
    # Combine all time data to get an average density
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']
    print(dbl_xlink_dset.shape[0])
    fil0_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 0]))
    fil1_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 1]))
    dbl_2D_distr, xedges, yedges = np.histogram2d(
        fil0_lambdas, fil1_lambdas, fil_bins)
    print(dbl_2D_distr)
    ax.set_aspect('equal')
    cf = ax.pcolormesh(fil_bins, fil_bins,
                       dbl_2D_distr.T / dbl_xlink_dset.shape[0])
    fig.colorbar(cf, ax=ax)
