#!/usr/bin/env python

"""@package docstring
File: ot-fix_graphs.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path

import numpy as np
# import powerlaw
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

from .otrap_sty import otrap_sty
from .helpers import get_suff, get_ot_force, collect_contiguous_intervals


# Unit conversions
SQRT_PI = np.sqrt(np.pi)
sec = .0357  # sec
um = .025  # um
pN = .1644  # pN
nm = 25.  # nm
nM = 106382.978723404  # nM for molecular simulations
nM_fp = 1147.604948473  # nM for PDE and ME simulations

#########################
#  Single-seed figures  #
#########################


def graph_fixed_OT_assays(opts):
    with plt.style.context(otrap_sty):
        plt.style.use(otrap_sty)
        if not opts.data_dir.exists():
            opts.data_dir.mkdir()
        if opts.run_type == 'single_seed':
            graph_seed_smoothed_interval(opts.input, opts.data_dir)
            print("Make single seed graphs.")
        elif opts.run_type == 'multi_seed':
            graph_multi_seed_force_intervals(opts.input, opts.data_dir)
            print("Make multi seed graphs")
        elif opts.run_type == 'param_scan':
            graph_param_scan_ot_force(opts)
            param_dirs = sorted(Path(opts.input).glob(opts.param + '*'),
                                key=get_suff)
            for pdir in param_dirs:
                graph_multi_seed_force_intervals(
                    pdir, opts.data_dir,
                    'multi_seed_force_intervals_' + pdir.name)
        else:
            raise IOError("No run type or incorrect run type was given.")


def graph_seed_smoothed_interval(data_path, data_dir):
    """TODO: Docstring for graph_seed_smoothed_interval.

    @param ot_path TODO
    @return: TODO

    """
    try:
        # data_path = ot_path /
        # '20-10-16_CG_OT_scan10.20_stat_rigid_zl-kin5_Edep0.0-1.0/simulations/Edep0.5/s0/ot_test_data.h5'
        h5_data = h5py.File(data_path, 'r')

        time_arr, ot_force = get_ot_force(
            h5_data, species_name='filament-trap')
        dt = time_arr[1] - time_arr[0]
        ot_force_smooth = savgol_filter(-ot_force[:, 1], 45, 3)
        (ot_force_deriv,
         _,
         pos_lengths,
         neg_lengths) = collect_contiguous_intervals(ot_force_smooth, dt)

        fig, axarr = plt.subplots(2, 3, figsize=(18, 10))

        axarr[0, 0].set_title('Raw')
        axarr[0, 0].plot(time_arr * sec, -ot_force[:, 1] * pN)
        # Smooth force using Savitsky-Golay filter with 3rd-order polynomial
        # and window of 45 points (.45 sec)
        axarr[0, 1].set_title('Smoothed \n (Savitsky-Golay p=3, w=.45 sec)')
        axarr[0, 1].plot(time_arr * sec, ot_force_smooth * pN)

        axarr[0, 2].set_title('Derivative of smoothed force')
        axarr[0, 2].plot(time_arr * sec, ot_force_deriv * pN / sec)

        pos_lengths *= dt * sec
        neg_lengths *= dt * sec

        axarr[1, 0].set_title('Increasing force intervals')
        axarr[1, 0].hist(pos_lengths)
        axarr[1, 1].set_title('Decreasing force intervals')
        axarr[1, 1].hist(neg_lengths)

        for ax in axarr.flatten()[:-1]:
            ax.set_xlabel('Time $t$ (sec)')

        axarr[0, 0].set_ylabel('Force $F_x$ (pN)')
        axarr[0, 2].set_ylabel('$dF_x/dt$ (pN/sec)')
        axarr[1, 0].set_ylabel("Number of intervals")

        plt.tight_layout()
        fig.savefig(data_dir / 'smoothed_interval.png')
    except BaseException:
        raise
    finally:
        h5_data.close()

########################
#  Multi-seed figures  #
########################


def graph_multi_seed_force_intervals(param_path, data_dir,
                                     save_name='multi_seed_force_intervals'):
    """TODO: Docstring for graph_param_scan_force_interval.

    @param arg1 TODO
    @return: TODO

    """
    # param_path = ot_path / \
    #     '20-10-16_CG_OT_scan10.20_stat_rigid_zl-kin5_Edep0.0-1.0/simulations/Edep0'
    seed_paths = param_path.glob('s*')

    pos_length_arr = None
    neg_length_arr = None

    fig, axarr = plt.subplots(1, 3, figsize=(18, 5))
    for s in seed_paths:
        try:
            h5_data = h5py.File(s / 'ot_test_data.h5', 'r')
            time, ot_force = get_ot_force(h5_data)

            dt = time[1] - time[0]
            ot_force_smooth = savgol_filter(-ot_force[:, 1], 45, 3)
            axarr[0].plot(time * sec, ot_force_smooth * pN)

            _, _, pos_lengths, neg_lengths = collect_contiguous_intervals(
                ot_force_smooth, dt)
            if pos_length_arr is None:
                pos_length_arr = pos_lengths
            else:
                pos_length_arr = np.concatenate(
                    (pos_length_arr, pos_lengths), axis=None)
            if neg_length_arr is None:
                neg_length_arr = neg_lengths
            else:
                neg_length_arr = np.concatenate(
                    (neg_length_arr, neg_lengths), axis=None)
        except BaseException:
            raise
        finally:
            h5_data.close()

    pos_length_arr *= dt * sec
    neg_length_arr *= dt * sec

    axarr[1].hist(pos_length_arr, bins=20, log=True)
    axarr[2].hist(neg_length_arr, bins=20, log=True)

    # axarr[1].set_title('Increasing force intervals')
    # axarr[2].set_title('Decreasing force intervals')
    axarr[0].set_ylabel("Force (pN)")
    axarr[1].set_ylabel("Number of intervals")
    for ax in axarr:
        ax.set_xlabel("Time (sec)")

    plt.tight_layout()
    fig.savefig(data_dir / (save_name + '.png'))

########################
#  Param scan figures  #
########################


def param_scan_force_plot(axarr, param_dirs, param_name=''):
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for param_dir, color in zip(param_dirs, color_cycle):
        seed_paths = param_dir.glob('s*')
        label_flag = True
        ot_force_arr = []
        for sd in seed_paths:
            try:
                h5_data = h5py.File(sd / 'ot_test_data.h5', 'r')
                time_arr, ot_force = get_ot_force(h5_data)
                ot_force_arr += [ot_force[:, 1]]
                if label_flag:
                    label = "{} = {}".format(
                        param_name, get_suff(param_dir))
                    axarr[1].plot(time_arr * sec, -ot_force[:, 1]
                                  * pN, c=color, label=label)
                    label_flag = False
                    print(label)
                else:
                    pass
            except BaseException:
                raise
                # ax.plot(time*sec, -ot_force[:,1]*pN, c=color)
            finally:
                h5_data.close()
        ot_force_arr = np.asarray(ot_force_arr)
        axarr[0].plot(time_arr * sec,
                      - (ot_force_arr.mean(axis=0)) * pN,
                      c=color)

    axarr[0].set_xlabel('Time $t$ (sec)')
    axarr[1].set_xlabel('Time $t$ (sec)')
    axarr[0].set_ylabel('Force $F_x$ (pN)')
    axarr[1].legend(loc='center left', bbox_to_anchor=(1.05, .5))


def graph_param_scan_ot_force(opts):
    """TODO: Docstring for graph_param_scan_ot_force.

    @param opts TODO
    @return: TODO

    """
    if opts.param is None:
        raise IOError("No parameter was given to graph parameter scan.")

    param_dirs = sorted(Path(opts.input).glob(opts.param + '*'), key=get_suff)
    fig, axarr = plt.subplots(1, 2, figsize=(14, 5))
    param_scan_force_plot(axarr, param_dirs, opts.param)
    plt.tight_layout()
    fig.savefig(opts.data_dir / 'param_scan_force.png')


# End #########################################
