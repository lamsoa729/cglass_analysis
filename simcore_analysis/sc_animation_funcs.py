#!/usr/bin/env python

"""@package docstring
File: sc_animation_funcs.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, MovieWriter
from pathlib import Path

from .sc_graphs import sc_graph_all_data_2d


def make_sc_animation(sd_data, writer='ffmpeg',
                      save_path=Path('./'), save_name=None):
    """!Make animation of time slices
    @return: TODO

    """
    fig = plt.figure(constrained_layout=True, figsize=(10, 10))
    graph_stl = {
        "axes.titlesize": 18,
        "axes.labelsize": 15,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "font.size": 15,
        "font.sans-serif": 'Helvetica',
        "text.usetex": False,
        'mathtext.fontset': 'cm',
    }
    with plt.style.context(graph_stl):
        plt.style.use(graph_stl)
        # fig, axarr = plt.subplots(1, 2, figsize=(10, 4))
        gs = fig.add_gridspec(2, 2)
        axarr = np.asarray([fig.add_subplot(gs[0, 0]),
                            fig.add_subplot(gs[0, 1]),
                            fig.add_subplot(gs[1, 0]),
                            fig.add_subplot(gs[1, 1]),
                            ])
        fig.suptitle(' ')
        nframes = sd_data.time.size
        frame_list = range(0, nframes, int(nframes / 200))
        print("  Number of frames =", nframes)
        # print(" Fig dpi =", fig.dpi)
        t0 = time.time()
        anim = FuncAnimation(
            fig,
            sc_graph_all_data_2d,
            frames=frame_list,
            fargs=(fig, axarr, sd_data),
            interval=50,
            blit=True)

    anim.save('{}_anim.mp4'.format(sd_data.run_name), writer=writer)
    t1 = time.time()
    print("Movie saved in: ", t1 - t0)
    return anim


def make_sc_animation_min(sd_data, writer='ffmpeg',
                          save_path=Path('./'), save_name=None):
    """!Make animation of time slices
    @return: TODO

    """
    graph_stl = {
        "axes.titlesize": 18,
        "axes.labelsize": 15,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "font.size": 18,
        "font.sans-serif": 'Helvetica',
        "text.usetex": False,
        'mathtext.fontset': 'cm',
    }
    # with plt.style.context(graph_stl):
    plt.style.use(graph_stl)
    # fig, axarr = plt.subplots(1, 2, figsize=(10, 4))
    fig = plt.figure(constrained_layout=True, figsize=(10, 5))
    gs = fig.add_gridspec(1, 2)
    axarr = np.asarray([fig.add_subplot(gs[0]),
                        fig.add_subplot(gs[1]),
                        ])
    fig.suptitle(' ')
    nframes = sd_data.time.size
    frame_list = range(0, nframes, int(nframes / 100))
    print("  Number of frames =", nframes)
    t0 = time.time()
    anim = FuncAnimation(
        fig,
        sc_graph_all_data_2d,
        frames=frame_list,
        fargs=(fig, axarr, sd_data),
        interval=50,
        blit=True)

    if save_name:
        anim.save(save_path / '{}_min.mp4'.format(save_name), writer=writer)
    else:
        anim.save(save_path / '{}_min.mp4'.format(sd_data.run_name),
                  writer=writer)
    t1 = time.time()
    print("Movie saved in: ", t1 - t0)
    return anim
