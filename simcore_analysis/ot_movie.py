#!/usr/bin/env python

import sys
import re

from pathlib import Path
# Analysis
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from math import *
# from spindle_unit_dict import SpindleUnitDict
# Image stuff
import cv2
"""@package docstring
File: ot_movie.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""


def make_generic_movie(work_dir, frame_dir='images/', name='ot_movie.mp4',
                       fps=60.0,):

    fps = float(fps)
    # uc = SpindleUnitDict()

    work_dir = Path(work_dir)
    # Init text writing
    # font = cv2.FONT_HERSHEY_SIMPLEX

    mov_path = work_dir / name
    if mov_path.exists():
        mov_path.unlink()

    # Get necessary directory path names
    frame_dir = work_dir / frame_dir

    # Format for video
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')

    # Make list of all the frames in the frame directory
    frame_list = [
        f for f in frame_dir.glob('*.bmp')]
    frame_num_list = [int(re.findall(r'\d+', str(f))[0]) for f in frame_list]
    # Make sure frame list is properly sorted based on frame number

    vid = None
    size = None
    fps = 60.0  # must be a float

    for f, n in zip(frame_list, frame_num_list):
        img = cv2.imread(str(frame_dir / f))
        # Find the size of these frames if not defined
        if size is None:
            size = img.shape[1], img.shape[0]
        # Resize image to correct size if size was defined
        elif size[0] != img.shape[1] and size[1] != img.shape[0]:
            img = cv2.resize(img, size)

        # Initialize video writer at the beginning
        if vid is None:
            vid = cv2.VideoWriter(
                str(mov_path),
                1900,
                fourcc=fourcc,
                fps=fps,
                frameSize=list(size))
        print("-- Frame number: {} --".format(n))
        # num += 1
        vid.write(img)
    return vid


##########################################
if __name__ == "__main__":
    if len(sys.argv) == 1:
        make_generic_movie(Path.cwd())
    elif len(sys.argv) == 2:
        make_generic_movie(Path.cwd(), sys.argv[1])
