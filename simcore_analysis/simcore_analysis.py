#!/usr/bin/env python

"""@package docstring
File: simcore_analysis.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import sys
import h5py
import numpy as np
from .sc_parse_data import collect_data
from .sc_analyze_data import (analyze_xlink_moments, analyze_avg_xlink_distr,
                              analyze_xlink_force)


def run_analysis(h5_data):
    if 'analysis' in h5_data:
        del h5_data['analysis']  # Start clean
    analysis_grp = h5_data.create_group('analysis')
    # analyze xlinks
    analyze_xlink_moments(h5_data)
    analyze_avg_xlink_distr(h5_data)
    analyze_xlink_force(h5_data)
    # analyze filaments (maybe)


def main(args):
    """!TODO: Docstring for main.

    Parameters
    ----------
    h5_file : TODO

    Returns
    -------
    TODO

    """
    try:
        h5_data = h5py.File(args[1], 'w')
        collect_data(h5_data, args[2], args[3], args[4])
        run_analysis(h5_data)
    except BaseException:
        print("Analysis failed")
        raise
    finally:
        h5_data.close()


##########################################
if __name__ == "__main__":
    main(sys.argv)
