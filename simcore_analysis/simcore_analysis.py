#!/usr/bin/env python

"""@package docstring
File: simcore_analysis.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import sys
import h5py
import yaml
import numpy as np
from .sc_parse_data import collect_data
from .sc_analyze_data import (analyze_xlink_moments, analyze_avg_xlink_distr,
                              analyze_xlink_force, analyze_singly_bound_xlinks)


def run_analysis(h5_data):
    if 'analysis' in h5_data:
        del h5_data['analysis']  # Start clean
    analysis_grp = h5_data.create_group('analysis')
    # analyze xlinks
    analyze_xlink_moments(h5_data)
    analyze_avg_xlink_distr(h5_data)
    analyze_singly_bound_xlinks(h5_data)
    analyze_xlink_force(h5_data)
    # analyze filaments (maybe)


def main(param_file=None):
    """!TODO: Docstring for main.

    Parameters
    ----------
    param_file: Yaml Parameter file

    Returns
    -------
    TODO

    """
    try:
        if param_file is None:
            param_file = sys.argv[1]
        with open(param_file, 'r') as pf:
            p_dict = yaml.safe_load(pf)
        print(p_dict)
        run_name = p_dict['run_name']
        xl_name = p_dict['crosslink'][0]['name']
        fil_name = p_dict['rigid_filament'][0]['name']
        h5_data = h5py.File(run_name + '_data.h5', 'w')
        collect_data(h5_data,
                     run_name + '_params.yaml',
                     run_name + '_crosslink_' + xl_name + '.spec',
                     run_name + '_rigid_filament_' + fil_name + '.posit')
        run_analysis(h5_data)
    except BaseException:
        print("Analysis failed")
        raise
    finally:
        h5_data.close()


##########################################
if __name__ == "__main__":
    main()
