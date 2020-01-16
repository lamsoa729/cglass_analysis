#!/usr/bin/env python

"""@package docstring
File: sc_analyze_param_scan.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path
import numpy as np
import yaml
import h5py


def get_param_from_dict(h5_data, param, spec=None):
    """!TODO: Docstring for get_param_from_dict.

    @param param: TODO
    @return: TODO

    """
    if spec is None:
        return h5_data.attrs[param]
    else:
        param_dict = yaml.safe_load(h5_data.attrs['param_dict'])
        # A little cludgy because you only check the parameters of the first
        # in the list in a given species.
        return param_dict[spec][0][param]


def collect_param_h5_files(dir_path, param, spec=None):
    """ Spider through directory structure to collect and put h5 files in a list"""
    h5_data_lst = sorted([h5py.File(hf, 'r+')
                          for hf in dir_path.glob('*/*.h5')],
                         lambda x: get_param_from_dict(x, param, spec))
    return h5_data_lst


def analyze_param_scan(sim_dir_path, param, spec=None):
    """!TODO: Docstring for analyze_run.

    @param sim_dir_path: TODO
    @param param: TODO
    @return: TODO

    """
    h5_data_lst = collect_param_h5_files(sim_dir_path, param, spec=None)
    pass


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
