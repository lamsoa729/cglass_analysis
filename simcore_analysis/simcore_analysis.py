#!/usr/bin/env python

"""@package docstring
File: simcore_analysis.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import argparse
from pathlib import Path
import h5py
import yaml

from .sc_parse_data import collect_data
from .sc_analyze_seed import (analyze_seed)
from .sc_analyze_seed_scan import analyze_seed_scan, collect_seed_h5_files
from .sc_analyze_param_scan import collect_param_h5_files
from .sc_analyze_run import analyze_run


def parse_args():
    parser = argparse.ArgumentParser(prog='simcore_analysis.py',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input", default=None,
                        help="Input for Simcore Analysis functions.")
    parser.add_argument("-s", "--seed", action="store_true", default=False,
                        help=("Run a single seed analysis. "
                              "Input is a parameter file name."))
    parser.add_argument("-S", "--seed_scan", action="store_true", default=False,
                        help=("Collect hdf5 data files from a directory tree "
                              "and create a file based on statistics of seed runs. "
                              "Input is the directory that holds seed directories."))
    parser.add_argument("-R", "--run", action="store_true", default=False,
                        help=("Create hdf5 files for both seed and parameter "
                              "scans. All seeds must be analyzed first."
                              "Input is parameter you wish to analyze."))
    parser.add_argument("--spec", type=str, default='',
                        help=(" Specify if parameter used in param scan or "
                              "full run anlaysis is a specific species "
                              "parameter. e.g. 'crosslink' "))
    opts = parser.parse_args()
    return opts


def run_seed_scan_analysis(param_dir_path):
    """!TODO: Docstring for prep_seed_scan_analysis.

    @param param_dir_path: TODO
    @return: TODO

    """
    try:
        if not isinstance(param_dir_path, Path):
            param_dir_path = Path(param_dir_path)
        # Get name of param_dir to make file name later
        name = param_dir_path.name
        file_path = param_dir_path / '{}.h5'.format(name)
        print(file_path)
        if file_path.exists():
            file_path.unlink()
        h5_out = h5py.File(file_path, 'w')
        # Collect seeds to analyze
        h5_data_lst = collect_seed_h5_files(param_dir_path)
        # analyze seeds
        analyze_seed_scan(h5_out, h5_data_lst)

    except BaseException:
        print("Analysis failed")
        raise
    finally:
        h5_out.close()
        for h5d in h5_data_lst:
            h5d.close()


def run_full_tree_analysis(param, spec=None):
    """!Run analysis to collect seed data files and combine into seed scan files
    and full run data files.

    @param sim_dir_path: path to simulation directory
    @return: TODO

    """
    try:
        sim_dir_path = Path('simulations')

        # Run seed scan analysis to consolidate first level of tree
        for pdirs in sim_dir_path.glob('*/'):
            run_seed_scan_analysis(pdirs)

        # Collect seed scan data files
        h5_data_lst = collect_param_h5_files(sim_dir_path, param, spec)

        # TODO: Run analyis on collected data files <15-01-20, ARL> #
        # analyze_run(h5_data_lst, param, spec)

        # name = Path.cwd().name
        # file_path = '{}.h5'.format(name)
        # h5_out = h5py.File(file_path, 'w')
    except BaseException:
        print("Analysis failed")
        raise
    finally:
        # h5_out.close()
        for h5d in h5_data_lst:
            h5d.close()


def run_seed_analysis(param_file=None):
    """!TODO: Docstring for prep_seed_analysis.

    @param param_file: TODO
    @return: TODO

    """
    try:
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
        analyze_seed(h5_data)
    except BaseException:
        print("Analysis failed")
        raise
    finally:
        h5_data.close()


def main():
    """!Main function for simcore_analysis

    Parameters
    ----------
    param_file: Yaml Parameter file

    Returns
    -------
    TODO

    """
    opts = parse_args()
    if opts.seed:
        run_seed_analysis(opts.input)
    elif opts.seed_scan:
        run_seed_scan_analysis(opts.input)
    elif opts.run:
        if opts.spec != '':
            run_full_tree_analysis(opts.input, opts.spec)
        else:
            run_full_tree_analysis(opts.input)

    else:
        raise IOError('No valid analysis type was given.')


##########################################
if __name__ == "__main__":
    main()
