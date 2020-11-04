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
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

from .sc_parse_data import collect_data, get_cpu_time_from_log
from .sc_analyze_seed import (analyze_seed)
from .sc_analyze_seed_scan import analyze_seed_scan, collect_seed_h5_files
from .sc_analyze_param_scan import collect_param_h5_files
from .sc_analyze_run import analyze_run
from .sc_seed_data import SeedData
from .sc_animation_funcs import make_sc_animation_min
from .ot_fix_graphs import graph_fixed_OT_assays


def parse_args():
    parser = argparse.ArgumentParser(
        prog='simcore_analysis.py',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("input", default=None,
                        help="Path used in Simcore Analysis functions.")

    parser.add_argument("-a ", "--analysis", type=str, default=None,
                        help=" Specify analysis type to determine if data will"
                        " be overwritten. Options include "
                        "(overwrite, analyze(default), or load.")
    parser.add_argument("-m", "--movie", action="store_true", default=False,
                        help=("Create an animation from a seed."))
    parser.add_argument("-g", "--graph", action="store_true", default=False,
                        help=("Create graph of a seed's end state."))

    parser.add_argument(
        "-T", "--run_type", type=str, default="single_seed",
        help=(
            "simcore_analysis can analyze multiple simulations and aggregate "
            "data according to various schemes. The 'run_type' argument "
            "specifies how to collect the data from nested data directories.\n"
            "Options:"
            "\tsingle_seed\n"
            "\tmulti_seed\n"
            "\tparam_scan\n"
            "\tparam_single_scan\n"
        ))
    parser.add_argument(
        "-A", "--assay_type", type=str, default="fixed-OT",
        help=(
            "Different experimental assays require different analysis. "
            "The 'assay_type' specifies what analysis to run on simulations."
        ))
    parser.add_argument(
        "-P", "--param", type=str, default=None,
        help=("Parameter to use in analysis or graphing functions."))
    parser.add_argument("--spec", type=str, default='',
                        help=("Specify if parameter used in param scan or "
                              "full run anlaysis is a specific species "
                              "parameter. e.g. 'crosslink' "))
    opts = parser.parse_args()

    # Post parsing changes to options
    opts.input = Path(opts.input)
    return opts


def run_full_tree_analysis(param, spec=None, analysis_type='analyze'):
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


def run_seed_scan_analysis(param_dir_path, analysis_type='analyze'):
    """!TODO: Docstring for prep_seed_scan_analysis.

    @param param_dir_path: TODO
    @return: TODO

    """
    if not isinstance(param_dir_path, Path):
        param_dir_path = Path(param_dir_path)
    # Get name of param_dir to make file name later
    name = param_dir_path.name
    h5_file = param_dir_path / '{}.h5'.format(name)
    print(h5_file)
    if not h5_file.exists() and analysis_type == 'load':
        print("!!! {} does not exist when trying to load !!!")
    if h5_file.exists():  # always delete this file. Can only analyze if
        h5_file.unlink()
    try:
        h5_out = h5py.File(h5_file, 'a')
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


def run_seed_analysis(param_file=None, analysis_type='analyze'):
    """!TODO: Docstring for prep_seed_analysis.

    @param param_file: TODO
    @return: TODO

    """
    if analysis_type is None:
        print("Not running analysis.")
        return

    with open(param_file, 'r') as pf:
        p_dict = yaml.safe_load(pf)
    run_name = p_dict['run_name']
    h5_file = Path(run_name + '_data.h5')

    if not h5_file.exists() and analysis_type == 'load':
        print("ANALYSIS: !!! {} does not exist when trying to load !!!".format(
            h5_file))
        return
    if h5_file.exists() and analysis_type == 'overwrite':
        print("ANALYSIS: Overwriting current h5 data")
        h5_file.unlink()

    try:
        h5_data = h5py.File(h5_file, 'a')
        if analysis_type != 'load' and ('xl_data' not in h5_data
                                        or 'filament_data' not in h5_data):
            print("ANALYSIS: Collecting data")
            collect_data(h5_data, run_name + '_params.yaml')
        print("ANALYSIS: Analyzing data")
        analyze_seed(h5_data)
        # Get run time statistics if they exist
        time_anal_flag = p_dict.get('time_analysis', False)
        if time_anal_flag:
            try:
                cpu_time = get_cpu_time_from_log(Path(run_name + '.log'))
                h5_data['analysis'].attrs['cpu_time'] = cpu_time
            except BaseException:
                print("ANALYSIS: !!! Could not collect time analysis !!!")
    except BaseException:
        print("ANALYSIS: failed")
        raise
    finally:
        h5_data.close()


def make_animation(param_file):
    """!TODO: Docstring for run_make_animation.
    @return: TODO

    """
    # try:
    sd_data = SeedData(param_file)
    Writer = FFMpegWriter
    writer = Writer(fps=25, metadata=dict(artist='Me'), bitrate=1800)
    make_sc_animation_min(sd_data, writer)


def run_analysis(opts):
    if opts.run_type == 'single_seed':
        run_seed_analysis(opts.input, opts.analysis)
        # graph_single_seed(opts.input, opts.graph)
    elif opts.run_type == 'multi_seed':
        run_seed_scan_analysis(opts.input, opts.analysis)
        # graph_multi_seed(opts.input, opts.graph)
    elif opts.run_type == 'param_scan':
        if opts.spec != '':
            run_full_tree_analysis(opts.input, opts.spec, opts.analysis)
        else:
            run_full_tree_analysis(opts.input, analysis_type=opts.analysis)
    else:
        raise IOError('No valid analysis type was given.')


def make_graphs(opts):
    """!TODO: Docstring for run_make_animation.
    @return: TODO

    """
    if opts.assay_type == 'fixed-OT':
        graph_fixed_OT_assays(opts)


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
    if opts.analysis:
        run_analysis(opts)

    if opts.graph:
        make_graphs(opts)


##########################################
if __name__ == "__main__":
    main()
