#!/usr/bin/env python

"""@package docstring
File: sc_seed_scan.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
from pathlib import Path
import numpy as np
import yaml
import h5py


def collect_seed_h5_files(dir_path):
    """ Spider through directory structure to collect and put h5 files in a list"""
    h5_data_lst = []
    for hf in dir_path.glob('[!.]*/*.h5'):
        h5d = h5py.File(hf, 'r+')
        h5_data_lst += [h5d]
        if not 'seed' in h5d.attrs:
            print("!!! {} does not have seed attribute.", hf)
        else:
            h5_data_lst += [h5d]
    h5_data_lst = sorted(h5_data_lst, key=lambda x: x.attrs['seed'])

    # h5_data_lst = sorted([h5py.File(hf,
    # 'r+') for hf in dir_path.glob('[!.]*/*.h5')],
    # key=lambda x: x.attrs['seed'])
    return h5_data_lst


def analyze_seed_scan(h5_out, h5_data_lst):
    """!TODO: Docstring for analyze_seed_scan_data.

    @param h5_output: TODO
    @param h5_file_lst: TODO
    @return: TODO

    """
    # Copy over params to h5 file
    h5_out.attrs['param_file'] = h5_data_lst[0].attrs['param_file']
    h5_out.attrs['n_seeds'] = len(h5_data_lst)
    # Create crosslink data file
    xl_grp = h5_out.create_group('xl_data')
    fil_grp = h5_out.create_group('filament_data')
    # Copy over params for crosslinkers
    for key, val in h5_data_lst[0]['xl_data'].attrs.items():
        xl_grp.attrs[key] = val
    for key, val in h5_data_lst[0]['filament_data'].attrs.items():
        fil_grp.attrs[key] = val
    h5_out.create_dataset('time', data=h5_data_lst[0]['xl_data/time'][...])

    # Analyze crosslink values
    analyze_avg_moments(xl_grp, h5_data_lst)
    analyze_avg_dbl_distr(xl_grp, h5_data_lst)
    analyze_avg_sgl_num(xl_grp, h5_data_lst)
    analyze_avg_sgl_distr(xl_grp, h5_data_lst)

    # analyze forces and torques
    analyze_avg_forces(h5_out, h5_data_lst)
    analyze_avg_work(h5_out, h5_data_lst)

    # Analyze filament positions
    analyze_avg_fil_dist(fil_grp, h5_data_lst)
    analyze_avg_fil_ang(fil_grp, h5_data_lst)


def analyze_avg_moments(xl_grp, h5_data_lst):
    """!TODO: Docstring for analyze_avg_moments.

    @param h5_out: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    # Collect data to analyze by combining all h5 data arrays in a 2D array
    mu00_arr = np.asarray([h5d['analysis/xl_zeroth_moment'][...]
                           for h5d in h5_data_lst])
    mu10_arr = np.asarray([h5d['analysis/xl_first_moments'][:, 0][...]
                           for h5d in h5_data_lst])
    mu01_arr = np.asarray([h5d['analysis/xl_first_moments'][:, 1][...]
                           for h5d in h5_data_lst])
    mu11_arr = np.asarray([h5d['analysis/xl_second_moments'][:, 0][...]
                           for h5d in h5_data_lst])
    mu20_arr = np.asarray([h5d['analysis/xl_second_moments'][:, 1][...]
                           for h5d in h5_data_lst])
    mu02_arr = np.asarray([h5d['analysis/xl_second_moments'][:, 2][...]
                           for h5d in h5_data_lst])

    xl_grp.create_dataset('zeroth_moment_mean', data=mu00_arr.mean(axis=0))
    xl_grp.create_dataset('zeroth_moment_std', data=mu00_arr.std(axis=0))

    first_mom_mean_arr = np.vstack((mu10_arr.mean(axis=0),
                                    mu01_arr.mean(axis=0)))
    first_mom_std_arr = np.vstack((mu10_arr.std(axis=0),
                                   mu01_arr.std(axis=0)))
    xl_grp.create_dataset(
        'first_moments_mean', data=first_mom_mean_arr.T)
    xl_grp.create_dataset(
        'first_moments_std', data=first_mom_std_arr.T)

    second_mom_mean_arr = np.vstack((mu11_arr.mean(axis=0),
                                     mu20_arr.mean(axis=0),
                                     mu02_arr.mean(axis=0)))
    second_mom_std_arr = np.vstack((mu11_arr.std(axis=0),
                                    mu20_arr.std(axis=0),
                                    mu02_arr.std(axis=0)))
    xl_grp.create_dataset(
        'second_moments_mean', data=second_mom_mean_arr.T)
    xl_grp.create_dataset(
        'second_moments_std', data=second_mom_std_arr.T)


def analyze_avg_dbl_distr(xl_grp, h5_data_lst):
    """!TODO: Docstring for analyze_avg_moments.

    @param h5_out: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    xl_grp.attrs['sgl_bin_edges'] = h5_data_lst[0]['analysis/singly_bound_distr'].attrs['bin_edges']
    fil0_sgl_avg_distr = np.asarray([h5d['analysis/singly_bound_distr'][0][...]
                                     for h5d in h5_data_lst])
    fil1_sgl_avg_distr = np.asarray([h5d['analysis/singly_bound_distr'][1][...]
                                     for h5d in h5_data_lst])
    sgl_avg_distr_mean = np.vstack((fil0_sgl_avg_distr.mean(axis=0),
                                    fil1_sgl_avg_distr.mean(axis=0)))
    sgl_avg_distr_std = np.vstack((fil0_sgl_avg_distr.std(axis=0),
                                   fil1_sgl_avg_distr.std(axis=0)))

    xl_grp.create_dataset('average_singly_bound_distr_mean',
                          data=sgl_avg_distr_mean)
    xl_grp.create_dataset('average_singly_bound_distr_std',
                          data=sgl_avg_distr_std)


def analyze_avg_sgl_distr(xl_grp, h5_data_lst):
    """!TODO: Docstring for analyze_avg_moments.

    @param h5_out: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    n_seeds = len(h5_data_lst)
    xl_grp.attrs['xedges'] = h5_data_lst[0]['analysis/average_doubly_bound_distr'].attrs['xedges']
    xl_grp.attrs['yedges'] = h5_data_lst[0]['analysis/average_doubly_bound_distr'].attrs['yedges']

    # Collect data to analyze by combining all h5 data arrays in a 2D array
    xl_dbl_distr = np.asarray([h5d['analysis/average_doubly_bound_distr'][...]
                               for h5d in h5_data_lst])
    xl_grp.create_dataset('average_doubly_bound_distr_mean',
                          data=xl_dbl_distr.mean(axis=0))
    xl_grp.create_dataset('average_doubly_bound_distr_std',
                          data=xl_dbl_distr.std(axis=0))


def analyze_avg_sgl_num(xl_grp, h5_data_lst):
    """!Analyze the number of singly bound motors or crosslinkers per step

    @param xl_grp: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    sgl_fil0_num_arr = np.asarray([h5d['analysis/singly_bound_number'][:, 0]
                                   for h5d in h5_data_lst])
    sgl_fil1_num_arr = np.asarray([h5d['analysis/singly_bound_number'][:, 1]
                                   for h5d in h5_data_lst])
    sgl_num_arr_mean = np.vstack((sgl_fil0_num_arr.mean(axis=0),
                                  sgl_fil1_num_arr.mean(axis=0)))
    sgl_num_arr_std = np.vstack((sgl_fil0_num_arr.std(axis=0),
                                 sgl_fil1_num_arr.std(axis=0)))
    xl_grp.create_dataset('singly_bound_number_mean',
                          data=sgl_num_arr_mean.T)
    xl_grp.create_dataset('singly_bound_number_std',
                          data=sgl_num_arr_std.T)


def analyze_avg_forces(h5_out, h5_data_lst):
    """!TODO: Docstring for analyze_avg_moments.

    @param h5_out: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    force_arr = np.asarray(
        [h5d['analysis/xl_forces'][...] for h5d in h5_data_lst])
    torque_arr = np.asarray(
        [h5d['analysis/xl_torques'][...] for h5d in h5_data_lst])
    h5_out.create_dataset('xl_forces_mean', data=force_arr.mean(axis=0))
    h5_out.create_dataset('xl_forces_std', data=force_arr.std(axis=0))
    h5_out.create_dataset('xl_torques_mean', data=torque_arr.mean(axis=0))
    h5_out.create_dataset('xl_torques_std', data=torque_arr.std(axis=0))


def analyze_avg_work(h5_out, h5_data_lst):
    """!TODO: Docstring for analyze_avg_moments.

    @param h5_out: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    lin_work_arr = np.asarray(
        [h5d['analysis/xl_linear_work'][...] for h5d in h5_data_lst])
    rot_work_arr = np.asarray(
        [h5d['analysis/xl_rotational_work'][...] for h5d in h5_data_lst])
    h5_out.create_dataset('xl_lin_work_mean', data=lin_work_arr.mean(axis=0))
    h5_out.create_dataset('xl_lin_work_std', data=lin_work_arr.std(axis=0))
    h5_out.create_dataset('xl_rot_work_mean', data=rot_work_arr.mean(axis=0))
    h5_out.create_dataset('xl_rot_work_std', data=rot_work_arr.std(axis=0))


def analyze_avg_fil_dist(fil_grp, h5_data_lst):
    """!Analyze the separation vectors between filament centers.

    @param fil_grp: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    r_ij_arr = []
    for h5d in h5_data_lst:
        fil_pos_dset = h5d['filament_data/filament_position']
        r_ij_arr += [fil_pos_dset[:, :, 1] - fil_pos_dset[:, :, 0]]

    r_ij_arr = np.asarray(r_ij_arr)

    fil_grp.create_dataset('fil_avg_sep_mean', data=r_ij_arr.mean(axis=0))
    fil_grp.create_dataset('fil_avg_sep_std', data=r_ij_arr.std(axis=0))


def analyze_avg_fil_ang(fil_grp, h5_data_lst):
    """!Analyze the separation vectors between filament centers.

    @param fil_grp: TODO
    @param h5_data_lst: TODO
    @return: TODO

    """
    uiuj_arr = np.zeros(
        (h5_data_lst[0]['filament_data/filament_orientation'].shape[0],
         len(h5_data_lst)))
    for i, h5d in enumerate(h5_data_lst):
        fil_orient_dset = h5d['filament_data/filament_orientation']
        uiuj_arr[:, i] = np.einsum('ij,ij->i', fil_orient_dset[:, :, 0],
                                   fil_orient_dset[:, :, 1])
    theta_arr = np.arccos(uiuj_arr)
    fil_grp.create_dataset('fil_avg_theta_mean', data=theta_arr.mean(axis=0))
    fil_grp.create_dataset('fil_avg_theta_std', data=theta_arr.std(axis=0))


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
