#!/usr/bin/env python

"""@package docstring
File: sc_analyze_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
import numpy as np
from pathlib import Path


def normalize(vec):
    """!TODO: Docstring for normalize.

    @param vec: TODO
    @return: TODO

    """
    norm = np.linalg.norm(vec, axis=-1)
    return np.divide(vec, norm[:, None],
                     out=np.zeros_like(vec), where=norm[:, None] != 0)


def flatten_dset(dset):
    """!Flattens a time series of positions of crosslinker heads on filament

    @param dset: TODO
    @return: TODO

    """
    return [item for sublist in dset for item in sublist]


def xl_zrl_force(r_i, r_j, u_i, u_j, s_i, s_j, ks):
    """!TODO: Docstring for xl_zrl_force.

    @param r_i: TODO
    @param r_j: TODO
    @param u_i: TODO
    @param u_j: TODO
    @param s_i: TODO
    @param s_j: TODO
    @param ks: TODO
    @return: TODO

    """
    return -ks * (r_j + (u_j * s_j) - r_i - (u_i * s_i))


def xl_zrl_stretch(r_i, r_j, u_i, u_j, s_i, s_j):
    return np.linalg.norm(r_j + (u_j * s_j) - r_i - (u_i * s_i))


def analyze_seed(h5_data):
    if 'analysis' in h5_data:
        del h5_data['analysis']  # Start clean
    h5_data.create_group('analysis')
    # analyze xlinks
    analyze_singly_bound_xlinks(h5_data)
    analyze_xlink_moments(h5_data)
    analyze_xlink_force(h5_data)
    analyze_xlink_work(h5_data)
    analyze_xlink_stretch_distr(h5_data)
    # if h5_data['filament'].attrs.get('stationary_flag', False):
    analyze_avg_xlink_distr(h5_data)
    # analyze filaments (maybe)


def analyze_xlink_moments(h5_data):
    anal_grp = h5_data['analysis']
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']

    dbl_num_arr = [xld_t.size for xld_t in dbl_xlink_dset[:, 0]]
    anal_grp.create_dataset('xl_zeroth_moment', data=dbl_num_arr)

    mu10_arr = [np.sum(xld_t) for xld_t in dbl_xlink_dset[:, 0]]
    mu01_arr = [np.sum(xld_t) for xld_t in dbl_xlink_dset[:, 1]]
    xl_first_mom_arr = np.vstack((mu10_arr, mu01_arr))
    anal_grp.create_dataset('xl_first_moments', data=xl_first_mom_arr.T)

    mu20_arr = [np.sum(np.power(xld_t, 2)) for xld_t in dbl_xlink_dset[:, 0]]
    mu02_arr = [np.sum(np.power(xld_t, 2)) for xld_t in dbl_xlink_dset[:, 1]]
    mu11_arr = [np.dot(xld_0, xld_1) for xld_0, xld_1 in
                zip(dbl_xlink_dset[:, 0], dbl_xlink_dset[:, 1])]
    xl_second_mom_arr = np.vstack((mu11_arr, mu20_arr, mu02_arr))
    anal_grp.create_dataset('xl_second_moments', data=xl_second_mom_arr.T)


def analyze_singly_bound_xlinks(h5_data):
    """!TODO: Docstring for analyze_singly_bound_xlink_num.

    @param h5_data: TODO
    @return: TODO

    """
    anal_grp = h5_data['analysis']
    xl_sgl_dset = h5_data['xl_data/singly_bound']

    xl_sgl_num_arr = np.zeros((xl_sgl_dset.shape[0], 2))
    for i, xl_arr in enumerate(xl_sgl_dset):
        xl_sgl_num_arr[i, 0] = xl_arr[0].size
        xl_sgl_num_arr[i, 1] = xl_arr[1].size
    anal_grp.create_dataset('singly_bound_number', data=xl_sgl_num_arr)

    half_l = h5_data['filament_data'].attrs['lengths'][0] * .5
    n_steps = h5_data.attrs['n_steps']
    n_spec = h5_data['xl_data'].attrs['n_spec']
    fil0_lambdas = np.asarray(flatten_dset(xl_sgl_dset[:, 0]))
    fil1_lambdas = np.asarray(flatten_dset(xl_sgl_dset[:, 1]))
    xl_fil0_avg_distr, bin_edges = np.histogram(fil0_lambdas, 50,
                                                range=[-half_l, half_l])
    xl_fil1_avg_distr, bin_edges = np.histogram(fil1_lambdas, 50,
                                                range=[-half_l, half_l])
    xl_sgl_avg_distr = np.stack(
        (xl_fil0_avg_distr, xl_fil1_avg_distr)).astype(float)
    xl_sgl_avg_distr *= float(n_spec / n_steps)
    xl_sgl_distr_dset = anal_grp.create_dataset('singly_bound_distr',
                                                data=xl_sgl_avg_distr)
    xl_sgl_distr_dset.attrs['bin_edges'] = bin_edges


def analyze_avg_xlink_distr(h5_data):
    """!TODO: Docstring for analyze_average_xlink_distr.

    @param h5_data: TODO
    @return: TODO

    """
    length = h5_data['filament_data'].attrs['lengths'][0]
    fil_bins = np.linspace(-.5 * length, .5 * length, 120)

    # Combine all time data to get an average density
    dbl_xlink_dset = h5_data['xl_data/doubly_bound']
    fil0_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 0]))
    fil1_lambdas = np.asarray(flatten_dset(dbl_xlink_dset[:, 1]))
    # print(fil1_lambdas)
    dbl_2D_distr, xedges, yedges = np.histogram2d(
        fil0_lambdas, fil1_lambdas, fil_bins)

    n_steps = h5_data.attrs['n_steps']
    n_spec = h5_data['xl_data'].attrs['n_spec']

    dbl_2D_distr *= float(n_spec / n_steps)

    xl_avg_distr_dset = h5_data['analysis'].create_dataset(
        'average_doubly_bound_distr', data=dbl_2D_distr)
    xl_avg_distr_dset.attrs['xedges'] = xedges
    xl_avg_distr_dset.attrs['yedges'] = yedges


def analyze_xlink_force(h5_data):
    """!Analyze the force on filament_j by crosslinkers attached to filament_i

    @param h5_data: TODO
    @return: TODO

    """
    xl_dbl_dset = h5_data['xl_data/doubly_bound']
    ks = h5_data['xl_data'].attrs['k_spring']
    fil_pos_dset = h5_data['filament_data/filament_position']
    fil_orient_dset = h5_data['filament_data/filament_orientation']
    u_i = fil_orient_dset[:, :, 0]
    u_j = fil_orient_dset[:, :, 1]
    r_i = fil_pos_dset[:, :, 0]
    r_j = fil_pos_dset[:, :, 1]
    xl_si = xl_dbl_dset[:, 0]
    xl_sj = xl_dbl_dset[:, 1]
    force_arr = np.zeros((len(u_i), 3))
    torque_arr = np.zeros((len(u_i), 2, 3))
    nframes = len(u_i)
    for i in range(nframes):
        for xl in range(len(xl_si[i])):
            force = xl_zrl_force(r_i[i], r_j[i], u_i[i], u_j[i],
                                 xl_si[i][xl], xl_sj[i][xl],
                                 ks)
            force_arr[i, :] += force
            torque_arr[i, 0, :] += np.cross(u_i[i] * xl_si[i][xl], -force)
            torque_arr[i, 1, :] += np.cross(u_j[i] * xl_sj[i][xl], force)

    h5_data['analysis'].create_dataset('xl_forces', data=force_arr)
    h5_data['analysis'].create_dataset('xl_torques', data=torque_arr)


def analyze_xlink_stretch_distr(h5_data):
    """!TODO: Docstring for analyze_xlink_stretch_distr.

    @param h5_data: TODO
    @return: TODO

    """

    xl_dbl_dset = h5_data['xl_data/doubly_bound']

    fil_pos_dset = h5_data['filament_data/filament_position']
    fil_orient_dset = h5_data['filament_data/filament_orientation']
    u_i = fil_orient_dset[:, :, 0]
    u_j = fil_orient_dset[:, :, 1]
    r_i = fil_pos_dset[:, :, 0]
    r_j = fil_pos_dset[:, :, 1]
    xl_si = xl_dbl_dset[:, 0]
    xl_sj = xl_dbl_dset[:, 1]
    nframes = len(u_i)
    print(nframes)
    stretch_frame_list = []
    max_h = 0
    for i in range(nframes):
        stretch_frame_list += [[]]
        for s_i, s_j in zip(xl_si[i], xl_sj[i]):
            stretch = xl_zrl_stretch(r_i[i], r_j[i],
                                     u_i[i], u_j[i],
                                     s_i, s_j)
            max_h = stretch if stretch > max_h else max_h
            stretch_frame_list[-1] += [stretch]

    step = .004
    fil_bins = np.arange(0, max_h + 2. * step, step)
    stretch_list_hist = np.zeros((nframes, len(fil_bins) - 1))
    for i, xl_list in enumerate(stretch_frame_list):
        stretch_list_hist[i] = np.histogram(xl_list, fil_bins)[0]

    stretch_dset = h5_data['analysis'].create_dataset('xl_stretch',
                                                      data=stretch_list_hist)
    try:
        stretch_dset.attrs['bin_edges'] = fil_bins
    except BaseException:
        pass

    h5_data['analysis'].create_dataset('xl_stretch_bin_edges', data=fil_bins)


def analyze_xlink_work(h5_data):
    """!TODO: Docstring for analyze_xlink_work.

    @param h5_data: TODO
    @return: TODO

    """
    if 'analysis/xl_forces' not in h5_data:
        analyze_xlink_force(h5_data)
    fil_pos_dset = h5_data['filament_data/filament_position']
    fil_orient_dset = h5_data['filament_data/filament_orientation']
    u_i = fil_orient_dset[:, :, 0]
    u_j = fil_orient_dset[:, :, 1]
    r_i = fil_pos_dset[:, :, 0]
    r_j = fil_pos_dset[:, :, 1]

    # Linear work calculations
    dr_i = np.zeros(r_i.shape)
    dr_i[1:] = r_i[1:] - r_i[:-1]
    f_i = -1. * h5_data['analysis/xl_forces'][...]
    dwl_i = np.zeros(r_i.shape[0])
    # Use trapezoid rule for numerical integration
    dwl_i[1:] = .5 * (np.einsum('ij,ij->i', dr_i[1:], f_i[:-1]) +
                      np.einsum('ij,ij->i', dr_i[1:], f_i[1:]))

    dr_j = np.zeros(r_j.shape)
    dr_j[1:] = r_j[1:] - r_j[:-1]
    f_j = h5_data['analysis/xl_forces'][...]
    dwl_j = np.zeros(r_j.shape[0])
    # Use trapezoid rule for numerical integration
    dwl_j[1:] = .5 * (np.einsum('ij,ij->i', dr_j[1:], f_j[:-1]) +
                      np.einsum('ij,ij->i', dr_j[1:], f_j[1:]))

    # Rotational work calculations
    dtheta_i_vec = np.zeros(u_i.shape)
    # Get the direction of small rotation
    dtheta_i_vec[1:] = normalize(np.cross(u_i[:-1], u_i[1:]))
    # Get amplitude of small rotation
    dtheta_i_vec[1:] *= np.arccos(
        np.clip(np.einsum('ij,ij->i', u_i[1:], u_i[:-1]), - 1., 1.))[:, None]
    tau_i = h5_data['analysis/xl_torques'][:, 0, :]
    dwr_i = np.zeros(u_i.shape[0])
    # Use trapezoid rule for numerical integration
    dwr_i[1:] = .5 * (np.einsum('ij,ij->i', dtheta_i_vec[1:], tau_i[:-1]) +
                      np.einsum('ij,ij->i', dtheta_i_vec[1:], tau_i[1:]))

    dtheta_j_vec = np.zeros(u_j.shape)
    # Get the direction of small rotation
    dtheta_j_vec[1:] = normalize(np.cross(u_j[:-1], u_j[1:]))
    # Get amplitude of small rotation
    dtheta_j_vec[1:] *= np.arccos(
        np.clip(np.einsum('ij,ij->i', u_j[1:], u_j[:-1]), -1., 1))[:, None]
    tau_j = h5_data['analysis/xl_torques'][:, 1, :]
    dwr_j = np.zeros(u_j.shape[0])
    # Use trapezoid rule for numerical integration
    dwr_j[1:] = .5 * (np.einsum('ij,ij->i', dtheta_j_vec[1:], tau_j[:-1]) +
                      np.einsum('ij,ij->i', dtheta_j_vec[1:], tau_j[1:]))
    xl_lin_work_dset = h5_data['analysis'].create_dataset(
        'xl_linear_work', data=np.stack((dwl_i, dwl_j), axis=-1),
        dtype=np.float32)
    xl_rot_work_dset = h5_data['analysis'].create_dataset(
        'xl_rotational_work', data=np.stack((dwr_i, dwr_j), axis=-1),
        dtype=np.float32)


#######
