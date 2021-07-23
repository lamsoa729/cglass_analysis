#!/usr/bin/env python

"""@package docstring
File: sc_parse_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import numpy as np
import yaml
import h5py
import re

HEADER_DT = np.dtype([('n_steps', np.int32),
                      ('n_posit', np.int32),
                      ('delta', np.double)])

ANCHOR_DT = np.dtype([('bound', np.bool),
                      ('active', np.bool),
                      ('static', np.bool),
                      ('pos', np.double, 3),
                      ('orient', np.double, 3),
                      ('lambda', np.double),
                      ('attached_id', np.int32),
                      ])

XLINK_DT = np.dtype([('doubly', np.bool),
                     ('diameter', np.double),
                     ('length', np.double),
                     ('pos', np.double, 3),
                     ('orient', np.double, 3),
                     ('anchors', ANCHOR_DT, 2),
                     ])

FIL_DT = np.dtype([('pos', np.double, 3),
                   ('spos', np.double, 3),
                   ('orient', np.double, 3),
                   ('diameter', np.double),
                   ('length', np.double),
                   ('mesh_id', np.int32),
                   ])

SITE_DT = np.dtype([('pos', np.double, 3)])

# FLEX_FIL_DT = np.dtype([()])

OTRAP_DT = np.dtype([('pos', np.double, 3),
                     ('spos', np.double, 3),
                     ('orient', np.double, 3),
                     ('diameter', np.double),
                     ('length', np.double),
                     ('bpos', np.double, 3),
                     ('bspos', np.double, 3),
                     ('attach_id', np.int32)])


def collect_data(h5_data, param_file_name):
    """!TODO: Docstring for collect_data.

    @param h5_data: TODO
    @param param_file_name: TODO
    @return: void, updates h5_data with information

    """
    init_data_file(h5_data, param_file_name)
    p_dict = yaml.safe_load(h5_data.attrs['param_file'])
    run_name = p_dict['run_name']

    if isinstance(p_dict['rigid_filament'], list):
        rg_fil_grp = h5_data.create_group('rigid_filament_data')
        for fil_p_dict in p_dict['rigid_filament']:
            get_rigid_filament_data(rg_fil_grp, run_name, fil_p_dict)
    if isinstance(p_dict['filament'], list):
        # fil_grp = h5_data.create_group('filament_data')
        for fil_p_dict in p_dict['filament']:
            print("WARNING: Flexible filament analysis not implemented yet.")
    if isinstance(p_dict['crosslink'], list):
        xl_grp = h5_data.create_group('crosslink_data')
        for xl_p_dict in p_dict['crosslink']:
            if xl_p_dict['concentration'] > 0:
                get_xlink_data(xl_grp, run_name, p_dict, xl_p_dict)
    if isinstance(p_dict['optical_trap'], list):
        ot_grp = h5_data.create_group('optical_trap_data')
        for ot_p_dict in p_dict['optical_trap']:
            get_optical_trap_data(ot_grp, run_name, ot_p_dict)
            # print("WARNING: Optical trap analysis not implemented yet.")


def init_data_file(h5_data, param_file_name):
    with open(param_file_name, 'r') as p_file:
        param_dict = yaml.safe_load(p_file)
        h5_data.attrs['param_file'] = yaml.dump(param_dict)
        for key, val in param_dict.items():
            if not isinstance(val, list) and not isinstance(val, dict):
                h5_data.attrs[key] = val


def get_xlink_data(h5_data, run_name, param_dict, xl_p_dict):
    # Get data from xlink file
    # FIXME: Make it so that it knows the actual lengths of the filaments
    # crosslinkers are attached to.
    half_length = param_dict['rigid_filament'][0]['length'] * .5

    xl_name = xl_p_dict['name']
    species_name = 'crosslink_' + xl_name
    xlink_spec_fname = run_name + '_' + species_name + '.spec'
    print("---- " + xlink_spec_fname + " header ----")

    with open(xlink_spec_fname, 'rb') as xlf:
        header = np.fromfile(xlf, HEADER_DT, count=1)[0]
        print(header)
        nframes = int(header[0] / header[1])
        xl_grp = h5_data.create_group(xl_name)
        xl_time_arr = np.arange(0, header[0], header[1]) * header[2]
        xl_grp.create_dataset('time', data=xl_time_arr)

        xl_grp.attrs['params'] = yaml.dump(xl_p_dict)

        # for key, val in xl_p_dict.items():
        # xl_grp.attrs[key] = val
        # Create a variable length data type
        xl_type = h5py.special_dtype(vlen=np.dtype(np.double))
        # Create a data set that can store all the frames of the doubly bound
        # motors
        xl_dbl_dset = xl_grp.create_dataset('doubly_bound',
                                            (nframes, 2,),
                                            dtype=xl_type)
        xl_sgl_dset = xl_grp.create_dataset('singly_bound',
                                            (nframes, 2,),
                                            dtype=xl_type)

        # Loop over number of frames which is the time / time step
        dbl_xlink_num = []
        sb_xlink_num = []
        for i in range(nframes):
            # Get number of crosslinks from file so we know how many to read in
            # next step
            xlink_num = np.fromfile(xlf, np.int32, count=1)[0]
            # Read number of xlink data so we can parse the entire frame and
            # store data
            xlinks = np.fromfile(xlf, XLINK_DT, count=xlink_num)
            sb_xlinks, db_xlinks = parse_xlink_frame(xlinks)
            # Load frame data for singly bound xlinks
            xl_sgl_dset[i, 0] = sb_xlinks[:][0]
            xl_sgl_dset[i, 1] = sb_xlinks[:][1]
            sb_xlink_num += [[len(sb_xlinks[:][0]), len(sb_xlinks[:][1])]]
            # Load frame data for doubly bound xlinks
            xl_dbl_dset[i, 0] = db_xlinks[:][0]
            xl_dbl_dset[i, 1] = db_xlinks[:][1]
            dbl_xlink_num += [len(db_xlinks[:][0])]
        # Subtract half the length of the filament from the lambda position so
        # zero corresponds to center of the filament.
        xl_sgl_dset[...] -= half_length
        xl_dbl_dset[...] -= half_length
        # xl_dbl_dset.attrs['dbl_xlink_num_arr'] = np.asarray(dbl_xlink_num)
        # xl_sgl_dset.attrs['sgl_xlink_num_arr'] = np.asarray(sb_xlink_num)
        dbl_num_dset = xl_grp.create_dataset(
            'singly_bound_num', data=sb_xlink_num)
        dbl_num_dset = xl_grp.create_dataset(
            'doubly_bound_num', data=dbl_xlink_num)


def get_rigid_filament_data(h5_data, run_name, fil_p_dict):
    """!Get data from rigid filament posit files. Includes lengths, mesh IDs,
    position, and orientations

    @param h5_data: hdf5 file to write too
    @param run_name: Name of CGLASS run
    @param fil_p_dict: Rigid filament parameter dictionary
    @return: void, changes h5_data to include rigid filament data.

    """

    fil_name = fil_p_dict['name']
    species_name = 'rigid_filament_' + fil_name
    fil_posit_fname = run_name + '_' + species_name + '.posit'
    print("---- " + fil_posit_fname + " header ----")
    fil_grp = h5_data.create_group(fil_name)
    for key, val in fil_p_dict.items():
        fil_grp.attrs[key] = val

    with open(fil_posit_fname, 'rb') as flf:
        header = np.fromfile(flf, HEADER_DT, count=1)[0]
        print(header)

        data_start = flf.tell()  # Save location of file right after header
        nframes = int(header[0] / header[1])  # Get number of frames to read

        # Get constant data that does not change. Lengths one day might change.
        fil_num = np.fromfile(flf, np.int32, count=1)[0]
        fils = np.fromfile(flf, FIL_DT, count=fil_num)

        lengths = [fil['length'] for fil in fils]
        mesh_ids = [fil['mesh_id'] for fil in fils]
        fil_grp.attrs['lengths'] = lengths
        fil_grp.attrs['mesh_ids'] = mesh_ids

        # Reset to beginning of data file.
        flf.seek(data_start)

        # Create data sets for storing filament data
        fil_time_arr = np.arange(0, header[0], header[1]) * header[2]
        fil_grp.create_dataset('time', data=fil_time_arr)
        fil_pos_dset = fil_grp.create_dataset(
            'position', (nframes, 3, fil_num,))
        fil_orient_dset = fil_grp.create_dataset(
            'orientation', (nframes, 3, fil_num,))

        for i in range(nframes):
            # Get number of filaments so we know how many to read in
            fil_num = np.fromfile(flf, np.int32, count=1)[0]
            # Read filament data so we can parse and store the entire frame
            fils = np.fromfile(flf, FIL_DT, count=fil_num)
            for fil in fils:
                mesh_id = fil['mesh_id']
                fil_pos_dset[i, :, mesh_ids.index(mesh_id)] = fil['pos']
                fil_orient_dset[i, :, mesh_ids.index(mesh_id)] = fil['orient']


def get_optical_trap_data(h5_data, run_name, ot_p_dict):
    """!Get data from optical trap spec files

    @param h5_data: hdf5 file to write too
    @param run_name: Name of CGLASS run
    @param ot_p_dict: Optical trap parameter dictionary
    @return: void, changes h5_data to include optical traps

    """
    ot_name = ot_p_dict['name']
    species_name = 'optical_trap_' + ot_name
    ot_spec_fname = run_name + '_' + species_name + '.spec'
    print("---- " + ot_spec_fname + " header ----")

    ot_grp = h5_data.create_group(ot_name)
    for key, val in ot_p_dict.items():
        ot_grp.attrs[key] = val

    with open(ot_spec_fname, 'rb') as otf:
        header = np.fromfile(otf, HEADER_DT, count=1)[0]
        print(header)
        nframes = int(header[0] / header[1])  # Get number of frames to read

        data_start = otf.tell()

        # Get constant data that does not change.
        ot_num = np.fromfile(otf, np.int32, count=1)[0]
        otraps = np.fromfile(otf, OTRAP_DT, count=ot_num)
        attach_ids = [ot['attach_id'] for ot in otraps]
        ot_grp.attrs['attach_ids'] = attach_ids

        otf.seek(data_start)
        # ot_num = 2  # TODO: Hard coded for the moment

        ot_time_arr = np.arange(0, header[0], header[1]) * header[2]

        ot_grp.create_dataset('time', data=ot_time_arr)
        ot_pos_dset = ot_grp.create_dataset(
            'trap_position', (nframes, 3, ot_num,))
        bead_pos_dset = ot_grp.create_dataset(
            'bead_position', (nframes, 3, ot_num,))

        print(nframes)
        for i in range(nframes):
            # Get number of optical traps to know how many to read in
            try:
                ot_num = np.fromfile(otf, np.int32, count=1)[0]
            except BaseException:
                print(" Could not get number of optical traps."
                      " Possibly hit end of file.")
                break
            # Read optical trap data so we can parse the entire frame and
            # store data
            otraps = np.fromfile(otf, OTRAP_DT, count=ot_num)
            # Store filament lengths for later data analysis
            for ot in otraps:
                attach_id = ot['attach_id']
                ot_pos_dset[i, :, attach_ids.index(attach_id)] = ot['pos']
                bead_pos_dset[i, :, attach_ids.index(attach_id)] = ot['bpos']


def get_cpu_time_from_log(log_file):
    """!TODO: Docstring for get_cpu_time_from_log.
    @return: TODO

    """
    if not log_file.exists():
        raise OSError

    pattern = re.compile(r'CPU Time: ([0-9]*\.[0-9]+)?$')
    with open(log_file, 'r') as lf:
        cpu_t = [re.findall(pattern, line)[0]
                 for line in lf if re.findall(pattern, line)][0]
    return float(cpu_t)


def parse_xlink_frame(xlink_data):
    sb_list = [[], []]
    db_list = [[], []]

    for xl in xlink_data:
        if xl['doubly']:
            if (xl['anchors'][0]['attached_id'] ==
                    xl['anchors'][1]['attached_id']):
                print("WARNING: Anchors are attached to same filament.")
            else:
                for anch in xl['anchors']:
                    if anch['attached_id'] < 0:
                        print("WARNING: ",
                              "Anchor not attached even though doubly bound.")
                    else:
                        # print("double bound",
                        #       xl['anchors'][0]['attached_id'],
                        #       xl['anchors'][1]['attached_id'])
                        db_list[anch['attached_id'] - 2] += [anch['lambda']]
        elif xl['anchors'][0]['bound']:
            # print("single bound",
            #       xl['anchors'][0]['attached_id'],
            #       xl['anchors'][1]['attached_id'])
            sb_list[xl['anchors'][0]['attached_id'] -
                    2] += [xl['anchors'][0]['lambda']]
    return sb_list, db_list


##########################################
