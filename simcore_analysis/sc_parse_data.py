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
SITE_DT = np.dtype(['pos', np.double, 3])

FLEX_FIL_DT = np.dtype([()])

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

    if isinstance(p_dict['rigid_filament']):
        for fil_p_dict in p_dict['rigid_filament']:
            fil_name = fil_p_dict['name']
            fil_file_name = run_name + '_rigid_filament_' + fil_name + '.posit'
            print(fil_file_name)
            get_rigid_filament_data(h5_data, fil_file_name, fil_p_dict)
    if isinstance(p_dict['filament']):
        for fil_p_dict in p_dict['filament']:
            fil_name = fil_p_dict['name']
            fil_file_name = run_name + '_filament_' + fil_name + '.spec'
            print(fil_file_name)
            print("WARNING: Flexible filament analysis not implemented yet.")
    if isinstance(p_dict['crosslink']):
        for xl_p_dict in p_dict['crosslink']:
            xl_name = xl_p_dict['name']
            xlink_file_name = run_name + '_crosslink_' + xl_name + '.spec'
            print(xlink_file_name)
            get_xlink_data(h5_data, xlink_file_name, p_dict, xl_p_dict)
    if isinstance(p_dict['optical_trap']):
        for ot_p_dict in p_dict['optical_trap']:
            ot_name = ot_p_dict['name']
            ot_file_name = run_name + '_optical_trap_' + ot_name + '.spec'
            print(ot_file_name)
            get_optical_trap_data(h5_data, ot_file_name, ot_p_dict)
            # print("WARNING: Optical trap analysis not implemented yet.")


def init_data_file(h5_data, param_file_name):
    with open(param_file_name, 'r') as p_file:
        param_dict = yaml.safe_load(p_file)
        h5_data.attrs['param_file'] = yaml.dump(param_dict)
        for key, val in param_dict.items():
            if not isinstance(val, list) and not isinstance(val, dict):
                h5_data.attrs[key] = val


def get_xlink_data(h5_data, xlink_spec_fname, param_dict, xl_p_dict):
    # Get data from xlink file
    half_length = param_dict['rigid_filament'][0]['length'] * .5

    with open(xlink_spec_fname, 'rb') as xlf:
        header = np.fromfile(xlf, HEADER_DT, count=1)[0]
        print(header)
        nframes = int(header[0] / header[1])
        xl_grp = h5_data.create_group('xl_data')
        xl_time_arr = np.arange(0, header[0], header[1]) * header[2]
        xl_grp.create_dataset('time', data=xl_time_arr)

        for key, val in xl_p_dict.items():
            xl_grp.attrs[key] = val
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
            # Load frame data for doubly bound xlinks
            xl_dbl_dset[i, 0] = db_xlinks[:][0]
            xl_dbl_dset[i, 1] = db_xlinks[:][1]
        # Subtract half the length of the filament from the lambda position so
        # zero corresponds to center of the filament.
        xl_sgl_dset[...] -= half_length
        xl_dbl_dset[...] -= half_length


def get_rigid_filament_data(h5_data, fil_posit_fname, fil_p_dict):

    with open(fil_posit_fname, 'rb') as flf:
        header = np.fromfile(flf, HEADER_DT, count=1)[0]
        # fil_num = 2  # Hard coded for two filaments
        data_start = flf.tell()
        fil_num = np.fromfile(flf, np.int32, count=1)[0]
        flf.seek(data_start)
        print(header)
        nframes = int(header[0] / header[1])  # Get number of frames to read

        # Setup of structures to store data
        fil_grp = h5_data.create_group('filament_data')
        fil_time_arr = np.arange(0, header[0], header[1]) * header[2]
        fil_grp.create_dataset('time', data=fil_time_arr)

        for key, val in fil_p_dict.items():
            fil_grp.attrs[key] = val
        # Create a data set that can store all the frames of the doubly bound
        # motors
        fil_pos_dset = fil_grp.create_dataset(
            'filament_position', (nframes, 3, fil_num,))
        fil_orient_dset = fil_grp.create_dataset(
            'filament_orientation', (nframes, 3, fil_num,))

        for i in range(int(header[0] / header[1])):
            # Get number of filaments from file so we know how many to read in
            # next step
            fil_num = np.fromfile(flf, np.int32, count=1)[0]
            # Read number of filament data so we can parse the entire frame and
            # store data
            fils = np.fromfile(flf, FIL_DT, count=fil_num)
            # Store filament lengths for later data analysis
            if 'lengths' not in fil_grp.attrs:
                lengths = [0.] * fil_num
                for fil in fils:
                    lengths[fil['mesh_id'] - 1] = fil['length']
                fil_grp.attrs['lengths'] = lengths
            for fil in fils:
                fil_pos_dset[i, :, fil['mesh_id'] - 1] = fil['pos']
                fil_orient_dset[i, :, fil['mesh_id'] - 1] = fil['orient']


def get_optical_trap_data(h5_data, ot_spec_fname, ot_p_dict):
    """!Get data from optical trap spec files

    @param h5_data: TODO
    @param ot_spec_fname: TODO
    @param ot_p_dict: TODO
    @return: TODO

    """
    with open(ot_spec_fname, 'rb') as otf:
        header = np.fromfile(otf, HEADER_DT, count=1)[0]
        data_start = otf.tell()
        ot_num = np.fromfile(otf, np.int32, count=1)[0]
        otf.seek(data_start)
        # ot_num = 2  # TODO: Hard coded for the moment
        print(header)
        nframes = int(header[0] / header[1])  # Get number of frames to read
        ot_grp = h5_data.create_group('optical_trap_data')
        ot_time_arr = np.arange(0, header[0], header[1]) * header[2]
        ot_grp.create_dataset('time', data=ot_time_arr)

        for key, val in ot_p_dict.items():
            ot_grp.attrs[key] = val

        ot_pos_dset = ot_grp.create_dataset(
            'optical_trap_position', (nframes, 3, ot_num,))
        bead_pos_dset = ot_grp.create_dataset(
            'bead_position', (nframes, 3, ot_num,))

        for i in range(int(header[0] / header[1])):
            # Get number of optical traps to know how many to read in
            ot_num = np.fromfile(otf, np.int32, count=1)[0]
            # Read optical trap data so we can parse the entire frame and
            # store data
            otraps = np.fromfile(otf, OTRAP_DT, count=ot_num)
            # Store filament lengths for later data analysis
            for ot in otraps:
                ot_pos_dset[i, :, ot['attach_id'] - 1] = ot['pos']
                bead_pos_dset[i, :, ot['attach_id'] - 1] = ot['bpos']


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
                        db_list[anch['attached_id'] - 1] += [anch['lambda']]
        elif xl['anchors'][0]['bound']:
            sb_list[xl['anchors'][0]['attached_id'] -
                    1] += [xl['anchors'][0]['lambda']]
    return sb_list, db_list


##########################################
