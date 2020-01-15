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
from pathlib import Path

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


def collect_data(h5_data, param_file_name,
                 xlink_file_name, filament_file_name):
    """!TODO: Docstring for collect_data.

    @param h5_data: TODO
    @param param_file_name: TODO
    @param xlink_file_name: TODO
    @param filament_file_name: TODO
    @return: TODO

    """
    init_data_file(h5_data, param_file_name)
    get_xlink_data(h5_data, xlink_file_name)
    get_filament_data(h5_data, filament_file_name)


def parse_xlink_frame(xlink_data):
    sb_list = [[], []]
    db_list = [[], []]

    for xl in xlink_data:
        if xl['doubly']:
            if xl['anchors'][0]['attached_id'] == xl['anchors'][1]['attached_id']:
                print("WARNING: Anchors are attached to same filament.")
            else:
                for anch in xl['anchors']:
                    if anch['attached_id'] < 0:
                        print(
                            "WARNING: Anchor not attached even though doubly bound.")
                    else:
                        db_list[anch['attached_id'] - 1] += [anch['lambda']]
        elif xl['anchors'][0]['bound']:
            sb_list[xl['anchors'][0]['attached_id'] -
                    1] += [xl['anchors'][0]['lambda']]
    return sb_list, db_list


def init_data_file(h5_data, param_file_name):
    with open(param_file_name, 'r') as p_file:
        param_dict = yaml.safe_load(p_file)
        h5_data.attrs['param_file'] = yaml.dump(param_dict)
        for key, val in param_dict.items():
            if not isinstance(val, list) and not isinstance(val, dict):
                h5_data.attrs[key] = val
    # xl_time_arr = h5_data.attrs['delta'] * np.arange(
    #     h5_data.attrs['n_steps'],
    #     param_dict['species']['n_spec'])
    # h5_data.create_dataset('time', data=xl_time_arr)


def get_xlink_data(h5_data, xlink_spec_fname):
    # Get data from xlink file
    param_dict = yaml.safe_load(h5_data.attrs['param_file'])
    h5_data.attrs['param_dict'] = yaml.dump(param_dict)
    half_length = param_dict['rigid_filament'][0]['length'] * .5

    with open(xlink_spec_fname, 'rb') as xlf:
        header = np.fromfile(xlf, HEADER_DT, count=1)[0]
        print(header)
        nframes = int(header[0] / header[1])
        xl_grp = h5_data.create_group('xl_data')
        xl_time_arr = np.arange(0, header[0], header[1]) * header[2]
        xl_grp.create_dataset('time', data=xl_time_arr)

        for key, val in param_dict['crosslink'][0].items():
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


def get_filament_data(h5_data, fil_posit_fname):
    param_dict = yaml.safe_load(h5_data.attrs['param_file'])

    with open(fil_posit_fname, 'rb') as flf:
        header = np.fromfile(flf, HEADER_DT, count=1)[0]
        fil_num = 2  # TODO: Hard coded for the moment
        print(header)
        nframes = int(header[0] / header[1])
        fil_grp = h5_data.create_group('filament_data')
        fil_time_arr = np.arange(0, header[0], header[1]) * header[2]
        fil_grp.create_dataset('time', data=fil_time_arr)

        for key, val in param_dict['rigid_filament'][0].items():
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


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
