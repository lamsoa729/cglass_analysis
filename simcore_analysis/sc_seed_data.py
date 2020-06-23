#!/usr/bin/env python

"""@package docstring
File: seed_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
from pathlib import Path
import h5py
import yaml
import numpy as np
from .sc_graphs import sc_graph_all_data_2d


class SeedData():

    """!Docstring for SeedData. """

    def __init__(self, param_file, sd_path=Path('./')):
        """!Initialize with parameter file

        @param param_file: parameter file for simcore seed

        """

        self._param_file = sd_path / param_file
        with open(self._param_file, 'r') as pf:
            self.p_dict = yaml.safe_load(pf)
        self.run_name = self.p_dict['run_name']
        self._h5_file = sd_path / '{}_data.h5'.format(self.run_name)
        self.h5_data = self.load()
        self.parse_data()
        self.init_flag = True

    def parse_data(self):
        """!TODO: Docstring for _set_data.
        @return: TODO

        """
        fil_grp = self.h5_data['filament_data']
        self.time = fil_grp['time'][:]

        self.r_i_arr = fil_grp['filament_position'][:, :, 0]
        self.r_j_arr = fil_grp['filament_position'][:, :, 1]
        self.u_i_arr = fil_grp['filament_orientation'][:, :, 0]
        self.u_j_arr = fil_grp['filament_orientation'][:, :, 1]

        self.xl_dbl_distr_arr, self.fil_bins = self.analyze_xlink_distr()
        self.xl_dbl_distr_max = np.amax(self.xl_dbl_distr_arr)

    def analyze_xlink_distr(self, bin_num=120):

        length = self.h5_data['filament_data'].attrs['lengths'][0]
        fil_bins = np.linspace(-.5 * length, .5 * length, bin_num)

        # print(fil_bins)
        # Combine all time data to get an average density
        dbl_xlink = self.h5_data['xl_data/doubly_bound'][...]
        # print(dbl_xlink_dset.shape[0])
        dbl_2D_distr = np.zeros((dbl_xlink.shape[0], bin_num - 1, bin_num - 1))
        for n in range(dbl_xlink.shape[0]):
            fil_i_lambdas = dbl_xlink[n, 0]
            fil_j_lambdas = dbl_xlink[n, 1]
            dbl_2D_distr[n], xedges, yedges = np.histogram2d(
                fil_i_lambdas, fil_j_lambdas, fil_bins)
        # print(dbl_2D_distr)
        # ax.set_aspect('equal')
        # cb = ax.pcolormesh(fil_bins, fil_bins,
            # dbl_2D_distr.T)
        return dbl_2D_distr, fil_bins

    def load(self):
        """!Load h5_data file
        @return: h5_data file

        """
        h5_data = h5py.File(self._h5_file, 'r+')
        # if analysis_type != 'load' and ('xl_data' not in h5_data
        #                                 or 'filament_data' not in h5_data):
        #     print("ANALYSIS: Collecting data")
        #     xl_name = p_dict['crosslink'][0]['name']
        #     fil_name = p_dict['rigid_filament'][0]['name']
        #     collect_data(h5_data,
        #                  run_name + '_params.yaml',
        #                  run_name + '_crosslink_' + xl_name + '.spec',
        #                  run_name + '_rigid_filament_' + fil_name + '.posit')
        # print("ANALYSIS: Analyzing data")
        # analyze_seed(h5_data)
        # # Get run time statistics if they exist
        # time_anal_flag = p_dict.get('time_analysis', False)
        # if time_anal_flag:
        #     try:
        #         cpu_time = get_cpu_time_from_log(Path(run_name + '.log'))
        #         h5_data['analysis'].attrs['cpu_time'] = cpu_time
        #     except BaseException:
        #         print("ANALYSIS: !!! Could not collect time analysis !!!")
        # except BaseException:
        #     print("ANALYSIS-MOVIE: failed")
        #     raise
        # finally:
        #     h5_data.close()
        return h5_data

    def animate(self, n, fig, axarr):
        gca_arts = sc_graph_all_data_2d(n, fig, axarr, self)
        return gca_arts

    def save(self):
        """!TODO: Docstring for save.
        @return: TODO

        """
        self.h5_data.close()


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
