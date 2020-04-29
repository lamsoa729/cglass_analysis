#!/usr/bin/env python

"""@package docstring
File: fp_steady_state.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
import numpy as np
from scipy.special import erf
from scipy.integrate import quad


def fp_steady_state_antipara(s_i, s_j, y, p_dict):
    L, ks, fs, ko, co, vo, beta = (p_dict['L'], p_dict['ks'], p_dict['fs'], p_dict['ko'],
                                   p_dict['co'], p_dict['vo'], p_dict['beta'])
    xi = s_i + s_j
    # xi = (s_i + s_j).flatten()
    # print(xi)
    lo = 2. * vo / ko
    ls = fs / ks
    a = beta * ks
    alpha = co * np.exp(-.5 * a * y * y) / lo

    def apara_region_1(z):
        pre_fact = alpha * np.sqrt(.5 * np.pi / a) * np.exp(.5 / (a * lo * lo)
                                                            - (z / lo))
        integral = (erf((a * lo * z - 1.) / (lo * np.sqrt(2. * a)))
                    - erf((-a * lo * L - 1.) / (lo * np.sqrt(2. * a))))
        return pre_fact * integral

    def apara_region_2(z):
        psi_0 = apara_region_1(0)
        pre_fact = alpha * ls * np.power(ls - z, (ls / lo) - 1.)
        def integrand(x): return np.power(
            ls - x, -ls / lo) * np.exp(-.5 * a * x * x)
        # print(z)
        integral = np.frompyfunc(lambda x: quad(integrand, 0, x)[0], 1, 1)
        return (pre_fact * integral(z)
                + np.power(1. - (z / ls), (ls / lo) - 1.) * psi_0)

    def apara_region_3(z):
        return co * np.exp(-.5 * a * (z * z + y * y))

    sol = np.piecewise(xi, [xi <= 0, (0 < xi) & (xi <= ls), ls < xi],
                       [apara_region_1, apara_region_2, apara_region_3])
    return sol


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
