# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 08:41:26 2016

@author: casasorozco
"""

from __future__ import division
import numpy as np
import scipy.misc as misc
import matplotlib.pyplot as plt


def etha(time, time_s, n_param, tau):
    """ Define a function to compute etha: smoothing function to replace
    step function.
    """
    n_sum = np.arange(n_param)
    etha_time = np.zeros(len(time))
    for ind, time in enumerate(time):
        sum_ind = 1/misc.factorial(n_sum) * (n_param * (time - t_s)/tau)**n_sum
        sum_term = sum(sum_ind)
        etha_time[ind] = 1 - np.exp(- n_param * (time - t_s)/tau) * sum_term

    return etha_time


if __name__ == '__main__':
    plt.style.use('~/Dropbox/phd_degree/general/python/daniel_latex.mplstyle')

    # ---------------------------------- Replicating Fig. 2
    fig_order, axis_order = plt.subplots(figsize=(3.3, 2.06))

    t_span = np.linspace(0, .05, 1000)
    t_s = 0  # seconds (this is important, take the time of event as zero!!)
    n_reg = 2
    tau_reg = 1e-2

    for order in range(1, 5):
        etha_eval = etha(t_span, t_s, order, tau_reg)
        axis_order.plot(t_span, etha_eval, label='$n = %s$' % order)

    axis_order.text(.02, .04, r'$\tau = %s$' % tau_reg)
    axis_order.set_xlim(xmax=0.05)
    axis_order.legend(frameon=False)

    plt.savefig('img/n_var.pdf', bbox_inches='tight')

    # ---------------------------------- Replicating Fig. 1
    fig_tau, axis_tau = plt.subplots(figsize=(3.3, 2.06))

    lines = ['-k', '--k', '-.k', ':k']
    tau_array = np.array([.001, .01, .05, .1])

    for i, tau in enumerate(tau_array):
        etha_fun_2 = etha(6*t_span, t_s, n_reg, tau)
        axis_tau.plot(6*t_span, etha_fun_2, lines[i],
                      label=r'$\tau = %s$' % tau, lw=1.5)

    axis_tau.legend(loc='best', frameon=0)
    axis_tau.set_xlim(xmax=6*t_span[-1])
    axis_tau.set_ylim()
    axis_tau.set_xlabel('time')
    axis_tau.set_ylabel('$\eta(t - t_s, n, $' r'$\tau)$')
    axis_tau.text(0.24, 0.8, '$n = %s$' % n_reg)

    plt.savefig('img/tau_var.pdf', bbox_inches='tight')
