# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 08:41:26 2016

@author: casasorozco
"""

from __future__ import division
import numpy as np
import scipy.misc as misc
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Define a function to compute etha: smoothing function for replacing step function
def etha(t, ts, n, tau):
    n_sum = np.arange(n)
    etha_time = np.zeros(len(t))
    for ind, time in enumerate(t):
        sum_ind        = 1/misc.factorial(n_sum) * (n * (time - t_s)/tau)**n_sum
        sum_term       = sum(sum_ind)
        etha_time[ind] = 1 - np.exp(- n * (time - t_s)/tau) * sum_term
    
    return etha_time
    
# --------------------------------------------------------------------- Figures
    
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['text.size'] =14
rcParams['text.usetex'] = True

# ---------------------------------- Replicating Fig. 2
plt.figure()
t_span  = np.linspace(0, .05, 1000)
t_s     = 0 # seconds (this is important, take the time of event as zero!!)
n_reg   = 3
tau_reg = 1e-2
    
for order in range(1, 5):
    etha_fun = etha(t_span, t_s, order, tau_reg)
    plt.plot(t_span, etha_fun, label='n = %s' % order)
    
plt.legend(loc=4)
plt.text(.02, .04, r'$\tau = %s$' %tau_reg)

# ---------------------------------- Replicating Fig. 1
plt.figure()

lines = ['-k', '--k', '-.k', ':k']

for i, tau_it in enumerate( np.array([.001, .01, .05, .1]) ):
    etha_fun_2 = etha(6*t_span, t_s, n_reg, tau_it)
    plt.plot(6*t_span, etha_fun_2, lines[i], label='tau = %s' % tau_it, lw=1.5)

plt.legend(loc='best', frameon=0)
plt.xlim(xmax=0.3)
plt.ylim(ymax=1.05)
plt.xlabel('t')
plt.ylabel(r'\eta')
plt.text(0.25, 0.18, 'n = %s' % n_reg)
plt.text(.15, .4, r'$\eta(t - t_s, n, \tau) = 1 - exp \left( -n \frac{(t - t_s)}{\tau} \right) \sum_{j = 0}^{n-1}{ \frac{1}{j!} \left( n \frac{(t-t_s)}{\tau} \right)}$',
         bbox=dict(facecolor='w', alpha=1))