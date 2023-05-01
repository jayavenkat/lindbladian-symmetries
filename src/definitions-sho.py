#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:21:55 2023

@author: jaya
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 16:58:40 2023

@author: jaya

"""

# Import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import qutip as qt
from matplotlib import cm

# Set the figure size in inches
def set_figure_size(width, height, ax=None):
    if not ax:
        ax = plt.gca()

    left = ax.figure.subplotpars.left
    right = ax.figure.subplotpars.right
    top = ax.figure.subplotpars.top
    bottom = ax.figure.subplotpars.bottom
    fig_width = float(width) / (right - left)
    fig_height = float(height) / (top - bottom)
    ax.figure.set_size_inches(fig_width, fig_height)

# Set the plotting settings
def set_plot_settings(pres_type):
    if pres_type == 'talk':
        s = 16
    elif pres_type == 'paper':
        s = 10
    mpl.rc('font',family='sans-serif')
    mpl.rc('font',size=s)
    mpl.rc('font',size=s)
    mpl.rcParams['font.family'] = 'CMU Sans Serif'
    mpl.rcParams['axes.formatter.useoffset'] = False
    mpl.rcParams['ytick.right'] = False
    mpl.rcParams['xtick.top'] = False
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['legend.fontsize'] = s
    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['axes.labelsize'] = s
    mpl.rcParams['xtick.labelsize'] = s
    mpl.rcParams['ytick.labelsize'] = s
    mpl.rcParams['xtick.major.pad']=  2 #3.5
    mpl.rcParams['ytick.major.pad']=  2 #3.5
    mpl.rcParams['axes.labelpad'] = 1 #4.0
    mpl.rcParams['legend.handlelength'] = 1.0#2.0
    mpl.rcParams['legend.handletextpad'] = 0.4# 0.8
    mpl.rcParams['legend.columnspacing'] = 1.2# 2.0,
    mpl.rcParams['lines.markersize'] = 4.0

# Construct the Hamiltonian of the Delta-Kerr-cat
def get_hamiltonian(n_fock=100, frequency=0):
    a = qt.destroy(n_fock)
    a_dag = a.dag()
    H = frequency * a_dag * a
    return H

# # Calculate the energy splittings of the system
# def get_splittings(n_fock=100, frequency=0, n_splitting=3):
#     H = get_hamiltonian(n_fock, frequency)
#     splittings = []
#     h_eigs = H.eigenenergies()[::-1]
#     for i in range(2, 2 * (n_splitting + 1), 2):
#         splittings.append(h_eigs[i + 1] - h_eigs[i])
#     return np.array(splittings)

# Create collapse operators for the system
def get_c_ops(n_fock=100, kappa_loss=1e-2, kappa_gain=1e-4):
    a = qt.destroy(n_fock)
    a_dag = a.dag()
    c_ops = [np.sqrt(kappa_loss) * a, np.sqrt(kappa_gain) * a_dag]
    return c_ops

# Construct the Lindbladian of the system
def get_lindbladian(n_fock=100,frequency=0, kappa_loss=1e-2, kappa_gain=1e-4):
    a = qt.destroy(n_fock)
    a_dag = a.dag()
    H = get_hamiltonian(n_fock, frequency)
    c_ops = get_c_ops(n_fock, kappa_loss, kappa_gain)
    lindbladian = qt.liouvillian(H, c_ops)
    return lindbladian

# Obtain sorted eigenstates of the system according to energy and parity
def get_sorted_eigenstates(n_fock=100, frequency=0):
    if n_fock % 2 == 1:
        raise ValueError('n_fock cannot be odd.')

    minus_H = -get_hamiltonian(n_fock, frequency)
    minus_parity_op = -(1j * np.pi * qt.num(n_fock)).expm()

    eigenvalues, eigenstates = qt.simdiag([minus_parity_op, minus_H], evals=True)

    # reorder the eigenstates and eigenvalues according to decreasing energy, and then decreasing parity.
    parities = -eigenvalues[0].reshape(2, n_fock//2).T.reshape(1, -1)
    energies = -eigenvalues[1].reshape(2, n_fock//2).T.reshape(1, -1)
    
    # the phase of the eigenstates are not determined
    plus_eigenstates = eigenstates[:n_fock//2]
    minus_eigenstates = eigenstates[n_fock//2:]
    reordered_eigenstates = [state for pair in zip(plus_eigenstates, minus_eigenstates) for state in pair]
    return parities, energies, reordered_eigenstates