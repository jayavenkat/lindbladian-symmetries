#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 09:31:13 2023

@author: jaya
"""

import numpy as np
import matplotlib.pyplot as plt
from qutip import *
#from qutip import basis, sigmaz, sigmap, sigmam, mesolve, Options
from scipy.constants import hbar
from definitions_sho import *

# path to save plots
save_path = r'/Users/jaya/My Drive/PhD Work/Papers/Iachello/lindbladian-symmetries/data/notes/lindbladian-symmetries/figures'
set_plot_settings('paper')
# time domain simulations
# system parameters
n_fock = 7               # Number of levels in the oscillator's Hilbert space
omega = 1.0                    # Oscillator frequency (arbitrary units)
kappa = 0.1           
timesteps = 1000               # Number of time steps for the simulation
t_final = 200                   # Final time for the simulation (arbitrary units)

# bath parameters
kappa = 0.1                    # Damping rate (arbitrary units)
n_th = 0.0 #* 0                  # thermal photon occupation number
kappa_loss = kappa * (1 + n_th) 
kappa_gain = kappa * n_th

# Define the Hamiltonian, collapse operator and initial state
a = destroy(n_fock)            # Annihilation operator
a_dag = a.dag()
h =  omega * a_dag * a       # Hamiltonian
c_ops = [np.sqrt(kappa_loss) * a, np.sqrt(kappa_gain) * a_dag]
psi0 = basis(n_fock, 2)        # Initial state (5 quanta)

# Time vector for the simulation
t = np.linspace(0, t_final, timesteps)

# Options for the solver
opts = Options(store_states=True, nsteps=5000)

# Solve the master equation
result = mesolve(h, psi0, t, c_ops, [a_dag*a], options=opts)

# Plot the results
plt.close('all')
fig, ax = plt.subplots()
ax.plot(kappa*t, result.expect[0], label = r'numerical$, n_{\mathrm{th}} = $' + str(n_th))
ax.plot(kappa*t, result.expect[0][0]*np.exp(-kappa*t), linestyle = 'dashed', \
        dashes = (2, 2), label = r'analytical$, \langle n \rangle = \langle n_0 \rangle e^{- \kappa t} $')
#ax.plot(kappa*t, result.expect[0][0]*np.exp(-kappa*t), 'k', label = 'time-domain analytical')
ax.semilogy()
plt.legend()
plt.xlabel(r'$\kappa t$')
plt.ylabel('Occupation number')
xticks = ax.get_xticks()
ax.set_xticklabels(replace_minus_with_hyphen(xticks))
yticks = ax.get_yticks()
ax.set_yticklabels(replace_minus_with_hyphen(yticks))
plt.title('Decay of a damped harmonic oscillator')
plt.show()
set_figure_size(3, 2, ax)

title = r'/average_photon_number_DSHO'+str(n_th)+'.png'
plt.savefig(save_path + title, dpi = 900, bbox_inches = 'tight')

#%%

# diagonalize the Lindbladian

#l = np.zeros((n_fock*n_fock, n_fock*n_fock))
l = get_lindbladian(n_fock,omega, kappa_loss, kappa_gain)
l_eigs = (l.eigenenergies()[::-1])
fig, ax = plt.subplots()

ax.scatter(np.real(l_eigs), np.imag(l_eigs), color =  'k', facecolor = 'none')
for i, eig in enumerate(l_eigs):
    if -np.real(eig) == 0.05:
        ax.scatter(np.real(eig), np.imag(eig), color =  'k', facecolor = 'red')
ax.set_xlim(-0.7,0.04)
ax.set_ylim(-6.2,6.2)

ax.set_xlabel(r'$\mathrm{Re}(\mathrm{eigs}(L))$')
ax.set_ylabel(r'$\mathrm{Im}(\mathrm{eigs}(L))$')
plt.suptitle('eigenvalues of the Lindbladian in complex plane')

xticks = list(np.linspace(-0.7, 0., 8))
ax.set_xticks(xticks)
yticks = list(np.arange(-6, 6+1, 2))
ax.set_yticks(yticks)

ax.set_xticklabels(replace_minus_with_hyphen(np.round(xticks, 1)))
ax.set_yticklabels(replace_minus_with_hyphen(np.round(yticks, 1)))
#save_path = r'/Users/jaya/My Drive/PhD Work/Papers/Iachello/lindbladian-symmetries/data/notes/lindbladian-symmetries/figures'
title = r'/Lindblad_eigenvals_DSHO'+str(n_th)+'.png'
plt.savefig(save_path + title, dpi = 900, bbox_inches = 'tight')
# set_plot_settings('paper')
set_figure_size(4, 2, ax)