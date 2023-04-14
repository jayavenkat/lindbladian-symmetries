#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 17:45:37 2023

@author: jaya
"""

import numpy as np
from definitions import *


n_fock = 50
kerr = 1
detunings = np.arange(-1*kerr, 12.1*kerr, 0.1*kerr)
squeezing = 0
etas = detunings/kerr

# splittings_detunings = np.zeros((len(detunings), n_fock))
h_eigs = np.zeros((len(detunings), n_fock))

# construct hamiltonian and eigenvalues
for i, detuning in enumerate(detunings):
    h = get_hamiltonian(n_fock, detuning, kerr, squeezing)
    h_eigs[i] = -h.eigenenergies()[::-1]
    

#%%

fig, ax = plt.subplots()
for i in range(12):
        ax.plot(etas, (h_eigs[:, i] - h_eigs[:, 0])/kerr, 'black')
ax.set_xlim(-1, 12)
ax.set_ylim(-0.5, 50)
ax.set_xlabel(r'$\eta = \Delta/K$')
ax.set_ylabel(r'$(E_n - E_0)/K$')
set_figure_size(3, 2.2, ax)
save_path = r'/Users/jaya/My Drive/PhD Work/Papers/Iachello/lindbladian-symmetries/data/notes/lindbladian-symmetries/figures'
title = r'/spectrum_delta'+str(squeezing)+'.png'
plt.savefig(save_path + title, dpi = 900, bbox_inches = 'tight')

#%%

# construct lindbladian and eigenvalues
n_fock = 20
kerr = 1
detunings = np.arange(-1*kerr, 12.1*kerr, 1*kerr)
squeezing = 0
etas = detunings/kerr

# splittings_detunings = np.zeros((len(detunings), n_fock))
l = np.zeros((n_fock*n_fock, n_fock*n_fock))
l_eigs = np.zeros((len(detunings), n_fock*n_fock))
lifetime = np.zeros(len(detunings))
#%%
# bath parameters
kappa = (1/40) * 1/(1/1000) #(kHz)
kerr_exp = 316.83*2*np.pi #kHz
n_th = 0
kappa_loss = kappa * (1 + n_th) / kerr_exp
kappa_gain = kappa * n_th / kerr_exp

# construct lindbladian and eigenvalues
for i, detuning in enumerate(detunings):
    l = get_lindbladian(n_fock, detuning, kerr, squeezing)
    l_eigs[i] = (l.eigenenergies()[::-1])
    sorted_idx = np.argsort(-np.real(l_eigs[i])) # sort according to negative real part
    l_eigs[i] = l_eigs[i, sorted_idx]
    lifetime[i] = -np.real(l_eigs[i, 1])

    
 #%%
print (detunings[3])
l_eigs_min = np.zeros(len(detunings))

fig, ax = plt.subplots()
# for i, detuning in enumerate(detunings):
#     l_eigs_min[i] = np.min(ls_eigs[i])
ax.plot(detunings, 1/(kerr*lifetime))

    