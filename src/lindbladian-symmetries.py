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
squeezing = 2
etas = detunings/kerr

# splittings_detunings = np.zeros((len(detunings), n_fock))
h_eigs = np.zeros((len(detunings), n_fock))

for i, detuning in enumerate(detunings):
    h = get_hamiltonian(n_fock, detuning, kerr, squeezing)
    h_eigs[i] = -h.eigenenergies()[::-1]

#%%

fig, ax = plt.subplots()
for i in range(12):
        ax.plot(etas, (h_eigs[:, i] - h_eigs[:, 0])/kerr, 'black')
ax.set_xlim(-1, 12)
ax.set_ylim(0, 50)
ax.set_xlabel(r'$\eta = \Delta/K$')
ax.set_ylabel(r'$(E_n - E_0)/K$')
set_figure_size(3, 2.2, ax)
save_path = r'/Users/jaya/My Drive/PhD Work/Papers/Iachello/lindbladian-symmetries/data/figures'
title = r'/spectrum_delta'+str(squeezing)+'.png'
plt.savefig(save_path + title, dpi = 900, bbox_inches = 'tight')