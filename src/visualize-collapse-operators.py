#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:24:23 2023

@author: jaya
"""
from definitions import *
from qutip import *

n_fock = 20
alpha  = 4

# visualizing coherent state under single photon loss

#define density matrix
rho_coherent = coherent_dm(n_fock, np.sqrt(alpha))

# plot photon number distribution of a coherent state
fig, axes = plt.subplots(1, 2, figsize=(12,3))
bar0 = axes[0].bar(np.arange(0, n_fock)-.5, rho_coherent.diag())
axes[0].set_xticks(list(np.arange(0, n_fock-1, 1)))
lbl0 = axes[0].set_title("Coherent state")

# define annihilation and creation operators
a = qt.destroy(n_fock)
a_dag = qt.create(n_fock)

# convert coherent state density matrix to vector
rho_coh_asvec = operator_to_vector(rho_coherent)
# rho_coh_asvec = np.reshape(rho_coherent, n_fock**2)

# act loss operator over coherent state 
rho_coh_vec_w_loss = lindblad_dissipator(a) * rho_coh_asvec

# coherent state vector under single photon loss to density matrix
# rho_coh_w_loss = qt.Qobj(np.reshape(rho_coh_vec_w_loss, ((n_fock,n_fock))))
rho_coh_w_loss = vector_to_operator(rho_coh_vec_w_loss)
rho_coh_w_loss = rho_coh_w_loss/rho_coh_w_loss.norm()

bar0 = axes[1].bar(np.arange(0, n_fock)-.5, rho_coh_w_loss.diag())
axes[1].set_xticks(list(np.arange(0, n_fock-1, 1)))
lbl0 = axes[1].set_title("Coherent state under single photon loss")

#%%
# visualizing increment of coherent state Wigner under photon loss / unit time
xvec = np.linspace(-alpha,alpha,200)

# coherent state wigner function
W_coh = wigner(rho_coherent, xvec, xvec)
nrm_W_coh = mpl.colors.Normalize(-W_coh.max(), W_coh.max())

# coherent state wigner function under single photon loss
W_coh_w_loss = wigner(rho_coh_w_loss, xvec, xvec)
nrm_W_coh_w_loss = mpl.colors.Normalize(-W_coh_w_loss.max(), W_coh_w_loss.max())
# axes.contourf(xvec, xvec, W_coh_w_loss, 100, cmap=cm.RdBu, norm=nrm_W_coh)

# plot the results

fig, axes = plt.subplots(1, 2, figsize=(12,3))

# coherent state (without and with loss)

cont0 = axes[0].contourf(xvec, xvec, W_coh, 100, cmap=cm.RdBu, \
                         norm=nrm_W_coh)
# coherent state under single photon loss
cont0 = axes[1].contourf(xvec, xvec, W_coh_w_loss, 100, cmap = cm.RdBu, \
                         norm = nrm_W_coh_w_loss)
    
set_figure_size(5, 2, axes[0])

#%%
