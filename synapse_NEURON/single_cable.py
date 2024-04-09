#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:54:45 2022

@author: spiros
"""

from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt

h.load_file("stdrun.hoc")

dend1 = h.Section(name='dend1')
dend1.L = 100
dend1.diam = 2
dend1.nseg = 11
dend1.insert('pas')

dend2 = h.Section(name='dend2')
dend2.L = 100
dend2.diam = 2
dend2.nseg = 11
dend2.insert('pas')

pre_spikes = [100, 150, 200, 250, 300, 350]  # list with presynaptic spikes (ms)
inputs = h.Vector(pre_spikes)
# inputs = h.Vector([100])

# VecStim
vecstim = h.VecStim()
vecstim.play(inputs)

# Synapse - plastic
syn1 = h.AmpaNmdaP(dend1(0.5))
syn1.eampa = 0
syn1.tau1_ampa = 2
syn1.tau2_ampa = 10
syn1.gbar_ampa = .001

syn1.enmda = 0
syn1.tau1_nmda = 8
syn1.tau2_nmda = 100
syn1.gbar_nmda = .001

# parameters of STP -- Depression
stp = 'facilitation'
if stp == 'depression':
    syn1.U0 = .5
    syn1.tauD = 100
    syn1.tauF = 50
elif stp == 'facilitation':
    syn1.U0 = .2
    syn1.tauD = 100
    syn1.tauF = 750


# NetCon - plastic
vs1 = h.NetCon(vecstim, syn1)
vs1.delay = 1.0  # delay in ms
vs1.weight[0] = 1

# Synapse - non-plastic
syn2 = h.AmpaNmda(dend2(0.5))
syn2.eampa = 0
syn2.tau1_ampa = 2
syn2.tau2_ampa = 10
syn2.gbar_ampa = .001

syn2.enmda = 0
syn2.tau1_nmda = 8
syn2.tau2_nmda = 50
syn2.gbar_nmda = .001

# NetCon - plastic
vs2 = h.NetCon(vecstim, syn2)
vs2.delay = 1.0  # delay in ms
vs2.weight[0] = 1

# Record vectors
v_vec1 = h.Vector().record(dend1(0.5)._ref_v)
v_vec2 = h.Vector().record(dend2(0.5)._ref_v)
t_vec = h.Vector().record(h._ref_t)

r_vec = h.Vector().record(syn1._ref_r)
u_vec = h.Vector().record(syn1._ref_u)


# Run simulation
h.finitialize(-70 * mV)
h.continuerun(1000 * ms)

v_vec1_norm = (v_vec1 - v_vec1.min())/(v_vec1.max() - v_vec1.min())
v_vec2_norm = (v_vec2 - v_vec2.min())/(v_vec2.max() - v_vec2.min())

# Plot the results
plt.figure()
plt.plot(t_vec, v_vec1_norm, label='plastic synapse')
plt.plot(t_vec, v_vec2_norm, label='non-plastic synapse')
plt.xlabel('time (ms)')
plt.ylabel('normalized voltage (a.u.)')
plt.legend()
plt.show()

# Plot the results
plt.figure()
plt.plot(t_vec, r_vec, label='resources')
plt.plot(t_vec, u_vec, label='release prob.')
plt.xlabel('time (ms)')
plt.ylabel('prob. (a.u.)')
plt.legend()
plt.show()