#!/usr/bin/env python
# coding: utf-8

import numpy as np
import rebound as rb
import reboundx as rx
import matplotlib.pyplot as plt

tup_num = 50
a_hist = np.linspace(1, 5, tup_num)
e_hist = np.linspace(0, 0.8, tup_num)

a_string = []
e_string = []

for num in a_hist:
    a_string.append(format(num, ".3f"))

for num in e_hist:
    e_string.append(format(num, ".3f"))

hist_sim_arch_directory = "/mnt/raid-cita/ksmith/cope/"
hist_xarch = "/mnt/raid-cita/ksmith/cope/xarchive_single_Qs10.bin"

def planet_hist(sim_arch_directory, xarch):
    last_snapshot = 63
    first_snapshot = 0

    ep_final = []
    ap_fmi = [] # final minus initial 

    for num in range(tup_num): # tup_num = 50
        eb_i = e_hist[num]
        ap_i = a_hist[num]
        # mu = 0.5
        a_c = 1.6 + 5.1*eb_i + (- 2.22*(eb_i**2)) + 4.12*0.5 + (- 4.27*eb_i*0.5) + (- 5.09*(0.5**2)) + 4.61*(eb_i**2)*(0.5**2)
    
        archive = hist_sim_arch_directory + f"sim_archive_Q10.0_eb{e_string[num]}_ap{a_string[num]}.npy"
        if ap_i <= a_c:
            sim = None
            sim = rx.SimulationArchive(archive, rebxfilename=xarch)
        
            for planet in range(2,17):
                if len(sim) == last_snapshot: # last snap that there should be for all full length simulations 
                    a_initial = sim[final_snapshot][0].particles[planet].a
                    a_final = sim[last_snapshot][0].particles[planet].a
                    a_final_minus_initial = a_final - a_initial
                    ap_fmi.append(a_final_minus_initial)
                
                    e_final = sim[last_snapshot][0].particles[planet].e
                    ep_final.append(e_final)
    np.save("/mnt/raid-cita/ksmith/planet_hist_ep_final.npy", ep_final)
    np.save("/mnt/raid-cita/ksmith/planet_hist_ap_fmi.npy", ap_fmi)

planet_hist(hist_sim_arch_directory, hist_xarch)
