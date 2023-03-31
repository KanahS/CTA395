#!/usr/bin/env python
# coding: utf-8

import numpy as np
import rebound as rb
import reboundx as rx
import matplotlib.pyplot as plt


e_val = 0.702 # INITIAL PLANETARY SEMI-MAJOR AXIS UNDER EVALUATION
a_val = 2.878 # INITIAL BINARY ECCENTRICITY UNDER EVALUATION

high_xarch = "/mnt/raid-cita/ksmith/cope_high-timestep/xarchive_single_Qs10.bin"

def planet_evolution(archive, extras):

    sim = None
    sim = rx.SimulationArchive(archive, rebxfilename=extras) #rb.SimulationArchive(archive)

    planet_semi = []
    planet_ecc = []
    sim_time = []

    for planet in range(2,17):
        a = []
        e = []
        t = []
        #for snap in range(len(sim)):
        for snap in range(2):    
            try:
                base = sim[snap][0]
                time = base.t

                plan_a = base.particles[planet].a
                plan_e = base.particles[planet].e

                a.append(plan_a)
                e.append(plan_e)
                t.append(time)
                #print("Particle Found")

            except: #ParticleNotFound: # move onto the next iteration/planet
                continue  # DO THE LIST INDEX ERROR
                #print(f"In snap {snap+1} Planet {planet-2} Not Found")

        planet_semi.append(a)
        planet_ecc.append(e)
        sim_time.append(t)
    
    np.save(high_dir_cope+f"sim_length_is_{len(sim)}", len(sim))
    np.save(high_dir_cope+"sim_time.npy", sim_time)
    np.save(high_dir_cope+"planet_semi.npy", planet_semi)
    np.save(high_dir_cope+"planet_ecc.npy", planet_ecc)

    #return sim_time, planet_semi, planet_ecc 


high_sim_cste = f"/mnt/raid-cita/ksmith/cste_high-timestep/raw_surv_time_Q10.0_eb{e_val}_ap{a_val}.npy"

# planetary stuff
#plan = planet_evolution(high_sim_cste, high_xarch)
# np.save(high_sim_planets+"high_sim_planets.npy", plan)
planet_evolution(high_sim_cste, high_xarch)

