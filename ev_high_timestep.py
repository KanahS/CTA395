#!/usr/bin/env python
# coding: utf-8

import numpy as np
import rebound as rb
import reboundx as rx
import matplotlib.pyplot as plt

high_dir_cope = "/mnt/raid-cita/ksmith/high_bin_plan_files/"

# BINARIES
def system_evolution(archive, extras):
    """ 
    Pulling data from the called upon simulation archive using REBOUND and REBOUNDx to collect evolution of
    the secondary binary in the simulation. 
    """ # UPDATE THIS DESCRIPTION
    sim = None
    sim = rx.SimulationArchive(archive, rebxfilename=extras)
    sim_time = []
    
    # SECONDARY BINARY
    b2_ecc = []
    b2_semi = []
    b2_period = []
    
    for snap in range(len(sim)):
        base = sim[snap][0]
        time = base.t
        
        # SECONDARY BINARY        
        b2_a = base.particles["Binary 2"].a
        b2_e = base.particles["Binary 2"].e
        b2_P = base.particles["Binary 2"].P
        
        b2_ecc.append(b2_e)
        b2_semi.append(b2_a)
        b2_period.append(b2_P)
        sim_time.append(time)
        
        #for key in base.particles.keys():
        #    print(key)
    
    np.save(high_dir_cope+"binary_sim_time.npy", sim_time)
    np.save(high_dir_cope+"binary_ecc.npy", b2_ecc)
    np.save(high_dir_cope+"binary_semi.npy", b2_semi)
    np.save(high_dir_cope+"binary_period.npy", b2_period)
    
    #return sim_time, b2_ecc, b2_semi, b2_period


e_val = 0.702 # INITIAL PLANETARY SEMI-MAJOR AXIS UNDER EVALUATION
a_val = 2.878 # INITIAL BINARY ECCENTRICITY UNDER EVALUATION


###### TURN THIS INTO A FUNCTION THAT DOES EITHER BINARY OR PLANETS
high_sim_arch = f"/mnt/raid-cita/ksmith/cope_high-timestep/sim_archive_Q10.0_eb{e_val}_ap{a_val}.bin"
high_xarch = "/mnt/raid-cita/ksmith/cope_high-timestep/xarchive_single_Qs10.bin"

# binary stuff
#high_sys_evolution = system_evolution(high_sim_arch, high_xarch)
system_evolution(high_sim_arch, high_xarch)
#high_sys_binaries = "/mnt/raid-cita/ksmith/"
#np.save(high_sys_binaries+"high_sim_planets.npy.npy", high_sys_evolution) 

# PLANETS

def planet_evolution(archive, extras):
    """
    This one includes the planet hashes
    """
    
    sim = None
    sim = rx.SimulationArchive(archive, rebxfilename=extras) #rb.SimulationArchive(archive)
    
    planet_semi = []
    planet_ecc = []
    sim_time = []
    
    for planet in range(2,17):
        a = []
        e = []
        t = []
        for snap in range(len(sim)):            
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
     
    np.save(high_dir_cope+"sim_time.npy", sim_time)
    np.save(high_dir_cope+"planet_semi.npy", planet_semi)
    np.save(high_dir_cope+"planet_ecc.npy", planet_ecc)

    #return sim_time, planet_semi, planet_ecc 


high_sim_cste = f"/mnt/raid-cita/ksmith/cste_high-timestep/raw_surv_time_Q10.0_eb{e_val}_ap{a_val}.npy"

# planetary stuff
#plan = planet_evolution(high_sim_cste, high_xarch)
#np.save(high_sim_planets+"high_sim_planets.npy", plan)
planet_evolution(high_sim_cste, high_xarch)
