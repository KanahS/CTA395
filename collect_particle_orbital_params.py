#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx

particle_hash_names = ["Binary 2", "Planet 1", "Planet 2", "Planet 3", "Planet 4", "Planet 5", 
                       "Planet 6", "Planet 7", "Planet 8", "Planet 9", "Planet 10", "Planet 11", "Planet 12", 
                       "Planet 13", "Planet 14", "Planet 15"] 

def collect_particle_orbital_params(archive):
    """

    """
    # initializing dictionary to collect evolution of all particle orbital elements
    particle_evolution_dict = {"Binary 2": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]},
                            "Planet 1": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 2": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 3": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 4": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]},
                            "Planet 5": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 6": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 7": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 8": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 9": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 10": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]},
                            "Planet 11": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 12": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 13": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 14": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]}, 
                            "Planet 15": {"semi-major_axis":[], 
                                          "eccentricity":[], 
                                          "time":[]},
                            
                            }
    
    sim_archive = rb.SimulationArchive(archive) # using to get number of snapshots
    snap_num = np.arange(len(sim_archive)) # number of snapshots to iterate trhough
    
    for snap in range(len(snap_num)):
        sim = None 
        sim = rb.Simulation(archive, snapshot=snap) 

        for particle_name in particle_hash_names:                            
            # try and except statement, to see if the particle is still in the simulation
            try:            
                # access the particle's orbital elements
                semi = sim.particles[particle_name].a
                ecc = sim.particles[particle_name].e
                time = sim.t
                
                # add ecc, semi, and time to dictionary
                particle_evolution_dict[particle_name]["eccentricity"].append(ecc)
                particle_evolution_dict[particle_name]["semi-major_axis"].append(semi)
                particle_evolution_dict[particle_name]["time"].append(time)
                
            except:
                continue 
                
    return particle_evolution_dict

high_time_step_collection = collect_particle_orbital_params("/mnt/raid-cita/ksmith/cope_high-timestep/sim_archive_Q10.0_eb0.702_ap2.878.bin")
dictionary_directory = "/mnt/raid-cita/ksmith/sim_dictionaries/"

# SAVE THE DIRECTIONARY
np.save(dictionary_directory+"high_time_step_dict_Q10_eb0.702_ap2.878.npy", high_time_step_collection)

np.save("/mnt/raid-cita/ksmith/DICT_DONE.npy", "DICTIONARY COLLECTION DONE")