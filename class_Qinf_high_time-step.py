#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx
import astropy.table
import astropy.units as u
from astropy.table import Table
from astropy.table import Column

tup_num = 50   # eb = 0.702 and ap = 2.878
e_b = [0.702]
a_p = [2.878]
Np = 15

Qex = []
tidal_lag = [np.inf]
for x in tidal_lag:
    Q = 10**x
    Qex.append(Q)

tup_list = []

for Q in Qex:
    for e in e_b:
        for a in a_p:
            tup_list.append((Q,e,a,Np))

Nq = len(Qex)
Ne = len(e_b)
Na = len(a_p)

def survival(initial):
    Q, eb, ap, Np = initial[0], initial[1], initial[2], initial[3]
    sim = rb.Simulation()
    sim.integrator = "whfast"

    mu = 0.5
    m1 = 1
    m2 = abs((m1*mu)/(1-mu))

    sim.add(m=m1, hash="Binary 1")
    sim.add(m=m2, a=1, e= eb, hash="Binary 2")

    R_star = 0.1*(1-0.5)**(1/3)
    sim.particles[0].r = R_star
    sim.particles[1].r = R_star

    #initializing Np massless planets
    list1 = np.arange(1,16)
    for i in list1:
        f_plan = np.random.rand()*2*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan, hash = f"Planet {i}") # MOST RECENTLY UPDATED FEB 16TH 2023 to iterate over both loops

    p_hash_list = [p.hash for p in sim.particles[1:]]
    p_hash_bin = ["Binary 2"]
    p_hash_planet = [f"Planet {i}" for i in range(1, Np+1)]
    p_hash_names = p_hash_bin + p_hash_planet

    a0_planet = [ap for i in range(16)]
    a0_bin = [1.0]
    e0_planet = [0. for i in range(16)]
    e0_bin = [eb]
    a0_data = a0_bin + a0_planet
    e0_data = e0_bin + e0_planet


    #array to keep track of survival times
    sim.move_to_com()

    # Adding Tidal Elements
    rebx = rx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    ps = sim.particles

    k2 = ps[0].params["tctl_k2"] = 0.035
    nb = ps[1].n

    if type(Q) != int:
        tau = ps[0].params["tctl_tau"] = 0
    else:
        tau = ps[0].params["tctl_tau"] = 3/(2*Q*k2*nb)

    # SIMULATION ARCHIVE ###filename_orbit = r"eb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_Q{:.1f}_tau{:.4f}.npy".format(eb,ap,Np,tup_num,Q,tau) # COPE
    #directory_orbit = "/mnt/raid-cita/ksmith/cope/" # cope
    #directory_orbit_high = "/mnt/raid-cita/ksmith/cope_high-timestep/" # cope_high-timestep
    #rebx.save(directory_orbit_high+f"test_reg_time-step_xarchive_single_Qs{Q}.bin")# rebx archive # COPE  # TURN BACK ON IF RUNNING ANOTHER SIM FOR DIFFERENT Q VALUES OR OTHER CHANGES
    #filename_orbit = r"test_reg_time-step_sim_archive_Q{:.1f}_eb{:.3f}_ap{:.3f}.bin".format(Q,eb,ap)# reb archive # cope
    #sim.automateSimulationArchive(directory_orbit_high+filename_orbit, interval=1e-3, deletefile=True) # COPE # for the old stuff do interval = 1e-3, for the new stuff do interval = 1e-2

    #*****************

    t = Table()

    #integrate
    #N_times = int(1000)
    N_orbit = (1e4)*2*np.pi  #(1e4)*2*np.pi for the regular stuff # for the new stuff, do (1e3)*2*np.pi
    T_orbit = 2.*np.pi*N_orbit
    N_times = 10*int(N_orbit)
    #         ***********
    #times = np.linspace(0,N_orbit,N_times)
    times = np.linspace(0, T_orbit, N_times)

    #array for survival times
    surv = []
    nbs = []
    ebs = []# SAVE ME, ONLY FOR HIGH E's
    #surv = np.zeros(Np) # CHANGE

    at = [[] for i in range(16)]
    et = [[] for i in range(16)]

    times_rec = []


    for i, time in enumerate(times):
        nb = ps[1].n
        nbs.append(nb)
        r_eb = ps[1].e
        ebs.append(r_eb)
        N_re = (1+(15./2.)*r_eb**2+(45./8.)*r_eb**4+(5./16.)*r_eb**6)/(1-r_eb**2)**6
        Omega_re = (1+3*r_eb**2+(3./8.)*r_eb**4)/(1-r_eb**2)**(9./2.)
        ps[0].params["Omega"] = N_re/Omega_re*nb

        sim.integrate(time, exact_finish_time=0)

        for num in reversed(range(2, sim.N)):
            p = sim.particles[num]
            o = p.calculate_orbit()
            p0 = sim.particles[0]
            p1 = sim.particles[1]
            d = np.sqrt(p.x**2 + p.y**2)
            d0 = np.sqrt((p.x - p0.x)**2 + (p.y - p0.y)**2)
            d1 = np.sqrt((p.x - p1.x)**2 + (p.y - p1.y)**2)

            if d > 20 or d0 < 0.25 or d1 < 0.25 or o.e > 1:
                #surv[num-2] = time  # CHANGE
                sim.remove(num)
                surv.append(time)

        p_hash_surv = [p.hash.value for p in sim.particles]
        for k,ph in enumerate(p_hash_list):
            if ph.value in p_hash_surv:
                at[k].append(sim.particles[ph].a)
                et[k].append(sim.particles[ph].e)
            else:
                at[k].append(0.)
                et[k].append(0.)

        if sim.N <= 2:
            break
    N_ejected = len(surv) ########################################################
    surv = np.pad(surv, Np - N_ejected)[Np - N_ejected:]
    surv[(surv==0)] = time

    t.add_column(Column(name="times", data=times))
    for k in range(16):
        t.add_column(Column(name=p_hash_names[k]+ " sma", data=at[k]))
        t.add_column(Column(name=p_hash_names[k]+ " ecc", data=et[k]))


    # SAVING RAW SURVIVAL TIMES
    #directory_surv = "/mnt/raid-cita/ksmith/cste_Q1/" #cste
    #directory_high = "/mnt/raid-cita/ksmith/cste_high-timestep/" # cste_high-timestep
    #file_surv = r"test_reg_time-step_raw_surv_time_Q{:.1f}_eb{:.3f}_ap{:.3f}.npy".format(Q,eb,ap) #cste
    #np.save(directory_high+file_surv, surv) #cste # CHANGE FROM BEING suvr OR high

    t.write("/mnt/raid-cita/ksmith/cope_high-timestep/class_Qinf_high_timestep.fits", overwrite=True)

    # SAVING BINARY ECCENTRICITIES W/ ebs, should only end up w/ 50 files as ap changes and eb remains the same, so the overwrite each other, don't bother trying to selectively have only 1 file for 1 eb
    ## DONT THINK I NEED THESE, CAN'T REMEBER (march 20th, 2023)
    #file_ebs = r"binary_evolution-initial_ebs{:.3f}_Q{Q}.npy".format(eb,Q) # cope
    #np.save("/mnt/raid-cita/ksmith/cope/"+file_ebs,ebs) # cope

    return np.mean(surv)

pool = rb.InterruptiblePool(processes=32)
mapping = pool.map(func= survival, iterable= tup_list)

# SAVING MAP
#directory_map = "/mnt/raid-cita/ksmith/all_Q_colour_maps/"# CSTE
#npy_surv = f"map_tup{tup_num}plan{Np}_Q{Qex[:]}.npy" # CSTE
#np.save(directory_map+npy_surv, mapping) # CSTE

# FINISH MESSAGE
directory_test = '/mnt/raid-cita/ksmith/'# DONE
completed = 'The simulation finished!'# DONE
print(completed)
name2 = "Qinf_high-time-step_DONE"
np.save(directory_test+name2, completed) # DONE
