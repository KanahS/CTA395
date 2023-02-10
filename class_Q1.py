#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx

tup_num = 50   
e_b = np.linspace(0, 0.8, tup_num)
a_p = np.linspace(1, 5, tup_num)
Np = 15

Qex = []
tidal_lag = [1]
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
    for i in range(Np):
        f_plan = np.random.rand()*2*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan, hash = f"Planet {i}") # MOST RECENTLY UPDATED FEB 9TH 2023
        
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
    
    directory_orbit = "/mnt/raid-cita/ksmith/cope/" # COPE
    rebx.save(directory_orbit+f"xarchive_single_Qs{Q}.bin")# rebx archive # COPE  # TURN BACK ON IF RUNNING ANOTHER SIM FOR DIFFERENT Q VALUES OR OTHER CHANGES
    #filename_orbit = r"eb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_Q{:.1f}_tau{:.4f}.npy".format(eb,ap,Np,tup_num,Q,tau) # COPE
    filename_orbit = r"sim_archive_Q{:.1f}_eb{:.3f}_ap{:.3f}.bin".format(Q,eb,ap)# reb archive # cope 
    sim.automateSimulationArchive(directory_orbit+filename_orbit, interval=1e3, deletefile=True) # COPE

    
    #integrate
    N_times = int(10000) 
    N_orbit = (1e4)*2*np.pi
    times = np.linspace(0,N_orbit,N_times)

    #array for survival times
    surv = []
    nbs = []
    ebs = []# SAVE ME, ONLY FOR HIGH E's
    #surv = np.zeros(Np) # CHANGE

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
        if sim.N <= 2:
            break
    N_ejected = len(surv) ########################################################
    surv = np.pad(surv, Np - N_ejected)[Np - N_ejected:]
    surv[(surv==0)] = time
   
    # Saving raw survival times
    #directory_surv = "/mnt/raid-cita/ksmith/cste/" #CSTE
    #file_surv = "raw_surv_time_Q{:.1f}_eb{:.3f}_ap{:.3f}.npy".format(Q,eb,ap) #CSTE
    #np.savetxt(directory_surv+file_surv, surv) #CSTE

    # SAVING BINARY ECCENTRICITIES W/ ebs, should only end up w/ 50 files as ap changes and eb remains the same, so the overwrite each other, don't bother trying to selectively have only 1 file for 1 eb
    file_ebs = "binary_evolution_initial_ebs{:.3f}.npy".format(eb) # COPE
    np.savetxt(directory_orbit,file_ebs,ebs)
    
    return np.mean(surv)
   
pool = rb.InterruptiblePool(processes=32)
mapping = pool.map(func= survival, iterable= tup_list)

# SAVING MAP
#directory_surv = "/mnt/raid-cita/ksmith/cste/"# CSTE
#npy_surv = f"map_tup{tup_num}plan{Np}_Q{Qex[:]}.npy" # CSTE
#np.savetxt(directory_surv+npy_surv, mapping) # CSTE

# FINISH MESSAGE
directory_test = '/mnt/raid-cita/ksmith/'# DONE
completed = 'The simulation finished!'# DONE
name = 'map_Q1_DONE' # DONE
np.save(directory_test+name, completed) # DONE
