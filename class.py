#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx

tup_num = 50   
e_b = np.linspace(0, 0.7, tup_num)
a_p = np.linspace(1, 5, tup_num)
Np = 15

Qex = []
tidal_lag = [None,2,3,4]
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
    for p in range(Np):
        f_plan = np.random.rand()*2.*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan)
    
    #array to keep track of survival times
    sim.move_to_com()

    # Adding Tidal Elements
    rebx = rx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)
    
    ps = sim.particles
    
    k2 = ps[0].params["tctl_k1"] = 0.035
    nb = ps[1].n
    
    if type(Q) != int:    
        tau = ps[0].params["tctl_tau"] = 0 
    else:
        tau = ps[0].params["tctl_tau"] = 3/(2*Q*k2*nb)
        
    #directory_orbit = "/mnt/raid-cita/ksmith/COPE_SUNNY_CONTROL/" # CONTROL
    #directory_orbit = '/mnt/raid-cita/ksmith/COPE_SUNNY_FINE_LONG/' # TIDES
    #rebx.save(directory_orbit+"xarchive_SUNNY_CONTROL.bin") # CONTROL
    #rebx.save(directory_orbit+r"xarchive_Q{Q}_SUNNY_FINE_LONG.bin") # TURN BACK ON IF RUNNING ANOTHER SIM FOR DIFFERENT Q VALUES OR OTHER CHANGES
    #filename_orbit = r"COPE_SUNNY_CONTROL_eb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_tau{:.1f}.bin".format(eb,ap,Np,tup_num,tau) # CONTROL
    #filename_orbit = r"COPE_SUNNY_FINE_LONG_eb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_Q{:.1f}_tau{:.4f}.bin".format(eb,ap,Np,tup_num,Q,tau) # TIDES
    #sim.automateSimulationArchive(directory_orbit+filename_orbit, interval=1e3, deletefile=True)

    
    #integrate
    N_times = int(10000) 
    N_orbit = (1e4)*2*np.pi
    times = np.linspace(0,N_orbit,N_times)

    #array for survival times
    surv = np.zeros(Np)

    for i, time in enumerate(times):

        nb = ps[1].n
        r_eb = ps[1].e
        
        N_re = (1+(15./2.)*r_eb**2+(45./8.)*r_eb**4+(5./16.)*r_eb**6)/(1-r_eb**2)**6
        Omega_re = (1+3*r_eb**2+(3./8.)*r_eb**4)/(1-r_eb**2)**(9./2.)
        ps[0].params["Omega"] = N_re/Omega_re*nb

        sim.integrate(time, exact_finish_time=0)

        for num in reversed(range(2, sim.N)):
            p = sim.particles[num] 

            if (p.x**2 + p.y**2) > (100)**2 or p.e > 1:
                surv[num-2] = time
               # print(f'removing planet {num}')
                sim.remove(num)

            if sim.N==2:
                break
    surv[(surv==0)] = time
    
    #print(f'simulation finished, {len(sim.particles)-2} planets remaining')
   
    # Saving raw survival times
    #directory_surv = '/mnt/raid-cita/ksmith/CSTE_SUNNY_FINE_LONG/' # TIDES
    #directory_surv = "/mnt/raid-cita/ksmith/CSTE_SUNNY_CONTROL/" # CONTROL    
    #np.savetxt(directory_surv+f'SUNNY_LONG_raw_survival_time_eb{eb}_ap{ap}_Q{Q}.npy', surv) # TIDES
    #np.savetxt(directory_surv+f"SUNNY_CONTROL_raw_survival_times_eb{eb}_ap{ap}.npy", surv) # CONTROL    
    return np.mean(surv)
   
pool = rb.InterruptiblePool()
mapping = pool.map(func= survival, iterable= tup_list)

directory_surv = "/mnt/raid-cita/ksmith/CSTE_SUNNY_FINE_LONG/"
#npy_surv = f"CSTE_SUNNY_CONTROL_map_tup{tup_num}plan{Np}.npy"
npy_surv = f"CSTE_SUNNY_LONG_map_tup{tup_num}plan{Np}_Qs{Qex}.npy" # TIDES
#bin_surv = f'CSTE_SUNNY_map_tup{tup_num}plan{Np}_Qi{Qex[0]}_Qf{Qex[-1]}.bin'
#txt_surv = f'CSTE_SUNNY_map_tup{tup_num}plan{Np}_Qi{Qex[0]}_Qf{Qex[-1]}.txt'

np.savetxt(directory_surv+npy_surv, mapping)
#np.savetxt(directory_surv+bin_surv, mapping)
#np.savetxt(directory_surv+txt_surv, mapping)

directory_test = '/mnt/raid-cita/ksmith/'
completed = 'The simulation finished!'
#name = "SUNNY_CONTROL_DONE"
name = 'SUNNY_FINE_DONE' # TIDES
np.save(directory_test+name, completed)
