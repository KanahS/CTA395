import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx
import astropy.table
import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy.table import Column
import astropy.units as u
import numpy as np
from astropy.io import fits

import math
pi = math.pi

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


# set axes labels to publication quality
import matplotlib as mpl

import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.ticker as ticker

tend = 1000

t1 = Table.read('/mnt/raid-cita/ksmith/cste_high-timestep/class_Q1_high_timestep.fits')
#tinf = Table.read('class_Qinf_high_timestep.fits')

planet_list = [f'Planet {i}' for i in range(1,16)]
binary_list = ['Binary 2']
name_list = binary_list + planet_list
sma_list = [name+' sma' for name in name_list]
ecc_list = [name+' ecc' for name in name_list]

at_1 = [t1[sma] for sma in sma_list]
#at_inf = [tinf[sma] for sma in sma_list]
et_1 = [t1[ecc] for ecc in ecc_list]
#et_inf = [tinf[ecc] for ecc in ecc_list]
time_1 = t1['times']
#time_inf = tinf['times']

fs_tick = 10
fs_lbl = 12
mpl.rcParams['xtick.labelsize']=fs_tick
mpl.rcParams['ytick.labelsize']=fs_tick
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lw_main = 0.2
lw_thick = 2.0
agval = 0.5
tlen = 5
twid = 1
dots = 0.1
fs_lgn = 5

fig, ax = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row', figsize=(7,7))

for i in range(16):
  ax[0,0].plot(time_1/2./np.pi, at_1[i], label=name_list[i])
#ax[0,0].set_xlabel(r'Orbits')
ax[0,0].set_ylabel(r'Semi-major Axis ($a_{\rm b}$)')
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].set_xlim(1,tend)
ax[0,0].legend(fontsize=fs_lgn)

for i in range(16):
  ax[1,0].plot(time_1/2./np.pi, et_1[i], label=name_list[i])
ax[1,0].set_xlabel(r'Orbits')
ax[1,0].set_ylabel(r'Eccentricity')
ax[1,0].set_xscale('log')
ax[1,0].set_xlim(1,tend)
ax[1,0].legend(fontsize=fs_lgn)

#for i in range(16):
#  ax[0,1].plot(time_inf/2./np.pi, at_inf[i], label=name_list[i])
#ax[0,0].set_xlabel(r'Orbits')
#ax[0,1].set_ylabel(r'Semi-major Axis ($a_{\rm b}$)')
ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[0,1].set_xlim(1,tend)
ax[0,1].legend(fontsize=fs_lgn)

#for i in range(16):
#  ax[1,1].plot(time_inf/2./np.pi, et_inf[i], label=name_list[i])
ax[1,1].set_xlabel(r'Orbits')
#ax[1,1].set_ylabel(r'Eccentricity')
ax[1,1].set_xscale('log')
ax[1,1].set_xlim(1,tend)
ax[1,1].legend(fontsize=fs_lgn)

plt.tight_layout()
plt.savefig("/mnt/raid-cita/ksmith/data_read_Aug4_2023.png", dpi=250)
plt.show()
