#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import pickle as pickle 
import bokeh 
from skimage.transform import rescale, resize

spec_kepler = np.load("/mnt/raid-cita/ksmith/spectrum_kepler.npy") 

def create_echelle_diagram(power, freq, dnu, numax, epsilon, min_interval, max_interval):
    """
    Create the echelle diagram corresponding to the given PSD contained in 'power'
    arguments:
        power = PSD associated with a star
        freq = corresponding frequencies
        dnu =  large frequency seperation
        numax = frequency of max aplitude
        epsilon = epsilon parameter
        min_interval = lower boundary (in number of dnu) of the desired frequency interval (default = 3.)
        max_interval = upper boundary (in number of dnu) of the desired frequency interval (default = 3.)
    """
    freq_l0 = int(numax/dnu)*dnu+(epsilon*dnu)     #frequency of the mode ell = 0, the closest to numax (and larger than numax)
    fbin = np.mean(np.diff(freq))     #size of a frequency bin (interval between 2 consecutive frequencies)
    
    #select the part of the PSD with the oscillations
    power_osc = power[np.where((freq>=(freq_l0-min_interval*dnu)) & (freq <= (freq_l0+max_interval*dnu)))]
    freq_osc = freq[np.where((freq>=(freq_l0-min_interval*dnu)) & (freq <= (freq_l0+max_interval*dnu)))]

    #number of frequency bins in dnu
    number_of_bin_for_dnu=int(dnu/fbin)
    #number of slices of the spectra
    sizeh=int(len(freq_osc)/number_of_bin_for_dnu)

    # makes sure that freqs_osc and power_osc have a length that is an integer number of slices
    freq_osc=freq_osc[:sizeh*number_of_bin_for_dnu]
    power_osc=power_osc[:sizeh*number_of_bin_for_dnu]

    # Echelle diagram
    # Stack the slices in a 2D array
    Z = np.reshape(np.transpose(power_osc), (sizeh, number_of_bin_for_dnu))

    return Z


for num in range(len(spec_kepler)):
    fig = plt.figure(figsize = (15,10))
    
    KIC_ech = spec_kepler[num][0]
    delnu_ech = spec_kepler[num][4]
    numax_ech = spec_kepler[num][3]
    freq_ech = spec_kepler[num][6] # spec_x , change to ech_kepler[num][1] for full PSD, ech_kepler[num][6] for sliced by 5
    power_ech = spec_kepler[num][7] # spec_y , change to ech_kepler[num][2] for full PSD, ech_kepler[num][7] for sliced by 5
    epsilon_ech = spec_kepler[num][5]
    
    ech = create_echelle_diagram(power_ech, freq_ech, delnu_ech, numax_ech, epsilon_ech, min_interval = 3.3, max_interval = 2.7)
    rs = resize(ech, (128, 128))
    
    plt.pcolor(rs, cmap='plasma', vmin=np.min(rs), vmax=np.max(rs))
    plt.axis("off")

    path_save = "/mnt/raid-cita/ksmith/runs/echelle_plasma/"
    img_name = "KIC_{}.png".format(KIC_ech)
    plt.savefig(path_save+img_name, dpi = 200)
    plt.close(fig)
    
