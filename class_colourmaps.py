import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx
from mpl_toolkits.axes_grid1 import make_axes_locatable

tup_num = 50
e_b = np.linspace(0, 0.8, tup_num)
a_p = np.linspace(1, 5, tup_num)

Np = 15

Qex = []
tidal_lag = [1,2,3,4,np.inf]
for x in tidal_lag:
    Q = 10**x
    Qex.append(Q)

tup_list = []

for Q in Qex:
    for e in e_b:
        for a in a_p:
            tup_list.append((Q,e,a,Np))

Ne = len(e_b)
Na = len(a_p)
Nq = len(Qex)

# THESE ARE THE COLOURMAPS GENERATED FOR EACH Q VALUE DURING THE SIMULATION RUNS 

#sunnyQ1 = np.load("/mnt/raid-cita/ksmith/all_Q_colour_maps/map_tup50plan15_Q[10].npy")
#sunnyQ2 = np.load("/mnt/raid-cita/ksmith/all_Q_colour_maps/map_tup50plan15_Q[100].npy")
#sunnyQ3 = np.load("/mnt/raid-cita/ksmith/all_Q_colour_maps/map_tup50plan15_Q[1000].npy")
#sunnyQ4 = np.load("/mnt/raid-cita/ksmith/all_Q_colour_maps/map_tup50plan15_Q[10000].npy")
#sunnyQinf = np.load("/mnt/raid-cita/ksmith/all_Q_colour_maps/map_tup50plan15_Q[inf].npy")

sunnyplot = np.concatenate((sunnyQ1,sunnyQ2,sunnyQ3,sunnyQ4,sunnyQinf))
# len(sunnyplot) = 12500

def stability(sunny):
    """
    Making colourplots
    """
    fig, ax = plt.subplots(1, Nq, figsize=(40,8), sharey=True)#, constrained_layout=True)

    ax = ax.ravel()
    SurvTimeAll = np.reshape(sunny, [Nq,Ne,Na])
    SurvTimeArr = [SurvTimeAll[i,:,:] for i in range(Nq)]
    colours = plt.colormaps["summer"] # cividis summer

    for i in range(Nq):                                  # divide here
        pcm = ax[i].pcolormesh(e_b, a_p, SurvTimeArr[i].T, vmin=0, vmax=1e4*2*np.pi, cmap=colours, shading='auto')

        a_b = 2.278 + 3.824*e_b - 1.71*(e_b**2) #blue
        a_c = 1.6 + 5.1*e_b + (- 2.22*(e_b**2)) + 4.12*0.5 + (- 4.27*e_b*0.5) + (- 5.09*(0.5**2)) + 4.61*(e_b**2)*(0.5**2) #green

        ax[i].plot(e_b, a_c, "o", color="crimson",label="$a_{c}(e_{b},\mu)$ = 1.6 + 5.1$e_{b}$ - 2.22$e_{b}^{2}$ + 4.12$\mu$ - 4.27$e_{b}$$\mu$ - 5.09$\mu^{2}$ + 4.61$e_{b}^{2}$$\mu^{2}$, $\mu$ = 0.5")
        ax[i].plot(e_b, a_c, color='crimson')

        ax[i].plot(e_b, a_b, "o", color="cornflowerblue",label="$a_{c}(e_{b})$ $\simeq$ 2.278 + 3.824$e_{b}$ - 1.71${e_{b}}^{2}$")
        ax[i].plot(e_b, a_b, color='cornflowerblue')

        ax[i].annotate(xy=(0.3,2), xytext=(0.4,2), text="$Unstable$", fontsize=30, bbox= dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=2, alpha=0.7))
        ax[i].annotate(xy=(0.1,4), xytext=(0.1,4), text="$Stable$", fontsize=30, bbox= dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=2, alpha=0.7))


        if Qex[i] == np.inf:
            ax[i].set_title(r'Q $\to$ $\infty$, No Tides', size=25)
        else:
            ax[i].set_title('Q={:.1e}'.format(Qex[i]), size=25)
        #ax[i].set_xlim(0.0,0.7)
        #ax[i].set_ylim(1,5)
        ax[i].tick_params(axis="both",labelsize=20)
        #ax[i].legend(prop={"size":20},loc="lower center")

    ax[0].set_ylabel('Initial Planetary Semi-Major Axis, $a_{p}$ [AU]',fontsize=30)
    fig.suptitle("Stability Criterion with Varying Tidal Force", fontsize=32)
    fig.supxlabel("Initial Binary Eccentricity, $e_{b}$", fontsize=30)
    lines, labels = ax[-1].get_legend_handles_labels()
    fig.legend(lines, labels, loc="lower center", prop={"size":30}, bbox_to_anchor=(0,-0.3,1,1))

    # creating invisible axis to plot colorbar
    axis = fig.add_axes([0,-0.5,1,1])
    axis.set_visible(False)

    cbar = plt.colorbar(pcm, ax=axis, location="bottom")
    #cbar.values([0,1e4*2*np.pi])
    cbar.set_label(label= 'Simulation Ejection Time [code units]', size=25) # CHANGE THIS NAME and change numbers on axis to be more easily undestood in terms of rebound time values
    cbar.ax.tick_params(labelsize=25)

    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(hspace=0.3)
    # MARTIN & fITZMAURICE fig 2 overplots resonances

    #plt.savefig("/home/ksmith/Downloads/AST425_colourmap.pdf", bbox_inches="tight")
    plt.show()

stability(sunnyplot)
