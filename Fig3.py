
#######################
#
#  The purpose of this figure is to show the inverse relationship between CO mixing ratio and transit time from convection
#
#######################

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import glob
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import scale_hgt as sh
import matplotlib
from trz_plots import read_OBSvar_both
#from ACCLIP_Plotter import get_ACCLIP_merge_data

def moving_window_smooth(data, width=1):
    data_copy = np.copy(data)
    smoothed = np.zeros(data_copy.shape)
    for ix in range(data_copy.shape[0]):
        for iy in range(data_copy.shape[1]):
            smoothed[ix,iy] = np.nanmean(data_copy[ix-width:ix+width,iy-width:iy+width])
    return smoothed    
    
def get_pix_data(OBSSAVEDIR, OBSvar_str, prs_bnds, src_region, max_days, int_choice):

    var, dep = read_OBSvar_both(OBSSAVEDIR, '', OBSvar_str, prs_bnds, src_region, max_days, int_choice)

    pixels = np.histogram2d(np.asarray(dep), np.asarray(var), bins = [50,40], range = [[-10,40],[0,400]])
    pixx = np.asarray(pixels[1])
    pixy = np.asarray(pixels[2])
    pixdata = np.transpose(pixels[0])
    pixdata = pixdata*100./len(dep)  # Normalize the fraction of all trajectories    (CHOOSE) 
    #pixdata = pixdata*100./np.max(pixdata)  # Normalize the numbers relative to all encounters      (CHOOSE)
    pixdata[pixdata == 0.0] = pixdata[pixdata == 0.0]-0.01
    
    print(np.min(pixdata), np.max(pixdata))
    
    return pixx, pixy, pixdata   

def tau_adiabat(ax, pixy, tau, conc_inits, col=[0,0,0], style='--'):
    for conc_init in conc_inits:
        #ax.plot(conc_init*np.exp(-pixy/tau), pixy, color=col, linestyle = style, linewidth=0.8)
        ax.plot(pixy, conc_init*np.exp(-pixy/tau), color=col, linestyle = style, linewidth=0.8)

if __name__ == '__main__':
    
    OBSvar_str = 'CO'
    vert_choice = 'z'
    prs_bnds = [0,500]
    src_regions = [[-1,360,-90,90],[65,95,18,35],[95,135,30,50]]
    max_days = 30
    int_choice = 'finalcld'
    OBSSAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
    PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'


    pixx, pixy, pixdata_glob = get_pix_data(OBSSAVEDIR, OBSvar_str, prs_bnds, src_regions[0], max_days, int_choice)
    pixx, pixy, pixdata_SA   = get_pix_data(OBSSAVEDIR, OBSvar_str, prs_bnds, src_regions[1], max_days, int_choice)
    pixx, pixy, pixdata_EA   = get_pix_data(OBSSAVEDIR, OBSvar_str, prs_bnds, src_regions[2], max_days, int_choice)

    midx = (pixx[1:] + pixx[0:-1])/2.
    midy = (pixy[1:] + pixy[0:-1])/2.

    fig = plt.figure(figsize = (5,4))
    ax = fig.add_subplot(1,1,1)

    cmap = plt.get_cmap('Greys')
    cmap.set_under('white')

    #levels = MaxNLocator(nbins=5).tick_values(0.0, 20)
    #levels = [0,2,10,20,30,40,50]
    levels = [0,0.02,0.1,0.25,0.4,0.55,0.7]
    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)

    pixplt = ax.pcolormesh(pixx, pixy, pixdata_glob, cmap = cmap, norm = norm)

    lev1 = 0.02   # 2
    lev2 = 0.30   # 30
    
    ax.contour(midx, midy, moving_window_smooth(pixdata_EA, width=1), levels=[lev1], colors=['red'], linewidths=1)
    ax.contour(midx, midy, moving_window_smooth(pixdata_EA, width=1), levels=[lev2], colors=['red'], linewidths=3)  
    ax.contour(midx, midy, moving_window_smooth(pixdata_SA, width=1), levels=[lev1], colors=['orange'], linewidths=1)  
    ax.contour(midx, midy, moving_window_smooth(pixdata_SA, width=1), levels=[lev2], colors=['orange'], linewidths=3)  
    
     
    # Add some "CO loss adiabats" to the plot
    if OBSvar_str == 'CO': 
        tau_adiabat(ax, np.arange(0,32,1), 28.4, [80,150,220], col='blue')  # 48.2 days uses the reaction rate from JPL (2.4E-13)
        #tau_adiabat(ax, np.arange(0,32,1), 60., np.arange(50,310,50), col='black', style='--')  # Assume 60 d CO lifetime in the ASM
    #    tau_adiabat(ax, np.arange(0,32,1), 60., [100], col='black', style='--')  # Assume 60 d CO lifetime in the ASM
    #    tau_adiabat(ax, np.arange(0,32,1), 10., [150], col='black', style=':')  # Assume 60 d CO lifetime in the ASM
    #    
    #if var_str == 'Prop-Eth': clock_curve(ax, np.arange(0,32,1), 3.65E-8, 7.11E-9, 0.4, 1E6, col=[0,0,0], style='--')   # NOTE: REACTION RATES CONVERTED TO  cm3/molec/day

        
    cbar = plt.colorbar(pixplt,fraction=0.045, pad=0.04, extend = 'max')
    cbar.set_label('Fraction of aircraft sampling (%)')

    ax.set_ylim(0,320)
    ax.set_xlim(0,30)    

    ax.set_ylabel('Observed ' + OBSvar_str + ' (ppbv)')
    ax.set_xlabel('Transit time to convection (days)')
    
    plt.figtext(0.3,0.8,'All Sampling', fontsize=12, fontweight='bold', color='gray')
    plt.figtext(0.3,0.75,'East Asia Convection', fontsize=12, fontweight='bold', color='red')
    plt.figtext(0.3,0.7,'South Asia Convection', fontsize=12, fontweight='bold', color='orange')

    plt.tight_layout()
    plt.savefig(PLOTDIR + 'Fig3.png', dpi=300)
    plt.savefig(PLOTDIR + 'Fig3.pdf', dpi=300)
    plt.close('all')

