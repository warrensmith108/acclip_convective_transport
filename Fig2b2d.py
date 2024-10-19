
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cf
import glob
import netCDF4 as nc4
import jdcal
import scale_hgt as sh
from plot_airplane_traj3d import get_processed_data

PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'  
WINDS = 'GFS_kin'
prs_bnds = [0,500]
max_days = 30

EA_box = [95,135,30,50]
SA_box = [65,95,18,35]

init_lon_GV, init_lat, init_prs, init_jd, CC_dep_g_GV, CC_lon, CC_lat, CC_pct_global_GV \
    = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')
init_lon_WB, init_lat, init_prs, init_jd, CC_dep_g_WB, CC_lon, CC_lat, CC_pct_global_WB \
    = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')

CC_dep_g = CC_dep_g_GV + CC_dep_g_WB
CC_pct_g = len(np.where(np.asarray(CC_dep_g) > -900)[0])*100/(len(init_lon_GV)+len(init_lon_WB))

print(CC_pct_g)

init_lon_GV, init_lat, init_prs, init_jd, CC_dep_EA_GV, CC_lon, CC_lat, CC_pct_EA_GV \
    = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, EA_box, max_days, 'finalcld', '', '')
init_lon_WB, init_lat, init_prs, init_jd, CC_dep_EA_WB, CC_lon, CC_lat, CC_pct_EA_WB \
    = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, EA_box, max_days, 'finalcld', '', '')

CC_dep_EA = CC_dep_EA_GV + CC_dep_EA_WB
CC_pct_EA = len(np.where(np.asarray(CC_dep_EA) > -900)[0])*100/(len(init_lon_GV)+len(init_lon_WB))

print(CC_pct_EA)

init_lon_GV, init_lat, init_prs, init_jd, CC_dep_SA_GV, CC_lon, CC_lat, CC_pct_SA_GV \
    = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, SA_box, max_days, 'finalcld', '', '')
init_lon_WB, init_lat, init_prs, init_jd, CC_dep_SA_WB, CC_lon, CC_lat, CC_pct_SA_WB \
    = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, SA_box, max_days, 'finalcld', '', '')

CC_dep_SA = CC_dep_SA_GV + CC_dep_SA_WB
CC_pct_SA = len(np.where(np.asarray(CC_dep_SA) > -900)[0])*100/(len(init_lon_GV)+len(init_lon_WB))

print(CC_pct_SA)


CC_hist_EA, binz = np.histogram(CC_dep_EA, bins = 30, range = (0,30), density = True)
CC_hist_EA = CC_hist_EA*CC_pct_EA  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
hist_GV = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = CC_hist_EA, histtype = 'step', label='East Asia: ' + str(round(CC_pct_EA,1)) + '%', linewidth=3)

CC_hist_SA, binz = np.histogram(CC_dep_SA, bins = 30, range = (0,30), density = True)
CC_hist_SA = CC_hist_SA*CC_pct_SA  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
hist_GV = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = CC_hist_SA, histtype = 'step', label='South Asia: ' + str(round(CC_pct_SA,1)) + '%', linewidth=3)

plt.xlabel('Transit time (days)', fontsize = 16)
plt.ylabel('Fraction of aircraft sampling (%)', fontsize = 16)
plt.xlim(0,30)
plt.legend(fontsize=16)
plt.figtext(0.6,0.73,'(Global: ' + str(round(CC_pct_g,1)) + '%)', fontsize=16)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.tight_layout()
if WINDS == 'ERA5_kin': plt.savefig(PLOTDIR + 'Both_AllRFs_' + WINDS + '_Fig2d.png', dpi=300)
if WINDS == 'GFS_kin': plt.savefig(PLOTDIR + 'Both_AllRFs_' + WINDS + '_Fig2b.png', dpi=300)
plt.close('all')



