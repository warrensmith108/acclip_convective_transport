


import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import glob
import netCDF4 as nc4
import jdcal
import scale_hgt as sh
from OBSmatch import get_OBSvar_in_srcregions, get_OBSvar_info, remove_missing_pts


SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'
src_strs = ['[105, 135, 40, 50]','[65, 95, 18, 35]']
prs_bnds = [0,500]
WINDS = 'ERA5_kin'
max_days = 30

# We're putting these all on one window so we might as well get that going now
fig = plt.figure(figsize=(12,4.5))
ax1 = plt.subplot(1,3,1)
ax2 = plt.subplot(1,3,2)
ax3 = plt.subplot(1,3,3)


# Read August 6 CO
OBSvar_traj_EA_A6, OBSvar_traj_SA_A6, CC_dep_EA_A6, CC_dep_SA_A6, pct_EA_A6, pct_SA_A6 = get_OBSvar_in_srcregions(SAVEDIR, ['20220806'], ['20220806'], 'CO', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
#OBSvar_traj_EA_A6, OBSvar_traj_SA_A6, CC_dep_EA_A6, CC_dep_SA_A6, pct_EA_A6, pct_SA_A6 = get_OBSvar_in_srcregions(SAVEDIR, ['20220806'], ['20220806'], 'CO', WINDS, prs_bnds, src_strs, max_days, 'finalcld')

# Read August 7 CO
OBSvar_traj_EA_A7, CC_dep_EA_A7, pct_EA_A7, init_prs_EA_A7, pct_EA_A7 = get_OBSvar_info(['20220807'], SAVEDIR, 'CO', 'GV', WINDS, prs_bnds, 'finalcld', src_strs[0], max_days, '', '')
OBSvar_traj_SA_A7, CC_dep_SA_A7, pct_SA_A7, init_prs_SA_A7, pct_SA_A7 = get_OBSvar_info(['20220807'], SAVEDIR, 'CO', 'GV', WINDS, prs_bnds, 'finalcld', src_strs[1], max_days, '', '') 

OBSvar_traj_EA_A7 = np.asarray(remove_missing_pts(OBSvar_traj_EA_A7, CC_dep_EA_A7))
OBSvar_traj_SA_A7 = np.asarray(remove_missing_pts(OBSvar_traj_SA_A7, CC_dep_SA_A7))

# Read Aug 6 and 7 prs
prs_traj_EA, prs_traj_SA, CC_dep_EA, CC_dep_SA, prs_pct_EA, prs_pct_SA = get_OBSvar_in_srcregions(SAVEDIR, ['20220806','20220807'], ['20220806'], 'prs', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
#prs_traj_EA, prs_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = get_OBSvar_in_srcregions(SAVEDIR, ['20220806','20220807'], ['20220806'], 'prs', WINDS, prs_bnds, src_strs, max_days, 'finalcld')


#####################  CO PANEL B

OBSdist_EA_A6, binz = np.histogram(OBSvar_traj_EA_A6, bins = 100, range = (0,500), density = True)
OBSdist_EA_A6 = OBSdist_EA_A6*pct_EA_A6*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio 
hist_EA = ax2.hist(binz[:-1], binz, edgecolor = (79/255., 38/255., 131/255.), weights = OBSdist_EA_A6, histtype = 'step', label='Northeast China Aug 6, ' + str(round(np.mean(OBSvar_traj_EA_A6[OBSvar_traj_EA_A6>-900]),1)) + ' ppbv', linewidth=3)

OBSdist_EA_A7, binz = np.histogram(OBSvar_traj_EA_A7, bins = 50, range = (0,500), density = True)
OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio 
hist_EA = ax2.hist(binz[:-1], binz, edgecolor = (79/255., 38/255., 131/255.), weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', \
    label='Northeast Asia Aug 7, ' + str(round(np.mean(OBSvar_traj_EA_A7[OBSvar_traj_EA_A7>-900]),3)) + ' ppbv', linewidth=3)
    
OBSdist_SA_A6, binz = np.histogram(OBSvar_traj_SA_A6, bins = 100, range = (0,500), density = True)
OBSdist_SA_A6 = OBSdist_SA_A6*pct_SA_A6*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio
hist_SA = ax2.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A6, histtype = 'step', label='South Asia Aug 6, ' + str(round(np.mean(OBSvar_traj_SA_A6[OBSvar_traj_SA_A6>-900]),1)) + ' ppbv', linewidth=3)

OBSdist_SA_A7, binz = np.histogram(OBSvar_traj_SA_A7, bins = 50, range = (0,500), density = True)
OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio
hist_SA = ax2.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7, ' + str(round(np.mean(OBSvar_traj_SA_A7[OBSvar_traj_SA_A7>-900]),3)) + ' ppbv', linewidth=3)

ax2.legend(loc='upper left', fontsize=8)
ax2.set_xlim(0,250)
ax2.set_ylim(0,8)
#ax2.set_title('CO Distributions for August 6-7, 2022', fontsize=10)
ax2.set_xlabel('CO (ppbv)', fontsize=10)
#plt.xlabel('Toluene/Benzene Ratio', fontsize=12)
ax2.set_ylabel("Fraction of airborne sampling (%)", fontsize=10)
#plt.ylabel('Relative Frequency', fontsize=16)
#ax2.set_xtickparams(fontsize=14)
#ax2.set_yticks(fontsize=14)

#plt.savefig(PLOTDIR + 'Fig8b.png', dpi=300)
#plt.close('all')


#######################   TTD PANEL A

OBSdist_EA_A6, binz = np.histogram(CC_dep_EA_A6, bins = 30, range = (0,30), density = True)
OBSdist_EA_A6 = OBSdist_EA_A6*pct_EA_A6
hist_EA = ax1.hist(binz[:-1], binz, edgecolor = (79/255., 38/255., 131/255.), weights = OBSdist_EA_A6, histtype = 'step', label='Northeast China Aug 6', linewidth=3)

OBSdist_EA_A7, binz = np.histogram(CC_dep_EA_A7, bins = 30, range = (0,30), density = True)
OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7
hist_EA = ax1.hist(binz[:-1], binz, edgecolor = (79/255., 38/255., 131/255.), weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', label='Northeast China Aug 7', linewidth=3)

OBSdist_SA_A6, binz = np.histogram(CC_dep_SA_A6, bins = 30, range = (0,30), density = True)
OBSdist_SA_A6 = OBSdist_SA_A6*pct_SA_A6
hist_SA = ax1.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A6, histtype = 'step', label='South Asia Aug 6', linewidth=3)

OBSdist_SA_A7, binz = np.histogram(CC_dep_SA_A7, bins = 30, range = (0,30), density = True)
OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7
hist_SA = ax1.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7', linewidth=3)

ax1.legend(fontsize=8)
ax1.set_xlim(0,30)
#ax1.set_title('Transit Time Distributions for 2022 ACCLIP Sampling', fontsize=12)
#ax1.set_title('Transit Time Distributions for Aug 6-7, 2022', fontsize=10)
ax1.set_xlabel('Transit time (days)', fontsize=10)
ax1.set_ylabel("Fraction of airborne sampling (%)", fontsize=10)
#plt.ylabel('Relative Frequency', fontsize=16)
#ax1.set_xticks(fontsize=14)
#ax1.set_yticks(fontsize=14)

#plt.savefig(PLOTDIR + 'Fig8.png', dpi=300)
#plt.close('all')


########################  INIT PRS PANEL C

# These numbers come from the calc_LZRH_from_pklfile.py file (currently on modeling1)
LZRH_lowq = 132.06
LZRH_med = 142.17
LZRH_upq = 187.75

init_alt_EA = sh.prs2alt(np.asarray(prs_traj_EA), 7.6, 1013.)
init_alt_SA = sh.prs2alt(np.asarray(prs_traj_SA), 7.6, 1013.)
#init_alt_EA_A7 = sh.prs2alt(np.asarray(init_prs_EA_A7), 7.6, 1013.)
#init_alt_SA_A7 = sh.prs2alt(np.asarray(init_prs_SA_A7), 7.6, 1013.)

prsdist_EA, binz = np.histogram(init_alt_EA, bins = 20, range = (0,20), density = True)
prsdist_EA_A6 = prsdist_EA*prs_pct_EA
hist_EA = ax3.hist(binz[:-1], binz, edgecolor = (79/255., 38/255., 131/255.), weights = prsdist_EA_A6, histtype = 'step', label='Northeast China Aug 6-7', linewidth=3, orientation='horizontal')

#OBSdist_EA_A7, binz = np.histogram(init_alt_EA_A7, bins = 20, range = (0,20), density = True)
#OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7
#hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'salmon', weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', label='Northeast China Aug 7', linewidth=3, orientation='horizontal')

prsdist_SA, binz = np.histogram(init_alt_SA, bins = 20, range = (0,20), density = True)
prsdist_SA = prsdist_SA*prs_pct_SA
hist_SA = ax3.hist(binz[:-1], binz, edgecolor = 'orange', weights = prsdist_SA, histtype = 'step', label='South Asia Aug 6-7', linewidth=3, orientation='horizontal')

#OBSdist_SA_A7, binz = np.histogram(init_alt_SA_A7, bins = 20, range = (0,20), density = True)
#OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7
#hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7', linewidth=3, orientation='horizontal')

ax3.fill_between([0,100],[sh.prs2alt(LZRH_lowq, 7.6, 1013.),sh.prs2alt(LZRH_lowq, 7.6, 1013.)],[sh.prs2alt(LZRH_upq, 7.6, 1013.),sh.prs2alt(LZRH_upq, 7.6, 1013.)], color=[0.9,0.9,0.9])
ax3.plot([0,100],[sh.prs2alt(LZRH_med, 7.6, 1013.),sh.prs2alt(LZRH_med, 7.6, 1013.)],color='k')

ax3.legend(fontsize=8)
ax3.set_xlim(0,15.7)
ax3.set_ylim(5,21)
#ax3.set_title('Aircraft Measurement Altitudes for Aug 6-7, 2022', fontsize=10)
ax3.set_ylabel('Aircraft Altitude (km)', fontsize=10)
#plt.xlabel('Toluene/Benzene Ratio', fontsize=12)
ax3.set_xlabel("Fraction of airborne sampling (%)", fontsize=10)
#plt.ylabel('Relative Frequency', fontsize=16)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)




plt.figtext(0.27, 0.71, '(a)', fontweight='bold', fontsize=20)
plt.figtext(0.40, 0.71, '(b)', fontweight='bold', fontsize=20)
plt.figtext(0.94, 0.77, '(c)', fontweight='bold', fontsize=20)

plt.tight_layout()

plt.savefig(PLOTDIR + 'Fig8.png', dpi=300)
plt.savefig(PLOTDIR + 'Fig8.pdf', dpi=300)
plt.close('all')





########################  DEMO OF INIT_THETA INSTEAD FOR PANEL C!!!

#th_traj_EA, th_traj_SA, CC_dep_EA, CC_dep_SA, prs_pct_EA, prs_pct_SA = get_OBSvar_in_srcregions(SAVEDIR, ['20220806','20220807'], ['20220806'], 'theta', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
#
#print(np.min(th_traj_SA), np.max(th_traj_SA))
#
##thdist_EA, binz = np.histogram(th_traj_EA, bins = 20, range = (0,20), density = True)
##thdist_SA, binz = np.histogram(th_traj_SA, bins = 20, range = (0,20), density = True)
#
#dist_SA, binz = np.histogram(th_traj_SA, bins = 24, range = (320,440), density = True)
#dist_SA = dist_SA*prs_pct_SA
#hist_SA = ax3.hist(binz[:-1], binz, edgecolor = 'orange', weights = dist_SA, histtype = 'step', label='South Asia Aug 6-7', linewidth=3, orientation='horizontal')
#
#dist_EA, binz = np.histogram(th_traj_EA, bins = 24, range = (320,440), density = True)
#dist_EA = dist_EA*prs_pct_EA
#hist_EA = ax3.hist(binz[:-1], binz, edgecolor = 'purple', weights = dist_EA, histtype = 'step', label='Northeast China Aug 6-7', linewidth=3, orientation='horizontal')

