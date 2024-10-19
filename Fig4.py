
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
from OBSmatch import replace_detection_limit, get_OBSvar_info, remove_missing_pts

def get_OBSvar_in_srcregions(SAVEDIR, GV_flts, WB_flts, OBSvar_str, WINDS, prs_bnds, src_strs, max_days, int_choice, do_WB=1):    

    
    if do_WB: 
        
        
        
        OBSvar_traj_EA_raw_GV, CC_dep_EA_GV, int_lon, int_lat, pct_EA_GV = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[0], max_days, '', '')
        OBSvar_traj_SA_raw_GV, CC_dep_SA_GV, int_lon, int_lat, pct_SA_GV = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[1], max_days, '', '')   
        OBSvar_traj_EA_raw_WB, CC_dep_EA_WB, int_lon, int_lat, pct_EA_WB = get_OBSvar_info(WB_flts, SAVEDIR, OBSvar_str, 'WB57', WINDS, prs_bnds, int_choice, src_strs[0], max_days, '', '')
        OBSvar_traj_SA_raw_WB, CC_dep_SA_WB, int_lon, int_lat, pct_SA_WB = get_OBSvar_info(WB_flts, SAVEDIR, OBSvar_str, 'WB57', WINDS, prs_bnds, int_choice, src_strs[1], max_days, '', '') 
        
        CC_dep_EA = np.asarray(list(CC_dep_EA_GV) + list(CC_dep_EA_WB))
        CC_dep_SA = np.asarray(list(CC_dep_SA_GV) + list(CC_dep_SA_WB))
        
        pct_EA = len(CC_dep_EA[CC_dep_EA > -900])*100./(len(CC_dep_EA_GV) + len(CC_dep_EA_WB))
        pct_SA = len(CC_dep_SA[CC_dep_SA > -900])*100./(len(CC_dep_SA_GV) + len(CC_dep_SA_WB))      
    
        OBSvar_traj_EA = remove_missing_pts(OBSvar_traj_EA_raw_GV + OBSvar_traj_EA_raw_WB, CC_dep_EA)
        OBSvar_traj_SA = remove_missing_pts(OBSvar_traj_SA_raw_GV + OBSvar_traj_SA_raw_WB, CC_dep_SA)

        
    else: 
        
        OBSvar_traj_EA_raw, CC_dep_EA, int_lon, int_lat, pct_EA = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[0], max_days, '', '')
        OBSvar_traj_SA_raw, CC_dep_SA, int_lon, int_lat, pct_SA = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[1], max_days, '', '')   
        
        CC_dep_EA = np.asarray(CC_dep_EA)
        CC_dep_SA = np.asarray(CC_dep_SA)  

        OBSvar_traj_EA = remove_missing_pts(OBSvar_traj_EA_raw, CC_dep_EA) 
        OBSvar_traj_SA = remove_missing_pts(OBSvar_traj_SA_raw, CC_dep_SA)  

    return OBSvar_traj_EA, OBSvar_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA
    
def plot_EA_SA_hist(ax, var_EA, var_SA, pct_EA, pct_SA, bins, rng, unit):
    
    OBSdist_EA, binz = np.histogram(var_EA, bins = bins, range = rng, density = True)
    OBSdist_EA = OBSdist_EA*pct_EA*(rng[1]-rng[0])/bins  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_EA = ax.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA, histtype = 'step', label='East Asia, ' + str(round(np.mean(var_EA[var_EA>-900]),1)) + ' ' + unit, linewidth=3)

    OBSdist_SA, binz = np.histogram(var_SA, bins = bins, range = rng, density = True)
    OBSdist_SA = OBSdist_SA*pct_SA*(rng[1]-rng[0])/bins  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage 
    hist_SA = ax.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA, histtype = 'step', label='South Asia, ' + str(round(np.mean(var_SA[var_SA>-900]),1)) + ' ' + unit, linewidth=3)    
    
    ax.legend()
    

WINDS = 'ERA5_kin'
SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'
#GV_FLT_DATES = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
GV_FLT_DATES = ['20220730','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV TOGA dates
WB_FLT_DATES = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']  # WB57 dates
#FLT_DATES = ['20220819']
prs_bnds = [0,500]  # The pressure surface below which we will not count trajectories
max_days = 30
src_strs = ['[95, 135, 30, 50]','[65, 95, 18, 35]']

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)

CO_EA, CO_SA, CC_dep_EA, CC_dep_SA, pct_EA_CO, pct_SA_CO = \
    get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, 'CO', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
DCM_EA, DCM_SA, CC_dep_EA, CC_dep_SA, pct_EA_DCM, pct_SA_DCM = \
    get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, 'CH2CL2', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
SO2_EA, SO2_SA, CC_dep_EA, CC_dep_SA, pct_EA_SO2, pct_SA_SO2 = \
    get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, 'SO2', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
Tol_EA_raw, Tol_SA_raw, CC_dep_EA, CC_dep_SA, pct_EA_Tol, pct_SA_Tol = \
    get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, 'Toluene', WINDS, prs_bnds, src_strs, max_days, 'finalcld', do_WB=0)

Tol_EA = replace_detection_limit(Tol_EA_raw, -888., 0.05)
Tol_SA = replace_detection_limit(Tol_SA_raw, -888., 0.05)

print(pct_EA_CO)
print(pct_EA_DCM)
print(pct_EA_SO2)
print(pct_EA_Tol)
    
plot_EA_SA_hist(ax1, CO_EA, CO_SA, pct_EA_CO, pct_SA_CO, 32, (0,320), 'ppbv')
plot_EA_SA_hist(ax2, DCM_EA, DCM_SA, pct_EA_CO, pct_SA_CO, 35, (0,700), 'pptv')
plot_EA_SA_hist(ax3, SO2_EA, SO2_SA, pct_EA_CO, pct_SA_CO, 60, (0,300), 'pptv')
plot_EA_SA_hist(ax4, Tol_EA, Tol_SA, pct_EA_Tol, pct_SA_Tol, 50, (0,100), 'pptv')  # The percentage is different for toluene because it only is measured on the GV

ax1.set_xlabel('CO (ppbv)')
ax1.set_ylabel('Fraction of airborne sampling (%)')

ax2.set_xlabel('CH$_2$Cl$_2$ (pptv)')
ax2.set_ylabel('Fraction of airborne sampling (%)')

ax3.set_xlabel('SO$_2$ (pptv)')
ax3.set_ylabel('Fraction of airborne sampling (%)')

ax4.set_xlabel('Toluene (pptv)')
ax4.set_ylabel('Fraction of airborne sampling (%)')

ax3.set_yscale('log')
ax4.set_yscale('log')

ax1.set_xlim(0,320)
ax2.set_xlim(0,650)
ax3.set_xlim(0,260)
ax4.set_xlim(0,100)

ax3.set_ylim(0.01,20)
ax4.set_ylim(0.01,20)

plt.figtext(0.42,0.77,'(a)',fontsize=20,fontweight='bold')
plt.figtext(0.83,0.77,'(b)',fontsize=20,fontweight='bold')
plt.figtext(0.42,0.35,'(c)',fontsize=20,fontweight='bold')
plt.figtext(0.83,0.35,'(d)',fontsize=20,fontweight='bold')

plt.savefig(PLOTDIR + 'Fig4.png', dpi=300)
plt.savefig(PLOTDIR + 'Fig4.pdf', dpi=300)
plt.close('all')

#if OBSvar_str == 'Toluene':



