
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

def get_OBSvar_info(FLT_DATES, SAVEDIR, OBSvar_str, AIRPLANE, WINDS, prs_bnds, int_choice, src_str, max_days, inst_str, med_str):

    OBSvar_traj = []
    dep      = []
    int_lon  = []
    int_lat  = []
    init_prs = []

    for FLT_DATE in FLT_DATES:    
        # Read the file 
        OBS_SAVEFILE = SAVEDIR + OBSvar_str + '_OBSmatch_' + int_choice + '_' + AIRPLANE + '_' + WINDS + inst_str + '_' + FLT_DATE + '_' + str(prs_bnds) + 'hPa_' + str(src_str) + '_' + str(max_days) + 'd' + med_str + '.pkl'
        with open(OBS_SAVEFILE, 'rb') as f:
            curr_OBSvar_traj, curr_dep, curr_int_lon, curr_int_lat, pct = pkl.load(f)
        OBSvar_traj.extend(curr_OBSvar_traj)        
        #BL_dep.extend(curr_BL_dep)
        dep.extend(curr_dep)
        int_lon.extend(curr_int_lon)
        int_lat.extend(curr_int_lat)
        
        #BLCC_dep.extend(curr_BLCC_dep)
        #init_prs.extend(curr_init_prs)

    return OBSvar_traj, np.asarray(dep), np.asarray(int_lon), np.asarray(int_lat), pct
    
def remove_missing_pts(OBSvar_traj_raw, CC_dep):
    OBSvar_traj = []
    for pt in range(len(CC_dep)):
        if CC_dep[pt] > -900 and OBSvar_traj_raw[pt] > -900: OBSvar_traj.append(OBSvar_traj_raw[pt])        
    return np.asarray(OBSvar_traj)

def remove_missing_pts_TT(OBSvar_traj_raw, CC_dep):    # We need a separate point remover for the tracer plots because we don't care about restricting the observations to positive values
    OBSvar_traj = []
    for pt in range(len(CC_dep)):
        if CC_dep[pt] > -900: OBSvar_traj.append(OBSvar_traj_raw[pt])        
    return np.asarray(OBSvar_traj)

def get_OBSvar_in_srcregions(SAVEDIR, GV_flts, WB_flts, OBSvar_str, WINDS, prs_bnds, src_strs, max_days, int_choice):    

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
    
    #init_prs_EA = remove_missing_pts(init_prs_EA_GV + init_prs_EA_WB, CC_dep_EA)
    #init_prs_SA = remove_missing_pts(init_prs_SA_GV + init_prs_SA_WB, CC_dep_SA)

    return OBSvar_traj_EA, OBSvar_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA

def replace_detection_limit(OBSvar_traj_raw, detlim_flag, detlim_replace):
    OBSvar_traj = []
    for pt in range(len(OBSvar_traj_raw)):
        if int(OBSvar_traj_raw[pt]) == detlim_flag: OBSvar_traj.append(detlim_replace)
        else: OBSvar_traj.append(OBSvar_traj_raw[pt])              
    return np.asarray(OBSvar_traj)
        
def get_OBSvar_in_srcregions_GV(SAVEDIR, GV_flts, OBSvar_str, WINDS, prs_bnds, src_strs, max_days, int_choice):    

    OBSvar_traj_EA_raw, CC_dep_EA, int_lon, int_lat, pct_EA = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[0], max_days, '', '')
    OBSvar_traj_SA_raw, CC_dep_SA, int_lon, int_lat, pct_SA = get_OBSvar_info(GV_flts, SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, int_choice, src_strs[1], max_days, '', '')   
    
    CC_dep_EA = np.asarray(CC_dep_EA)
    CC_dep_SA = np.asarray(CC_dep_SA)  


    
    if OBSvar_str == 'Toluene':
        OBSvar_traj_EA_raw = replace_detection_limit(OBSvar_traj_EA_raw, -888., 0.05)
        OBSvar_traj_SA_raw = replace_detection_limit(OBSvar_traj_SA_raw, -888., 0.05)
        
    OBSvar_traj_EA = remove_missing_pts(OBSvar_traj_EA_raw, CC_dep_EA)   # TT removed
    OBSvar_traj_SA = remove_missing_pts(OBSvar_traj_SA_raw, CC_dep_SA)   # TT removed
    
    #init_prs_EA = remove_missing_pts(init_prs_EA, CC_dep_EA)
    #init_prs_SA = remove_missing_pts(init_prs_SA, CC_dep_SA)

    #print(np.min(OBSvar_traj_EA), np.max(OBSvar_traj_EA))

    return OBSvar_traj_EA, OBSvar_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA
        
if __name__ == '__main__':

    AIRPLANE = 'GV'
    WINDS = 'ERA5_kin'
    SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
    #GV_FLT_DATES = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
    GV_FLT_DATES = ['20220730','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV TOGA dates
    WB_FLT_DATES = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']  # WB57 dates
    #FLT_DATES = ['20220819']
    prs_bnds = [0,500]  # The pressure surface below which we will not count trajectories
    max_days = 30
    src_strs = ['[95, 135, 30, 50]','[65, 95, 18, 35]']   # [95, 130, 30, 45]  '[105, 135, 40, 50]'  '[-1, 360, -90, 90]',
    
    #OBSvar_str = 'CO'
    #OBSvar_lab = 'CO'
    #OBSvar_unit = 'ppbv'
    #hist_rng = (0,320)
    #bins = 50
    #rng = (0,500)
    
    #OBSvar_str = 'SO2'
    #OBSvar_lab = 'SO$_2$'
    #OBSvar_unit = 'pptv'    
    #hist_rng = (0,250)
    #bins=100
    #rng = (0,500)
    
    OBSvar_str = 'CH2CL2'  # This is the variable that will be "matched" to the trajectories
    OBSvar_lab = 'CH$_2$Cl$_2$'
    OBSvar_unit = 'pptv'
    hist_rng = (0,700)
    bins=35
    rng = (0,700)
    
    #OBSvar_str = 'Toluene'  # This is the variable that will be "matched" to the trajectories
    #OBSvar_lab = 'Toluene'
    #OBSvar_unit = 'pptv'
    #hist_rng = (0,100)
    #bins=50
    #rng = (0,100)

    if len(GV_FLT_DATES) > 1: PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/' + AIRPLANE + '_AllRFs_' + WINDS + '/'
    else: PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/' + AIRPLANE + '_' + FLT_DATES[0] + '_' + WINDS + '/'

    # For different source regions, make distributions (these plots are used for 2023 AGU)


        
    ############################################
    #  Below are some plots for the full campaign
    ############################################
        
    if OBSvar_str == 'Toluene': 
        OBSvar_traj_EA, OBSvar_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = \
            get_OBSvar_in_srcregions_GV(SAVEDIR, GV_FLT_DATES, OBSvar_str, WINDS, prs_bnds, src_strs, max_days, 'finalcld')           
    else:
        OBSvar_traj_EA, OBSvar_traj_SA, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = \
            get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, OBSvar_str, WINDS, prs_bnds, src_strs, max_days, 'finalcld')    
        
    OBSdist_EA, binz = np.histogram(OBSvar_traj_EA, bins = bins, range = rng, density = True)
    OBSdist_EA = OBSdist_EA*pct_EA*(rng[1]-rng[0])/bins  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA, histtype = 'step', label='East Asia, ' + str(round(np.mean(OBSvar_traj_EA[OBSvar_traj_EA>-900]),1)) + ' ' + OBSvar_unit, linewidth=3)

    OBSdist_SA, binz = np.histogram(OBSvar_traj_SA, bins = bins, range = rng, density = True)
    OBSdist_SA = OBSdist_SA*pct_SA*(rng[1]-rng[0])/bins  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage 
    hist_SA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA, histtype = 'step', label='South Asia, ' + str(round(np.mean(OBSvar_traj_SA[OBSvar_traj_SA>-900]),1)) + ' ' + OBSvar_unit, linewidth=3)
    
    plt.legend(loc='upper right', fontsize=10)
    plt.title(OBSvar_lab + ' Distributions for 2022 ACCLIP Sampling', fontsize=16)
    plt.xlabel(OBSvar_lab + ' (' + OBSvar_unit + ')', fontsize=14)
    #plt.xlabel('Toluene/Benzene Ratio', fontsize=12)
    plt.ylabel("Fraction of aircraft sampling (%)", fontsize=14)
    #plt.ylabel('Relative Frequency', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    if OBSvar_str == 'SO2': 
        plt.yscale('log')
        plt.ylim(0.01,10)
        
    if OBSvar_str == 'Toluene':
        plt.yscale('log')
        plt.ylim(0.01,30)
    

        
    plt.tight_layout()
    plt.savefig('/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/' + OBSvar_str + '_dists.png', dpi=300)
    plt.close('all')




    ###############  Do a tracer relationship plot for the two regions during August 6

    MeForm_traj_EA_A6, MeForm_traj_SA_A6, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = get_OBSvar_in_srcregions_GV(SAVEDIR, ['20220806','20220807'], 'MeForm', WINDS, prs_bnds, ['[105, 135, 40, 50]','[65, 95, 18, 35]'], max_days, 'finalcld')
    #MeForm_A6 = np.asarray(list(MeForm_traj_EA_A6) + list(MeForm_traj_SA_A6))
    CO_traj_EA_A6,  CO_traj_SA_A6,  CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = get_OBSvar_in_srcregions_GV(SAVEDIR, ['20220806','20220807'], 'CO_A', WINDS, prs_bnds, ['[105, 135, 40, 50]','[65, 95, 18, 35]'], max_days, 'finalcld')
    #CO_A6 = np.asarray(list(CO_traj_EA_A6) + list(CO_traj_SA_A6))

    MeForm_traj_EA_all, MeForm_traj_SA_all, CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = get_OBSvar_in_srcregions_GV(SAVEDIR, GV_FLT_DATES, 'MeForm', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
    MeForm_all = np.asarray(list(MeForm_traj_EA_all) + list(MeForm_traj_SA_all))
    CO_traj_EA_all,  CO_traj_SA_all,  CC_dep_EA, CC_dep_SA, pct_EA, pct_SA = get_OBSvar_in_srcregions_GV(SAVEDIR, GV_FLT_DATES, 'CO_A', WINDS, prs_bnds, src_strs, max_days, 'finalcld')
    CO_all = np.asarray(list(CO_traj_EA_all) + list(CO_traj_SA_all))
    
    plt.scatter(CO_all, MeForm_all, s=2, marker='o', color='k', label='All ACCLIP sampling')
    #plt.scatter(CO_A6, MeForm_A6, s=3, marker='o', color='red', label='August 6, 2022 only')
    plt.scatter(CO_traj_EA_A6, MeForm_traj_EA_A6, s=3, marker='o', color='red', label='Northeast Asia, Aug 6-7 2022 only')
    
    
    plt.xlim(0,350)
    plt.ylim(0,300)
    
    plt.title('Methyl Formate vs CO for 2022 ACCLIP Sampling')
    plt.xlabel('CO (ppbv)')
    plt.ylabel('Methyl Formate (pptv)')
    
    plt.legend(loc='upper left')
    
    plt.savefig(PLOTDIR + 'MeForm_vs_CO.png',dpi=300)
    plt.close('all')
    



        
        
        

        
    ############################################
    #  Below are a lot of plots related to the Aug6/7 case study
    ############################################                   
        
    # Read August 6
    OBSvar_traj_EA_A6, OBSvar_traj_SA_A6, CC_dep_EA_A6, CC_dep_SA_A6, init_prs_EA_A6, init_prs_SA_A6, pct_EA_A6, pct_SA_A6 = get_OBSvar_in_srcregions(SAVEDIR, ['20220806'], ['20220806'], OBSvar_str, WINDS, prs_bnds, src_strs, max_days, 'finalcld')
    OBSvar_traj_EA_A6, OBSvar_traj_SA_A6, CC_dep_EA_A6, CC_dep_SA_A6, init_prs_EA_A6, init_prs_SA_A6, pct_EA_A6, pct_SA_A6 = get_OBSvar_in_srcregions(SAVEDIR, ['20220806'], ['20220806'], OBSvar_str, WINDS, prs_bnds, src_strs, max_days, 'finalcld')
    #OBSvar_traj_EA_A6, OBSvar_traj_SA_A6, CC_dep_EA_A6, CC_dep_SA_A6, init_prs_EA_A6, init_prs_SA_A6, pct_EA_A6, pct_SA_A6 = get_OBSvar_in_srcregions(SAVEDIR, GV_FLT_DATES, WB_FLT_DATES, OBSvar_str, WINDS, prs_bnds, src_strs, max_days)

    # Read August 7 
    OBSvar_traj_EA_A7, CC_dep_EA_A7, pct_EA_A7, init_prs_EA_A7 = get_OBSvar_info(['20220807'], SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, 'CC', src_strs[0], max_days)
    OBSvar_traj_SA_A7, CC_dep_SA_A7, pct_SA_A7, init_prs_SA_A7 = get_OBSvar_info(['20220807'], SAVEDIR, OBSvar_str, 'GV', WINDS, prs_bnds, 'CC', src_strs[1], max_days) 

    OBSvar_traj_EA_A7 = np.asarray(remove_missing_pts(OBSvar_traj_EA_A7, CC_dep_EA_A7))
    OBSvar_traj_SA_A7 = np.asarray(remove_missing_pts(OBSvar_traj_SA_A7, CC_dep_SA_A7))
    
    OBSdist_EA_A6, binz = np.histogram(OBSvar_traj_EA_A6, bins = 100, range = (0,500), density = True)
    OBSdist_EA_A6 = OBSdist_EA_A6*pct_EA_A6*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio 
    hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A6, histtype = 'step', label='East Asia, ' + str(round(np.mean(OBSvar_traj_EA_A6[OBSvar_traj_EA_A6>-900]),1)) + ' ' + OBSvar_unit, linewidth=3)

    OBSdist_EA_A7, binz = np.histogram(OBSvar_traj_EA_A7, bins = 50, range = (0,500), density = True)
    OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio 
    #hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', \
    #    label='Northeast Asia Aug 7, ' + str(round(np.mean(OBSvar_traj_EA_A7[OBSvar_traj_EA_A7>-900]),3)) + ' ' + OBSvar_unit, linewidth=3)
        
    OBSdist_SA_A6, binz = np.histogram(OBSvar_traj_SA_A6, bins = 100, range = (0,500), density = True)
    OBSdist_SA_A6 = OBSdist_SA_A6*pct_SA_A6*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio
    hist_SA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A6, histtype = 'step', label='South Asia, ' + str(round(np.mean(OBSvar_traj_SA_A6[OBSvar_traj_SA_A6>-900]),1)) + ' ' + OBSvar_unit, linewidth=3)

    OBSdist_SA_A7, binz = np.histogram(OBSvar_traj_SA_A7, bins = 50, range = (0,500), density = True)
    #OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7*10.  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage.  10 necessary because of CO bin size ratio
    #hist_SA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7, ' + str(round(np.mean(OBSvar_traj_SA_A7[OBSvar_traj_SA_A7>-900]),3)) + ' ' + OBSvar_unit, linewidth=3)
    
    plt.legend(loc='upper right', fontsize=10)
    plt.xlim(0,200)
    plt.ylim(0.1,20)
    plt.title(OBSvar_lab + ' Distributions for 2022 ACCLIP Sampling', fontsize=16)
    plt.xlabel(OBSvar_lab + ' (' + OBSvar_unit + ')', fontsize=14)
    #plt.xlabel('Toluene/Benzene Ratio', fontsize=12)
    plt.ylabel("Fraction of aircraft's trajectories (%)", fontsize=14)
    #plt.ylabel('Relative Frequency', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.yscale('log')

    plt.tight_layout()
    plt.savefig(PLOTDIR + OBSvar_str + '_dists.png', dpi=300)
    plt.close('all')
    
    
    #  Use the same data as above to make a TTD for the two regions too

    OBSdist_EA_A6, binz = np.histogram(CC_dep_EA_A6, bins = 30, range = (0,30), density = True)
    OBSdist_EA_A6 = OBSdist_EA_A6*pct_EA_A6
    hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A6, histtype = 'step', label='Northeast Asia Aug 6', linewidth=3)

    OBSdist_EA_A7, binz = np.histogram(CC_dep_EA_A7, bins = 30, range = (0,30), density = True)
    OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7
    #hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', label='Northeast Asia Aug 7', linewidth=3)
    
    OBSdist_SA_A6, binz = np.histogram(CC_dep_SA_A6, bins = 30, range = (0,30), density = True)
    OBSdist_SA_A6 = OBSdist_SA_A6*pct_SA_A6
    hist_SA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A6, histtype = 'step', label='South Asia Aug 6', linewidth=3)

    OBSdist_SA_A7, binz = np.histogram(CC_dep_SA_A7, bins = 30, range = (0,30), density = True)
    OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7
    #hist_SA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7', linewidth=3)
    
    plt.legend(fontsize=10)
    plt.xlim(0,30)
    plt.title('Transit Time Distributions for 2022 ACCLIP Sampling', fontsize=16)
    #plt.title('Transit Time Distributions for Aug 6 2022', fontsize=16)
    plt.xlabel('Transit time (days)', fontsize=14)
    plt.ylabel("Fraction of aircraft's trajectories (%)", fontsize=14)
    #plt.ylabel('Relative Frequency', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.tight_layout()
    plt.savefig(PLOTDIR + 'TTDs.png', dpi=300)
    plt.close('all')
    
    
    # Use this same data to separate the aircraft sampling altitudes
    
    init_alt_EA_A6 = sh.prs2alt(np.asarray(init_prs_EA_A6), 7.6, 1013.)
    init_alt_SA_A6 = sh.prs2alt(np.asarray(init_prs_SA_A6), 7.6, 1013.)
    init_alt_EA_A7 = sh.prs2alt(np.asarray(init_prs_EA_A7), 7.6, 1013.)
    init_alt_SA_A7 = sh.prs2alt(np.asarray(init_prs_SA_A7), 7.6, 1013.)
    
    OBSdist_EA_A6, binz = np.histogram(init_alt_EA_A6, bins = 20, range = (0,20), density = True)
    OBSdist_EA_A6 = OBSdist_EA_A6*pct_EA_A6
    hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A6, histtype = 'step', label='Northeast Asia Aug 6', linewidth=3, orientation='horizontal')

    OBSdist_EA_A7, binz = np.histogram(init_alt_EA_A7, bins = 20, range = (0,20), density = True)
    OBSdist_EA_A7 = OBSdist_EA_A7*pct_EA_A7
    #hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'red', weights = OBSdist_EA_A7, histtype = 'step', linestyle = 'dashed', label='Northeast Asia Aug 7', linewidth=3, orientation='horizontal')
    
    OBSdist_SA_A6, binz = np.histogram(init_alt_SA_A6, bins = 20, range = (0,20), density = True)
    OBSdist_SA_A6 = OBSdist_SA_A6*pct_SA_A6
    hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A6, histtype = 'step', label='South Asia Aug 6', linewidth=3, orientation='horizontal')

    OBSdist_SA_A7, binz = np.histogram(init_alt_SA_A7, bins = 20, range = (0,20), density = True)
    OBSdist_SA_A7 = OBSdist_SA_A7*pct_SA_A7
    #hist_EA = plt.hist(binz[:-1], binz, edgecolor = 'orange', weights = OBSdist_SA_A7, histtype = 'step', linestyle = 'dashed', label='South Asia Aug 7', linewidth=3, orientation='horizontal')

    plt.legend(fontsize=10)
    plt.ylim(5,21)
    plt.title('Aircraft Measurement Altitudes for 2022-08-06', fontsize=16)
    plt.ylabel('Aircraft Altitude (km)', fontsize=14)
    #plt.xlabel('Toluene/Benzene Ratio', fontsize=12)
    plt.xlabel("Fraction of aircraft's trajectories (%)", fontsize=14)
    #plt.ylabel('Relative Frequency', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.tight_layout()
    plt.savefig(PLOTDIR + 'init_alt.png', dpi=300)
    plt.close('all')
    
    

    
    

    #for src_str in src_strs:
    #
    #    for int_type in ['CC','BL']:
    #    
    #        ###### Make a pixel plot aggregating all the data
    #
    #        OBSvar_traj, dep, pct = get_OBSvar_info(FLT_DATES, SAVEDIR, OBSvar_str, AIRPLANE, WINDS, prs_bnds, int_type, src_str, max_days)
    #
    #        pixels = np.histogram2d(np.asarray(OBSvar_traj), np.asarray(dep), bins = [40,60], range = [[0,400],[0,60]])
    #
    #        pixx = pixels[1]
    #        pixy = pixels[2]
    #        pixdata = np.transpose(pixels[0])
    #        pixdata = pixdata*pct/len(dep)  # Normalize the fraction of all trajectories (need pct, not 100, to account for all the traj which don't reach an influence)
    #        #pixdata = pixdata/np.max(pixdata)  # Normalize the numbers relative to all encounters    
    #
    #        pixdata[pixdata == 0.0] = pixdata[pixdata == 0.0]-0.01
    #
    #        cmap = plt.get_cmap('YlOrRd')
    #        cmap.set_under([0.8,0.8,0.8])
    #
    #        levels = MaxNLocator(nbins=5).tick_values(0.0, 0.5)
    #        norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    #
    #        fig = plt.figure(figsize = (10,8))
    #        ax = fig.add_subplot(1,1,1)
    #        pixplt = ax.pcolormesh(pixx, pixy, pixdata, cmap = cmap, norm = norm)
    #
    #        cbar = plt.colorbar(pixplt,fraction=0.045, pad=0.04, extend = 'max')
    #        #cbar.ax.set_ylabel('Normalized Frequency of Intercept', fontsize = 16)
    #        cbar.ax.set_ylabel("Fraction of aircraft's trajectories (%)", fontsize = 16)
    #        cbar.ax.tick_params(labelsize = 16)
    #
    #        #plt.scatter(OBSvar_traj[::100], CC_dep[::100], s=0.1)
    #
    #        ax.set_xlim(0,320)
    #        ax.set_ylim(0,30)
    #        #ax.set_title('ACCLIP ' + AIRPLANE + ' RFs ' + src_str, fontsize=16)
    #        ax.set_xlabel(OBSvar_str + OBSvar_unit, fontsize=20)
    #        if int_type == 'CC': ax.set_ylabel('Transit time to convective influence (days)', fontsize=20)        
    #        if int_type == 'BL': ax.set_ylabel('Transit time to top of boundary layer (days)', fontsize=20)        
    #        ax.tick_params(labelsize = 20)
    #
    #        plt.savefig(PLOTDIR + OBSvar_str + '_vs_' + int_type + 'dep_' + src_str + '.png')
    #        plt.close('all')


        ####### Compute TTDs for the different bins of CO

        #OBS_binss = [[0,50],[50,100],[100,150],[150,200],[200,250],[250,300]]
        #cols = [[0.9,0.9,0.9],[0.6,0.6,0.6],[0.3,0.3,0.3],[0.0,0.0,0.0],'cyan','blue','red']

        #for OBS in range(len(OBS_binss)):

        #    OBS_bins = OBS_binss[OBS]
        #    #OBSvar_traj_bin = []
        #    CC_dep_bin = []
        #    for pt in range(len(OBSvar_traj)):
        #        curr_OBSvar_traj = OBSvar_traj[pt]
        #        if curr_OBSvar_traj >= OBS_bins[0] and curr_OBSvar_traj <= OBS_bins[1] and CC_dep[pt] > -900:
        #            #OBSvar_traj_bin.append(curr_OBSvar_traj)
        #            CC_dep_bin.append(CC_dep[pt])
            
        #    plt.hist(CC_dep_bin, bins=60, range=[0,60], density=True, histtype='step', label='CO ' + str(OBS_bins) + ' ppbv, Mean = ' + str(round(np.mean(CC_dep_bin),1)) + 'd', color = cols[OBS])

        #plt.xlim(0,40)    
        #plt.xlabel('Transit time to convective encounter (days)', fontsize=12)
        #plt.ylabel('Relative frequency', fontsize=12)
        #plt.legend()
        #plt.savefig(PLOTDIR + 'CCdep_CObins_' + src_str + '.png')
        #plt.close('all')




