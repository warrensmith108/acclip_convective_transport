
####################################################################
#
#  This script is intended to read "processed" output and plot chemical tracer relationships
#
####################################################################

from OBSmatch import get_OBSvar_info
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os

def remove_missing_pts_2vars(OBSvar1_raw, OBSvar2_raw, CC_dep, int_lon, int_lat):
    OBSvar1 = []
    OBSvar2 = []
        
    for pt in range(len(CC_dep)):
        if CC_dep[pt] > -900 and OBSvar1_raw[pt] > -900 and OBSvar2_raw[pt] > -900:           
            OBSvar1.append(OBSvar1_raw[pt])        
            OBSvar2.append(OBSvar2_raw[pt])        
    return np.asarray(OBSvar1), np.asarray(OBSvar2)


def read_OBSvars(SAVEDIR, flt_choice, OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, src_str, max_days, int_choice):    
      
    OBSvar1_traj_raw, CC_dep1, int_lon1, int_lat1, pct = get_OBSvar_info(flt_choice, SAVEDIR, OBSvar1_str, 'GV', WINDS, prs_bnds, int_choice, src_str, max_days, '', '')
    OBSvar2_traj_raw, CC_dep2, int_lon2, int_lat2, pct = get_OBSvar_info(flt_choice, SAVEDIR, OBSvar2_str, 'GV', WINDS, prs_bnds, int_choice, src_str, max_days, '', '') 
    
    print(len(OBSvar1_traj_raw))
    print(len(OBSvar2_traj_raw))
    print(len(CC_dep1))
    print(len(CC_dep2))
    
    #CC_dep_EA = np.asarray(CC_dep_EA_GV + CC_dep_EA_WB)
    #CC_dep_SA = np.asarray(CC_dep_SA_GV + CC_dep_SA_WB)
    
    #pct_EA = len(CC_dep_EA[CC_dep_EA > -900])*100./(len(CC_dep_EA_GV) + len(CC_dep_EA_WB))
    #pct_SA = len(CC_dep_SA[CC_dep_SA > -900])*100./(len(CC_dep_SA_GV) + len(CC_dep_SA_WB))      
    
    OBSvar1_traj, OBSvar2_traj = remove_missing_pts_2vars(OBSvar1_traj_raw, OBSvar2_traj_raw, CC_dep1, int_lon1, int_lat1)   
    
    #OBSvar1_traj = remove_missing_pts(OBSvar1_traj_raw, CC_dep1)
    #OBSvar2_traj = remove_missing_pts(OBSvar2_traj_raw, CC_dep2)
    
    #init_prs_EA = remove_missing_pts(init_prs_EA_GV + init_prs_EA_WB, CC_dep_EA)
    #init_prs_SA = remove_missing_pts(init_prs_SA_GV + init_prs_SA_WB, CC_dep_SA)

    return OBSvar1_traj, OBSvar2_traj, CC_dep1, int_lon1, int_lat1
    
    
    

AIRPLANE = 'GV'
WINDS = 'ERA5_kin'
int_choice = 'finalcld'
prs_bnds = [0,500]
max_days = 30
SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/chemical_relationships/'

#OBSvar_labs = ['CO_P','SO2','N2O','BC','CH2CL2','CHCl3','MeForm','Toluene','Benzene','Propane','Ethane','CFC11','CFC12','CFC113','HCFC22','HCFC141b','DMS','CH3Cl','CH3Br','ClBenzene','CCl4','iButane','iPentane','MTBE','CH3CCl3','Ethanol']
OBSvar_strs = ['CO','Toluene','Ethanol','MeForm','Benzene','ClBenzene','iButane','iPentane','nButane','nPentane','MTBE']
OBSvar_labs = ['CO','Toluene','Ethanol','Methyl Formate','Benzene','Chlorobenzene','iButane','iPentane','nButane','nPentane','MTBE']
units       = ['ppbv','pptv','pptv', 'pptv', 'pptv', 'pptv', 'pptv', 'pptv', 'pptv', 'pptv', 'pptv'] 
 
#for v1 in range(len(OBSvar_strs)):
for v1 in [0]:  # To just do CO on the x-axis for now 
 
    OBSvar1_str = OBSvar_strs[v1]
    OBSvar1_lab = OBSvar_labs[v1]
    unit1 = units[v1]

    for v2 in range(len(OBSvar_strs)):
        
        OBSvar2_str = OBSvar_strs[v2]
        OBSvar2_lab = OBSvar_labs[v2]
        unit2 = units[v2]
        
        if v2 > v1: 
                               
            print(OBSvar1_str, OBSvar2_str)
            
            OBSvar1_glob, OBSvar2_glob, dep, int_lon, int_lat = read_OBSvars(SAVEDIR, ['AllRFs'], OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[-1, 360, -90, 90]', max_days, int_choice)
            OBSvar1_NEC,  OBSvar2_NEC, dep, int_lon, int_lat  = read_OBSvars(SAVEDIR, ['20220806','20220807'], OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[105, 135, 40, 50]', max_days, int_choice)
            OBSvar1_A19,  OBSvar2_A19, dep_A19, int_lon_A19, int_lat_A19  = read_OBSvars(SAVEDIR, ['20220819'], OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[95, 135, 30, 50]', max_days, int_choice)
            OBSvar1_A15,  OBSvar2_A15, dep_A15, int_lon_A15, int_lat_A15  = read_OBSvars(SAVEDIR, ['20220815'], OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[95, 135, 30, 50]', max_days, int_choice)
            
            plt.scatter(OBSvar1_glob, OBSvar2_glob, s=2, marker='o', color='gray', label='All ACCLIP sampling')            
            plt.scatter(OBSvar1_A15, OBSvar2_A15, s=5, marker='o', color='blue', label='E Asia, August 15 only')
            plt.scatter(OBSvar1_A19, OBSvar2_A19, s=5, marker='o', color='black', label='E Asia, August 19 only')
            plt.scatter(OBSvar1_NEC, OBSvar2_NEC, s=5, marker='o', color='red', label='NE China, Aug 6-7 only')

            #OBSvar1_EA,   OBSvar2_EA   = read_OBSvars(SAVEDIR, OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[95, 135, 30, 50]', max_days, int_choice)
            #OBSvar1_SA,   OBSvar2_SA   = read_OBSvars(SAVEDIR, OBSvar1_str, OBSvar2_str, WINDS, prs_bnds, '[65, 95, 18, 35]', max_days, int_choice)  
            
            #plt.scatter(OBSvar1_SA, OBSvar2_SA, s=3, marker='o', color='orange', label='South Asia')
            #plt.scatter(OBSvar1_EA, OBSvar2_EA, s=3, marker='o', color='red', label='East Asia')
            
            plt.xlim(0,np.max(OBSvar1_glob)*1.05)
            plt.ylim(0,np.max(OBSvar2_glob)*1.05)
            
            plt.title(OBSvar2_lab + ' vs ' + OBSvar1_lab)
            
            plt.xlabel(OBSvar1_lab + ' (' + unit1 + ')')
            plt.ylabel(OBSvar2_lab + ' (' + unit2 + ')')
            
            plt.legend()
            
            plt.savefig(PLOTDIR + OBSvar2_str + 'vs' + OBSvar1_str + '.png', dpi=300)
            plt.savefig(PLOTDIR + OBSvar2_str + 'vs' + OBSvar1_str + '.pdf', dpi=300)
            plt.close('all')
            
            
            if OBSvar1_str == 'Benzene' and OBSvar2_str == 'Ethanol':
                
                A15_int_lon = []
                A15_int_lat = []
                A15_dep = []
                
                for pt in range(len(OBSvar1_A15)):
                    
                    if OBSvar1_A15[pt] <= 140 and OBSvar2_A15[pt] > 1900:   # This is an attempt to capture the upper branch of this chemical relationship

                        A15_int_lon.append(int_lon_A15[pt])
                        A15_int_lat.append(int_lat_A15[pt])
                        A15_dep.append(dep_A15[pt])
                        
                plt.hist(A15_int_lon, range=[100,150], bins=50)
                plt.savefig(PLOTDIR + OBSvar2_str + 'vs' + OBSvar1_str + 'Aug15_int_lons.png')
                plt.close('all')

                plt.hist(A15_int_lat, range=[30,50], bins=30)
                plt.savefig(PLOTDIR + OBSvar2_str + 'vs' + OBSvar1_str + 'Aug15_int_lats.png')
                plt.close('all')

                plt.hist(A15_dep, range=[0,5], bins=20)
                plt.savefig(PLOTDIR + OBSvar2_str + 'vs' + OBSvar1_str + 'Aug15_dep.png')
                plt.close('all')                

            

                                                