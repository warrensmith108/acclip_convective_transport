
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

AIRPLANE = 'GV'
WINDS = 'ERA5_kin'
RAWDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/'
PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'
OBSSAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'
if AIRPLANE == 'GV': OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_60SEC_Merge/R1.0/icartt/'       # Note the 1.0!  These files have the tropopause in them
#if AIRPLANE == 'GV': OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_TOGA_Merge/R1.0/icartt/'      
if AIRPLANE == 'WB57': OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/WB57/NCAR_60SEC_Merge/R2/icartt/' 

if AIRPLANE == 'GV': FLT_DATES = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
#if AIRPLANE == 'GV': FLT_DATES = ['20220730','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  FOR TOGA!!!!!
if AIRPLANE == 'WB57': FLT_DATES = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']  # WB57 Osan Flights
#FLT_DATES = ['20220714','20220716','20220718','20220721_RA_1','20220721_RA_2','20220724','20220725','20220727','20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901','20220909','20220912','20220913','20220914']  # WB57 Flights
#OBS_DATES = ['20220714','20220716','20220718','20220721'     ,'20220721'     ,'20220724','20220725','20220727','20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901','20220909','20220912','20220913','20220914']  # WB57 Flights

prs_bndss = [[0,500]]  # The pressure surfaces between which we will count trajectories  [[100,150],[150,200]]

max_days = 30

int_choice = 'finalcld'
            
src_regions = [[-1,360,-90,90],[65,95,18,35],[105,135,40,50],[95,135,30,50]]#,[125,150,10,30]] 

inst_str = ''  # either '_INST' or '', where INST is something like TOGA, AWAS, etc

do_all = 1  # This determines whether we will also save all the flights together (set to 0 when just processing a single flight)

med_str = ''   # This string will be appended to the save files, either '_medians' or ''
 
do_OBS = 1   # This determines whether the code will stitch observations to the same grid as the trajectories

# List all the variables that will be "matched" to the trajectories, do this en masse now
if do_OBS:
    if AIRPLANE == 'WB57':
        #OBSvar_strs = ['CO_COLD2_ppbv','SO2_LIF','CH2Cl2_WAS']#'N2O','BC_mass_90_550_nm_std','CHCl3_WAS','Propane_WAS','Ethane_WAS',\
        #               #'CFC11_WAS','CFC12_WAS','CFC113_WAS','HCFC22_WAS','HCFC141b_WAS','DMS_WAS','CH3Cl_WAS','CH3Br_WAS','C6H5Cl_WAS','CCl4_WAS','iButane_WAS','iPentane_WAS','G_ALT_MMS']
        #OBSvar_labs = ['CO','SO2','CH2CL2','N2O','BC','CHCl3','Propane','Ethane','CFC11','CFC12','CFC113','HCFC22','HCFC141b','DMS','CH3Cl','CH3Br','ClBenzene','CCl4','iButane','iPentane','z']
        OBSvar_strs = ['POT_MMS']#'P_MMS']
        OBSvar_labs = ['theta']#'prs']
    if AIRPLANE == 'GV':                                           
        #OBSvar_strs = ['CO','SO2_GTCIMS','CH2Cl2_AWAS','Toluene_TOGA','C2H5OH_TOGA']#,'N2O_AERODYNE','BC_mass_90_550_nm_std_SP2','CHCl3_AWAS','MeFormate_TOGA','Benzene_TOGA',\
        #               #'Propane_AWAS','Ethane_AWAS','CFC11_AWAS','CFC12_AWAS','CFC113_AWAS','HCFC22_AWAS','HCFC141b_AWAS','DMS_AWAS','CH3Cl_AWAS','CH3Br_AWAS','C6H5Cl_AWAS',\
        #               #'CCl4_AWAS','iButane_AWAS','iPentane_AWAS','nButane_AWAS','nPentane_AWAS','MTBE_TOGA','CH3CCl3_TOGA','GGALT','PSXC']
        #OBSvar_labs = ['CO','SO2','CH2CL2','Toluene','Ethanol']#,'N2O','BC','CHCl3','MeForm','Benzene','Propane','Ethane','CFC11','CFC12','CFC113','HCFC22','HCFC141b','DMS','CH3Cl','CH3Br',\
        #               #'ClBenzene','CCl4','iButane','iPentane','MTBE','CH3CCl3','z','nButane','nPentane','prs']
        OBSvar_strs = ['THETA']
        OBSvar_labs = ['theta']        
else:
    OBSvar_strs = ['']
    OBSvar_labs = ['']

if int_choice == 'tropth': unit = 'K'
else: unit = 'hPa'
            
for prs_bnds in prs_bndss: 

    for src_region in src_regions:
        
        for obs in range(len(OBSvar_strs)):
            
            OBSvar_str = OBSvar_strs[obs]
            OBSvar_lab = OBSvar_labs[obs]           

            # We may wish to save all the flight information, for the specific source regions
            if do_all:    
                init_lon_all = []
                init_lat_all = []
                init_prs_all = []
                init_jd_all  = []
                #BL_dep_all = []
                #BL_lon_all = []
                #BL_lat_all = []
                int_dep_all = []  # Now generalized
                int_lon_all = []
                int_lat_all = []
                #CC_prs_all = []
                #CC_temp_all = []
                #OBSvar_traj_BL_all = []
                OBSvar_traj_int_all = []
                

            for fff in range(len(FLT_DATES)):

                FLT_DATE = FLT_DATES[fff]
                #OBS_DATE = OBS_DATES[fff]
            
                print('***', FLT_DATE, '***', src_region, '***', int_choice, '***', OBSvar_str)

                RAWFILE = RAWDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + inst_str + '_' + int_choice + '.pkl'                 
                
                init_lon_raw, init_lat_raw, init_prs_raw, init_jd_raw, int_dep_raw, int_lon_raw, int_lat_raw = pkl.load(open(RAWFILE, 'rb'))
              
                # If we are doing medians of each cluster, subset those now        
                if len(med_str) > 2: 
                    init_lon_raw = np.median(init_lon_raw, axis=1)
                    init_lat_raw = np.median(init_lat_raw, axis=1)
                    init_prs_raw = np.median(init_prs_raw, axis=1)
                    init_jd_raw  = np.median(init_jd_raw, axis=1)
                    int_dep_raw   = np.median(int_dep_raw, axis=1)
                    int_lon_raw   = np.median(int_lon_raw, axis=1)
                    int_lat_raw   = np.median(int_lat_raw, axis=1)      
                # Otherwise, just flatten the arrays
                else: 
                    init_lon_raw = init_lon_raw.flatten()
                    init_lat_raw = init_lat_raw.flatten()
                    init_prs_raw = init_prs_raw.flatten()
                    init_jd_raw  = init_jd_raw.flatten()
                    int_dep_raw   = int_dep_raw.flatten()
                    int_lon_raw   = int_lon_raw.flatten()
                    int_lat_raw   = int_lat_raw.flatten()

                
                # Create blank fields to be filled with subsetted values
                init_lon = []
                init_lat = []
                init_prs = [] 
                init_jd = []
                int_dep = []
                int_lon = []
                int_lat = []
                #BL_init_jd = []
                #CC_prs = []
                #CC_temp = []

                for pt in range(len(init_lon_raw)):
                    
                    curr_prs = init_prs_raw[pt]
                    curr_init_lat = init_lat_raw[pt]
                    
                    if curr_prs > prs_bnds[0] and curr_prs < prs_bnds[1] and curr_init_lat > -90:
                    
                        init_lon.append(init_lon_raw[pt])
                        init_lat.append(init_lat_raw[pt])
                        init_prs.append(init_prs_raw[pt])
                        init_jd.append(init_jd_raw[pt])
                        
                        if int_dep_raw[pt] <= max_days and int_lon_raw[pt] >= src_region[0] and int_lon_raw[pt] <= src_region[1] and int_lat_raw[pt] >= src_region[2] and int_lat_raw[pt] <= src_region[3]:
                            int_dep.append(int_dep_raw[pt])
                            int_lon.append(int_lon_raw[pt])
                            int_lat.append(int_lat_raw[pt])
                            #BL_init_jd.append(init_jd_raw[pt])  
                        else:
                            int_dep.append(-99999)
                            int_lon.append(-99999)
                            int_lat.append(-99999)
                            #BL_init_jd.append(-99999)
                            
                                    
                int_good = len(np.where(np.asarray(int_lon) > -900)[0])
                if len(init_lon) > 0: int_pct = int_good*100./len(init_lon)
                else: int_pct = -999           
                
                PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + inst_str + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd_' + int_choice + med_str + '_processed.pkl'   
                with open(PROCFILE, 'wb') as f:
                    pkl.dump((init_lon, init_lat, init_prs, init_jd, int_dep, int_lon, int_lat, int_pct), f)
                    
                    
                if do_OBS:

                    curr_year = int(FLT_DATE[0:4])
                    curr_month = int(FLT_DATE[4:6])
                    curr_day = int(FLT_DATE[6:8])
                    
                    ###########  Compare against ACCLIP Obs
                    
                    AIRFILE = glob.glob(OBSROOT + '*' + FLT_DATE[0:8] + '*.ict')[0]     
                    AIRheaders = np.loadtxt(AIRFILE, delimiter=',', usecols = (0), dtype='str')
                    AIRdata = np.loadtxt(AIRFILE, delimiter=',', skiprows=int(AIRheaders[0]))  # Read the file
                    
                    # Now that we calculate the species in a mass-produced way, we can move the ratio calc to the analysis and plotting script
                    #def ratio_calc(var1_str, var2_str):
                    #    OBSvar1 = np.asarray(AIRdata[:,int(list(AIRheaders).index(var1_str))-11])
                    #    OBSvar2 = np.asarray(AIRdata[:,int(list(AIRheaders).index(var2_str))-11])
                    #    OBSvar = []
                    #    for pt in range(len(OBSvar1)):
                    #        if OBSvar1[pt] >= 0 and OBSvar2[pt] >= 0: OBSvar.append(OBSvar1[pt]/OBSvar2[pt])
                    #        else: OBSvar.append(-999)                    
                    #    return np.asarray(OBSvar)                
                    
                    #if AIRPLANE == 'GV':                         
                    #    OBSsec = np.asarray(AIRdata[:,int(list(AIRheaders).index('Mid_Time_UTC'))-11])
                    #    if OBSvar_str == 'CO': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('CO_AERODYNE'))-11])
                    #    if OBSvar_str == 'SO2': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('SO2_GTCIMS'))-11])
                    #    if OBSvar_str == 'N2O': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('N2O_AERODYNE'))-11])
                    #    if OBSvar_str == 'BC': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('BC_mass_90_550_nm_std_SP2'))-11])
                    #    if OBSvar_str == 'CH2Cl2': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('CH2Cl2_AWAS'))-11])
                    #    if OBSvar_str == 'CHCl3': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('CHCl3_AWAS'))-11])
                    #    if OBSvar_str == 'MeForm': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('MeFormate_TOGA'))-11])
                    #    if OBSvar_str == 'Toluene': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('Toluene_TOGA'))-11])
                    #    if OBSvar_str == 'Tol-Benz': OBSvar = ratio_calc('Toluene_TOGA', 'Benzene_TOGA')
                    #    if OBSvar_str == 'Prop-Eth': OBSvar = ratio_calc('Propane_AWAS', 'Ethane_AWAS')
                    #    if OBSvar_str == 'z': OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index('GGALT'))-11])*1.e-3
                    #    if OBSvar_str == 'troprelz': OBSvar = (np.asarray(AIRdata[:,int(list(AIRheaders).index('GGALT'))-11])*1.e-3) - np.asarray(AIRdata[:,int(list(AIRheaders).index('wmo_1st_z_era5'))-11])                        

                    OBSsec = np.asarray(AIRdata[:,int(list(AIRheaders).index('Mid_Time_UTC'))-11])

                    # We need an ugly bypass for GV CO to use Picarro if Aerodyne is missing
                    if AIRPLANE == 'GV' and OBSvar_str == 'CO':
                        CO_A = np.asarray(AIRdata[:,int(list(AIRheaders).index('CO_AERODYNE'))-11])
                        CO_P = np.asarray(AIRdata[:,int(list(AIRheaders).index('CO_PICARRO'))-11])
                        OBSvar = []
                        for pt in range(len(CO_A)):
                            if CO_A[pt] > -900: OBSvar.append(CO_A[pt])
                            else: OBSvar.append(CO_P[pt])  # If Picarro is also missing, just fill with the missing value
                        
                    else: OBSvar = np.asarray(AIRdata[:,int(list(AIRheaders).index(OBSvar_str))-11])  # For every normal specie

                    #  Convert the AIRsec field to a JD
                    OBS_jd = []
                    for t in range(len(AIRdata)):
                        curr_OBSsec = OBSsec[t]
                        curr_jd = sum(jdcal.gcal2jd(curr_year, curr_month, curr_day)) + curr_OBSsec/86400.
                        OBS_jd.append(curr_jd)
                    OBS_jd = np.asarray(OBS_jd)
                        
                    # At each trajectory point, save the corresponding observation
                    # BL departure
                    OBSvar_traj_int = []
                    for pt in range(len(int_lon)):
                        curr_jd = init_jd[pt]
                        OBSvar_traj_int.append(OBSvar[np.argmin(np.abs(curr_jd-OBS_jd))])

                    # Save the file   (create an error if running with median)
                    OBS_SAVEFILE = OBSSAVEDIR + OBSvar_lab + '_OBSmatch_' + int_choice + '_' + AIRPLANE + '_' + WINDS + inst_str + '_' + FLT_DATE + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd' + med_str + '.pkl'
                    with open(OBS_SAVEFILE, 'wb') as f:
                        pkl.dump((OBSvar_traj_int, int_dep, int_lon, int_lat, int_pct), f)

                if do_all:
                    init_lon_all.extend(init_lon)
                    init_lat_all.extend(init_lat)
                    init_prs_all.extend(init_prs)
                    init_jd_all.extend(init_jd)
                    int_dep_all.extend(int_dep)
                    int_lon_all.extend(int_lon)
                    int_lat_all.extend(int_lat)
                    if do_OBS:
                        OBSvar_traj_int_all.extend(OBSvar_traj_int)
                    
            if do_all: 

                int_good = len(np.where(np.asarray(int_lon_all) > -900)[0])
                if len(init_lon_all) > 0: int_pct_all = int_good*100./len(init_lon_all)          
                else: int_pct_all = -999            
               
                PROCFILE = PROCDIR + 'AllRFs_' + AIRPLANE + '_' + WINDS + inst_str + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd_' + int_choice + '' + med_str + '_processed.pkl'   
                with open(PROCFILE, 'wb') as f:
                    pkl.dump((init_lon_all, init_lat_all, init_prs_all, init_jd_all, int_dep_all, int_lon_all, int_lat_all, int_pct_all), f)

                if do_OBS:

                    OBS_SAVEFILE = OBSSAVEDIR + OBSvar_lab + '_OBSmatch_' + int_choice + '_' + AIRPLANE + '_' + WINDS + inst_str + '_AllRFs_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd' + med_str + '.pkl'
                    with open(OBS_SAVEFILE, 'wb') as f:
                        pkl.dump((OBSvar_traj_int_all, int_dep_all, int_lon_all, int_lat_all, int_pct_all), f)
                        

