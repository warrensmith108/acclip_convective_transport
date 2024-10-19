
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

def get_processed_data(PROCDIR, FLT_DATES, AIRPLANE, WINDS, prs_bnds, src_region, max_days, int_choice, inst_str, med_str, unit='hPa'):

    init_lon = []
    init_lat = []
    init_prs = []
    init_jd = []
    int_dep = []
    int_lon = []
    int_lat = []
    int_pct = []

    for FLT_DATE in FLT_DATES:
        
        PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + inst_str + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd_' + int_choice + '' + med_str + '_processed.pkl'

        #if len(INST) > 1: 
        #    if do_median: PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd' + INST + '_' + int_type + '_medians_processed.pkl'     
        #    else:  PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd' + INST + '_' + int_type + '_processed.pkl'      
        #elif do_median: PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd_' + int_type + '_medians_processed.pkl'  
        #else: PROCFILE = PROCDIR + FLT_DATE + '_' + AIRPLANE + '_' + WINDS + '_' + str(prs_bnds) + unit + '_' + str(src_region) + '_' + str(max_days) + 'd_' + int_type + '_processed.pkl'      
        
        with open(PROCFILE, 'rb') as f:
            curr_init_lon, curr_init_lat, curr_init_prs, curr_init_jd, curr_int_dep, curr_int_lon, curr_int_lat, curr_int_pct = pkl.load(f)
            
        init_lon = init_lon + curr_init_lon
        init_lat = init_lat + curr_init_lat
        init_prs = init_prs + curr_init_prs
        init_jd = init_jd + curr_init_jd
        int_dep = int_dep + curr_int_dep
        int_lon = int_lon + curr_int_lon
        int_lat = int_lat + curr_int_lat
        curr_int_pct = [curr_int_pct]        
        int_pct = int_pct + curr_int_pct  # This is sloppy... but we don't have time for rational solutions
    
    return init_lon, init_lat, init_prs, init_jd, int_dep, int_lon, int_lat, int_pct
                
def moving_window_smooth(data, width=1):
    data_copy = np.copy(data)
    smoothed = np.zeros(data_copy.shape)
    for ix in range(data_copy.shape[0]):
        for iy in range(data_copy.shape[1]):
            smoothed[ix,iy] = np.nanmean(data_copy[ix-width:ix+width,iy-width:iy+width])
    return smoothed   
    
def plot_dep_dist(dep, pct, label, FLT_DATE, prs_bnds, src_str, col):

    plt.hist(dep, range = (0,30), bins = 30, density = True, color = col)
    plt.title(FLT_DATE)
    plt.xlabel(label + ' Time (d)')
    plt.ylabel('Relative Frequency')

    plt.savefig(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + label + '_dep_hist_' + FLT_DATE + '_' + str(prs_bnds) + '_' + str(src_str) + '.png', dpi=300)
    plt.close('all')
    
    #plt.figtext(0.6,0.70, '30-day pct: ' + str(round(pct, 2)) + '%', color='purple')
    #plt.figtext(0.6,0.75, 'Source: ' + src_str, color='purple')
    #plt.figtext(0.6,0.8, 'Prs Bounds: ' + str(prs_bnds) + ' hPa', color='purple')

def plot_prs_dist(prs, pct, label, FLT_DATE, prs_bnds, src_str, col):

    plt.hist(prs, range = (0,1000), bins = 50, density = True, color = col, orientation='horizontal')
    plt.title(FLT_DATE)
    
    plt.ylim(1000,0)
    plt.ylabel(label + ' Pressure (hPa)')
    plt.xlabel('Relative Frequency')
    
    plt.savefig(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + label + '_prs_hist_' + FLT_DATE + '_' + str(prs_bnds) + '_' + str(src_str) + '.png', dpi=300)
    plt.close('all')
    
##########  Intercept Plots

def draw_box(ax, left, right, bot, top, col, thk):
    ax.plot([left, left], [bot, top], color=col, linewidth=thk)
    ax.plot([right, right], [bot, top], color=col, linewidth=thk)
    ax.plot([left, right], [bot, bot], color=col, linewidth=thk)
    ax.plot([left, right], [top, top], color=col, linewidth=thk)

def plot_pix_dist(enc_lon, enc_lat, init_lon, init_lat, FLT_DATE, pct, label, prs_bnds, prslon, prslat, PBLT, GPH_plot, src_region, src_str, AIRPLANE, WINDS, cld_choice, WB_lon=[-999], WB_lat=[-999], INST=''):

    BLpixels = np.histogram2d(enc_lon, enc_lat, bins = [360,180], range = [[0,360],[-90,90]])
    pixlon = BLpixels[1]
    pixlat = BLpixels[2]
    pixdata = np.transpose(BLpixels[0])
    #pixdata = pixdata/np.sum(pixdata)  # Normalize the numbers relative to all encounters
    pixdata = pixdata*100./len(enc_lon)  # Normalize the fraction of all trajectories (pct replaced with 100, because processed files now have missing points included)
   
    pixdata[pixdata == 0.0] = pixdata[pixdata == 0.0]-0.01  # Set the zeros to below zero so they show up gray

    #cmap = plt.get_cmap('Reds')
    cmap = ListedColormap(np.asarray([[0.8,0.8,0.8],[0.6,0.6,0.6],[1,1,0],[100/256,60/256,3/256],[250/256,150/256,0]]))
    #cmap.set_under([0.8,0.8,0.8])
    cmap.set_under([1,1,1])
    cmap.set_over([1,0,0])

    #levels = np.arange(0.02,0.21,0.02)
    levels = [0.01,0.04,0.07,0.10,0.20]  # Note changed from 0.01
    norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
   
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
        
    pixplt = ax.pcolormesh(pixlon, pixlat, pixdata, cmap = cmap, norm = norm)
    
    ax.coastlines()
    ax.add_feature(cf.BORDERS)
    
    if AIRPLANE == 'GV' or AIRPLANE == 'Both': init_col = 'cyan'
    if AIRPLANE == 'WB57': init_col = 'blue'

    plt.scatter(np.asarray(init_lon), np.asarray(init_lat), s=0.01, color=init_col)
    
    if AIRPLANE == 'Both': plt.scatter(np.asarray(WB_lon), np.asarray(WB_lat), s=0.01, color='blue')
    
    #if FLT_DATE == '20220730' or FLT_DATE == '20220830': ax.set_title(AIRPLANE + ' Aircraft: ' + FLT_DATE[0:4] + '-' + FLT_DATE[4:6] + '-31', fontsize=16)
    #else: ax.set_title(AIRPLANE + ' Aircraft: ' + FLT_DATE[0:4] + '-' + FLT_DATE[4:6] + '-' + FLT_DATE[6:], fontsize=16)
    
    #plt.figtext(0.17,0.27,'30-day pct: ' + str(round(pct, 2)) + '%', fontsize=10)
    #plt.figtext(0.20,0.33,'Prs Bounds: ' + str(prs_bnds) + 'hPa')     
    #plt.figtext(0.20,0.27,'Max Days: ' + str(max_days) + 'd')    
    
    #plt.figtext(0.35,0.26,str(round(pct, 2)) + '% influence within ' + str(max_days) + 'd', fontsize=24, backgroundcolor='white', fontweight='bold')
    #if AIRPLANE == 'GV': plt.figtext(0.12,0.79,'NSF/NCAR GV', fontsize=24, backgroundcolor='white', fontweight='bold', color=init_col)
    #if AIRPLANE == 'WB57': plt.figtext(0.12,0.79,'NASA WB-57', fontsize=24, backgroundcolor='white', fontweight='bold', color=init_col)
     
    # Tibetan Plateau
    #PBLT_s = moving_window_smooth(PBLT, width = 2)
    #ax.contour(prslon, prslat, PBLT_s, [609.], colors=['gray'], linewidths = 1)
    
    GPH_plot = moving_window_smooth(GPH_plot, width = 2)
    ax.contour(prslon, prslat, GPH_plot, [1.0E7, 2.5E7], colors=['black'], linewidths = 2)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,  color='gray', alpha=0.5, linestyle='--', xlocs=[30,60,90,120,150,179.999], ylocs=[-10,0,10,20,30,40,50,60])
    deg = u"\N{DEGREE SIGN}"
    ax.set_xticks([60,90,120,150])
    ax.set_xticklabels(['60' + deg,'90' + deg,'120' + deg,'150' + deg], fontsize=16)#, fontweight='bold')
    ax.set_yticks([0,20,40,60])
    ax.set_yticklabels(['0' + deg,'20' + deg,'40' + deg,'60' + deg], fontsize=16)#, fontweight='bold')
    #ax.set_yticks([-10,0,10,20,30,40,50,60])
    #ax.set_yticklabels(['-10' + deg,'0' + deg,'10' + deg,'20' + deg,'30' + deg,'40' + deg,'50' + deg,'60' + deg], fontsize=16, fontweight='bold')
    
    cbaxes = fig.add_axes([0.82,0.31,0.015,0.39])
    cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'vertical')
    cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'vertical')
    cbar.ax.tick_params(labelsize=14)
    #cbar.ax.set_ylabel('Normalized Frequency of Intercept', fontsize = 16)
    #cbar.ax.set_xlabel("Fraction of aircraft's trajectories (%)", fontsize = 16)
    cbar.ax.set_ylabel("Fraction of aircraft sampling (%)", fontsize = 14)#, fontweight='bold')
    
    for l in cbar.ax.yaxis.get_ticklabels():
        #l.set_weight("bold")
        l.set_fontsize(16)
    
    ax.set_position([0.08,0.15,0.72,0.72])
    
    ax.set_xlim(50,160)
    ax.set_ylim(5,55)    

    #draw_box(ax, src_region[0], src_region[1], src_region[2], src_region[3], 'red', 1.5)
    #if src_str == 'Global':
    #draw_box(ax,  30, 150, 15, 45, 'black', 1.5)     
    #draw_box(ax, 105, 135, 40, 50, 'red', 3)      # Northeast China (Aug 6-7)   
    draw_box(ax, 95,  135, 30, 50, 'red', 2)      # EASF
    draw_box(ax,  65,  95, 18, 35, 'orange', 2)   # S Asia
    
    try: os.mkdir(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/')
    except: zzzz = -999
                
    print(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + 'Intercept_Map_' + label + '_' + str(prs_bnds) + '_' + str(src_str) + '_finalcld' + INST + '.png')                
    plt.savefig(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + 'Intercept_Map_' + label + '_' + str(prs_bnds) + '_' + str(src_str) + '_finalcld' + INST + '.png', dpi=300)
    plt.close('all')
                

if __name__ == '__main__':

    SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/'
    PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'
    PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'  
    OBSSAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'    
    #src_regions = [[-1,360,-90,90],[95,130,30,45],[65,95,18,35]]


    ######  A TTD plot for the relative contributions of E and S Asia

    #init_lon, init_lat, init_prs, BL_dep_EA, BL_lon_EA, BL_lat_EA, BL_pct_EA, CC_dep_EA, CC_lon_EA, CC_lat_EA, CC_pct_EA = get_processed_data(PROCDIR, ['AllRFs'], 'GV', 'ERA5', [0,250], [95,130,30,45], max_days, cld_choice)
    #init_lon, init_lat, init_prs, BL_dep_SA, BL_lon_SA, BL_lat_SA, BL_pct_SA, CC_dep_SA, CC_lon_SA, CC_lat_SA, CC_pct_SA = get_processed_data(PROCDIR, ['AllRFs'], 'GV', 'ERA5', [0,250], [65,95,18,35], max_days, cld_choice)
    #
    #plt.hist(CC_dep_SA, range = (0,30), bins = 30, density = False, color = 'orange', histtype='step', linewidth=3, label='South Asia Monsoon')
    #plt.hist(CC_dep_EA, range = (0,30), bins = 30, density = False, color = 'red', histtype='step', linewidth=3, label='East Asia Monsoon')
    #
    #plt.title('Transit Times for GV Flights')
    #plt.xlabel('Transit time (days)')
    #plt.ylabel('Number of Trajectories')
    #plt.xlim(0,10)
    #plt.legend()
    #
    #plt.savefig(PLOTDIR + 'disttest.png')
    #plt.close('all')


    ##### A map plot

    prs_bnds = [0,500]  # The pressure surfaces between which we will count trajectories
    AIRPLANE = 'WB57'
    if AIRPLANE == 'GV':   flts = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
    #if AIRPLANE == 'GV':   flts = ['20220806_30m','20220807_30m']
    if AIRPLANE == 'WB57': flts = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']
    WINDS = 'GFS_kin'
    max_days = 30
    int_choice = 'finalcld'
    do_med = 0    # Do medians?
    inst_str = '' # either '_INST' or '', where INST is something like TOGA, AWAS, etc

    print('Reading in analysis information...')
    
    PRSFILE = '/home/wsmith/traject/src/ACCLIP_season_input/GFS/TRAJ_WINDS_GFS_2022043000_2022090500_0.0.nc'
    prs_data = nc4.Dataset(PRSFILE)
    PBLT = np.mean(prs_data.variables['PSFC'][:]*0.87, axis=0)
    prslon = prs_data.variables['Longitude'][:]
    prslat = prs_data.variables['Latitude'][:]
    
    SF_FILES = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas1.202208*')    
    SFdata = nc4.Dataset(SF_FILES[0])
    SF_prs = SFdata.variables['level0'][:]
    SF_prs_index = list(SF_prs).index(150.)    
    SF_2022 = np.squeeze(SFdata.variables['STRM_L100'][:,SF_prs_index,:,:])
    
    for SF_FILE in SF_FILES[1:]:
        curr_SFdata = nc4.Dataset(SF_FILE)
        curr_SF = np.squeeze(curr_SFdata.variables['STRM_L100'][:,SF_prs_index,:,:])
        SF_2022 = np.append(SF_2022, curr_SF, axis=0)
    SF_2022 = np.flip(SF_2022, axis=1)  # This is necessary because the lat fields are flipped in the GFS / CFS versions used here        
    SF_mean = np.squeeze(np.mean(SF_2022, axis=0))

        
    ####### Old version: get the GFS GPH at 100 hPa
    #GPH = prs_data.variables['GPH'][::sub,:,:,:]
    #GPHprs = prs_data.variables['Altitude'][:]
    #GPHprs_index = list(GPHprs).index(100.)
    #GPHjul = prs_data.variables['Julian_day'][::sub]
    #GPHsec = prs_data.variables['Seconds'][::sub]
    #GPHjd = GPHjul+(GPHsec/86400.)-0.5        
    #GPH_tind0 = list(GPHjd).index(np.sum(jdcal.gcal2jd(2022,8,1)))                
    #GPH_tind1 = list(GPHjd).index(np.sum(jdcal.gcal2jd(2022,8,31)))                
    #GPH_plot = np.squeeze(np.mean(GPH[GPH_tind0:GPH_tind1,GPHprs_index,:,:], axis=0))
    #######
         
    # Pixel map plots for individual flights
    #for fltx in flts: 
    #    
    #    flt = fltx[0:8]
    #    print(flt)
    #
    #    # Grab the current day's streamfunction
    #    SF_FILE = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas*' + flt + '*')[0]    
    #    SFdata = nc4.Dataset(SF_FILE)
    #    curr_SF = np.flip(np.squeeze(np.mean(SFdata.variables['STRM_L100'][:,SF_prs_index,:,:], axis=0)), axis=0)
    #
    #    init_lon, init_lat, init_prs, init_jd, BL_dep, BL_lon, BL_lat, BL_pct, CC_dep, CC_lon, CC_lat, CC_prs, CC_temp, CC_pct \
    #        = get_processed_data(PROCDIR, [fltx], AIRPLANE, WINDS, prs_bnds, [-1,360,-90,90], max_days, cld_choice, INST=inst_str)
    #    plot_pix_dist(CC_lon, CC_lat, init_lon, init_lat, fltx, CC_pct, 'CC', prs_bnds, prslon, prslat, PBLT, curr_SF, 'Global', 'Global', AIRPLANE, WINDS, cld_choice, INST=inst_str) 
    ##    plot_pix_dist(BL_lon, BL_lat, init_lon, init_lat, flt, BL_pct, 'BL', prs_bnds, prslon, prslat, PBLT, curr_SF, 'Global', 'Global', AIRPLANE, WINDS, cld_choice) 
    ##    plot_dep_dist(CC_dep, CC_pct, 'CC', flt, prs_bnds, 'Global', 'k')
    ##    plot_prs_dist(CC_prs, CC_pct, 'CC', flt, prs_bnds, 'Global', 'k')

    # For all the flights at once
    #init_lon, init_lat, init_prs, BL_dep, BL_lon, BL_lat, BL_pct, CC_dep, CC_lon, CC_lat, CC_pct = get_processed_data(PROCDIR, ['AllRFs'], AIRPLANE, WINDS, prs_bnds, [-1,360,-90,90], max_days, cld_choice)
    #plot_pix_dist(CC_lon, CC_lat, init_lon, init_lat, 'AllRFs', CC_pct, 'CC', prs_bnds, prslon, prslat, PBLT, GPH_plot, 'Global', 'Global', AIRPLANE, WINDS, cld_choice)
    #plot_pix_dist(BL_lon, BL_lat, init_lon, init_lat, 'AllRFs', BL_pct, 'BL', prs_bnds, prslon, prslat, PBLT, GPH_plot, 'Global', 'Global', AIRPLANE, WINDS, cld_choice)




    # For all flights from both planes at the same time
    
    RF_choice = 'AllRFs'    
    if RF_choice != 'AllRFs':
        SF_FILE = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas*' + RF_choice + '*')[0]    
        SFdata = nc4.Dataset(SF_FILE)
        SF_mean = np.flip(np.squeeze(np.mean(SFdata.variables['STRM_L100'][:,SF_prs_index,:,:], axis=0)), axis=0)    
    
    init_lon_GV, init_lat_GV, init_prs_GV, init_jd, CC_dep_GV, CC_lon_GV, CC_lat_GV, CC_pct_GV \
        = get_processed_data(PROCDIR, [RF_choice], 'GV', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '', unit='hPa')
    init_lon_WB, init_lat_WB, init_prs_WB, init_jd, CC_dep_WB, CC_lon_WB, CC_lat_WB, CC_pct_WB \
        = get_processed_data(PROCDIR, [RF_choice], 'WB57', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '', unit='hPa')
        
    # Read these in to fill in gaps for the flight tracks
    init_lon_GV_all, init_lat_GV_all, aaa, bbb, ccc, ddd, eee, fff = get_processed_data(PROCDIR, [RF_choice], 'GV', WINDS, [0,500], [-1,360,-90,90], max_days, 'finalcld', '', '')
    init_lon_WB_all, init_lat_WB_all, aaa, bbb, ccc, ddd, eee, fff = get_processed_data(PROCDIR, [RF_choice], 'WB57', WINDS, [0,500], [-1,360,-90,90], max_days, 'finalcld', '', '')
    
    # Plot the planes separately
    plot_pix_dist(CC_lon_GV, CC_lat_GV, init_lon_GV_all, init_lat_GV_all, RF_choice, CC_pct_GV, 'CC', prs_bnds, prslon, prslat, PBLT, SF_mean, 'Global', 'Global', 'GV', WINDS, 'finalcld')
    plot_pix_dist(CC_lon_WB, CC_lat_WB, init_lon_WB_all, init_lat_WB_all, RF_choice, CC_pct_WB, 'CC', prs_bnds, prslon, prslat, PBLT, SF_mean, 'Global', 'Global', 'WB57', WINDS, 'finalcld')
    
    # Plot the planes together
    CC_lon = np.asarray(CC_lon_GV + CC_lon_WB)
    CC_lat = np.asarray(CC_lat_GV + CC_lat_WB)
    CC_pct = len(CC_lon[CC_lon > -900])*100./(len(init_lon_GV) + len(init_lat_GV))
                    
    plot_pix_dist(CC_lon, CC_lat, init_lon_GV_all, init_lat_GV_all, RF_choice, CC_pct, int_choice, prs_bnds, prslon, prslat, PBLT, SF_mean, 'Global', 'Global', 'Both', WINDS, int_choice, WB_lon = init_lon_WB_all, WB_lat = init_lat_WB_all)

    # Plot the planes in vertical altitude space, colored by convective influence time
    
    #with open(OBSSAVEDIR + 'CO_OBSmatch_CC_GV_ERA5_AllRFs_' + str(prs_bnds) + '_[-1, 360, -90, 90]_30d.pkl', 'rb') as f:
    #    OBSvar_traj_CC_GV, CC_dep, CC_lon, CC_lat, CC_pct = pkl.load(f)
    #with open(OBSSAVEDIR + 'CO_OBSmatch_CC_WB57_ERA5_AllRFs_' + str(prs_bnds) + '_[-1, 360, -90, 90]_30d.pkl', 'rb') as f:
    #    OBSvar_traj_CC_WB, CC_dep, CC_lon, CC_lat, CC_pct = pkl.load(f)     
    #
    #CC_dep = np.asarray(CC_dep_GV + CC_dep_WB)
    #OBSvar = np.asarray(OBSvar_traj_CC_GV + OBSvar_traj_CC_WB)
    #init_lon = np.asarray(init_lon_GV + init_lon_WB)
    #init_prs = np.asarray(init_prs_GV + init_prs_WB)
    #init_alt = sh.prs2alt(init_prs, 8., 1013.)
    
    def vert_track_clr(init_lon, init_alt, OBSvar, OBSvar_str, OBSvar_unit, levels, cbar_rev = 0):
    
        if cbar_rev: 
            cmap = plt.get_cmap('YlOrRd_r')
            cmap.set_over('white')
            cmap.set_under('black')
            ext = 'min'
        else: 
            cmap = plt.get_cmap('YlOrRd')
            cmap.set_over('black')
            cmap.set_under('white')
            ext = 'both'
        norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
        
        plt.scatter(init_lon[OBSvar > -900], init_alt[OBSvar > -900], s=1, c=OBSvar[OBSvar > -900], cmap=cmap, norm=norm)
        plt.xlabel('Longitude (E)')
        plt.ylabel('Altitude (km)')
        plt.xlim(120,160)
        
        cbar = plt.colorbar(extend=ext)
        cbar.set_label(OBSvar_str + ' ' + OBSvar_unit)
        
        plt.savefig(PLOTDIR + 'Vertical_Tracks_' + OBSvar_str + '.png', dpi=300)
        plt.close('all')
        
    #vert_track_clr(init_lon, init_alt, OBSvar, 'CO', '(ppbv)', np.arange(20,171,30))    
    #vert_track_clr(init_lon, init_alt, CC_dep, 'Convective_Influence_Time', '(days)', [1,3,7,14,21,30], cbar_rev=1) 

    ##### A distribution plot separated by airplane 

    init_lon, init_lat, init_prs, init_jd, CC_dep_GV, CC_lon, CC_lat, CC_pct_GV \
        = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')
    init_lon, init_lat, init_prs, init_jd, CC_dep_WB, CC_lon, CC_lat, CC_pct_WB \
        = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')

    CC_hist_GV, binz = np.histogram(CC_dep_GV, bins = 30, range = (0,30), density = True)
    CC_hist_GV = CC_hist_GV*CC_pct_GV  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_GV = plt.hist(binz[:-1], binz, edgecolor = 'cyan', weights = CC_hist_GV, histtype = 'step', label='NSF/NCAR GV, ' + str(round(CC_pct_GV,1)) + '%', linewidth=3)

    CC_hist_WB, binz = np.histogram(CC_dep_WB, bins = 30, range = (0,30), density = True)
    CC_hist_WB = CC_hist_WB*CC_pct_WB  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_WB = plt.hist(binz[:-1], binz, edgecolor = 'blue', weights = CC_hist_WB, histtype = 'step', label='NASA WB-57, ' + str(round(CC_pct_WB,1)) + '%', linewidth=3)
    
    plt.xlabel('Transit time (days)', fontsize = 16)
    plt.ylabel('Fraction of Trajectories (%)', fontsize = 16)
    plt.xlim(0,30)
    plt.legend(fontsize=16)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()
    plt.savefig(PLOTDIR + 'Both_AllRFs_' + WINDS + '/Dist_Both_' + str(prs_bnds) + '_finalcld.png', dpi=300)
    plt.close('all')
    
    ###### Overplotted TTDs for both planes for various pressure bins
    
    def make_both_TTD_hist(prs_bnds, color):
    
        init_lon_GV, init_lat, init_prs, init_jd, CC_dep_GV, CC_lon, CC_lat, CC_pct_GV \
            = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')   
        init_lon_WB, init_lat, init_prs, init_jd, CC_dep_WB, CC_lon, CC_lat, CC_pct_WB \
            = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '')  
                       
        CC_dep = np.asarray(CC_dep_GV + CC_dep_WB)
        CC_pct = len(CC_dep[CC_dep > -900])*100./(len(init_lon_GV)+len(init_lon_WB))
            
        CC_hist, binz = np.histogram(CC_dep, bins = 30, range = (0,30), density = True)
        CC_hist = CC_hist*CC_pct  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
        hist_GV = plt.hist(binz[:-1], binz, edgecolor = color, weights = CC_hist, histtype = 'step', label=str(prs_bnds) + ' hPa, ' + str(round(CC_pct, 1)) + '%', linewidth=3)        
        
    make_both_TTD_hist( [50,100], [0.50,0,0])
    make_both_TTD_hist([100,150], [0.65,0,0])
    make_both_TTD_hist([150,200], [0.80,0,0])
    make_both_TTD_hist([200,250], [0.95,0,0])
    
    plt.xlabel('Transit time (days)', fontsize = 16)
    plt.ylabel('Fraction of Trajectories (%)', fontsize = 16)
    plt.xlim(0,30)
    plt.legend(fontsize=16)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()
    
    plt.legend()
    plt.savefig(PLOTDIR + 'BothPlane_TTD_hists.png', dpi=300)
    plt.close('all')
    


    ###### Overplotted CO dists for both planes for various pressure bins
    
    def make_both_CO_hist(prs_bnds, color):
    
        with open(OBSSAVEDIR + 'CO_OBSmatch_CC_GV_ERA5_AllRFs_' + str(prs_bnds) + '_[-1, 360, -90, 90]_30d.pkl', 'rb') as f:
            OBSvar_traj_CC_GV, CC_dep, CC_lon, CC_lat, CC_pct = pkl.load(f)
        with open(OBSSAVEDIR + 'CO_OBSmatch_CC_WB57_ERA5_AllRFs_' + str(prs_bnds) + '_[-1, 360, -90, 90]_30d.pkl', 'rb') as f:
            OBSvar_traj_CC_WB, CC_dep, CC_lon, CC_lat, CC_pct = pkl.load(f)
            
        OBSvar_traj = np.asarray(OBSvar_traj_CC_GV + OBSvar_traj_CC_WB)
        OBSvar_mean = str(round(np.mean(OBSvar_traj[OBSvar_traj > -900]), 1))
        plt.hist(OBSvar_traj, edgecolor = color, histtype = 'step', range=(0,350), bins=35, label=str(prs_bnds) + ' hPa, ' + OBSvar_mean + ' ppbv', linewidth=3, density=True) 
                                  
    make_both_CO_hist( [50,100], [0.50,0,0])
    make_both_CO_hist([100,150], [0.65,0,0])
    make_both_CO_hist([150,200], [0.80,0,0])
    make_both_CO_hist([200,250], [0.95,0,0])
    
    plt.xlabel('CO (ppbv)', fontsize = 16)
    plt.ylabel('Relative Frequency', fontsize = 16)
    plt.xlim(0,300)
    plt.legend(fontsize=16)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.tight_layout()
    
    plt.legend()
    plt.savefig(PLOTDIR + 'BothPlane_CO_hists.png', dpi=300)
    plt.close('all')    
    



    


