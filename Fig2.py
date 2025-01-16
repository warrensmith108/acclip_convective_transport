
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
    
##########  Intercept Plots

def draw_box(ax, left, right, bot, top, col, thk):
    ax.plot([left, left], [bot, top], color=col, linewidth=thk)
    ax.plot([right, right], [bot, top], color=col, linewidth=thk)
    ax.plot([left, right], [bot, bot], color=col, linewidth=thk)
    ax.plot([left, right], [top, top], color=col, linewidth=thk)

def plot_pix_dist(ax, enc_lon, enc_lat, init_lon, init_lat, FLT_DATE, pct, label, prs_bnds, prslon, prslat, GPH_plot, src_region, src_str, AIRPLANE, WINDS, cld_choice, WB_lon=[-999], WB_lat=[-999], INST=''):

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
        
    pixplt = ax.pcolormesh(pixlon, pixlat, pixdata, cmap = cmap, norm = norm)
    
    ax.coastlines()
    ax.add_feature(cf.BORDERS)
    
    if AIRPLANE == 'GV' or AIRPLANE == 'Both': init_col = 'cyan'
    if AIRPLANE == 'WB57': init_col = 'blue'

    #ax.plot(init_lon, init_lat, color=init_col, linewidth=0.2)
    #ax.plot(WB_lon, WB_lat, color='blue', linewidth=0.2)

    ax.scatter(np.asarray(init_lon[::10]), np.asarray(init_lat[::10]), s=0.01, color=init_col)                   # Subset by 20 to save file space
    if AIRPLANE == 'Both': ax.scatter(np.asarray(WB_lon[::10]), np.asarray(WB_lat[::10]), s=0.01, color='blue')
    
    GPH_plot = moving_window_smooth(GPH_plot, width = 2)
    ax.contour(prslon, prslat, GPH_plot, [1.0E7, 2.5E7], colors=['black'], linewidths = 2)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,  color='gray', alpha=0.5, linestyle='--', xlocs=[30,60,90,120,150,179.999], ylocs=[-10,0,10,20,30,40,50,60])
    deg = u"\N{DEGREE SIGN}"
    ax.set_xticks([60,90,120,150])
    ax.set_xticklabels(['60' + deg,'90' + deg,'120' + deg,'150' + deg], fontsize=16)#, fontweight='bold')
    ax.set_yticks([0,20,40,60])
    ax.set_yticklabels(['0' + deg,'20' + deg,'40' + deg,'60' + deg], fontsize=16)#, fontweight='bold')
     
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
                
    #print(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + 'Intercept_Map_' + label + '_' + str(prs_bnds) + '_' + str(src_str) + '_finalcld' + INST + '.png')                
    #plt.savefig(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/' + 'Intercept_Map_' + label + '_' + str(prs_bnds) + '_' + str(src_str) + '_finalcld' + INST + '.png', dpi=300)
    #plt.close('all')
    
    return pixplt



SAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/'
PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'  
OBSSAVEDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/OBSmatch/'  


##### A map plot

def plot_map_axes(ax, WINDS):

    prs_bnds = [0,500]  # The pressure surfaces between which we will count trajectories
    #if AIRPLANE == 'GV':   flts = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
    #if AIRPLANE == 'GV':   flts = ['20220806_30m','20220807_30m']
    #if AIRPLANE == 'WB57': flts = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']
    #WINDS = 'GFS_kin'
    max_days = 30
    #int_choice = 'finalcld'
    #do_med = 0    # Do medians?
    #inst_str = '' # either '_INST' or '', where INST is something like TOGA, AWAS, etc

    init_lon_GV, init_lat_GV, init_prs_GV, init_jd, CC_dep_GV, CC_lon_GV, CC_lat_GV, CC_pct_GV \
        = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '', unit='hPa')
    init_lon_WB, init_lat_WB, init_prs_WB, init_jd, CC_dep_WB, CC_lon_WB, CC_lat_WB, CC_pct_WB \
        = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, [-1,360,-90,90], max_days, 'finalcld', '', '', unit='hPa')
       
    SF_FILES = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas1.202208*')    
    SFdata = nc4.Dataset(SF_FILES[0])
    SF_lon = SFdata.variables['lon'][:]
    SF_lat = SFdata.variables['lat'][:]    
    SF_prs = SFdata.variables['level0'][:]
    SF_prs_index = list(SF_prs).index(150.)    
    SF_2022 = np.squeeze(SFdata.variables['STRM_L100'][:,SF_prs_index,:,:])

    for SF_FILE in SF_FILES[1:]:
        curr_SFdata = nc4.Dataset(SF_FILE)
        curr_SF = np.squeeze(curr_SFdata.variables['STRM_L100'][:,SF_prs_index,:,:])
        SF_2022 = np.append(SF_2022, curr_SF, axis=0)
    #SF_2022 = np.flip(SF_2022, axis=1)  # This is necessary because the lat fields are flipped in the GFS / CFS versions used here        
    SF_mean = np.squeeze(np.mean(SF_2022, axis=0))   

    # Plot the planes together
    CC_lon = np.asarray(CC_lon_GV + CC_lon_WB)
    CC_lat = np.asarray(CC_lat_GV + CC_lat_WB)
    CC_pct = len(CC_lon[CC_lon > -900])*100./(len(init_lon_GV) + len(init_lat_GV))
    
    # Read these in to fill in gaps for the flight tracks
    init_lon_GV_all, init_lat_GV_all, aaa, bbb, ccc, ddd, eee, fff = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, [0,500], [-1,360,-90,90], max_days, 'finalcld', '', '')
    init_lon_WB_all, init_lat_WB_all, aaa, bbb, ccc, ddd, eee, fff = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, [0,500], [-1,360,-90,90], max_days, 'finalcld', '', '')    
                    
    pixplt = plot_pix_dist(ax, CC_lon, CC_lat, init_lon_GV_all, init_lat_GV_all, 'AllRFs', CC_pct, 'finalcld', prs_bnds, SF_lon, SF_lat, SF_mean, 'Global', 'Global', 'Both', WINDS, 'finalcld', WB_lon = init_lon_WB_all, WB_lat = init_lat_WB_all)

    return pixplt
    
    
def plot_dist_axes(ax, WINDS):
    
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

    init_lon_GV, init_lat, init_prs, init_jd, CC_dep_EA_GV, CC_lon, CC_lat, CC_pct_EA_GV \
        = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, EA_box, max_days, 'finalcld', '', '')
    init_lon_WB, init_lat, init_prs, init_jd, CC_dep_EA_WB, CC_lon, CC_lat, CC_pct_EA_WB \
        = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, EA_box, max_days, 'finalcld', '', '')

    CC_dep_EA = CC_dep_EA_GV + CC_dep_EA_WB
    CC_pct_EA = len(np.where(np.asarray(CC_dep_EA) > -900)[0])*100/(len(init_lon_GV)+len(init_lon_WB))

    init_lon_GV, init_lat, init_prs, init_jd, CC_dep_SA_GV, CC_lon, CC_lat, CC_pct_SA_GV \
        = get_processed_data(PROCDIR, ['AllRFs'], 'GV', WINDS, prs_bnds, SA_box, max_days, 'finalcld', '', '')
    init_lon_WB, init_lat, init_prs, init_jd, CC_dep_SA_WB, CC_lon, CC_lat, CC_pct_SA_WB \
        = get_processed_data(PROCDIR, ['AllRFs'], 'WB57', WINDS, prs_bnds, SA_box, max_days, 'finalcld', '', '')

    CC_dep_SA = CC_dep_SA_GV + CC_dep_SA_WB
    CC_pct_SA = len(np.where(np.asarray(CC_dep_SA) > -900)[0])*100/(len(init_lon_GV)+len(init_lon_WB))

    CC_hist_EA, binz = np.histogram(CC_dep_EA, bins = 30, range = (0,30), density = True)
    CC_hist_EA = CC_hist_EA*CC_pct_EA  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_GV = ax.hist(binz[:-1], binz, edgecolor = 'red', weights = CC_hist_EA, histtype = 'step', label='East Asia: ' + str(round(CC_pct_EA,1)) + '%', linewidth=3)

    CC_hist_SA, binz = np.histogram(CC_dep_SA, bins = 30, range = (0,30), density = True)
    CC_hist_SA = CC_hist_SA*CC_pct_SA  # Normalize so that the integral is the encounter percentage, don't use 100 to maintain percentage
    hist_GV = ax.hist(binz[:-1], binz, edgecolor = 'orange', weights = CC_hist_SA, histtype = 'step', label='South Asia: ' + str(round(CC_pct_SA,1)) + '%', linewidth=3)

    ax.set_xlabel('Transit time (days)', fontsize = 16)
    ax.set_ylabel('Fraction of aircraft sampling (%)', fontsize = 16)
    ax.set_xlim(0,30)
    ax.legend(fontsize=16)
    #ax.figtext(0.6,0.73,'(Global: ' + str(round(CC_pct_g,1)) + '%)', fontsize=16)

    ax.tick_params(axis='both', which='major', labelsize=16)
    
    return CC_pct_g
    

#################  Commence Plotting

if __name__ == '__main__':
    
    fig = plt.figure(figsize=(16,10))  

    ax1 = fig.add_subplot(2,2,1,projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3,projection=ccrs.PlateCarree())
    ax4 = fig.add_subplot(2,2,4)

    pixplt = plot_map_axes(ax1, 'GFS_kin')
    GFS_dep = plot_dist_axes(ax2, 'GFS_kin')
    pixplt = plot_map_axes(ax3, 'ERA5_kin')
    ERA5_dep = plot_dist_axes(ax4, 'ERA5_kin')

    print(GFS_dep)
    print(ERA5_dep)



    ax1.set_position([0.09,0.46,0.50,0.50])
    ax3.set_position([0.09,0.00,0.50,0.50])

    ax2.set_position([0.715,0.53,0.27,0.36])
    ax4.set_position([0.715,0.07,0.27,0.36])

    cbaxes = fig.add_axes([0.60,0.24,0.015,0.50])
    cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'vertical')
    #cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'vertical')
    cbar.ax.tick_params(labelsize=16)
    #cbar.ax.set_ylabel('Normalized Frequency of Intercept', fontsize = 16)
    #cbar.ax.set_xlabel("Fraction of aircraft's trajectories (%)", fontsize = 16)
    cbar.ax.set_ylabel("Fraction of aircraft sampling (%)", fontsize = 16 )#, fontweight='bold')

    #for l in cbar.ax.yaxis.get_ticklabels():
    #    #l.set_weight("bold")
    #    l.set_fontsize(16)

    plt.figtext(0.84,0.77,'(Global: ' + str(round(GFS_dep,1)) + '%)',fontsize=16)
    plt.figtext(0.84,0.31,'(Global: ' + str(round(ERA5_dep,1)) + '%)',fontsize=16)
        
    plt.figtext(0.015,0.175,'ERA5',fontsize=48,fontweight='bold',rotation=90)
    plt.figtext(0.015,0.655,'GFS',fontsize=48,fontweight='bold',rotation=90)

    plt.figtext(0.53,0.55,'(a)',fontsize=36,fontweight='bold')
    plt.figtext(0.53,0.09,'(c)',fontsize=36,fontweight='bold')

    plt.figtext(0.92,0.65,'(b)',fontsize=36,fontweight='bold')
    plt.figtext(0.92,0.15,'(d)',fontsize=36,fontweight='bold')

    plt.savefig(PLOTDIR + 'Fig2.png', dpi=300)
    plt.savefig(PLOTDIR + 'Fig2.pdf')
    plt.close('all')







