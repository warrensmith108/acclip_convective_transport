
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
from Fig2 import get_processed_data, moving_window_smooth, draw_box
import matplotlib

def get_ACCLIP_merge_data(OBSROOT, FLT_DATES, var_str):
    
    OBSvar = []

    for FLT_DATE in FLT_DATES:
    
        AIRFILE = glob.glob(OBSROOT + '*' + FLT_DATE + '*.ict')[0]     
        AIRheaders = np.loadtxt(AIRFILE, delimiter=',', usecols = (0), dtype='str')
        AIRdata = np.loadtxt(AIRFILE, delimiter=',', skiprows=int(AIRheaders[0]))  # Read the file
        OBSvar.extend(AIRdata[:,int(list(AIRheaders).index(var_str))-11])
        
    return np.asarray(OBSvar)
   
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


fig = plt.figure(figsize=(10,5))  

ax1 = fig.add_subplot(2,2,1,projection=ccrs.PlateCarree())
ax2 = fig.add_subplot(2,2,2)

   
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'  
PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'


# Grab the current day's streamfunction
SF_FILES = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas*')    
SFdata = nc4.Dataset(SF_FILES[0])
    
SF_FILE = glob.glob('/UTLS/Model_Output/ACCLIP_2022/CFSV2_SF/2022/cdas*20220806*')[0]    
SFdata = nc4.Dataset(SF_FILE)
SF_prs = SFdata.variables['level0'][:]
SF_lon = SFdata.variables['lon'][:]
SF_lat = SFdata.variables['lat'][:]
SF_prs_index = list(SF_prs).index(150.)
curr_SF = np.squeeze(np.mean(SFdata.variables['STRM_L100'][:,SF_prs_index,:,:], axis=0))

# Grab the flight tracks and intercept information
init_lon_GV, init_lat_GV, init_prs, init_jd, int_dep, int_lon_GV, int_lat_GV, int_pct \
    = get_processed_data(PROCDIR, ['20220806','20220807'], 'GV', 'ERA5_kin', [0,500], [-1,360,-90,90], 30, 'finalcld', '', '', unit='hPa')
init_lon_WB, init_lat_WB, init_prs, init_jd, int_dep, int_lon_WB, int_lat_WB, int_pct \
    = get_processed_data(PROCDIR, ['20220806'], 'WB57', 'ERA5_kin', [0,500], [-1,360,-90,90], 30, 'finalcld', '', '', unit='hPa')

int_lon = int_lon_GV + int_lon_WB
int_lat = int_lat_GV + int_lat_WB



#################  PANEL A - MAP

pixels = np.histogram2d(int_lon, int_lat, bins = [360,180], range = [[0,360],[-90,90]])
pixlon = pixels[1]
pixlat = pixels[2]
pixdata = np.transpose(pixels[0])
#pixdata = pixdata/np.sum(pixdata)  # Normalize the numbers relative to all encounters
pixdata = pixdata*100./len(int_lon)  # Normalize the fraction of all trajectories (pct replaced with 100, because processed files now have missing points included)

pixdata[pixdata == 0.0] = pixdata[pixdata == 0.0]-0.01  # Set the zeros to below zero so they show up gray

#cmap = plt.get_cmap('Reds')
cmap = ListedColormap(np.asarray([[0.8,0.8,0.8],[0.6,0.6,0.6],[1,1,0],[100/256,60/256,3/256],[250/256,150/256,0]]))
#cmap.set_under([0.8,0.8,0.8])
cmap.set_under([1,1,1])
cmap.set_over([1,0,0])

#levels = np.arange(0.02,0.21,0.02)
levels = [0.01,0.04,0.07,0.10,0.20]  # Note changed from 0.01
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    
pixplt = ax1.pcolormesh(pixlon, pixlat, pixdata, cmap = cmap, norm = norm)

ax1.coastlines()
ax1.add_feature(cf.BORDERS)

ax1.scatter(np.asarray(init_lon_WB), np.asarray(init_lat_WB), s=0.01, color='blue')
ax1.scatter(np.asarray(init_lon_GV), np.asarray(init_lat_GV), s=0.01, color='cyan')

GPH_plot = moving_window_smooth(curr_SF, width = 2)
ax1.contour(SF_lon, SF_lat, GPH_plot, [1.0E7, 2.5E7], colors=['black'], linewidths = 2)

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,  color='gray', alpha=0.5, linestyle='--', xlocs=[30,60,90,120,150,179.999], ylocs=[-10,0,10,20,30,40,50,60])
deg = u"\N{DEGREE SIGN}"
ax1.set_xticks([60,90,120,150])
ax1.set_xticklabels(['60' + deg,'90' + deg,'120' + deg,'150' + deg], fontsize=14)#, fontweight='bold')
ax1.set_yticks([0,20,40,60])
ax1.set_yticklabels(['0' + deg,'20' + deg,'40' + deg,'60' + deg], fontsize=14)#, fontweight='bold')
#ax.set_yticks([-10,0,10,20,30,40,50,60])
#ax.set_yticklabels(['-10' + deg,'0' + deg,'10' + deg,'20' + deg,'30' + deg,'40' + deg,'50' + deg,'60' + deg], fontsize=16, fontweight='bold')

cbaxes = fig.add_axes([0.11,0.12,0.40,0.02])
cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'horizontal')
cbar.ax.tick_params(labelsize=14)
#cbar.ax.set_ylabel('Normalized Frequency of Intercept', fontsize = 16)
#cbar.ax.set_xlabel("Fraction of aircraft's trajectories (%)", fontsize = 16)
cbar.ax.set_xlabel("Fraction of aircraft sampling (%)", fontsize = 14)#, fontweight='bold')

for l in cbar.ax.yaxis.get_ticklabels():
    #l.set_weight("bold")
    l.set_fontsize(14)

ax1.set_xlim(50,160)
ax1.set_ylim(5,55)    
    
draw_box(ax1, 105, 135, 40, 50, (79/255., 38/255., 131/255.), 2)      # Northeast China (Aug 6-7)      # Vikings colors
draw_box(ax1,  65,  95, 18, 35, 'orange', 2)   # S Asia
    

########################  PANEL B - DISTRIBUTION    


    
# Define paths and such
GV_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_60SEC_Merge/R1/icartt/'      # Instrument/AERODYNE/R1.1/'
WB_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/WB57/NCAR_60SEC_Merge/R1/icartt/'

GV_FLT_DATES = ['20220807','20220806']   # These are reversed on purpose so that the more polluted Aug 6 dots show up on top
WB_FLT_DATES = ['20220806']

WB_lon = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_LONG_MMS')
WB_lat = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_LAT_MMS')
WB_alt = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_ALT_MMS')*1.e-3
WB_prs = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'P_MMS')
WB_th  = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'POT_MMS')
WB_CO  = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'CO_COLD2_ppbv')

GV_lon = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGLON')
GV_lat = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGLAT')
GV_alt = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGALT')*1.e-3
GV_prs = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'PSXC')
GV_th  = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'THETA')
GV_CO_P  = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'CO_PICARRO')
GV_CO_A  = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'CO_AERODYNE')

GV_CO = []  # We want to merge the CO field together and use Picarro where Aerodyne is missing
for pt in range(len(GV_CO_A)):
    if GV_CO_A[pt] > -900: GV_CO.append(GV_CO_A[pt])
    else: GV_CO.append(GV_CO_P[pt])

# Read in the tropopause information for the vertical sections
with open('/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/trop_prs_file.pkl', 'rb') as f:
    ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_jd, ERA5_trop_prs = pkl.load(f) 
    
trop_jd_ind0  = list(ERA5_trop_jd).index(2459797.5)  # Aug 6 0Z
trop_jd_ind1  = list(ERA5_trop_jd).index(2459798.0)  # Aug 6 12Z

trop_lat_ind0 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-40.))
trop_lat_ind1 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-45.))

mean_ERA5_trop = np.mean(ERA5_trop_prs[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))
std_ERA5_trop  = np.std(ERA5_trop_prs[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))

mean_ERA5_trop = smooth(mean_ERA5_trop,5)

#print(mean_ERA5_trop)
#print(np.shape(mean_ERA5_trop))

#fig = plt.figure()
#ax  = fig.add_subplot(111)
ax2b = ax2.twinx()

# Set up a colorbar
cmap = plt.get_cmap('YlOrRd')
cmap.set_over('black')
cmap.set_under('white')
#levels = [70,100,150,200,300,400]
levels = np.arange(20,181,40)
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    
#ax.fill_between(ERA5_trop_lon, mean_ERA5_trop-std_ERA5_trop, mean_ERA5_trop+std_ERA5_trop, color=[0.7,0.7,0.7])
ax2.plot(ERA5_trop_lon, mean_ERA5_trop, color='gray', linestyle='-', linewidth=2)
ax2b.plot(np.zeros(10), sh.prs2alt(np.arange(100,1001,100), 7.6, 1013.))

ax2.scatter(WB_lon, WB_prs, s=20, c=WB_CO, cmap=cmap, norm=norm, label='NASA WB-57', linewidths=0.1, zorder=990)
sc = ax2.scatter(GV_lon, GV_prs, s=20, c=GV_CO, cmap=cmap, norm=norm, label='NSF NCAR GV', linewidths=0.1, zorder=999)   # edgecolors='k',

#sc = ax.scatter(WB_lon, WB_prs, s=2, c='blue')
#ax.scatter(GV_lon, GV_prs, s=2, c='cyan')

cax = fig.add_axes([0.65,0.12,0.28,0.02])
cbar = plt.colorbar(sc, orientation='horizontal', extend='both', cax=cax)
cbar.ax.set_xlabel('CO (ppbv)', fontsize = 14)#, fontweight='bold')

#plt.ylim(1000,40)
#ax.set_yscale('log')
plt.xticks(fontsize=14)
#plt.yticks([10,12,14,16,18], fontsize=12)

ax2.set_xlabel('Longitude (\N{DEGREE SIGN}E)', fontsize=14)
ax2b.set_ylabel('Altitude (km)', fontsize=14)
ax2.set_ylabel('Pressure (hPa)', fontsize=14)

#ax.set_xlim(127,154)
#ax.set_ylim(9,19.5)
ax2.set_xlim(125,153)

#ax2.set_yticks([0,5,10,15,20])

ax2.set_yscale('log')
ax2.set_yticks([300,200,100,60])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.get_yaxis().get_major_formatter().labelOnlyBase = False

ax2.set_ylim(300,60)
ax2b.set_ylim(sh.prs2alt(300,7.6,1013.),sh.prs2alt(60,7.6,1013.))
ax2b.set_yticks([10,12,14,16,18,20])

ax2.tick_params(axis='both', which='major', labelsize=14)
ax2b.tick_params(axis='both', which='major', labelsize=14)



######################  FINALIZE

ax1.set_position([0.06,0.25,0.50,0.50])
ax2.set_position([0.65,0.275,0.28,0.445])

plt.figtext(0.07,0.29,'(a)',fontsize=20,fontweight='bold')
plt.figtext(0.87,0.64,'(b)',fontsize=20,fontweight='bold')
                  
plt.savefig(PLOTDIR + 'Fig7.png', dpi=300)
plt.savefig(PLOTDIR + 'Fig7.pdf', dpi=300)
plt.close('all')


