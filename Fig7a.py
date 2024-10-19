
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
from plot_airplane_traj3d import get_processed_data, moving_window_smooth, draw_box

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

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    
pixplt = ax.pcolormesh(pixlon, pixlat, pixdata, cmap = cmap, norm = norm)

ax.coastlines()
ax.add_feature(cf.BORDERS)

plt.scatter(np.asarray(init_lon_WB), np.asarray(init_lat_WB), s=0.01, color='blue')
plt.scatter(np.asarray(init_lon_GV), np.asarray(init_lat_GV), s=0.01, color='cyan')

#plt.figtext(0.17,0.27,'30-day pct: ' + str(round(pct, 2)) + '%', fontsize=10)
#plt.figtext(0.20,0.33,'Prs Bounds: ' + str(prs_bnds) + 'hPa')     
#plt.figtext(0.20,0.27,'Max Days: ' + str(max_days) + 'd')    

#plt.figtext(0.35,0.26,str(round(pct, 2)) + '% influence within ' + str(max_days) + 'd', fontsize=24, backgroundcolor='white', fontweight='bold')
#if AIRPLANE == 'GV': plt.figtext(0.12,0.79,'NSF/NCAR GV', fontsize=24, backgroundcolor='white', fontweight='bold', color=init_col)
#if AIRPLANE == 'WB57': plt.figtext(0.12,0.79,'NASA WB-57', fontsize=24, backgroundcolor='white', fontweight='bold', color=init_col)
 
# Tibetan Plateau
#PBLT_s = moving_window_smooth(PBLT, width = 2)
#ax.contour(prslon, prslat, PBLT_s, [609.], colors=['gray'], linewidths = 1)

GPH_plot = moving_window_smooth(curr_SF, width = 2)
ax.contour(SF_lon, SF_lat, GPH_plot, [1.0E7, 2.5E7], colors=['black'], linewidths = 2)

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,  color='gray', alpha=0.5, linestyle='--', xlocs=[30,60,90,120,150,179.999], ylocs=[-10,0,10,20,30,40,50,60])
deg = u"\N{DEGREE SIGN}"
ax.set_xticks([60,90,120,150])
ax.set_xticklabels(['60' + deg,'90' + deg,'120' + deg,'150' + deg], fontsize=16)#, fontweight='bold')
ax.set_yticks([0,20,40,60])
ax.set_yticklabels(['0' + deg,'20' + deg,'40' + deg,'60' + deg], fontsize=16)#, fontweight='bold')
#ax.set_yticks([-10,0,10,20,30,40,50,60])
#ax.set_yticklabels(['-10' + deg,'0' + deg,'10' + deg,'20' + deg,'30' + deg,'40' + deg,'50' + deg,'60' + deg], fontsize=16, fontweight='bold')

cbaxes = fig.add_axes([0.2,0.18,0.6,0.02])
cbar = plt.colorbar(pixplt, cax=cbaxes, extend = 'both', orientation = 'horizontal')
cbar.ax.tick_params(labelsize=14)
#cbar.ax.set_ylabel('Normalized Frequency of Intercept', fontsize = 16)
#cbar.ax.set_xlabel("Fraction of aircraft's trajectories (%)", fontsize = 16)
cbar.ax.set_xlabel("Fraction of aircraft sampling (%)", fontsize = 14)#, fontweight='bold')

for l in cbar.ax.yaxis.get_ticklabels():
    #l.set_weight("bold")
    l.set_fontsize(16)

#ax.set_position([0.08,0.15,0.72,0.72])

ax.set_xlim(50,160)
ax.set_ylim(5,55)    

#draw_box(ax, src_region[0], src_region[1], src_region[2], src_region[3], 'red', 1.5)
#if src_str == 'Global':
#draw_box(ax,  30, 150, 15, 45, 'black', 1.5)     
draw_box(ax, 105, 135, 40, 50, (79/255., 38/255., 131/255.), 3)      # Northeast China (Aug 6-7)      # Vikings colors
#draw_box(ax, 95,  135, 30, 50, 'red', 2)      # EASF
draw_box(ax,  65,  95, 18, 35, 'orange', 2)   # S Asia

try: os.mkdir(PLOTDIR + AIRPLANE + '_' + FLT_DATE + '_' + WINDS + '/')
except: zzzz = -999
                           
plt.savefig(PLOTDIR + 'Fig7a.png', dpi=300)
plt.close('all')

