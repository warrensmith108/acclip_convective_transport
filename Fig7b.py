
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import glob
import netCDF4 as nc4
import jdcal
from ACCLIP_Plotter import get_ACCLIP_merge_data
import scale_hgt as sh

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
# Define paths and such
GV_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_60SEC_Merge/R1/icartt/'      # Instrument/AERODYNE/R1.1/'
WB_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/WB57/NCAR_60SEC_Merge/R1/icartt/'
PLOTDIR = '/home/wsmith/ACCLIP/plots/'

#GV_FLT_DATES = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
#WB_FLT_DATES = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']  # WB57 Osan Flights 

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

fig = plt.figure()
ax  = fig.add_subplot(111)
ax2 = ax.twinx()

# Set up a colorbar
cmap = plt.get_cmap('YlOrRd')
cmap.set_over('black')
cmap.set_under('white')
#levels = [70,100,150,200,300,400]
levels = np.arange(20,181,40)
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    
#ax.fill_between(ERA5_trop_lon, mean_ERA5_trop-std_ERA5_trop, mean_ERA5_trop+std_ERA5_trop, color=[0.7,0.7,0.7])
ax.plot(ERA5_trop_lon, mean_ERA5_trop, color='gray', linestyle='-', linewidth=2)
ax2.plot(np.zeros(10), sh.prs2alt(np.arange(100,1001,100), 7.6, 1013.))

ax.scatter(WB_lon, WB_prs, s=20, c=WB_CO, cmap=cmap, norm=norm, label='NASA WB-57', linewidths=0.1, zorder=990)
sc = ax.scatter(GV_lon, GV_prs, s=20, c=GV_CO, cmap=cmap, norm=norm, label='NSF NCAR GV', linewidths=0.1, zorder=999)   # edgecolors='k',

#sc = ax.scatter(WB_lon, WB_prs, s=2, c='blue')
#ax.scatter(GV_lon, GV_prs, s=2, c='cyan')

cax = fig.add_axes([0.2,0.1,0.6,0.03])
cbar = plt.colorbar(sc, orientation='horizontal', extend='both', label='CO (ppbv)', cax=cax)

#plt.ylim(1000,40)
#ax.set_yscale('log')
plt.xticks(fontsize=12)
#plt.yticks([10,12,14,16,18], fontsize=12)

ax.set_xlabel('Longitude (E)', fontsize=12)
ax2.set_ylabel('Altitude (km)', fontsize=12)
ax.set_ylabel('Pressure (hPa)', fontsize=12)

#ax.set_xlim(127,154)
#ax.set_ylim(9,19.5)
ax.set_xlim(125,153)

#ax2.set_yticks([0,5,10,15,20])

ax.set_yscale('log')
ax.set_yticks([300,200,100,60])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().get_major_formatter().labelOnlyBase = False

ax.set_ylim(300,60)
ax2.set_ylim(sh.prs2alt(300,7.6,1013.),sh.prs2alt(60,7.6,1013.))

ax.set_position([0.2,0.25,0.6,0.6])

#plt.legend(loc='lower right')

#plt.savefig(PLOTDIR + 'Vertical_FltTracks.png', dpi=300)
#plt.tight_layout()
plt.savefig(PLOTDIR + 'Figure7b.png', dpi=300)
plt.close('all')


