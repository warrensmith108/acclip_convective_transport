
import numpy as np
import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import netCDF4 as nc4
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import gridspec as gspec
from scale_hgt import alt2prs, prs2alt
import jdcal
import pickle as pkl
#from plot_airplane_traj3d import get_processed_data
from ACCLIP_Plotter import get_ACCLIP_merge_data
import math
import scale_hgt as sh
import matplotlib


# Define file paths
cld_choice = 'final'
z_choice = 'alt'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'

date_rng = ['2022080100','2022083118']  #######  FOR BEST RESULTS, USE ALL OF AUGUST

GV_FLT_DATES = ['20220730','20220804','20220806','20220807','20220812','20220815','20220816','20220819','20220822','20220823','20220825','20220826','20220829','20220830'] #  GV dates
WB_FLT_DATES = ['20220802','20220804','20220806','20220812','20220813','20220815','20220816','20220819','20220821','20220823','20220825','20220826','20220829','20220831','20220901']  # WB57 Osan Flights 

GV_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_60SEC_Merge/R0/icartt/'      # Instrument/AERODYNE/R1.1/'
WB_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/WB57/NCAR_60SEC_Merge/R1/icartt/'

##################  Retrieve flight tracks

WB_lon = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_LONG_MMS')
WB_lat = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_LAT_MMS')
WB_alt = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'G_ALT_MMS')*1.e-3
WB_prs = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'P_MMS')
WB_th  = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'POT_MMS')
#WB_CO  = get_ACCLIP_merge_data(WB_OBSROOT, WB_FLT_DATES, 'CO_COLD2_ppbv')

GV_lon = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGLON')
GV_lat = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGLAT')
GV_alt = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'GGALT')*1.e-3
GV_prs = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'PSXC')
GV_th  = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'THETA')
#GV_CO  = get_ACCLIP_merge_data(GV_OBSROOT, GV_FLT_DATES, 'CO_PICARRO')

################## Retrieve cloud information

cld_file = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/cloud_data_file_' + cld_choice + '.pkl'
with open(cld_file, 'rb') as f: cloud_jd, cloudlon, cloudlat, cloudalt, cloudtheta, trop_alt, trop_theta = pkl.load(f)

# Convert date range to JD
cloud_jd = np.asarray(cloud_jd)
init_jd = sum(jdcal.gcal2jd(int(date_rng[0][0:4]), int(date_rng[0][4:6]), int(date_rng[0][6:8]))) + int(date_rng[0][8:10])/24.
end_jd  = sum(jdcal.gcal2jd(int(date_rng[1][0:4]), int(date_rng[1][4:6]), int(date_rng[1][6:8]))) + int(date_rng[1][8:10])/24.
tind0 = list(cloud_jd).index(init_jd)
tind1 = list(cloud_jd).index(end_jd)
NTC_sub = tind1 - tind0 + 1

if z_choice == 'alt':
    cloud_data = cloudalt
    cld_minhgt = 14.
    cbar_max   = 20
    z_unit = 'km'
    
cld_subset = cloud_data[tind0:tind1+1,:,:]
cld_frac   = np.zeros(np.shape(np.squeeze(cld_subset[0,:,:])))-999

# Compute the cloud fraction
for ln in range(len(cloudlon)):
    for lt in range(len(cloudlat)):
        #if rainlon[ln] <= plt_lonrng[1] and rainlon[ln] >= plt_lonrng[0] and rainlat[lt] <= plt_latrng[1] and rainlat[lt] >= plt_latrng[0]:        
        cld_frac[lt,ln] = len(np.where(cld_subset[:,lt,ln] >= cld_minhgt)[0])*100./NTC_sub

cld_frac[cld_frac == 0] = cld_frac[cld_frac == 0]-0.001  # Set 0 to slightly below zero

##################  Retrieve dynamical background

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
SF_2022_mean = np.squeeze(np.mean(SF_2022, axis=0))

##################  Read in the tropopause information for the vertical sections

with open('/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/trop_prs_file.pkl', 'rb') as f:
    ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_jd, ERA5_trop_prs = pkl.load(f) 
trop_jd_ind0  = list(ERA5_trop_jd).index(2459792.5)
trop_jd_ind1  = list(ERA5_trop_jd).index(2459823.5)
trop_lat_ind0 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-20.))
trop_lat_ind1 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-40.))

mean_ERA5_trop_prs = np.mean(ERA5_trop_prs[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))
std_ERA5_trop_prs  = np.std(ERA5_trop_prs[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))

with open('/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/trop_theta_file.pkl', 'rb') as f:
    ERA5_trop_lon, ERA5_trop_lat, ERA5_trop_jd, ERA5_trop_th = pkl.load(f) 
trop_jd_ind0  = list(ERA5_trop_jd).index(2459792.5)
trop_jd_ind1  = list(ERA5_trop_jd).index(2459823.5)
trop_lat_ind0 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-20.))
trop_lat_ind1 = np.argmin(np.abs(np.asarray(ERA5_trop_lat)-40.))

mean_ERA5_trop_th = np.mean(ERA5_trop_th[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))
std_ERA5_trop_th  = np.std(ERA5_trop_th[trop_jd_ind0:trop_jd_ind1,trop_lat_ind0:trop_lat_ind1,:], axis=(0,1))

#################  Commence Plotting

fig = plt.figure(figsize=(9,8))   # Or 12x5 for the horizontal configuration

#####  Panel a: Convective fraction

ax1 = fig.add_subplot(2,2,(1,2),projection=ccrs.PlateCarree())

cmap = plt.get_cmap('YlOrRd')
cmap.set_under([1,1,1])

levels = MaxNLocator(nbins=10).tick_values(0.0, cbar_max)
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    
pixplt = ax1.pcolormesh(cloudlon[::4], cloudlat[::4], cld_frac[::4,::4], cmap = cmap, norm = norm)  # ::4 done to save file space
ax1.coastlines()
ax1.add_feature(cf.BORDERS)
#ax1.set_title(cld_choice + ' clouds above ' + str(cld_minhgt) + z_unit + ' for ' + str(date_rng), fontsize = 10)
#ax1.set_title('Deep convective fraction during August 2022', fontsize = 12)
#ax1.set_position([0.05,0.2,0.6,0.4])
ax1.set_position([0.08,0.15,0.72,0.72])

ax1.contour(SF_lon, SF_lat, SF_2022_mean, [1.0E7, 2.5E7], colors=['black'], linewidths = 2)
    
ax1.scatter(GV_lon, GV_lat, s=0.05, color='cyan')
#ax1.scatter(WB_lon, WB_lat, s=0.05, color='blue')

for WBflt in WB_FLT_DATES:  # This should not be necessary....but here we are
    WB_lon_flt = get_ACCLIP_merge_data(WB_OBSROOT, [WBflt], 'G_LONG_MMS')
    WB_lat_flt = get_ACCLIP_merge_data(WB_OBSROOT, [WBflt], 'G_LAT_MMS')
    ax1.scatter(WB_lon_flt, WB_lat_flt, s=0.05, color='blue')
    
#for WBflt in WB_DATES:  # This should not be necessary....but here we are
#    WB_lon = get_ACCLIP_merge_data(WB_OBSROOT, [WBflt], 'G_LONG_MMS')
#    WB_lat = get_ACCLIP_merge_data(WB_OBSROOT, [WBflt], 'G_LAT_MMS')
#    ax1.scatter(WB_lon, WB_lat, s=0.05, color='blue')    

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,  color='gray', alpha=0.5, linestyle='--', xlocs=[30,60,90,120,150,179.999], ylocs=[0,20,40,60])
deg = u"\N{DEGREE SIGN}"
ax1.set_xticks([60,90,120,150])
ax1.set_xticklabels(['60' + deg,'90' + deg,'120' + deg,'150' + deg], fontsize=12)
ax1.set_yticks([0,20,40,60])
ax1.set_yticklabels(['0' + deg,'20' + deg,'40' + deg,'60' + deg], fontsize=12)

ax1.set_xlim(50,160)
ax1.set_ylim(0,50)



#####  Panel b: Vertical flight tracks in pressure and altitude

ax2 = fig.add_subplot(2,2,3)
ax2t = ax2.twinx()

# Set up a colorbar
cmap = plt.get_cmap('YlOrRd')
cmap.set_over('black')
#cmap.set_under('white')
#levels = [70,100,150,200,300,400]
levels = np.arange(20,201,20)
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
    
ax2.fill_between(ERA5_trop_lon, mean_ERA5_trop_prs-std_ERA5_trop_prs, mean_ERA5_trop_prs+std_ERA5_trop_prs, color=[0.7,0.7,0.7])
ax2.plot(ERA5_trop_lon, mean_ERA5_trop_prs, color='black', linestyle='-')
ax2t.plot(np.zeros(10), sh.prs2alt(np.arange(100,1001,100), 7.6, 1013.))

#ax.scatter(WB_lon, WB_alt, s=2, c=WB_CO, cmap=cmap, norm=norm, marker='s', label='NASA WB-57')
#sc = ax.scatter(GV_lon, GV_alt, s=2, c=GV_CO, cmap=cmap, norm=norm, marker='^', label='NSF NCAR GV')

sc = ax2.scatter(WB_lon, WB_prs, s=1, c='blue')
ax2.scatter(GV_lon, GV_prs, s=1, c='cyan')

#plt.scatter(WB_lon, WB_prs, s=1, color='blue')
#plt.scatter(GV_lon, GV_prs, s=1, color='cyan')

#cbar = plt.colorbar(sc, extend='both', label='CO (ppbv)')

#plt.ylim(1000,40)
#plt.xticks(fontsize=12)
#plt.yticks([10,12,14,16,18], fontsize=12)

ax2.set_xlabel('Longitude (E)', fontsize=12)
ax2t.set_ylabel('Altitude (km)', fontsize=12)
ax2.set_ylabel('Pressure (hPa)', fontsize=12)

ax2.set_xlim(120,160)

ax2t.set_yticks([0,5,10,15,20])

ax2.set_yscale('log')
ax2.set_yticks([1000,400,200,100])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.get_yaxis().get_major_formatter().labelOnlyBase = False

ax2.set_ylim(1000,60)
ax2t.set_ylim(0,sh.prs2alt(60,8.0,1013.))

#####  Panel c: Vertical flight tracks in theta

ax3 = fig.add_subplot(2,2,4)

ax3.fill_between(ERA5_trop_lon, mean_ERA5_trop_th-std_ERA5_trop_th, mean_ERA5_trop_th+std_ERA5_trop_th, color=[0.7,0.7,0.7])
ax3.plot(ERA5_trop_lon, mean_ERA5_trop_th, color='black', linestyle='-')

#ax.scatter(WB_lon, WB_alt, s=2, c=WB_CO, cmap=cmap, norm=norm, marker='s', label='NASA WB-57')
#sc = ax.scatter(GV_lon, GV_alt, s=2, c=GV_CO, cmap=cmap, norm=norm, marker='^', label='NSF NCAR GV')
 
sc = ax3.scatter(WB_lon, WB_th, s=1, c='blue')
ax3.scatter(GV_lon, GV_th, s=1, c='cyan')

#plt.scatter(WB_lon, WB_prs, s=1, color='blue')
#plt.scatter(GV_lon, GV_prs, s=1, color='cyan')

#cbar = plt.colorbar(sc, extend='both', label='CO (ppbv)')

#plt.ylim(1000,40)
#ax.set_yscale('log')
#plt.xticks(fontsize=12)
plt.yticks([300,350,400,450], fontsize=12)
ax3.set_xlabel('Longitude (E)', fontsize=12)
ax3.set_ylabel('Potential Temperature (K)', fontsize=12)

ax3.set_xlim(120,160)
ax3.set_ylim(300,460)

#####  Finalize 

plt.figtext(0.12,0.52,'(a)',fontsize=20, fontweight='bold')
plt.figtext(0.40,0.11,'(b)',fontsize=20, fontweight='bold')
plt.figtext(0.91,0.11,'(c)',fontsize=20, fontweight='bold')



#ax1.set_position([0.00,0.30,0.4,0.3])
#ax2.set_position([0.42,0.30,0.2,0.3])
#ax3.set_position([0.73,0.30,0.2,0.3])

cax = fig.add_axes([0.87,0.50,0.02,0.32])
cbar = plt.colorbar(pixplt, extend = 'max', orientation='vertical', cax=cax)
cbar.ax.set_ylabel('Aug 2022 Conv. Frac. >' + str(int(cld_minhgt)) + z_unit + ' (%)', fontsize = 12)
cbar.ax.tick_params(labelsize = 12)

plt.tight_layout()

ax1.set_position([0.11,0.45,0.73,0.4])

plt.savefig(PLOTDIR + 'Fig1.png', dpi=300)
plt.savefig(PLOTDIR + 'Fig1.pdf', dpi=300)
plt.close('all')



