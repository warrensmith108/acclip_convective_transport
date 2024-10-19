
import numpy as np
import numpy.ma as ma
import glob
import os
import math
import netCDF4 as nc4
import pickle as pkl
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import jdcal
from ACCLIP_Plotter import get_ACCLIP_merge_data
from plot_airplane_traj3d import get_processed_data
from Fig3 import moving_window_smooth as mws

def get_plot_fields(plot_datestr, plot_hr, ERA5_wind_lev, PROCDIR, GV_dates, WB_dates):

    # Calculate some necessary time variables 
    plot_jd = np.sum(jdcal.gcal2jd(2022, int(plot_datestr[0:2]), int(plot_datestr[2:]))) + (int(plot_hr[0:2])/24.) + (int(plot_hr[2:])/1440.)
    days_since_jan1 = math.floor(plot_jd - np.sum(jdcal.gcal2jd(2022, 1, 1))) + 1   # Adding 1 is necessary because the year begins from 1, not 0

    # Get the Tb field from the corresponding HIMAWARI file 
    HIMAWARI_FILE = '/home/wsmith/ACCLIP/HIMAWARI_Data/HIM8.LD.2022' + str(days_since_jan1) + '.' + plot_hr + '.04km.C14.nc'
    HIMnc = nc4.Dataset(HIMAWARI_FILE)
    HIMlon = HIMnc.variables['longitude'][:]
    HIMlat = HIMnc.variables['latitude'][:]
    HIMtb  = HIMnc.variables['TEMP'][:]
    HIMtb  = ma.getdata(HIMtb)/10.  # Unmask and scale the brightness temperature field

    # Get the wind vectors from ERA5
    ERA5nc = nc4.Dataset(ERA5_FILE)
    ERA5lon = ERA5nc.variables['Longitude'][:]
    ERA5lat = ERA5nc.variables['Latitude'][:]
    ERA5prs = ERA5nc.variables['Altitude'][:]
    ERA5jd  = ERA5nc.variables['Julian_day'][:] + (ERA5nc.variables['Seconds'][:])/86400. - 0.5   # Dumb TRAJ3D correction

    prs_ind0  = list(ERA5prs).index(ERA5_wind_lev[0])
    prs_ind1  = list(ERA5prs).index(ERA5_wind_lev[1])
    time_ind = list(ERA5jd).index(plot_jd)
    ERA5_u0 = ERA5nc.variables['u'][time_ind,prs_ind0,:,:]
    ERA5_v0 = ERA5nc.variables['v'][time_ind,prs_ind0,:,:]
    ERA5_u1 = ERA5nc.variables['u'][time_ind,prs_ind1,:,:]
    ERA5_v1 = ERA5nc.variables['v'][time_ind,prs_ind1,:,:]

    # Read the relevant ACCLIP flight tracks
    GV_lon = get_ACCLIP_merge_data(GV_OBSROOT, GV_dates, 'GGLON')
    GV_lat = get_ACCLIP_merge_data(GV_OBSROOT, GV_dates, 'GGLAT')
    WB_lon = get_ACCLIP_merge_data(WB_OBSROOT, WB_dates, 'G_LONG_MMS')
    WB_lat = get_ACCLIP_merge_data(WB_OBSROOT, WB_dates, 'G_LAT_MMS')

    #GV_Aug6_lon = get_ACCLIP_merge_data(GV_OBSROOT, ['20220819'], 'GGLON')
    #GV_Aug6_lat = get_ACCLIP_merge_data(GV_OBSROOT, ['20220819'], 'GGLAT')
    #WB_Aug6_lon = get_ACCLIP_merge_data(WB_OBSROOT, ['20220819'], 'G_LONG_MMS')
    #WB_Aug6_lat = get_ACCLIP_merge_data(WB_OBSROOT, ['20220819'], 'G_LAT_MMS')

    #GV_CO = get_ACCLIP_merge_data(GV_OBSROOT, ['20220819'], 'CO_AERODYNE')
    #WB_CO = get_ACCLIP_merge_data(WB_OBSROOT, ['20220819'], 'CO_COLD2_ppbv')


    # Read the cloud intercept information to line it up with the brightness temperature
    init_lon_GV, init_lat_GV, init_prs_GV, init_jd, dep_GV, int_lon_GV, int_lat_GV, pct \
            = get_processed_data(PROCDIR, GV_dates, 'GV', 'ERA5_kin', [0,500], [-1,360,-90,90], 30, 'finalcld', '', '', unit='hPa')
    init_lon_WB, init_lat_WB, init_prs_WB, init_jd, dep_WB, int_lon_WB, int_lat_WB, pct \
            = get_processed_data(PROCDIR, WB_dates, 'WB57', 'ERA5_kin', [0,500], [-1,360,-90,90], 30, 'finalcld', '', '', unit='hPa')

    int_lon = np.asarray(int_lon_GV + int_lon_WB)
    int_lat = np.asarray(int_lat_GV + int_lat_WB)

    pixels = np.histogram2d(int_lon, int_lat, bins = [360,180], range = [[0,360],[-90,90]])
    pixlon = (pixels[1][0:-1] + pixels[1][1:])/2.  # Note we need to average these to use contour, we need the pixel centers, not the edges
    pixlat = (pixels[2][0:-1] + pixels[2][1:])/2.
    pixdata = np.transpose(pixels[0])
    pixdata = pixdata*100./len(int_lon)  # Normalize the fraction of all trajectories (pct replaced with 100, because processed files now have missing points included)

    return HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata
    
    
def plot_HIMAWARI_overlay(ax, HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata, cmap, norm):

    pixplt = ax.pcolormesh(HIMlon[::10], HIMlat[::10], HIMtb[::10,::10], cmap = cmap, norm = norm)  # Subsetting by 10 is done to save PDF file size

    ax.coastlines()
    ax.add_feature(cf.BORDERS)

    ax.plot(GV_lon, GV_lat, color='cyan', linewidth=2)
    ax.plot(WB_lon, WB_lat, color='blue', linewidth=2)

    ax.contour(pixlon, pixlat, mws(pixdata, width=1), levels=[0.3], colors=['red'], linewidths=[3])

    #cmap2 = plt.get_cmap('YlOrRd')
    #cmap2.set_under('white')
    #cmap2.set_over('magenta')
    #levels2 = MaxNLocator(nbins=10).tick_values(0,200)
    #norm2 = BoundaryNorm(levels2, ncolors = cmap.N, clip = False)

    #ax.scatter(GV_Aug6_lon, GV_Aug6_lat, s=2, c=GV_CO, cmap=cmap2, norm=norm2)
    #sc = ax.scatter(WB_Aug6_lon, WB_Aug6_lat, s=2, c=WB_CO, cmap=cmap2, norm=norm2)

    #cbar2 = plt.colorbar(sc,fraction=0.015, pad=0.03, extend = 'both')
    #cbar2.ax.set_ylabel('Observed CO (ppbv)', fontsize = 8)
    #cbar2.ax.tick_params(labelsize = 10)
    #cbar2.ax.set_position([0.85,0.55,0.02,0.3])

    #ax.quiver(ERA5lon[::wind_sub], ERA5lat[::wind_sub], ERA5_u0[::wind_sub,::wind_sub], ERA5_v0[::wind_sub,::wind_sub], minlength = 0.1, scale = 500, width=0.005, color='blue')
    ax.quiver(ERA5lon[::wind_sub], ERA5lat[::wind_sub], ERA5_u0[::wind_sub,::wind_sub], ERA5_v0[::wind_sub,::wind_sub], minlength = 0.1, scale = 500, color='gray')      # Use this for the Aug 6 overview
    #ax.quiver(ERA5lon[::wind_sub], ERA5lat[::wind_sub], ERA5_u1[::wind_sub,::wind_sub], ERA5_v1[::wind_sub,::wind_sub], minlength = 0.1, scale = 500, color='magenta')      

    ax.set_xlim(100,155)
    ax.set_ylim(25,50)

    #ax.set_title('HIMAWARI Brightness Temperature for 2022' + plot_datestr + '_' + plot_hr[0:2] + 'Z', fontsize=12)

# File paths
PROCDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/vars/processed/'
PLOTDIR = '/home/wsmith/traject/TRAJ3D_ACCLIP_FLIGHTS/2022_RFs/plots/'
ERA5_FILE = '/data/wsmith/ERA5/ERA5_TRAJ3D_2022.nc'
GV_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/GV/NCAR_60SEC_Merge/R1/icartt/'      # Instrument/AERODYNE/R1.1/'
WB_OBSROOT = '/UTLS/Field_Data/ACCLIP_Data/Final_Aircraft_Data/WB57/NCAR_60SEC_Merge/R1/icartt/'

# User settings
plot_datestrs = ['0804','0815','0818']  # We assume 2022 because you know, ACCLIP
plot_hrs = ['1200','0000','1800']
ERA5_wind_lev = [150.,850.]  # hPa
wind_sub = 16


############## Make a dank plot

fig = plt.figure(figsize=(8,10))
ax1 = fig.add_subplot(3,1,1,projection=ccrs.PlateCarree())
ax2 = fig.add_subplot(3,1,2,projection=ccrs.PlateCarree())
ax3 = fig.add_subplot(3,1,3,projection=ccrs.PlateCarree())

cmap = plt.get_cmap('Greys_r')
cmap.set_under([0,0,0])
cmap.set_over([1,1,1])
levels = MaxNLocator(nbins=7).tick_values(210, 280)
norm = BoundaryNorm(levels, ncolors = cmap.N, clip = False)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    
HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata = \
    get_plot_fields(plot_datestrs[0], plot_hrs[0], ERA5_wind_lev, PROCDIR, ['20220806','20220807'], ['20220806'])
plot_HIMAWARI_overlay(ax1, HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata, cmap, norm)

HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata = \
    get_plot_fields(plot_datestrs[1], plot_hrs[1], ERA5_wind_lev, PROCDIR, ['20220815'], ['20220815'])
plot_HIMAWARI_overlay(ax2, HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata, cmap, norm)

HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata = \
    get_plot_fields(plot_datestrs[2], plot_hrs[2], ERA5_wind_lev, PROCDIR, ['20220819'], ['20220819'])
plot_HIMAWARI_overlay(ax3, HIMlon, HIMlat, HIMtb, ERA5lon, ERA5lat, ERA5_u0, ERA5_v0, ERA5_u1, ERA5_v1, GV_lon, GV_lat, WB_lon, WB_lat, pixlon, pixlat, pixdata, cmap, norm)
   
   

cbar = plt.colorbar(sm, extend = 'both')
cbar.ax.set_ylabel('HIMAWARI Brightness Temperature (K)', fontsize = 12)
#cbar.ax.set_ylabel('Brightness Temperature (K)', fontsize=8)
cbar.ax.tick_params(labelsize = 10)
cbar.ax.set_position([0.88,0.3,0.02,0.4])

ax3.set_position([0.2,0.08,0.62,0.28])  # Not sure why this axis is misbehaving

plt.figtext(0.50,0.665,'(a) Aug 6-7, 2022', fontsize=16, fontweight='bold', backgroundcolor='white', color='k')
plt.figtext(0.50,0.395,'(b) Aug 15, 2022', fontsize=16, fontweight='bold', backgroundcolor='white', color='k')
plt.figtext(0.50,0.12,'(c) Aug 19, 2022', fontsize=16, fontweight='bold', backgroundcolor='white', color='black')

#ax1.scatter([116.4], [39.9], s=20, color='yellow', zorder=999)
#ax2.scatter([116.4], [39.9], s=20, color='yellow', zorder=999)   # Beijing, but this doesn't help
#ax3.scatter([116.4], [39.9], s=20, color='yellow', zorder=999)

#ax.set_xlim(110,135)
#ax.set_ylim(24,44)

#ax.set_title('Overview for ACCLIP Flights on August 19, 2022')

plt.savefig(PLOTDIR + 'Fig5.png', dpi=300)
plt.savefig(PLOTDIR + 'Fig5.pdf', dpi=300)
plt.close('all')

