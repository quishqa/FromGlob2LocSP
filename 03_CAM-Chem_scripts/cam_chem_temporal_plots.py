#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 01:07:38 2020

@author: mgavidia
"""

import pickle
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import wrf_sp_eval.qualar_py as qr
import wrf_sp_eval.data_preparation as dp
import aas4wrf as a4w
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# Opening global model output
cam = xr.open_dataset("./cam_oct_14.nc")

# Opening wrfinput 
wrf = xr.open_dataset('wrfinput_d02')
xlat = wrf.XLAT.values[0, :, :]
xlon = wrf.XLONG.values[0, :, :]

# loading RMSP shp
reader = shpreader.Reader("masp_shp/masp_shp.shp")
rmsp = list(reader.geometries())
RMSP = cfeature.ShapelyFeature(rmsp, ccrs.PlateCarree())


# Selecting SPMA from CAM - Using WRF-Chem D02 reference q  
cam_sp = cam.sel(lat=slice(xlat.min(), xlat.max()),
                 lon=slice(xlon.min() % 360, xlon.max() % 360))

lon1d = cam_sp.lon.values
lat1d = cam_sp.lat.values

# Changing Longitud 0-360 to -180 - 180
lon1d3 = ((lon1d + 180) % 360) - 180

## Calculating cell borders
lon1d_b = a4w.cell_bound(lon1d)
lat1d_b = a4w.cell_bound(lat1d)
lon1d_b3 = ((lon1d_b + 180) % 360) - 180
lat1d_b3 = ((lat1d_b + 180) % 360) - 180

lon, lat = np.meshgrid(lon1d, lat1d)
lon_b, lat_b = np.meshgrid(lon1d_b, lat1d_b)


# Opening station loactions
aqs = pd.read_csv("cetesb2017_latlon.dat")
aqs['lon360'] = aqs.lon % 360


# Creating filter to select aqs inside each cell
filter_rmsp = ((aqs.lon360 > lon1d_b[0]) & (aqs.lon360 < lon1d_b[-1]) &
                (aqs.lat > lat1d_b[0]) & (aqs.lat < lat1d_b[-1]))

filter_A = ((aqs.lon360 > lon1d_b[0]) & (aqs.lon360 < lon1d_b[1]) & 
            (aqs.lat > lat1d_b[1]) & (aqs.lat < lat1d_b[2]))
filter_B = ((aqs.lon360 > lon1d_b[1]) & (aqs.lon360 < lon1d_b[2]) &
            (aqs.lat > lat1d_b[1]) & (aqs.lat < lat1d_b[2]))
filter_C = ((aqs.lon360 > lon1d_b[0]) & (aqs.lon360 < lon1d_b[1]) &
            (aqs.lat > lat1d_b[0]) & (aqs.lat < lat1d_b[1]))
filter_D = ((aqs.lon360 > lon1d_b[1]) & (aqs.lon360 < lon1d_b[2]) &
            (aqs.lat > lat1d_b[0]) & (aqs.lat < lat1d_b[1]))

# Filter AQS by cell
aqsA = aqs[filter_A]
aqsB = aqs[filter_B]
aqsC = aqs[filter_C]
aqsD = aqs[filter_D]
aqsRMSP = aqs[filter_rmsp]

print("AQS in CellA: ", len(aqsA.index))
print("AQS in CellB: ", len(aqsB.index))
print("AQS in CellC: ", len(aqsC.index))
print("AQS in CellD: ", len(aqsD.index))
print("AQS in RMSP: ", len(aqsRMSP.index))


# A plot with stations inside CAM-Chem cells
fig = plt.figure(figsize=(11,10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon1d_b.min(), lon1d_b.max(),
               lat1d_b.min(), lat1d_b.max()])
ax.coastlines(resolution='10m', zorder=150)
ax.add_feature(cfeature.STATES.with_scale('10m'), zorder=150)
ax.add_feature(RMSP, facecolor="none", edgecolor="black", zorder=150)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  color="darkgray", linewidth=1.5, zorder=50)
ax.scatter(aqsA.lon360, aqsA.lat, marker="o", color="black", zorder=200)
ax.scatter(aqsB.lon360, aqsB.lat, marker="o", color="red", zorder=200)
ax.scatter(aqsD.lon360, aqsD.lat, marker="o", color="blue",zorder=200)
gl.xlocator = mticker.FixedLocator(lon1d_b3)
gl.ylocator = mticker.FixedLocator(lat1d_b)
ax.set_yticks(np.concatenate((lat[:,0], lat_b[:,0])),
              crs=ccrs.PlateCarree())
ax.set_xticks(((np.concatenate((lon[0], lon_b[0])) + 180) % 360) - 180,
              crs=ccrs.PlateCarree())
ax.text(lon1d3[0], lat1d[1] - 0.075, "A",
        transform=ccrs.PlateCarree(), fontsize=22)
ax.text(lon1d3[1], lat1d[1] - 0.075, "B",
        transform=ccrs.PlateCarree(), fontsize=22)
ax.text(lon1d3[0], lat1d[0] - 0.075, "C",
        transform=ccrs.PlateCarree(), fontsize=22)
ax.text(lon1d3[1], lat1d[0] - 0.075, "D",
        transform=ccrs.PlateCarree(), fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.savefig("CAM-Chem-cell.pdf", dpi=300, bbox_inches='tight')

def plot_aqs_cell(aqs_df, ax='None'):
    if ax is None:
        ax.plt.gca()
    # ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lon1d_b.min(), lon1d_b.max(),
                   lat1d_b.min(), lat1d_b.max()])
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.STATES.with_scale('10m'))
    ax.add_feature(RMSP, facecolor="none", edgecolor="black")
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, alpha=0.9)
    ax.scatter(aqs_df.lon360, aqs_df.lat, marker="o", color="black")
    ax.set_title("Air quality stations")
    gl.xlocator = mticker.FixedLocator(lon1d_b3)
    gl.ylocator = mticker.FixedLocator(lat1d_b)




# Loading all data from CETESB aqs
a_dict = open("cetesb_2014_all_pol.pkl", "rb")
cetesb = pickle.load(a_dict)
a_dict.close()

cet_aqs_A = {aqs: cetesb[aqs] for aqs in aqsA.name}
cet_aqs_B = {aqs: cetesb[aqs] for aqs in aqsB.name}
cet_aqs_D = {aqs: cetesb[aqs] for aqs in aqsD.name}
cet_aqs_MASP = {aqs: cetesb[aqs] for aqs in aqsRMSP.name}



def stations_with_data(cet_df):
    means = {aqs:cet_df[aqs].mean().to_frame().T for aqs in list(cet_df)}
    means_df = pd.concat(means)
    return means_df.notna().sum()
                                                                    
aqs_dict = dict(A=stations_with_data(cet_aqs_A), 
                B=stations_with_data(cet_aqs_B),
                D=stations_with_data(cet_aqs_D))       

aqs_total_df = pd.DataFrame.from_dict(aqs_dict, orient="index")                                                            


def group_aqs_df(aqs_dic):
    aqs_merge = pd.concat(aqs_dic.values())
    aqs_group = aqs_merge.groupby(aqs_merge.index).mean()
    return aqs_group

cet_A = group_aqs_df(cet_aqs_A)['2014-10-06': '2014-10-12']
cet_B = group_aqs_df(cet_aqs_B)['2014-10-06': '2014-10-12']
cet_D = group_aqs_df(cet_aqs_D)['2014-10-06': '2014-10-12']
cet_MASP = group_aqs_df(cet_aqs_MASP)['2014-10-06': '2014-10-12']


# Function to mol/mol to ug/m3

def molmol_to_ugm3(df, pol, M):
    R = 8.3142 # J/K/mol
    pol_ugm3 = (df[pol] * 10**6 * # to ppm
                df.lev * 10 **2 * # to Pa
                M / R / df.tk)
    return pol_ugm3


# cam to datframe
cell_A_lon, cell_A_lat = lon1d[0], lat1d[1]
cell_B_lon, cell_B_lat = lon1d[1], lat1d[1]
cell_C_lon, cell_C_lat = lon1d[0], lat1d[0]
cell_D_lon, cell_D_lat = lon1d[1], lat1d[0]


def ds_to_df(cam_sp, cell_lat, cell_lon):
    df = pd.DataFrame({
        'date': cam_sp.time.values,
        'tk': cam_sp['T'].sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values,
        'O3': cam_sp.O3.sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values,
        'lev': cam_sp.lev.values[-1],
        'CO': cam_sp.CO.sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values,
        'NO': cam_sp.NO.sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values,
        'NO2': cam_sp.NO2.sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values,
        'PM25': cam_sp.PM25.sel(lat=cell_lat, lon=cell_lon, lev=1000,
                              method='nearest').values
        })
    df.set_index('date', inplace=True)
    df.index = df.index.tz_localize('UTC').tz_convert("America/Sao_Paulo")
    df['o3'] = molmol_to_ugm3(df, 'O3', 48)
    df['no'] = molmol_to_ugm3(df, 'NO', 30)
    df['no2'] = molmol_to_ugm3(df, 'NO2', 46)
    df['co'] = df.CO * 10**6
    df['pm25'] = df.PM25 * 10**9
    
    return df
    


def cam_sp_df(cam_sp):
    df = pd.DataFrame({
        'date': cam_sp.time.values,
        'tk': cam_sp['T'].sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values,
        'O3': cam_sp.O3.sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values,
        'lev': cam_sp.lev.values[-1],
        'CO': cam_sp.CO.sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values,
        'NO': cam_sp.NO.sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values,
        'NO2': cam_sp.NO2.sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values,
        'PM25': cam_sp.PM25.sel(lev=1000, method='nearest').mean(dim=['lat', 'lon']).values
        })
    df.set_index('date', inplace=True)
    df.index = df.index.tz_localize('UTC').tz_convert("America/Sao_Paulo")
    df['o3'] = molmol_to_ugm3(df, 'O3', 48)
    df['no'] = molmol_to_ugm3(df, 'NO', 30)
    df['no2'] = molmol_to_ugm3(df, 'NO2', 46)
    df['co'] = df.CO * 10**6
    df['pm25'] = df.PM25 * 10**9    
    return df


cam_D = ds_to_df(cam_sp, cell_D_lat, cell_D_lon)['2014-10-06': '2014-10-12']
cam_B = ds_to_df(cam_sp, cell_B_lat, cell_B_lon)['2014-10-06': '2014-10-12']
cam_A = ds_to_df(cam_sp, cell_A_lat, cell_A_lon)['2014-10-06': '2014-10-12']
cam_MASP = cam_sp_df(cam_sp)['2014-10-06':'2014-10-12']


# A plot 
# fig, ax = plt.subplots(1, 1)
# cet_A.o3.plot(label="CETESB", ax=ax, linewidth=1.5, marker='o')
# cam_A.O3.plot(label="CAM-Chem", ax =ax, linewidth=1.5, marker='d')
# ax.set_ylabel("$O_3 \; (\mu g / m^3)$")
# ax.set_xlabel('')
# ax.legend()


def pol_plot(cet, cam, pol, ylab, ax=None):
    if ax is None:
        ax=plt.gca()
    date_format = mdates.DateFormatter("%d")
    ax.plot(cet[pol], color='black', 
            label="CETESB", marker="o")
    ax.plot(cam[pol], color='red',
            label="CAM-Chem", marker="d")
    if pol != "co":
        ax.set_ylabel("$" + ylab + "\; (\mu g \;  m^{-3})$")
        ax.set_ylim([0, max(cam[pol].max(),
                    cet[pol].max()) + 10])
    else:
        ax.set_ylabel("$" + ylab + "\; (ppm)$")
        ax.set_ylim([0, max(cam[pol].max(),
                    cet[pol].max()) + 0.25])
    ax.set_xlabel("October 2014")
    ax.xaxis.set_major_formatter(date_format)
    ax.legend(frameon=False)


def hourly_prof(cet, cam, pol, ylab, ax=None):
    if ax is None:
        ax.plt.gca()
    ax.plot(cet[pol], color="black",
            label="CETESB", marker="o")
    ax.plot(cam[pol], color="red",
            label="CAM-Chem", marker="d")
    ax.set_xlabel("Hours (LT)")
    if pol != "co":
        ax.set_ylabel("$" + ylab + "\; (\mu g \;  m^{-3})$")
        ax.set_ylim([0, max(cam[pol].max(),
                    cet[pol].max()) + 10])
    else:
        ax.set_ylabel("$" + ylab + "\; (ppm)$")
        ax.set_ylim([0, max(cam[pol.upper()].max(),
                    cet[pol].max()) + 0.25])
    ax.legend(frameon=False)


def all_hourly_day(cet, cam):
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    axes = ax.flatten()
    pol_plot(cet, cam, "o3", "O_3", axes[0])
    pol_plot(cet, cam, "co", "CO", axes[1])
    pol_plot(cet, cam, "no", "NO", axes[2])
    pol_plot(cet, cam, "no2", "NO_2", axes[3])
    pol_plot(cet, cam, "pm25", "PM_{2.5}", axes[4])
    axes[5].axis("off")


    
def all_daily_day(cet, cam):
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    axes = ax.flatten()
    pol_plot(cet, cam, "o3", "O_3", axes[0])
    pol_plot(cet, cam, "co", "CO", axes[1])
    pol_plot(cet, cam, "no", "NO", axes[2])
    pol_plot(cet, cam, "no2", "NO_2", axes[3])
    pol_plot(cet, cam, "pm25", "PM_{2.5}", axes[4])
    axes[5].axis("off")

    
def all_hourly_prof(cet, cam):
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    axes = ax.flatten()
    hourly_prof(cet, cam, "o3", "O_3", axes[0])
    hourly_prof(cet, cam, "co", "CO", axes[1])
    hourly_prof(cet, cam, "no", "NO", axes[2])
    hourly_prof(cet, cam, "no2", "NO_2", axes[3])
    hourly_prof(cet, cam, "pm25", "PM_{2.5}", axes[4])
    axes[5].axis("off")


def daily_and_hour_mean(cet, cam):
    daily_cet = cet.groupby(cet.index.day).mean()
    daily_cam = cam.groupby(cam.index.day).mean()
    hourly_cet = cet.groupby(cet.index.hour).mean()
    hourly_cam = cam.groupby(cam.index.hour).mean()
    return daily_cet, hourly_cet, daily_cam, hourly_cam


day_cet_A, hour_cet_A, day_cam_A, hour_cam_A = daily_and_hour_mean(cet_A, cam_A)
day_cet_B, hour_cet_B, day_cam_B, hour_cam_B = daily_and_hour_mean(cet_B, cam_B)
day_cet_D, hour_cet_D, day_cam_D, hour_cam_D = daily_and_hour_mean(cet_D, cam_D)
day_cet_MASP, hour_cet_MASP, day_cam_MASP, hour_cam_MASP = daily_and_hour_mean(cet_MASP, cam_MASP)

# all_hourly_day(cet_A, cam_A)
# all_daily_day(day_cet_A, day_cam_A)
# all_hourly_prof(hour_cet_A, hour_cam_A)

# all_hourly_day(cet_B, cam_B)
# all_daily_day(day_cet_B, day_cam_B)
# all_hourly_prof(hour_cet_B, hour_cam_B)

# all_hourly_day(cet_D, cam_D)
# all_hourly_day(day_cet_D, day_cam_D)
# all_hourly_prof(hour_cet_D, hour_cam_D)


# # CAMMASP vs CET

# all_hourly_day(cet_MASP, cam_MASP)
# all_daily_day(day_cet_MASP, day_cam_MASP)
# all_hourly_prof(hour_cet_MASP, hour_cam_MASP)


def map_all_hourly_day(cet, cam, aqs_df):
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(321,  projection=ccrs.PlateCarree())
    ax2 = plt.subplot(322)
    ax3 = plt.subplot(323)
    ax4 = plt.subplot(324)
    ax5 = plt.subplot(325)
    ax6 = plt.subplot(326)
    plot_aqs_cell(aqs_df, ax1)
    pol_plot(cet, cam, "o3", "O_3", ax2)
    pol_plot(cet, cam, "no", "NO", ax3)
    pol_plot(cet, cam, "no2", "NO_2", ax4)
    pol_plot(cet, cam, "co", "CO", ax5)
    pol_plot(cet, cam, "pm25", "PM_{2.5}", ax6)



def map_all_daily_day(cet, cam, aqs_df):
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(321,  projection=ccrs.PlateCarree())
    ax2 = plt.subplot(322)
    ax3 = plt.subplot(323)
    ax4 = plt.subplot(324)
    ax5 = plt.subplot(325)
    ax6 = plt.subplot(326)
    plot_aqs_cell(aqs_df, ax1)
    pol_plot(cet, cam, "o3", "O_3", ax2)
    pol_plot(cet, cam, "no", "NO", ax3)
    pol_plot(cet, cam, "no2", "NO_2", ax4)
    pol_plot(cet, cam, "co", "CO", ax5)
    pol_plot(cet, cam, "pm25", "PM_{2.5}", ax6)

    
def map_all_hourly_prof(cet, cam, aqs_df):
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(321,  projection=ccrs.PlateCarree())
    ax2 = plt.subplot(322)
    ax3 = plt.subplot(323)
    ax4 = plt.subplot(324)
    ax5 = plt.subplot(325)
    ax6 = plt.subplot(326)
    plot_aqs_cell(aqs_df, ax1)
    hourly_prof(cet, cam, "o3", "O_3", ax2)
    hourly_prof(cet, cam, "no", "NO", ax3)
    hourly_prof(cet, cam, "no2", "NO_2", ax4)
    hourly_prof(cet, cam, "co", "CO", ax5)
    hourly_prof(cet, cam, "pm25", "PM_{2.5}", ax6)


map_all_hourly_day(cet_MASP, cam_MASP, aqsRMSP)
plt.savefig('masp_all_data.pdf', dpi = 300, bbox_inches='tight')
map_all_daily_day(day_cet_MASP, day_cam_MASP, aqsRMSP)
plt.savefig('masp_all_day.pdf', dpi = 300, bbox_inches='tight')
map_all_hourly_prof(hour_cet_MASP, hour_cam_MASP, aqsRMSP)
plt.savefig('masp_all_profile.pdf', dpi = 300, bbox_inches='tight')

map_all_hourly_day(cet_A, cam_A, aqsA)
plt.savefig('cell_A_all_data.pdf', dpi = 300, bbox_inches='tight')
map_all_daily_day(day_cet_A, day_cam_A, aqsA)
plt.savefig('cell_A_all_day.pdf', dpi = 300, bbox_inches='tight')
map_all_hourly_prof(hour_cet_A, hour_cam_A, aqsA)
plt.savefig('cell_A_all_profile.pdf', dpi = 300, bbox_inches='tight')


map_all_hourly_day(cet_B, cam_B, aqsB)
plt.savefig('cell_B_all_data.pdf', dpi = 300, bbox_inches='tight')
map_all_daily_day(day_cet_B, day_cam_B, aqsB)
plt.savefig('cell_B_all_day.pdf', dpi = 300, bbox_inches='tight')
map_all_hourly_prof(hour_cet_B, hour_cam_B, aqsB)
plt.savefig('cell_B_all_profile.pdf', dpi = 300,bbox_inches='tight')

map_all_hourly_day(cet_D, cam_D, aqsD)
plt.savefig('cell_D_all_data.pdf', dpi = 300, bbox_inches='tight')
map_all_daily_day(day_cet_D, day_cam_D, aqsD)
plt.savefig('cell_D_all_day.pdf', dpi = 300, bbox_inches='tight')
map_all_hourly_prof(hour_cet_D, hour_cam_D, aqsD)
plt.savefig('cell_D_all_profile.pdf', dpi = 300,bbox_inches='tight')


# Statistics 

import wrf_sp_eval.model_stats as ms

cam_MASP['name'] = 'CAM-Chem'
cet_MASP['name'] = 'CETESB'

vars_2_eval = ['o3', 'no', 'no2', 'co', 'pm25']
masp_stats = ms.some_vars_stats_per_station(cam_MASP, cet_MASP, 
                                            vars_2_eval, to_df=True)


def stats_calc(cam, cet, output_name):
    cam['name'] = 'CAM-Chem'
    cet['name'] = 'CETESB'
    vars_2_eval = ['o3', 'no', 'no2', 'co', 'pm25']
    stats = ms.some_vars_stats_per_station(cam, cet, 
                                           vars_2_eval, 
                                           to_df=True)
    output_file = output_name + '_stats.csv'
    stats.to_csv(output_file, sep=',', index_label="pol")
    return stats




st_MASP = stats_calc(cam_MASP, cet_MASP, "masp")    
st_A = stats_calc(cam_A, cet_A, "cell_A")    
st_B = stats_calc(cam_B, cet_B, "cell_B")    
st_D = stats_calc(cam_D, cet_D, "cell_D")    

print(st_MASP
      .drop(columns=["aqs"])
      .to_latex(float_format="%.2f"))
print(st_A
      .drop(columns=["aqs"])
      .to_latex(float_format="%.2f"))
print(st_D
      .drop(columns=["aqs"])
      .to_latex(float_format="%.2f"))
print(st_B
      .drop(columns=["aqs"])
      .to_latex(float_format="%.2f"))



import scipy.stats
import numpy as np

# From: http://onlinestatbook.com/2/estimation/correlation_ci.html
# https://medium.com/@shandou/how-to-compute-confidence-interval-for-pearsons-r-a-brief-guide-951445b9cb2d

o3_n = st_MASP.loc["o3"].N
o3_r = st_MASP.loc["o3"].R
alpha = 0.05
df = o3_n -2

t_critical = scipy.stats.t.ppf(1 - alpha /2, o3_n -2)

def t_person_test_for_R(r, n):
    t = r * np.sqrt(n - 2)/ np.sqrt(1 - r**2)
    return t

t_cal = t_person_test_for_R(o3_r, o3_n)





def r_pearson_significance(n, r, alpha, deg_free = 2):
    '''
    Calculate Pearson's R significance. With a two-tail
    test (non-directional).
    Based on:
    https://medium.com/@shandou/how-to-compute-confidence-interval-for-pearsons-r-a-brief-guide-951445b9cb2d

    Parameters
    ----------
    n : int
        sample size.
    r : float
        Pearson R.
    alpha : float
        test significant level.
    deg_free : int, optional
        degrees of freedom. The default is 2.

    Returns
    -------
    t_cal : float
        Calcualted t value.
    t_cri : float
        Critical t value.

    '''
    
    t_cri = scipy.stats.t.ppf(1 - alpha / 2.0, deg_free)
    t_cal = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
    if t_cal > t_cri:
        print("Significant linear relationship")
    else:
        print("No significant linear relationship")
    return (t_cal, t_cri)

# CI

alpha = 0.05 / 2 # two_tail test
z_critical = scipy.stats.norm.ppf(1 - alpha)
z_prime = 0.5 * np.log((1 + o3_r) / (1 - o3_r))

se = 1 / np.sqrt(o3_n - 3) # Sample standard error
CI_lower = z_prime - z_critical * se
CI_upper = z_prime + z_critical * se


def z_prime_to_r(z_prime):
    r = (np.exp(2 * z_prime) - 1) / (np.exp(2 * z_prime) + 1)
    return r


CI = (np.tanh(CI_lower), np.tanh(CI_upper))
def r_pearson_confidence_interval(n, r, alpha):
    '''
    Calculate Pearson's R confidence intervals, 
    using two-tail test.
    Based on:
    http://onlinestatbook.com/2/estimation/correlation_ci.html
    https://medium.com/@shandou/how-to-compute-confidence-interval-for-pearsons-r-a-brief-guide-951445b9cb2d

    Parameters
    ----------
    n : int
        sample size.
    r : float
        Pearson's R.
    alpha : float
        confidence level (e.g. if 95% then alpha = 0.05).

    Returns
    -------
    r_lower : float
        lower CI .
    r_upper : float
        upper CI.

    '''
    alph = 0.05 / 2.0 # two-tail test:
    z_critical = scipy.stats.norm.ppf(1 - alph)
    # r to z' by Fisher's z' transform:
    z_prime =0.5 * np.log((1 + r) / (1 - r))
    # Sample standard error:
    se = 1 / np.sqrt(n - 3)
    # Computing CI using z':
    ci_lower = z_prime - z_critical * se
    ci_upper = z_prime + z_critical * se
    # Converting z' back to r values:
    r_lower = np.tanh(ci_lower)
    r_upper = np.tanh(ci_upper)
    return (r_lower, r_upper)
    


#  Checking R significance
for n, r in zip(st_D.N.values, st_D.R.values):
    r_pearson_significance(n, r, 0.05)
    print(r_pearson_confidence_interval(n, r, 0.05))
    
