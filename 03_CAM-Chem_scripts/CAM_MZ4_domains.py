#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 18:13:11 2020

@author: mgavidia
"""

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


# Opening global model output
cam = xr.open_dataset("./cam_oct_14.nc")
mz4 = xr.open_dataset("./mz4_oct_14.nc", decode_times=False)
mz4.coords["time"] = cam.coords["time"]

# Opening wrfinput 
wrf = xr.open_dataset('wrfinput_d01')
wrf_d02 = xr.open_dataset('wrfinput_d02')
xlat = wrf.XLAT.values[0, :, :]
xlon = wrf.XLONG.values[0, :, :]

xlat_d02 = wrf_d02.XLAT.values[0, :, :]
xlon_d02 = wrf_d02.XLONG.values[0, :, :]

# loading RMSP shp
reader = shpreader.Reader("masp_shp/masp_shp.shp")
rmsp = list(reader.geometries())
RMSP = cfeature.ShapelyFeature(rmsp, ccrs.PlateCarree())


# Cutting for the domain 1
cam_sp = cam.sel(lat=slice(xlat.min(), xlat.max()),
                 lon=slice(xlon.min() % 360, xlon.max() % 360))
mz4_sp = mz4.sel(lat=slice(xlat.min(), xlat.max()),
                 lon=slice(xlon.min() % 360, xlon.max() % 360))

# Cutting for domain 2
cam_rmsp = cam.sel(lat=slice(xlat_d02.min(), xlat_d02.max()),
                 lon=slice(xlon_d02.min() % 360, xlon_d02.max() % 360))



# coordinates for doms1
lon1d_cam = cam_sp.lon.values
lat1d_cam = cam_sp.lat.values

lon1d_mz4 = mz4_sp.lon.values
lat1d_mz4 = mz4_sp.lat.values

lon1d3_cam = ((lon1d_cam + 180) % 360) - 180
lon1d3_mz4 = ((lon1d_mz4 + 180) % 360) - 180

## Calculating cell borders for doms 1
lon1d_cam_b = a4w.cell_bound(lon1d_cam)
lat1d_cam_b = a4w.cell_bound(lat1d_cam)

lon1d_mz4_b = a4w.cell_bound(lon1d_mz4)
lat1d_mz4_b = a4w.cell_bound(lat1d_mz4)


lon_cam, lat_cam = np.meshgrid(lon1d_cam, lat1d_cam)
lon_cam_b, lat_cam_b = np.meshgrid(lon1d_cam_b, lat1d_cam_b)

lon_mz4, lat_mz4 = np.meshgrid(lon1d_mz4, lat1d_mz4)
lon_mz4_b, lat_mz4_b = np.meshgrid(lon1d_mz4_b, lat1d_mz4_b)


# coordinates for dom1
lon1d_cam_rmsp = cam_rmsp.lon.values
lat1d_cam_rmsp = cam_rmsp.lat.values
lon1d3_cam_rmsp = ((lon1d_cam_rmsp + 180) % 360) - 180

lon1d_cam_rmsp_b = a4w.cell_bound(lon1d_cam_rmsp)
lat1d_cam_rmsp_b = a4w.cell_bound(lat1d_cam_rmsp)

lon_cam_rmsp, lat_cam_rmsp = np.meshgrid(lon1d_cam_rmsp, lat1d_cam_rmsp)
lon_cam_rmsp_b, lat_cam_rmsp_b = np.meshgrid(lon1d_cam_rmsp_b, lat1d_cam_rmsp_b)

# Getting data to plot
mz4_co_mean = (mz4_sp.CO_VMR_inst
               .mean(dim='time')
               .sel(lev=1000, method='nearest') * 10**6)
cam_co_mean = (cam_sp.CO
               .mean(dim='time')
               .sel(lev=1000, method='nearest') * 10**6)

cam_co_rmsp_mean = (cam_rmsp.CO
                    .mean(dim='time')
                    .sel(lev=1000, method='nearest') * 10**6)



xticks = np.round(np.arange(xlon.min(), xlon.max(), 2.0), 1)[1:]
yticks = np.round(np.arange(xlat.min(), xlat.max(), 1.5), 1)[1:]

xticks_d02 = np.round(np.arange(xlon_d02.min(), xlon_d02.max(), 1.5), 1)[1:]
yticks_d02 = np.round(np.arange(xlat_d02.min(), xlat_d02.max(), 1.5), 1)[1:]


# Plot comparing MZ4 and CAMS domains
fig, axes = plt.subplots(1, 2, figsize=[10, 4],
                         subplot_kw={'projection': ccrs.PlateCarree()})
mz4_co_mean.plot(ax=axes[0], cmap='viridis', vmin=0, vmax=0.45, 
                 cbar_kwargs={'shrink':0.5, 'label':'CO [ppm]'})
axes[0].set_title("MOZART-4")
cam_co_mean.plot(ax=axes[1], cmap='viridis', vmin=0, vmax=0.45,
                 cbar_kwargs={'shrink':0.5, 'label':'CO [ppm]'})
axes[1].set_title("CAM-Chem")
for ax in axes:
    ax.set_extent([xlon.min(), xlon.max(), xlat.min(), xlat.max()])
    ax.coastlines('10m')
    ax.add_feature(cfeature.STATES.with_scale("10m"))
    ax.add_feature(RMSP, facecolor="white", edgecolor="none", alpha=0.65)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
plt.savefig('mz4_cam_spatial_resolution.pdf', bbox_inches='tight',
            dpi=300)
plt.savefig('mz4_cam_spatial_resolution.png', bbox_inches='tight',
            dpi=300)



# # Opening station data
# aqs = pd.read_csv("cetesb2017_latlon.dat")
# aqs['lon360'] = aqs.lon % 360


# # Creating filter to select aqs inside each cell
# filter_rmsp_cam = ((aqs.lon360 > lon1d_cam_b[0]) & (aqs.lon360 < lon1d_cam_b[-1]) &
#                 (aqs.lat > lat1d_cam_b[0]) & (aqs.lat < lat1d_cam_b[-1]))
# filter_rmsp_mz4 = ((aqs.lon360 > lon1d_mz4_b[0]) & (aqs.lon360 < lon1d_mz4_b[-1]) &
#                 (aqs.lat > lat1d_mz4_b[0]) & (aqs.lat < lat1d_mz4_b[-1]))

# filter_A = ((aqs.lon360 > lon1d_cam_rmsp_b[0]) & (aqs.lon360 < lon1d_cam_rmsp_b[1]) & 
#             (aqs.lat > lat1d_cam_rmsp_b[1]) & (aqs.lat < lat1d_cam_rmsp_b[2]))
# filter_B = ((aqs.lon360 > lon1d_cam_rmsp_b[1]) & (aqs.lon360 < lon1d_cam_rmsp_b[2]) &
#             (aqs.lat > lat1d_cam_rmsp_b[1]) & (aqs.lat < lat1d_cam_rmsp_b[2]))
# filter_C = ((aqs.lon360 > lon1d_cam_rmsp_b[0]) & (aqs.lon360 < lon1d_cam_rmsp_b[1]) &
#             (aqs.lat > lat1d_cam_rmsp_b[0]) & (aqs.lat < lat1d_cam_rmsp_b[1]))
# filter_D = ((aqs.lon360 > lon1d_cam_rmsp_b[1]) & (aqs.lon360 < lon1d_cam_rmsp_b[2]) &
#             (aqs.lat > lat1d_cam_rmsp_b[0]) & (aqs.lat < lat1d_cam_rmsp_b[1]))

# # Filter AQS by cell
# aqsA = aqs[filter_A]
# aqsB = aqs[filter_B]
# aqsC = aqs[filter_C]
# aqsD = aqs[filter_D]
# aqsRMSP = aqs[filter_rmsp_cam]

# print("AQS in CellA: ", len(aqsA.index))
# print("AQS in CellB: ", len(aqsB.index))
# print("AQS in CellC: ", len(aqsC.index))
# print("AQS in CellD: ", len(aqsD.index))
# print("AQS in RMSP: ", len(aqsRMSP.index))

# import matplotlib.ticker as mticker

# fig = plt.figure(figsize=(11,10))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines(resolution='10m')
# ax.add_feature(cfeature.STATES.with_scale('10m'))
# ax.add_feature(RMSP, facecolor="none", edgecolor="black")
# # cam_co_rmsp_mean.plot(ax=ax, vmin=0, vmax=0.45,
# #                       cbar_kwargs={'shrink':0.6, 'label': 'ppb'})
# ax.set_extent([lon_cam_rmsp_b.min() - 0.25, lon_cam_rmsp_b.max() + 0.25,
#                lat_cam_rmsp_b.min() - 0.25, lat_cam_rmsp_b.max() + 0.25])
# # ax.scatter(lon_cam_rmsp_b, lat_cam_rmsp_b, s=100, marker="o", color='orange')
# ax.scatter(aqsD.lon360, aqsD.lat, marker="o", color="blue")
# ax.scatter(aqsB.lon360, aqsB.lat, marker="o", color="red")
# ax.scatter(aqsA.lon360, aqsA.lat, marker="o", color="black")
# # ax.set_yticks(np.concatenate((lat_cam_rmsp[:,0], lat_cam_rmsp_b[:,0])),
# #               crs=ccrs.PlateCarree())
# # ax.set_xticks(((np.concatenate((lon_cam_rmsp[0], lon_cam_rmsp_b[0])) + 180) % 360) - 180,
# #               crs=ccrs.PlateCarree())
# ax.set_ylabel("")
# ax.set_xlabel("")
# ax.set_title("")
# ax.text(lon1d3_cam_rmsp[0], lat1d_cam_rmsp[1], "A",
#         transform=ccrs.PlateCarree(), fontsize=22, horizontalalignment='center')
# ax.text(lon1d3_cam_rmsp[1], lat1d_cam_rmsp[1], "B",
#         transform=ccrs.PlateCarree(), fontsize=22, horizontalalignment='center')
# ax.text(lon1d3_cam_rmsp[0], lat1d_cam_rmsp[0], "C",
#         transform=ccrs.PlateCarree(), fontsize=22, horizontalalignment='center')
# ax.text(lon1d3_cam_rmsp[1], lat1d_cam_rmsp[0], "D",
#         transform=ccrs.PlateCarree(), fontsize=22, horizontalalignment='center')
# gl = ax.gridlines(crs=ccrs.PlateCarree(), color="black")
# gl.xlocator = mticker.FixedLocator(((lon_cam_rmsp_b[0] + 180) % 360 ) -180)
# gl.ylocator = mticker.FixedLocator(lat_cam_rmsp_b[:, 0])
