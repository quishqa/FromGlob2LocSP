#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:39:04 2020

@author: mgavidia
"""

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import cartopy.mpl.geoaxes
from shapely.geometry.polygon import LinearRing
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

wrf_sp_d01 = xr.open_dataset('./wrfinput_d01')
wrf_sp_d02 = xr.open_dataset('./wrfinput_d02')

# Loading coordinates
xlat1 = wrf_sp_d01.XLAT[0, :, :]
xlon1 = wrf_sp_d01.XLONG[0, :, :]
xlat2 = wrf_sp_d02.XLAT[0, :, :]
xlon2 = wrf_sp_d02.XLONG[0, : ,:]

# Defining d01 domains
d01_lon = [xlon1.min().values, xlon1.min().values,
           xlon1.max().values, xlon1.max().values]
d01_lat = [xlat1.min().values, xlat1.max().values,
           xlat1.max().values, xlat1.min().values]
d01_sqr = LinearRing(list(zip(d01_lon, d01_lat)))


# Defining d02 domains
d02_lon = [xlon2.min().values, xlon2.min().values,
           xlon2.max().values, xlon2.max().values]
d02_lat = [xlat2.min().values, xlat2.max().values,
           xlat2.max().values, xlat2.min().values]
d02_sqr = LinearRing(list(zip(d02_lon, d02_lat)))

# Loading shapefile
reader = shpreader.Reader('masp_shp/masp_shp.shp')
muns = list(reader.geometries())
MUNS = cfeature.ShapelyFeature(muns,ccrs.PlateCarree())

xticks = np.round(np.arange(xlon1.min(), xlon1.max(), 1.5), 2)[1:]
yticks = np.round(np.arange(xlat1.min(), xlat1.max(), 1.5), 2)[1:]

fig = plt.figure(figsize=(16, 9))
ax = plt.axes(projection = ccrs.PlateCarree())
cbax = ax.contourf(xlon1, xlat1, wrf_sp_d01.HGT.isel(Time=0),
                   transform=ccrs.PlateCarree(),
                   cmap="terrain")
ax.add_feature(cfeature.STATES.with_scale('10m'))
ax.coastlines('10m', zorder=15)
ax.add_feature(MUNS, facecolor='none', edgecolor='black', alpha=0.5)
ax.add_feature(cfeature.OCEAN, zorder=10)
ax.add_geometries([d02_sqr], ccrs.PlateCarree(), facecolor='none', 
                  edgecolor="black", linewidth=1.5, zorder=20)
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.set_yticklabels(yticks, fontsize=15)
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_xticklabels(xticks, fontsize=15)
ax.text(xlon1.min(), xlat1.max() + 0.075, "D01", fontsize=25, 
        transform=ccrs.PlateCarree())
ax.text(xlon2.min(), xlat2.max() + 0.075, "D02", fontsize=25,
        transform=ccrs.PlateCarree())
axins = inset_axes(ax, width="40%", height="40%", loc="lower right", 
                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(
                       map_projection=cartopy.crs.Orthographic(
                           central_longitude=-50.0)
                       ))
axins.coastlines()
axins.add_feature(cfeature.OCEAN)
axins.add_feature(cfeature.LAND)
axins.add_feature(cfeature.BORDERS)
axins.add_geometries([d01_sqr], ccrs.PlateCarree(),facecolor='red',
                     edgecolor='red', alpha=0.5)
cb = plt.colorbar(cbax, ax=ax, label = 'HGT (m)', shrink=0.7)
cb.ax.tick_params(labelsize=15)
cb.set_label(label="HGT (m)", fontsize=15)
plt.tick_params(labelsize=20)
plt.savefig('wrf_sp_domains.png', dpi=300,  bbox_inches="tight")
plt.savefig('wrf_sp_domains.pdf', dpi=300,  bbox_inches="tight")



