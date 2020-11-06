#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 23:55:11 2020

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

# Selecting SA
cam_sa = cam.sel(lat=slice(-60, 20),
                 lon=slice(-90 % 360, -30 % 360))

# CO mean
# Some examples
CO_mean_col = (cam_sa
               .CO
               .mean(dim="time")
               .mean(dim="lev"))


O3_mean = (cam_sa
           .O3
           .sel(lev=1000, method="nearest")
           .mean(dim="time"))





def molmol2ugm3(ds, pol, M):
    '''
    Transform mol/mol to ug/m3

    Parameters
    ----------
    ds : xarray dataset
        CAM-Chem output as xr dataset.
    pol : str
        Pollutant name.
    M : float
        Molecular weight.

    Returns
    -------
    pol_m : xarray dataarray
        Pollutant as ug/m3.

    '''
    pol_sfc = (ds[pol].sel(lev=1000, method="nearest"))
    t_sfc = (ds["T"].sel(lev=1000, method="nearest"))
    pol_m = (pol_sfc * 10**6 * 
             pol_sfc.lev.values * 10**2 * 
             48 / 8.3142 / t_sfc.values).mean(dim="time")
    return pol_m


o3_m = molmol2ugm3(cam_sa, "O3", 48)
no_m = molmol2ugm3(cam_sa, "NO", 30)
no2_m = molmol2ugm3(cam_sa, "NO2", 46)
co_m =  (cam_sa
            .CO
            .sel(lev=1000, method="nearest")
            .mean(dim="time")) * 10 ** 6 # to ppm
pm25_m = (cam_sa
             .PM25
             .sel(lev=1000, method="nearest")
             .mean(dim="time")) * 10**9 # to ug/m3



def plot_pol_SA_sfc(da, lab, title, ax = None):
    '''
    Plot mean surface concentration 

    Parameters
    ----------
    da : xarray dataarray
        Pollutant mean sfc.
    lab : str
        Colorbar label.
    title : str
        Plot title.
    ax : plt axes, optional
        ax to plot. The default is None.

    Returns
    -------
    None.

    '''
    if ax is None:
        ax = plt.gca()
        
    cbar_dict = {'shrink': 0.45,
                 'label': lab}
    da.plot(ax=ax,
            cbar_kwargs=cbar_dict,
            cmap="inferno_r")
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS.with_scale("10m"))
    ax.set_extent([cam_sa.lon.min(), cam_sa.lon.max(),
                   cam_sa.lat.min() + 1.8, cam_sa.lat.max() - 1.5])
    ax.set_title(title)
    
    

fig, axes = plt.subplots(3, 2, figsize=[6, 10],
                        subplot_kw={'projection': ccrs.PlateCarree()})
plot_pol_SA_sfc(o3_m, "$O_3 \; (\mu  g \; m^{-3})$",
                "Ozone", ax=axes[0, 0])
plot_pol_SA_sfc(no_m, "$NO \; (\mu\; g m^{-3})$",
                "Nitric oxide", ax=axes[0, 1])
plot_pol_SA_sfc(no2_m, "$NO_2 \; (\mu g \;m^{-3})$",
                "Nitrogen dioxide", ax=axes[1, 0])
plot_pol_SA_sfc(co_m, "$CO \; (ppm)$", 
                "Carbon monoxide", ax=axes[1, 1])
plot_pol_SA_sfc(pm25_m, "$PM_{2.5} \; (\mu g \; m^{-3})$",
                "Fine particles", ax=axes[2, 0])
axes[2, 1].axis("off")
plt.savefig("cam_SA_megan.pdf", dpi=300, bbox_inches="tight")


# Retriving data for são paulo

sp_lat, sp_lon = [-24.031, -46.250]

o3_sp = o3_m.sel(lat=sp_lat, lon=sp_lon % 360, method="nearest").values
no_sp = no_m.sel(lat=sp_lat, lon=sp_lon % 360, method="nearest").values
no2_sp = no2_m.sel(lat=sp_lat, lon=sp_lon % 360, method="nearest").values
co_sp = co_m.sel(lat=sp_lat, lon=sp_lon % 360, method="nearest").values
pm25_sp = pm25_m.sel(lat=sp_lat, lon=sp_lon % 360, method="nearest").values

print("O3: ", o3_sp)
print("NO: ", no_sp)
print("NO2: ", no2_sp)
print("CO: ", co_sp)
print("PM25: ", pm25_sp)
