#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:50:50 2020

@author: mgavidia
"""

import wrf as wrf
from netCDF4 import Dataset
import wrf_sp_eval.data_preparation as dp
import wrf_sp_eval.qualar_py as qr
import wrf_sp_eval.model_stats as ms
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
import pickle

wrf_file = "/scr2/mgavidia/wrf4/WRF/test/em_real/wrfout_d02_2014-10-04_00:00:00_cbc_megan"
# wrf_file = "/scr2/mgavidia/WRF4/WRF/test/em_real/wrfout_d02_2014-10-03_00:00:00_cbmz"
# wrf_file = "/scr2/mgavidia/wrf4/WRF/test/em_real/wrfout_d02_2014-10-03_00:00:00_no_cbc"

wrfout = Dataset(wrf_file)

# Extracting met variables
t2 = wrf.getvar(wrfout, "T2", timeidx=wrf.ALL_TIMES, method="cat")
rh2 = wrf.getvar(wrfout, "rh2", timeidx=wrf.ALL_TIMES, method="cat")
psfc = wrf.getvar(wrfout, "PSFC", timeidx=wrf.ALL_TIMES, method="cat")


# Extracting pollutants variables a4w
o3_wrf= wrf.getvar(wrfout, "o3", timeidx=wrf.ALL_TIMES, method="cat")
co_wrf = wrf.getvar(wrfout, "co", timeidx=wrf.ALL_TIMES, method="cat")
no_wrf = wrf.getvar(wrfout, "no", timeidx=wrf.ALL_TIMES, method="cat")
no2_wrf = wrf.getvar(wrfout, "no2", timeidx=wrf.ALL_TIMES, method="cat")



# Retrieving pollutants from surface
o3_sfc = o3_wrf.isel(bottom_top=0)
co_sfc = co_wrf.isel(bottom_top=0)
no_sfc = no_wrf.isel(bottom_top=0)
no2_sfc = no2_wrf.isel(bottom_top=0)



# Transform surface polutant from ppm to ug/m3
o3_u_wrf = dp.ppm_to_ugm3(o3_sfc, t2, psfc, 48)
no_u_wrf = dp.ppm_to_ugm3(no_sfc, t2, psfc, 30)
no2_u_wrf = dp.ppm_to_ugm3(no2_sfc, t2, psfc, 46)
nox_u_wrf = (no_sfc + no2_sfc) * 1000
nox_u_wrf = nox_u_wrf.rename('nox')


# Getting dates to download from wrfout t2 variable
start_date, end_date = dp.qualar_st_end_time(t2)

# Loading List of stations and use only the stations inside wrfout
# cetesb_dom = dp.stations_in_domains("./aqs_htap_vs_a4w.dat", wrfout, t2)
cetesb_dom = dp.stations_in_domains("./cetesb2017_latlon.dat", wrfout, t2)

# Loading CETESB Polutants 
a_dict = open("cetesb_pols_2014_all_pol.pkl", "rb")
cetesb_pol = pickle.load(a_dict)
a_dict.close()


wrf_pol =  dp.cetesb_from_wrf(cetesb_dom, 
                              (o3_u_wrf, no_u_wrf, no2_u_wrf, 
                               nox_u_wrf, co_sfc),
                               to_local=True)


cetesb_pol_dom = {aqs:cetesb_pol[aqs] for aqs in list(wrf_pol)}


for aqs in list(cetesb_pol_dom):
    cetesb_pol_dom[aqs]['name'] = aqs
 

WRF, CET =dp.model_eval_setup(wrf_pol, cetesb_pol_dom, date_start='2014-10-06')


wrf_df = pd.concat(WRF)
cet_df = pd.concat(CET)

pol_params = ["o3", "no", "no2", "co"]


# for param in met_params:
#     cet_df[param] = pd.to_numeric(cet_df[param])


def cetesb_aqs_pol_mean_df(CET, params, cetesb_dom):
    '''
    Calculates pollutant mean for all CETESB AQS measurements

    Parameters
    ----------
    CET : dict
        Contains a dataframe per aqs with hourly pollutant
        measurements.
    params : list
        list of columns names to transform to numeric.
    cetesb_dom : pandas dataframe
        data frame with AQS lat and lon.

    Returns
    -------
    cet_df_mean : pandas DataFrame
        pollutants means per aqs.

    '''
    cet_df = pd.concat(CET)
    for param in params:
        cet_df[param] = pd.to_numeric(cet_df[param])
        
    cet_df_mean = (cet_df.groupby('name')
                   .mean()
                   .join(cetesb_dom.set_index("name")))
    
    return (cet_df_mean)


def filter_pol_mean_aqs(cet_df_mean, param):
    """
    Select column (param) from mean aqs data frame
    and remove aqs with no measurements

    Parameters
    ----------
    cet_df_mean : pandas Data Frame
        Pollutnas means per aqs.
    param : str
        Pollutant abreviation.

    Returns
    -------
    param_mean_no_nan : pandas data frame
        Aqs with mean values.

    """
    param_mean_no_nan = (cet_df_mean
                       [[param, "lon", "lat"]]
                       .dropna()
                       .copy())
    return (param_mean_no_nan)
    


def get_wrf_mean_cet_mean_range(cet_aqs, wrf_mean, param):
    '''
    Retrieve min and max concetration value
    from cetesb mean measurements and wrf mean concentrations


    Parameters
    ----------
    cet_aqs : pandas DataFrame
        AQS with mean concentration.
    wrf_mean : xarray dataarray
        WRF mean values.
    param : str
        Parameter (column) name.

    Returns
    -------
    list
        [min, max].

    '''
    cet_range = [cet_aqs[param].min(), cet_aqs[param].max()]
    wrf_range = [wrf_mean.min().values.tolist(),
                 wrf_mean.max().values.tolist()]
    cet_wrf = cet_range + wrf_range
    min_range = np.floor(min(cet_wrf))
    max_range = np.ceil(max(cet_wrf))
    return ([min_range, max_range])



wrf_local_time = (rh2.Time
                  .to_index()
                  .tz_localize("UTC")
                  .tz_convert("America/Sao_Paulo"))

# rh2.Time.to_series().tz_localize("UTC").tz_convert("America/Sao_Paulo")


cet_df_mean = cetesb_aqs_pol_mean_df(CET, pol_params, cetesb_dom)
no_aqs = filter_pol_mean_aqs(cet_df_mean, "no")

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

#loading MASP shapefile
reader_masp=shpreader.Reader("masp_shp/masp_shp.shp")
masp = list(reader_masp.geometries())
MASP = cfeature.ShapelyFeature(masp, ccrs.PlateCarree())





def plot_spatial_wrf_cet_points(wrf_pol,cet_df_mean, pol, label_name,
                                ax=None):
    wrf_mean = (wrf_pol
                .sel(Time=slice("2014-10-06 03:00",
                                "2014-10-13 00:00"))
                .mean(dim="Time"))
    cet_aqs = filter_pol_mean_aqs(cet_df_mean, pol)
    
    pol_min, pol_max = get_wrf_mean_cet_mean_range(cet_aqs,
                                                   wrf_mean,
                                                   pol)
    norm = matplotlib.colors.Normalize(0, pol_max)
    if ax is None:
        ax = plt.gca()
        
    
    cs = ax.contourf(wrf_mean.XLONG, wrf_mean.XLAT, wrf_mean.values,
                 cmap="magma_r", norm=norm)
    ax.scatter(cet_aqs.lon.values, cet_aqs.lat.values, c=cet_aqs[pol].values,
               marker="o", edgecolor="k", linewidth=0.5, norm=norm,
               transform=ccrs.PlateCarree(), cmap="magma_r", zorder=100)
    ax.add_feature(MASP, facecolor="none", linewidth=0.5, edgecolor="k")
    ax.add_feature(cfeature.STATES.with_scale("10m"))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%",
                              axes_class=plt.Axes)
    cb= ax.get_figure().colorbar(cs, cax=cax, orientation="vertical")
    cb.set_label(label_name)
    



def filter_clean_stats_per_pol(wrf_stats, cetesb_dom, pol, stats):
    '''
    Create a data frame for a statistic indicator for a selected
    pollutant with the lat and lon of the AQS.

    Parameters
    ----------
    wrf_stats : pandas DataFrame
        DataFrame with the statistic indicator of a 
        selected pollutant per aqs.
    cetesb_dom : pandas DataFrame
        Cetesb aqs in domain.
    pol : str
        Pollutant name.
    stats : str
        atbreviation of statistic indicator.

    Returns
    -------
    wrf_stats_pol : pandas data frame
        statistic indicator per aqs for a selected
        pollutant.

    '''
    cetesb_aqs = cetesb_dom.copy().rename(columns={"name": "aqs"})
    wrf_aqs = (wrf_stats
               .loc[pol]
               .copy().
               merge(cetesb_aqs, on=["aqs"])
               [[stats, "lon", "lat"]]
               .dropna())
    return wrf_aqs
    


def plot_stats_cet_points(wrf_stats, cetesb_dom, pol, stats,
                          label_name, ax=None):
    pol_stats = filter_clean_stats_per_pol(wrf_stats, cetesb_dom,
                                           pol, stats)
    if stats == "R":
        max_val = 1.0
    else:
        max_val = np.ceil(pol_stats[stats].abs().max())
    
    if ax is None:
        ax = plt.gca()
    divider = make_axes_locatable(ax)
    ax.set_extent([t2.XLONG.min(), t2.XLONG.max(),
                   t2.XLAT.min(), t2.XLAT.max()])
    cs = ax.scatter(pol_stats.lon.values, pol_stats.lat.values,
                    c=pol_stats[stats].values, marker='o', edgecolor='k',
                    transform=ccrs.PlateCarree(), 
                    cmap="coolwarm",
                    linewidths=0.5,
                    zorder=500, 
                    vmin=-max_val, vmax=max_val)
    ax.add_feature(MASP, facecolor='none', linewidth=0.25, edgecolor='k')
    ax.add_feature(cfeature.STATES.with_scale('10m'))
    ax.coastlines("10m")
    cax = divider.append_axes("right", "5%", pad="3%", 
                              axes_class=plt.Axes)
    cb = ax.get_figure().colorbar(cs, cax=cax, orientation="vertical")
    cb.set_label(label_name, labelpad=-1)

wrf_stats = ms.all_aqs_some_vars(WRF, CET,
                                 pol_params)
o3_R = filter_clean_stats_per_pol(wrf_stats, cetesb_dom, 
                                  "o3", "R")
o3_MB = filter_clean_stats_per_pol(wrf_stats, cetesb_dom, 
                                  "o3", "MB")


# For d02 (20, 15)
fig, axes = plt.subplots(3, 4, figsize=(20, 15),
                          subplot_kw={'projection': ccrs.PlateCarree()})
plot_spatial_wrf_cet_points(o3_u_wrf, cet_df_mean, "o3", "$\mu g \; m^{-3}$",
                            ax=axes[0,0])
axes[0, 0].set_title("Ozone")
plot_spatial_wrf_cet_points(no_u_wrf, cet_df_mean, "no", "$\mu g \; m^{-3}$",
                            ax=axes[0,1])
axes[0, 1].set_title("Nitric oxide")
plot_spatial_wrf_cet_points(no2_u_wrf, cet_df_mean, "no2", "$\mu g \; m^{-3}$", 
                            ax=axes[0,2])
axes[0, 2].set_title("Nitrogen dioxide")
plot_spatial_wrf_cet_points(co_sfc, cet_df_mean, "co", "$ppm$", 
                            ax=axes[0,3])
axes[0, 3].set_title("Carbon monoxide")

plot_stats_cet_points(wrf_stats, cetesb_dom, "o3", "R", "", ax=axes[1, 0])
axes[1, 0].set_title("R for ozone")
plot_stats_cet_points(wrf_stats, cetesb_dom, "no", "R","", ax=axes[1, 1])
axes[1, 1].set_title("R for nitric oxide")
plot_stats_cet_points(wrf_stats, cetesb_dom, "no2", "R","",  ax=axes[1, 2])
axes[1, 2].set_title("R for nitrogen dioxide")
plot_stats_cet_points(wrf_stats, cetesb_dom, "co", "R", "",     ax=axes[1, 3])
axes[1, 3].set_title("R for carbon monoxide")

plot_stats_cet_points(wrf_stats, cetesb_dom, "o3", "MB", "$\mu g \; m^{-3}$", 
                      ax=axes[2, 0])
axes[2, 0].set_title("MB for ozone")
plot_stats_cet_points(wrf_stats, cetesb_dom, "no", "MB", "$\mu g \; m^{-3}$", 
                      ax=axes[2, 1])
axes[2, 1].set_title("MB for nitric oxide")
plot_stats_cet_points(wrf_stats, cetesb_dom, "no2", "MB", "$\mu g \; m^{-3}$",
                      ax=axes[2, 2])
axes[2, 2].set_title("MB for nitrogen dioxide")
plot_stats_cet_points(wrf_stats, cetesb_dom, "co", "MB", "$ppm$",
                      ax=axes[2, 3])
axes[2, 3].set_title("MB for carbon monoxide")
plt.subplots_adjust(wspace=0.25, 
                    hspace=-0.6) #0.6 for d02 and 0.7 for d01
plt.savefig("d02_spatial_pol_megan.pdf", dpi=300, bbox_inches="tight")