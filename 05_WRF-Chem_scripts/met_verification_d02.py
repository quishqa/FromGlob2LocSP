#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:52:11 2020

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

wrfout = Dataset(wrf_file)

# Extracting met variables
t2 = wrf.getvar(wrfout, "T2", timeidx=wrf.ALL_TIMES, method="cat")
rh2 = wrf.getvar(wrfout, "rh2", timeidx=wrf.ALL_TIMES, method="cat")
psfc = wrf.getvar(wrfout, "PSFC", timeidx=wrf.ALL_TIMES, method="cat")
wind = wrf.getvar(wrfout, "uvmet10_wspd_wdir", timeidx=wrf.ALL_TIMES,
                  method="cat")
ws = wind.sel(wspd_wdir="wspd")
wd = wind.sel(wspd_wdir="wdir")
u10 = wrf.getvar(wrfout, "U10", timeidx=wrf.ALL_TIMES,
                  method="cat")
v10 = wrf.getvar(wrfout, "V10", timeidx=wrf.ALL_TIMES,
                  method="cat")



# Getting dates to download from wrfout t2 variable
start_date, end_date = dp.qualar_st_end_time(t2)

# Loading List of stations and use only the stations inside wrfout
# cetesb_dom = dp.stations_in_domains("./aqs_htap_vs_a4w.dat", wrfout, t2)
cetesb_aqs = dp.stations_in_domains("./cetesb2017_latlon.dat", wrfout, t2)
cetesb_masp = dp.stations_in_domains("./aqs_met_masp.dat", wrfout, t2)
cetesb_rural = dp.stations_in_domains("./aqs_met_rural.dat", wrfout, t2)


# Wind components
def wd_degree_to_components(cet):
    # Removing wind with ws < 0.5, measurement
    # dominated by turbulence
    df = cet.copy()
    df.loc[df.ws <= 0.5, ["ws", "wd"]] = np.nan, np.nan
    if df.wd.dropna().empty:
        df["u10"] = np.nan
        df["v10"] = np.nan
    else:
        df["u10"] = - df.ws.values * np.sin(np.deg2rad(df.wd.values))
        df["v10"] = - df.ws.values * np.cos(np.deg2rad(df.wd.values))
    return(df)

def wd_components_to_degree(df):
    if df.u10.dropna().empty:
        df["wd2"] = np.nan
    else:
        r2d = 45.0 /np.arctan(1.0)
        df["wd2"] = np.arctan2(df.u10.values, df.v10.values) * r2d + 180
    return(df)


# Loading CETESB Meteorology
met_filename = "cetesb_met_2014.pkl"
with open(met_filename, 'rb') as a_dict:
    cetesb_met = pickle.load(a_dict)

# Adding AQS name in cetesb_pol dataframe
# And changing T in C to K
## TODO: add AQS name in qualR_py
for aqs in list(cetesb_met):
    cetesb_met[aqs]['name'] = aqs
    cetesb_met[aqs]['t2'] = cetesb_met[aqs].t2 + 273.15

# Calculating wind components for CETESB AQS
cetesb_met = {aqs:wd_degree_to_components(cetesb_met[aqs]) 
              for aqs in list(cetesb_met)}

# Extracting AQS from WRF's wrfout from all domain
wrf_dom = dp.cetesb_from_wrf(cetesb_aqs,
                             (t2, rh2, wind, u10, v10),
                             to_local=True)
wrf_masp = {aqs:wrf_dom[aqs] for aqs in cetesb_masp.name}
wrf_rural = {aqs:wrf_dom[aqs] for aqs in cetesb_rural.name}

# Extracting from cetesb_met, AQS for each case (dom, masp, and rural)
cetesb_dom = {aqs:cetesb_met[aqs] for aqs in cetesb_aqs.name}
cetesb_dom_masp = {aqs:cetesb_met[aqs] for aqs in cetesb_masp.name}
cetesb_dom_rural = {aqs:cetesb_met[aqs] for aqs in cetesb_rural.name}

# Creating model and observation pairs
# Removing first three days as spin-up time
WRF, CET = dp.model_eval_setup(wrf_dom, cetesb_dom,
                               date_start='2014-10-06')
WRF_masp, CET_masp = dp.model_eval_setup(wrf_masp, cetesb_dom_masp,
                                         date_start='2014-10-06')
WRF_rural, CET_rural = dp.model_eval_setup(wrf_rural, cetesb_dom_rural,
                                           date_start='2014-10-06')

# Function to plot
# wind plots



def plot_wind_vector(wrf, cet, ax=None, hide=False,
                     legend=True, title=False):
    cet_uv = wd_components_to_degree(cet)
    if ax is None:
        ax = plt.gca()
    formater = mdates.DateFormatter("%d")
    q = ax.quiver(cet.index, 0, cet_uv.u10, cet_uv.v10, 
              color="black", label="CETESB")
    ax.quiver(wrf.index, 0, wrf.u10, wrf.v10,
              color="red", label="WRF-Chem", alpha=0.65)
    ax.quiverkey(q, 0.15, 0.1, 5, r'$5 \frac{m}{s}$', labelpos='E',
             coordinates='axes')
    ax.yaxis.set_ticks([])
    ax.set_ylabel("Scale wind vector")
    if legend:
        ax.legend(frameon=False)
    if title:
        ax.set_title(wrf.name.unique()[0])
    if hide:
        ax.set_xticklabels(())
    else:
        ax.set_xlabel("October 2014")
        ax.xaxis.set_major_formatter(formater)

def hourly_mean_plot_err_bar(wrf, cet, var, ylab, ylim, ax=None,
                             save_fig=False, fmt=None, hide=False,
                             legend=True, title=False):
    wrf_var = wrf[var]
    cet_var = cet[var]
    
    df = pd.concat([wrf_var, cet_var], axis=1,
                   keys=["wrf", "cet"])
    df_wrf_h = (df.groupby(df.index.hour)
                .agg({"mean", "std"})
                ['wrf']
                .add_prefix("wrf_"))
    
    df_cet_h = (df.groupby(df.index.hour)
                .agg({"mean", "std"})
                ['cet']
                .add_prefix("cet_"))
    if cet_var.dropna().empty:
        df_cet_h["cet_mean"] = np.nan
        df_cet_h["cet_std"] = np.nan
    if ax is None:
        ax = plt.gca()
    
    ax.errorbar(df_cet_h.index, df_cet_h.cet_mean, yerr=df_cet_h.cet_std,
                color="black", label="CETESB.", fmt="-o", markersize=4)
    ax.errorbar(df_wrf_h.index, df_wrf_h.wrf_mean, yerr=df_wrf_h.wrf_std,
                color="red", label="WRF-Chem", fmt="-o", markersize=4)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylab)
    if hide:
        ax.set_xlabel("")
    else:
        ax.set_xlabel("Hour of day (LT)")
    if legend:
        ax.legend(frameon=False)
    if title:
        ax.set_title(wrf.name.unique()[0])


def plot_wrf_vs_cetesb(wrf, cet, var, ylab, ylim, ax=None,
                       save_fig=False, fmt=None, hide=False,
                       legend=True, title=False):
    if ax is None:
        ax =plt.gca()
    formater = mdates.DateFormatter("%d")
    ax.plot(cet[var], color="black", label='CETESB', linewidth=2.5)
    ax.plot(wrf[var], color="red", linewidth=2.5,
            label="WRF-Chem")
    ax.set_ylim(ylim)
    ax.set_ylabel(ylab)
    if hide:
        ax.set_xticklabels(())
    else:
        ax.set_xlabel("October 2014")
        ax.xaxis.set_major_formatter(formater)
    if legend:
        ax.legend(frameon=False)
    if title:
        ax.set_title(wrf.name.unique()[0])
    if save_fig:
        file_name = ("wrf_" + var + '_'+ wrf.name.unique()[0] 
                      + fmt)
        plt.savefig(file_name, bbox_inches="tight", dpi=300)
        plt.clf()


def aqs_with_four_pollutants(WRF, CET, filename, ylims=None):
    aqs = list(WRF)
    pol = ['t2', 'rh2', 'ws', 'co']
    ylabs = ["$O_3 \; (\mu g / m^3)$",
              "$NO \; (\mu g / m^3)$",
              "$NO_2\; (\mu g / m^3)$",
              "$CO \; (ppm)$"]
    if ylims==None:
        ylims = [ (0, 150), (0, 500), (0, 150), (0, 8)]
    
    
    fig, ax = plt.subplots(4, 4, figsize=(22, 10))
    for i in  range(4):
        for j in range(4):
            if i == 0:
                hourly_mean_plot_err_bar(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                 ylims[i], ax=ax[i, j], title=True, hide=True)
            elif i<= 2:
                hourly_mean_plot_err_bar(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                 ylims[i], ax=ax[i, j], legend=False, hide=True)
            else:
                hourly_mean_plot_err_bar(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                 ylims[i], ax=ax[i, j], legend=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")




def wrf_vs_cetesb_four_parameters(WRF, CET, filename, ylims=None):
    aqs = list(WRF)
    pol = ['t2', 'rh2', 'ws', 'wd']
    ylabs = ["T2 $(K)$",
              "RH2 $(\%)$",
              "WS10 $(m \;s^{-1})$",
              "Scaled wind vector"]
    if ylims==None:
        ylims = [ (270, 320), (0, 100), (0, 10)]    
    
    fig, ax = plt.subplots(len(pol), len(aqs), figsize=(22, 10))
    for i, _ in  enumerate(pol):
        for j, _ in enumerate(aqs):
            if i == 0:
                plot_wrf_vs_cetesb(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                 ylims[i], ax=ax[i, j], legend=True, 
                                 hide=True, title=True)
            elif i<= 2:
                plot_wrf_vs_cetesb(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                   ylims[i], ax=ax[i, j], legend=False, 
                                   hide=True)
            else:
                plot_wind_vector(WRF[aqs[j]], CET[aqs[j]], 
                                 ax=ax[i, j], legend=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    
    
wrf_vs_cetesb_four_parameters(WRF_masp, CET_masp, "met_masp_d02_megan.pdf")
# wrf_vs_cetesb_four_parameters(WRF_masp, CET_masp, "met_masp_d01_megan.pdf")
# wrf_vs_cetesb_four_parameters(WRF_rural, CET_rural, "met_rural_d01_megan.pdf")



met_global_eval = ms.global_stat_some_vars(WRF, CET, 
                                            ["t2", "rh2", "ws", "wd", "u10", "v10"])

print(met_global_eval.to_latex(float_format="%.2f",
                               na_rep="-"))


# pin_wrf = WRF_masp["Pinheiros"]
# pin_cet = CET_masp["Pinheiros"]


# def plot_wind_vector2(wrf, cet, ax=None, hide=False,
#                      legend=True, title=False):
#     cet_uv = wd_components_to_degree(cet)
#     if ax is None:
#         ax = plt.gca()
#     formater = mdates.DateFormatter("%d")
#     q = ax.quiver(cet.index, 0, cet_uv.u10, cet_uv.v10, 
#               color="black", label="CETESB")
#     ax.quiver(wrf.index, 0, wrf.u10, wrf.v10,
#               color="red", label="WRF-Chem", alpha=0.65)
#     ax.quiverkey(q, 0.15, 0.05, 5, r'$5 \frac{m}{s}$', labelpos='E',
#              coordinates='axes')
#     ax.yaxis.set_ticks([])
#     ax.set_ylabel("Scale wind vector")
#     if legend:
#         ax.legend(frameon=False)
#     if title:
#         ax.set_title(wrf.name.unique()[0])
#     if hide:
#         ax.set_xticklabels(())
#     else:
#         ax.set_xlabel("October 2014")
#         ax.xaxis.set_major_formatter(formater)

# fig, ax = plt.subplots()
# plot_wind_vector2(pin_wrf, pin_cet, ax=ax)
# fig, axes = plt.subplots(4, 1, figsize=(5, 15))
# plot_wrf_vs_cetesb(pin_wrf, pin_cet, "t2", "T2", (270, 310), ax=axes[0])
# plot_wrf_vs_cetesb(pin_wrf, pin_cet, "rh2", "RH2", (0, 100), ax=axes[1])
# plot_wrf_vs_cetesb(pin_wrf, pin_cet, "ws", "$m s^{-1}$", (0, 10), ax=axes[2])
# plot_wind_vector(pin_wrf, pin_cet, ax=axes[3])


