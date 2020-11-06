#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 23:00:56 2020

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

# Domain 1: here we can compare "rural" vs urban
wrf_file = "/scr2/mgavidia/wrf4/WRF/test/em_real/wrfout_d02_2014-10-04_00:00:00_cbc_megan"

# Loading wrfout
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
nox_u_wrf = (no_sfc + no2_sfc) * 1000 # oservation are in ppb
nox_u_wrf = nox_u_wrf.rename('nox')

# Getting dates to download from wrfout t2 variable
start_date, end_date = dp.qualar_st_end_time(t2)

# Loading List of stations and use only the stations inside wrfout
cetesb_dom = dp.stations_in_domains("./cetesb2017_latlon.dat", wrfout, t2)
cetesb_masp = dp.stations_in_domains("./aqs_masp.dat", wrfout, t2)
cetesb_rural = dp.stations_in_domains("./aqs_rural.dat", wrfout, t2)

# Loading CETESB Polutants 
## Previously downloaded using qualR_py
## cetesb_pol is a dictionary with AQS name as keys
## and pollutant measurements as values.
a_dict = open("cetesb_pols_2014_all_pol.pkl", "rb")
cetesb_pol = pickle.load(a_dict)
a_dict.close()

# Adding AQS name in cetesb_pol dataframe
## TODO: add AQS name in qualR_py
for aqs in list(cetesb_pol):
    cetesb_pol[aqs]['name'] = aqs
    
# Extracting AQS from WRF's wrfout
# From all domain
wrf_dom =  dp.cetesb_from_wrf(cetesb_dom, 
                              (o3_u_wrf, no_u_wrf, no2_u_wrf, 
                               nox_u_wrf, co_sfc),
                               to_local=True)
wrf_masp = dp.cetesb_from_wrf(cetesb_masp, 
                              (o3_u_wrf, no_u_wrf, no2_u_wrf, 
                               nox_u_wrf, co_sfc),
                               to_local=True)
wrf_rural = dp.cetesb_from_wrf(cetesb_rural, 
                              (o3_u_wrf, no_u_wrf, no2_u_wrf, 
                               nox_u_wrf, co_sfc),
                               to_local=True)

# Extraction from cetesb_pol, AQS inside domains
cetesb_dom = {aqs:cetesb_pol[aqs] for aqs in list(wrf_dom)}
cetesb_dom_masp = {aqs:cetesb_pol[aqs] for aqs in list(wrf_masp)}
cetesb_dom_rural = {aqs:cetesb_pol[aqs] for aqs in list(wrf_rural)}

# Creating model and observation pairs
# Removing first three days as spin-up time
WRF, CET = dp.model_eval_setup(wrf_dom, cetesb_dom, 
                               date_start='2014-10-06')
WRF_masp, CET_masp = dp.model_eval_setup(wrf_masp, cetesb_dom_masp,
                                         date_start='2014-10-06')
WRF_rural, CET_rural = dp.model_eval_setup(wrf_rural, cetesb_dom_rural,
                                           date_start='2014-10-06')


# Function to plot
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
    pol = ['o3', 'no', 'no2', 'co']
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


def aqs_with_four_pollutants_rural(WRF, CET, filename, ylims=None):
    aqs = list(WRF)
    pol = ['o3', 'no', 'no2']
    ylabs = ["$O_3 \; (\mu g / m^3)$",
              "$NO \; (\mu g / m^3)$",
              "$NO_2\; (\mu g / m^3)$"]
    if ylims==None:
        ylims = [ (0, 150), (0, 500), (0, 150), (0, 8)]
    
    
    fig, ax = plt.subplots(3, 4, figsize=(22, 10))
    for i in  range(3):
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





def wrf_vs_cetesb_four_pollutants(WRF, CET, filename, ylims=None):
    aqs = list(WRF)
    pol = ['o3', 'no', 'no2', 'co']
    ylabs = ["$O_3 \; (\mu g / m^3)$",
              "$NO \; (\mu g / m^3)$",
              "$NO_2\; (\mu g / m^3)$",
              "$CO \; (ppm)$"]
    if ylims==None:
        ylims = [ (0, 200), (0, 450), (0, 150), (0, 4)]    
    
    fig, ax = plt.subplots(len(pol), len(aqs), figsize=(22, 10))
    for i in  range(4):
        for j in range(4):
            if i == 0:
                plot_wrf_vs_cetesb(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                 ylims[i], ax=ax[i, j], legend=True, 
                                 hide=True, title=True)
            elif i<= 2:
                plot_wrf_vs_cetesb(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                   ylims[i], ax=ax[i, j], legend=False, 
                                   hide=True)
            else:
                plot_wrf_vs_cetesb(WRF[aqs[j]], CET[aqs[j]], pol[i], ylabs[i], 
                                   ylims[i], ax=ax[i, j], legend=False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")


masp_ylims = [(0, 200), (0, 250), (0, 150), (0, 2.5)]
rural_ylims = [(0, 200), (0, 100), (0, 100), (0, 2.5)]

aqs_with_four_pollutants(WRF_masp, CET_masp, 'pol_masp_d02_megan.pdf',
                         ylims=masp_ylims)

masp_ylims = [(0, 250), (0, 300), (0, 150), (0, 4)]
wrf_vs_cetesb_four_pollutants(WRF_masp, CET_masp, 'pol_masp_ts_d02_megan.pdf',
                         ylims=masp_ylims)





pol_global_eval = ms.global_stat_some_vars(WRF, CET, 
                                           ["o3", "no", "no2", "nox", "co"])

print(pol_global_eval.to_latex(float_format="%.2f",
                               na_rep="-"))
aqs_stats = ms.all_aqs_some_vars(WRF, CET, ["o3", "no", "no2", "nox", "co"])




aqs_stats.loc["o3"].dropna()
