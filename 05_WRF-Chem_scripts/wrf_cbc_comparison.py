#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 20:28:08 2020

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


wrf_file = "/scr2/mgavidia/WRF4/WRF/test/em_real/wrfout_d01_2014-10-04_00:00:00_megan"
wrf_file_cbc = "/scr2/mgavidia/wrf4/WRF/test/em_real/wrfout_d01_2014-10-04_00:00:00_cbc_megan"


wrfout = Dataset(wrf_file)
wrfout_cbc = Dataset(wrf_file_cbc)


# Extracting met variables
t2 = wrf.getvar(wrfout, "T2", timeidx=wrf.ALL_TIMES, method="cat")
rh2 = wrf.getvar(wrfout, "rh2", timeidx=wrf.ALL_TIMES, method="cat")
psfc = wrf.getvar(wrfout, "PSFC", timeidx=wrf.ALL_TIMES, method="cat")

t2_cbc = wrf.getvar(wrfout_cbc, "T2", timeidx=wrf.ALL_TIMES, method="cat")
rh2_cbc = wrf.getvar(wrfout_cbc, "rh2", timeidx=wrf.ALL_TIMES, method="cat")
psfc_cbc = wrf.getvar(wrfout_cbc, "PSFC", timeidx=wrf.ALL_TIMES, method="cat")


# Extracting pollutants variables wrf no cbc
o3_wrf= wrf.getvar(wrfout, "o3", timeidx=wrf.ALL_TIMES, method="cat")
co_wrf = wrf.getvar(wrfout, "co", timeidx=wrf.ALL_TIMES, method="cat")
no_wrf = wrf.getvar(wrfout, "no", timeidx=wrf.ALL_TIMES, method="cat")
no2_wrf = wrf.getvar(wrfout, "no2", timeidx=wrf.ALL_TIMES, method="cat")

# Extracting pollutants variables wrf cbc
o3_wrf_cbc = wrf.getvar(wrfout_cbc, "o3", timeidx=wrf.ALL_TIMES, method="cat")
co_wrf_cbc = wrf.getvar(wrfout_cbc, "co", timeidx=wrf.ALL_TIMES, method="cat")
no_wrf_cbc = wrf.getvar(wrfout_cbc, "no", timeidx=wrf.ALL_TIMES, method="cat")
no2_wrf_cbc = wrf.getvar(wrfout_cbc, "no2", timeidx=wrf.ALL_TIMES, method="cat")


# Retrieving pollutants from surface
o3_sfc = o3_wrf.isel(bottom_top=0)
co_sfc = co_wrf.isel(bottom_top=0)
no_sfc = no_wrf.isel(bottom_top=0)
no2_sfc = no2_wrf.isel(bottom_top=0)

o3_sfc_cbc = o3_wrf_cbc.isel(bottom_top=0)
co_sfc_cbc = co_wrf_cbc.isel(bottom_top=0)
no_sfc_cbc = no_wrf_cbc.isel(bottom_top=0)
no2_sfc_cbc = no2_wrf_cbc.isel(bottom_top=0)

# Transform surface polutant from ppm to ug/m3
o3_u_wrf = dp.ppm_to_ugm3(o3_sfc, t2, psfc, 48)
no_u_wrf = dp.ppm_to_ugm3(no_sfc, t2, psfc, 30)
no2_u_wrf = dp.ppm_to_ugm3(no2_sfc, t2, psfc, 46)
nox_u_wrf = (no_sfc + no2_sfc) * 1000
nox_u_wrf = nox_u_wrf.rename('nox')

# Transform surface polutant from ppm to ug/m3
o3_u_wrf_cbc = dp.ppm_to_ugm3(o3_sfc_cbc, t2, psfc, 48)
no_u_wrf_cbc = dp.ppm_to_ugm3(no_sfc_cbc, t2, psfc, 30)
no2_u_wrf_cbc = dp.ppm_to_ugm3(no2_sfc_cbc, t2, psfc, 46)
nox_u_wrf_cbc = (no_sfc_cbc + no2_sfc_cbc) * 1000
nox_u_wrf_cbc = nox_u_wrf_cbc.rename('nox')


# Getting dates to download from wrfout t2 variable
start_date, end_date = dp.qualar_st_end_time(t2)

# Loading List of stations and use only the stations inside wrfout
# cetesb_dom = dp.stations_in_domains("./cetesb2017_latlon.dat", wrfout, t2)
cetesb_dom = dp.stations_in_domains("./aqs_cbc.dat", wrfout, t2)
cetesb_masp = dp.stations_in_domains("./aqs_masp.dat", wrfout, t2)
cetesb_rural = dp.stations_in_domains("./aqs_rural.dat", wrfout, t2)


# Loading CETESB data
pol_filename = "cetesb_pols_2014_all_pol.pkl"
with open(pol_filename, 'rb') as a_dict:
    cetesb_pol = pickle.load(a_dict)




wrf_pol = dp.cetesb_from_wrf(cetesb_dom, 
                             (o3_u_wrf, no_u_wrf, no2_u_wrf, nox_u_wrf, co_sfc),
                             to_local=True)
cbc_pol = dp.cetesb_from_wrf(cetesb_dom, 
                             (o3_u_wrf_cbc, no_u_wrf_cbc, nox_u_wrf_cbc, 
                              no2_u_wrf_cbc, co_sfc_cbc),
                             to_local=True)

# Subsetting aqs in domains from main CETESB AQS
cetesb_pol_dom = {aqs:cetesb_pol[aqs] for aqs in list(wrf_pol)}

for aqs in list(cetesb_pol_dom):
    cetesb_pol_dom[aqs]['name'] = aqs



WRF, CET = dp.model_eval_setup(wrf_pol, cetesb_pol_dom, date_start="2014-10-06")
CBC, _ = dp.model_eval_setup(cbc_pol, cetesb_pol_dom, date_start="2014-10-06")


def plot_wrf_vs_cbc(wrf, cbc, cet, var, ylab, ylim, ax=None,
                    save_fig=False, fmt=None, hide=False,
                    legend=True, title=False):
    if ax is None:
        ax =plt.gca()
    formater = mdates.DateFormatter("%d")
    ax.plot(cet[var], color="black", label='Obs', linewidth=2.5)
    ax.plot(wrf[var], color="blue", linewidth=2.5,
            label="WRF-Chem/NALROM")
    ax.plot(cbc[var], color="red", linewidth=2.5,
            label="WRF-Chem/CAM-Chem")
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
        file_name = ("wrf_cbc_" + var + '_'+ wrf.name.unique()[0]
                      + fmt)
        plt.savefig(file_name, bbox_inches="tight", dpi=300)
        plt.clf()



def hourly_mean_plot(wrf, cbc, cet, var, ylab, ylim, ax=None, save_fig=False,
                      fmt=None, hide=False, legend=True, title=False):
    wrf_var = wrf[var]
    cbc_var = cbc[var]
    cet_var = cet[var]
    df = pd.concat([wrf_var, cbc_var, cet_var],
                    axis=1, keys=['wrf', 'cbc', 'obs'])
    df_h = df.groupby(df.index.hour).mean()
    if cet_var.dropna().empty:
        df_h['obs'] = np.nan
    if ax is None:
        ax = plt.gca()
    ax.plot(df_h.obs, color="black", label="Obs.", linewidth=2.5)
    ax.plot(df_h.wrf, color="blue", label="WRF-Chem/NALROM", linewidth=2.5)
    ax.plot(df_h.cbc, color="red", label="WRF-Chem/CAM-Chem", linewidth=2.5)
    ax.set_ylim(ylim)
    ax.grid(False)
    if hide:
        # ax.set_xticklabels(())
        ax.set_xlabel("")
    else:
        ax.set_xlabel("Hours (LT)")
    ax.set_ylabel(ylab)
    if legend:
        ax.legend(frameon=False)
    if title:
        ax.set_title(wrf.name.unique()[0])
        

aqs = cetesb_dom.name.values
pol = ['o3', 'no', 'no2', 'co']
ylims = [ (0, 250), (0, 250), (0, 110), (0, 4)]
ylabs = ["$O_3 \; (\mu g / m^3)$",
         "$NO \; (\mu g / m^3)$",
         "$NO_2\; (\mu g / m^3)$",
         "$CO \; (ppm)$"]



# fig, ax = plt.subplots(4, 4, figsize=(22, 10))
# for i in range(4):
#     for j in range(4):
#         if i == 0:
#             hourly_mean_plot(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                            ylabs[i], ylims[i], ax=ax[i, j], title=True, hide=True)
#         elif i <=2:
#             hourly_mean_plot(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                                ylabs[i], ylims[i], ax=ax[i, j], legend=False, hide=True)
#         else:
#             hourly_mean_plot(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                                ylabs[i], ylims[i], ax=ax[i, j], legend=False)
# plt.savefig('CBC_WRF_ts.png', dpi=300, bbox_inches="tight")


# fig, ax = plt.subplots(4, 4, figsize=(22, 10))
# for i in range(4):
#     for j in range(4):
#         if i == 0:
#             plot_wrf_vs_cbc(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                            ylabs[i], ylims[i], ax=ax[i, j], title=True, hide=True)
#         elif i <=2:
#             plot_wrf_vs_cbc(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                                ylabs[i], ylims[i], ax=ax[i, j], legend=False, hide=True)
#         else:
#             plot_wrf_vs_cbc(WRF[aqs[j]], CBC[aqs[j]], CET[aqs[j]], pol[i],
#                                ylabs[i], ylims[i], ax=ax[i, j], legend=False)

fig, ax = plt.subplots(1, 4, figsize=(16, 3))
for i, aqs in enumerate(cetesb_dom.name):
    if i == 0:
        hourly_mean_plot(WRF[aqs], CBC[aqs], CET[aqs], "o3",
                         "$O_3 \;  (\mu g \; m^ {-3})$", (0, 150), ax=ax[i],
                         title=True, legend=False)
    elif i <=2:
        hourly_mean_plot(WRF[aqs], CBC[aqs], CET[aqs], "o3",
                         "", (0, 150), ax=ax[i],
                         title=True, legend=False)
    else:
        hourly_mean_plot(WRF[aqs], CBC[aqs], CET[aqs], "o3",
                         "", (0, 150), ax=ax[i],
                         title=True, legend=True)        
plt.savefig("o3_cbc_dif_megan.pdf", dpi=300, bbox_inches="tight")



# Spatial difference
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

#loading MASP shapefile
reader_masp=shpreader.Reader("masp_shp/masp_shp.shp")
masp = list(reader_masp.geometries())
MASP = cfeature.ShapelyFeature(masp, ccrs.PlateCarree())


u10 = wrf.getvar(wrfout, "U10", timeidx=wrf.ALL_TIMES,
                  method="cat")
v10 = wrf.getvar(wrfout, "V10", timeidx=wrf.ALL_TIMES,
                  method="cat")


cbc_diff = ((o3_u_wrf - o3_u_wrf_cbc)
            .sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00"))
            .mean(dim="Time"))
u10_m = u10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00")).mean(dim="Time")
v10_m = v10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00")).mean(dim="Time")

val_limits = np.max([np.abs(cbc_diff.min()), cbc_diff.max()])

fig, ax = plt.subplots(subplot_kw = {'projection':ccrs.PlateCarree()})
cbc_diff.plot(x="XLONG", y="XLAT", ax =ax, vmin=-16, vmax=16,
              cmap="RdBu_r",
              cbar_kwargs={"shrink": 1, 
                           "label": "$\mu g \; m^{-3}$"})
ax.add_feature(cfeature.STATES.with_scale("10m"))
ax.add_feature(MASP, facecolor="none", linewidth=0.25, edgecolor="k")
ax.set_extent([cbc_diff.XLONG.min(), cbc_diff.XLONG.max(),
               cbc_diff.XLAT.min(), cbc_diff.XLAT.max()])
Q =ax.quiver(u10_m.XLONG.values[::5, ::5], u10_m.XLAT.values[::5, ::5],
          u10_m.values[::5, ::5], v10_m.values[::5, ::5])
ax.quiverkey(Q, 0.15, 0.15, 5, r'$5 \frac{m}{s}$', labelpos='E',
             coordinates='figure')

def spatial_diff_wrf_cbc(wrf, cbc, u10, v10, label, title, qkey=False, ax = None):
    cbc_diff = ((wrf - cbc)
                 .sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00"))
                 .mean(dim="Time"))
    u10_m = u10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00")).mean(dim="Time")
    v10_m = v10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00")).mean(dim="Time")
    
    sp_x, sp_y = (57, 49) # Sp coords in wrf grid
    sp_diff = cbc_diff.sel(south_north=sp_y, west_east=sp_x).values
    print(cbc_diff.name + " diff = ", str(sp_diff))

    if ax is None:
        ax = plt.gca()

    cbc_diff.plot(x="XLONG", y="XLAT", ax =ax,
                  cbar_kwargs={"shrink": 0.3, 
                               "label": label})
    ax.add_feature(cfeature.STATES.with_scale("10m"))
    ax.add_feature(MASP, facecolor="none", linewidth=0.25, edgecolor="k")
    ax.set_extent([cbc_diff.XLONG.min(), cbc_diff.XLONG.max(),
                   cbc_diff.XLAT.min(), cbc_diff.XLAT.max()])
    ax.set_title(title)
    Q = ax.quiver(u10_m.XLONG.values[::5, ::5], u10_m.XLAT.values[::5, ::5],
                  u10_m.values[::5, ::5], v10_m.values[::5, ::5])
    if qkey:
        ax.quiverkey(Q, 1.25, -0.25, 5, r'$5 \frac{m}{s}$', labelpos='E',
                     coordinates='axes')

fig, ax = plt.subplots(2, 2, figsize=(10, 10),
                       subplot_kw = {'projection':ccrs.PlateCarree()})
spatial_diff_wrf_cbc(o3_u_wrf, o3_u_wrf_cbc, u10, v10,"$\mu g \; m^{-3}$",
                     "$O_3$", qkey=True, ax =ax[0, 0])
spatial_diff_wrf_cbc(no_u_wrf, no_u_wrf_cbc, u10, v10,"$\mu g \; m^{-3}$",
                     "$NO$", ax =ax[0, 1])
spatial_diff_wrf_cbc(no2_u_wrf, no2_u_wrf_cbc, u10, v10,"$\mu g \; m^{-3}$",
                     "$NO_2$", ax =ax[1,0])
spatial_diff_wrf_cbc(co_sfc, co_sfc_cbc, u10, v10,"$ppm$",
                     "$CO$", ax =ax[1, 1])
plt.subplots_adjust(wspace=0.25, hspace=-0.5)
plt.savefig("cbc_effects_megan.pdf", dpi=300, bbox_inches="tight")



day = list(range(9, 22))
night = list(range(22,24)) + list(range(0, 9))


def spatial_diff_wrf_cbc_day(wrf, cbc, u10, v10, day, label, title, lims=None,
                             qkey=False, ax = None):
    cbc_diff_all = ((wrf - cbc)
                    .sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00")))                    
    u10_m_all = u10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00"))
    v10_m_all = v10.sel(Time=slice("2014-10-05 21:00", "2014-10-13 00:00"))
    
    cbc_diff = cbc_diff_all.sel(Time=np.in1d(cbc_diff_all["Time.hour"], day)).mean(dim="Time")
    u10_m = u10_m_all.sel(Time=np.in1d(u10_m_all["Time.hour"], day)).mean(dim="Time")
    v10_m = v10_m_all.sel(Time=np.in1d(v10_m_all["Time.hour"], day)).mean(dim="Time")
    
    sp_x, sp_y = (57, 49) # Sp coords in wrf grid
    sp_diff = cbc_diff.sel(south_north=sp_y, west_east=sp_x).values
    print(cbc_diff.name + " diff = ", str(sp_diff))
    
    if ax is None:
        ax = plt.gca()
    if lims is None:
        lims = np.max([np.abs(cbc_diff.min()), cbc_diff.max()])

    cbc_diff.plot(x="XLONG", y="XLAT", ax =ax, vmin = -lims, vmax= lims,
                  cmap="RdBu_r",
                  cbar_kwargs={"shrink": 0.3, 
                               "label": label})
    ax.add_feature(cfeature.STATES.with_scale("10m"))
    ax.add_feature(MASP, facecolor="none", linewidth=0.25, edgecolor="k")
    ax.set_extent([cbc_diff.XLONG.min(), cbc_diff.XLONG.max(),
                   cbc_diff.XLAT.min(), cbc_diff.XLAT.max()])
    ax.set_title(title)
    Q = ax.quiver(u10_m.XLONG.values[::5, ::5], u10_m.XLAT.values[::5, ::5],
                  u10_m.values[::5, ::5], v10_m.values[::5, ::5])
    if qkey:
        ax.quiverkey(Q, 1.25, -0.25, 5, r'$5 \frac{m}{s}$', labelpos='E',
                     coordinates='axes')

fig, axes = plt.subplots(4, 2, subplot_kw = {'projection':ccrs.PlateCarree()},
                         figsize=(10, 15))
axs = axes.flatten()
spatial_diff_wrf_cbc_day(o3_u_wrf, o3_u_wrf_cbc, u10, v10, day, "$\mu g \; m^{-3}$", "O$_3$ (daylight)",
                         ax=axs[0], lims=21)
spatial_diff_wrf_cbc_day(o3_u_wrf, o3_u_wrf_cbc, u10, v10, night,"$\mu g \; m^{-3}$", "O$_3$ (nighttime)",
                         ax=axs[1], lims=21)
spatial_diff_wrf_cbc_day(no_u_wrf, no_u_wrf_cbc, u10, v10, day, "$\mu g \; m^{-3}$", "NO (daylight)",
                         ax=axs[2], lims=6)
spatial_diff_wrf_cbc_day(no_u_wrf, no_u_wrf_cbc, u10, v10, night,"$\mu g \; m^{-3}$", "NO (nighttime)",
                         ax=axs[3],lims=6)
spatial_diff_wrf_cbc_day(no2_u_wrf, no2_u_wrf_cbc, u10, v10, day, "$\mu g \; m^{-3}$", "NO$_2$ (daylight)",
                         ax=axs[4], lims=8)
spatial_diff_wrf_cbc_day(no2_u_wrf, no2_u_wrf_cbc, u10, v10, night,"$\mu g \; m^{-3}$", "NO$_2$ (nighttime)",
                         ax=axs[5], lims=8)
spatial_diff_wrf_cbc_day(co_sfc, co_sfc_cbc, u10, v10, day, "$ppm$", "CO (daylight)",
                         ax=axs[6], lims=0.15)
spatial_diff_wrf_cbc_day(co_sfc, co_sfc_cbc, u10, v10, night,"$ppm$", "CO (nighttime)",
                         ax=axs[7], lims=0.15)
plt.subplots_adjust(wspace=0.25, hspace=-0.5)
plt.savefig("cbc_effects_day_night_megan.pdf", dpi=300, bbox_inches="tight")


# Getting SP coods in wrfout south_nort, west_east coordinanetes
sp_lat, sp_lon = [-23.5505, -46.6333]

import wrf as wrf
from netCDF4 import Dataset


x, y = wrf.ll_to_xy(wrfout, latitude=sp_lat, longitude=sp_lon)


