#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 02:05:02 2020

@author: mgavidia
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt


# Emission paths
# Vehicular emissions using PyChEmiss
veic00_path = "/scr2/mgavidia/python_stunts/PyChEmiss/wrfchemi_00z_d02"
veic12_path = "/scr2/mgavidia/python_stunts/PyChEmiss/wrfchemi_12z_d02"

# Industrial emission from ANTHRO emiss
ind00_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_ind/ind_wrfchemi_00z_d02"
ind12_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_ind/ind_wrfchemi_12z_d02"

# Residential emission from ANTHRO emiss
res00_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_res/res_wrfchemi_00z_d02"
res12_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_res/res_wrfchemi_12z_d02"

# Ships emission from SHIPS emiss
shp00_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_shp/shp_wrfchemi_00z_d02"
shp12_path = "/scr2/mgavidia/wrf_utils/ANTHRO/src/sp_shp/shp_wrfchemi_12z_d02"



# Loading emissions nc
veic00 = xr.open_dataset(veic00_path)
veic12 = xr.open_dataset(veic12_path)

ind00 = xr.open_dataset(ind00_path)
ind12 = xr.open_dataset(ind12_path)

res00 = xr.open_dataset(res00_path)
res12 = xr.open_dataset(res12_path)

shp00 = xr.open_dataset(shp00_path)
shp12 = xr.open_dataset(shp12_path)

# Emissions in WRF-Chem
emi_esp = ['E_CO', 'E_HCHO', 'E_C2H5OH', 'E_KET',
           'E_NH3', 'E_XYL', 'E_TOL', 'E_ISO', 'E_OLI',
           'E_OLT', 'E_OL2', 'E_HC8', 'E_HC5', 'E_ORA2',
           'E_ETH', 'E_ALD', 'E_CSL', 'E_SO2', 'E_HC3',
           'E_NO2', 'E_NO', 'E_CH3OH', 'E_PM25I', 'E_PM25J',
           'E_SO4I', 'E_SO4J', 'E_NO3I', 'E_NO3J', 'E_ORGI',
           'E_ORGJ', 'E_ECI', 'E_ECJ', 'E_SO4C', 'E_NO3C',
           'E_ORGC', 'E_ECC']

# Function to add
def htap_plus_local(local, htap, emi_esp):
    '''
    Function to add global emission reggrided from ANTHRO_EMISS
    to veicular emission from LAPAT preprocessor

    Parameters
    ----------
    local : xr.Dataset
        Veicular emissions from LAPAT preprocessor.
    htap : xr.Dataset
        Emission for HTAP Global emissions from ANTHRO_EMISS.
    emi_sp : List
        List of chem species in local.

    Returns
    -------
    A wrfchemi with added emissions.

    '''
    to_add = set(htap.data_vars).intersection(emi_esp)
    wrfchemi = local.copy()
    for emi in to_add:
        wrfchemi[emi] += htap[emi].values
    return wrfchemi


wrfchemi00zA = htap_plus_local(veic00, ind00, emi_esp)
wrfchemi12zA = htap_plus_local(veic12, ind12, emi_esp)

wrfchemi00zB = htap_plus_local(wrfchemi00zA, ind00, emi_esp)
wrfchemi12zB = htap_plus_local(wrfchemi12zA, ind12, emi_esp)

wrfchemi00z = htap_plus_local(wrfchemi00zB, shp00, emi_esp)
wrfchemi12z = htap_plus_local(wrfchemi12zB, shp12, emi_esp)


# If you want some plots uncomment this part


# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import cartopy.io.shapereader as shpreader
# import matplotlib.colors
#
# #loading MASP shapefile
# reader_masp=shpreader.Reader("masp_shp/masp_shp.shp")
# masp = list(reader_masp.geometries())
# MASP = cfeature.ShapelyFeature(masp, ccrs.PlateCarree())
#
# # Opening wrfinput
# wrfinput_path = "./wrfinput_d02"
# wrfinput = xr.open_dataset(wrfinput_path)
#
# xlat = wrfinput.XLAT[0, :, :]
# xlon = wrfinput.XLONG[0, :, :]
#
# import numpy as np
#
# lon = xlon.values[1, :]
# lat = xlat.values[:, 1]
# lons, lats = np.meshgrid(lon, lat)
#
#
# e_co_v12 = veic12.isel(Time=0, emissions_zdim=0).E_CO
# e_co_i12 = ind12.isel(Time=0, emissions_zdim_stag=0).E_CO
# e_co_r12 = res12.isel(Time=0, emissions_zdim_stag=0).E_CO
# e_co_t12 = wrfchemi12z.isel(Time=0, emissions_zdim=0).E_CO
#
#
# fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()},
#                          figsize=(8, 8))
# axs = axes.flatten()
# cv = axs[0].pcolormesh(lons, lats, e_co_v12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[0].set_title("Vehicular")
# ci = axs[1].pcolormesh(lons, lats, e_co_i12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[1].set_title("Industrial")
# cr = axs[2].pcolormesh(lons, lats, e_co_r12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[2].set_title("Residential")
# ct = axs[3].pcolormesh(lons, lats, e_co_t12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[3].set_title("Total")
# for ax in axs:
#     ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()],
#                   crs=ccrs.PlateCarree())
#     ax.coastlines('10m', color='white', linewidth=0.15)
#     ax.add_feature(MASP, facecolor='none', edgecolor='white',
#                    linewidth=0.15)
# ccs = [cv, ci, cr, ct]
# for i in range(4):
#     plt.colorbar(ccs[i], ax=axs[i], orientation='horizontal',
#                  shrink=0.7, pad=0.025)
# fig.subplots_adjust(hspace=0.1)
# plt.savefig('emissions_plot_CO.pdf', dpi=300, bbox_inches='tight')
# plt.savefig('emissions_plot_CO.png', dpi=300, bbox_inches='tight')
#
#
#
# e_no_v12 = veic12.isel(Time=0, emissions_zdim=0).E_NO
# e_no_i12 = ind12.isel(Time=0, emissions_zdim_stag=0).E_NO
# e_no_r12 = res12.isel(Time=0, emissions_zdim_stag=0).E_NO
# e_no_t12 = wrfchemi12z.isel(Time=0, emissions_zdim=0).E_NO
#
#
# fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()},
#                          figsize=(8, 8))
# axs = axes.flatten()
# cv = axs[0].pcolormesh(lons, lats, e_no_v12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[0].set_title("Vehicular")
# ci = axs[1].pcolormesh(lons, lats, e_no_i12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[1].set_title("Industrial")
# cr = axs[2].pcolormesh(lons, lats, e_no_r12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[2].set_title("Residential")
# ct = axs[3].pcolormesh(lons, lats, e_no_t12.values, transform=ccrs.PlateCarree(),
#                  alpha=0.95)
# axs[3].set_title("Total")
# for ax in axs:
#     ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()],
#                   crs=ccrs.PlateCarree())
#     ax.coastlines('10m', color='white', linewidth=0.15)
#     ax.add_feature(MASP, facecolor='none', edgecolor='white',
#                    linewidth=0.15)
# ccs = [cv, ci, cr, ct]
# for i in range(4):
#     plt.colorbar(ccs[i], ax=axs[i], orientation='horizontal',
#                  shrink=0.7, pad=0.025)
# fig.subplots_adjust(hspace=0.1)
# plt.savefig('emissions_plot_NO.pdf', dpi=300, bbox_inches='tight')
# plt.savefig('emissions_plot_NO.png', dpi=300, bbox_inches='tight')
#
#
#
#
#
# wrfchemi00z.to_netcdf("wrfchemi_00z_d02_simple_voc_cal_edgar")
# wrfchemi12z.to_netcdf("wrfchemi_12z_d02_simple_voc_cal_edgar")
#
#
