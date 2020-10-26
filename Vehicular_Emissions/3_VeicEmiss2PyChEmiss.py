#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:18:19 2020

This script creates the input to PyChEmiss
based on vehicular emissions spatially disaggregated
by steet length.

@author: mgavidia
"""

import pandas as pd
import numpy as np

osm_output = pd.read_csv('s3_test2.txt', 
                         names=['id', 'x', 'y', 'lonkm', 'mainkm', 'urban'],
                         delim_whitespace=True)

osm = osm_output[['id', 'x', 'y', 'lonkm']].copy()



total_vehicles = 49365478 # calculated in all domain

# From CETESB report
# E_CO_year = 417.026 # kTn/year
# E_HC_year = 85.216# kTn/year
# E_NOX_year = 171.384 * 0.195#* 0.2# KTn/year
# E_SO2_year = 5.574  # KTn/year
# E_PM_year = 4.848   # kTn/year

# MASP totals from CETESB report first cal
# E_CO_year = 160.57 # kTn/year
# E_HC_year = 37.87 * 3.5 /4 # kTn/year
# E_NOX_year = 55.77 / 4# KTn/year
# E_SO2_year = 4.23  # KTn/year
# E_PM_year = 1.54   # kTn/year

E_CO_year = 6837.5 * 0.55# kTn/year
E_HC_year = 2851.65 * 0.30  # kTn/year
E_NOX_year = 1134.08 * 0.30 # KTn/year
E_SO2_year = 72.3604 * 0.30 # KTn/year
E_PM_year = 144.312 # * 0.35 # kTn/year

E_NO_year = E_NOX_year * 0.9
E_NO2_year = E_NOX_year * 0.1

cell_area = 81.0 # km2

hrsplt_co = [0.019, 0.012, 0.008, 0.004, 0.003, 0.003,
             0.006, 0.017, 0.047, 0.074, 0.072, 0.064,
             0.055, 0.052, 0.051, 0.048, 0.052, 0.057,
             0.068, 0.087, 0.085, 0.057, 0.035, 0.024]
hrsplt_no = [0.019, 0.015, 0.012, 0.010, 0.008, 0.009,
             0.015, 0.030, 0.048, 0.053, 0.051, 0.061,
             0.064, 0.064, 0.061, 0.060, 0.060, 0.065,
             0.065, 0.066, 0.056, 0.044, 0.035, 0.027]


emitted_species = ["E_SO2","E_NO","E_ALD","E_HCHO","E_ORA2","E_NH3","E_HC3","E_HC5", 
                   "E_HC8","E_ETH","E_CO","E_CO2","E_OL2","E_OLT","E_OLI","E_TOL",   
                   "E_XYL","E_KET","E_CSL","E_ISO","E_NO2","E_CH3OH","E_C2H5OH",     
                   "E_PM25I","E_PM25J","E_SO4I","E_SO4J","E_NO3I","E_NO3J","E_ORGI", 
                   "E_ORGJ","E_ECI","E_ECJ","E_NAAI","E_NAAJ","E_SO4C","E_NO3C",     
                   "E_ORGC","E_ECC","E_ORGI_A","E_ORGJ_A","E_ORGI_BB","E_ORGJ_BB",   
                   "E_CO_A","E_CO_BB","E_PM_10"]

# E_ORA2 (benzoic acid)
MW = [64.00, 30.00, 44.05, 30.09, 122.12, 17.00, 42.66, 60.00,
      96.00, 30.07, 28.00, 44.00, 28.05, 56.00, 56.00, 92.00,
      104.00, 58.08, 108.14, 68.12, 46.00, 32.00, 46,
      1., 1., 1., 1., 1., 1., 1.,
      1., 1., 1., 1., 1., 1., 1., 
      1., 1., 1., 1., 1., 1., 
      1., 1., 1.]

# Dataframe with species and their molecular weight
CBMZ = pd.DataFrame({'emi_spe': emitted_species,
                     'mw': MW})

# VOC fractions from air quality forecast file
voc_split = pd.read_csv("./voc_frac.csv",
                        names=['emi_spe', 'tot', 'frac', 'frac_round', 'per'])
voc_split = voc_split.merge(CBMZ, on="emi_spe")
voc_split['total_emi'] = E_HC_year * voc_split.frac.values

# Getting the ratio between vehicles/road_length
veic_per_cell = total_vehicles / osm.lonkm.sum()



def from_kTn_year_to_mol_day(emiss_year, MW):
    '''
    Convert emission in kTn/year to mol/day

    Parameters
    ----------
    emiss_year : float
        Emission of specie.
    MW : float
        Molecular weigth.

    Returns
    -------
    emiss_mol_day : float
        emission in mol/day.

    '''
    emiss_mol_day = emiss_year * 10**9 / 365.0 / MW
    return emiss_mol_day



def emiss_per_cell(emiss_mol_day, total_vehicles, veic_per_cell):
    '''
    Calculate emissions per cell.At the end emission is distributed by 
    lonKm

    Parameters
    ----------
    emiss_mol_day : float
        Total emission in mol/day.
    total_vehicles : float
        Total number of vehicles.
    veic_per_cell : float
        Vehicles per cell.

    Returns
    -------
    emiss_per_cell : float
        DESCRIPTION.

    '''
    emiss_per_cell = emiss_mol_day /total_vehicles * veic_per_cell
    return emiss_per_cell


osm['veh'] = osm.lonkm.values * veic_per_cell
osm['E_CO'] = emiss_per_cell(from_kTn_year_to_mol_day(E_CO_year, 28), 
                             total_vehicles, osm.veh.values) /cell_area
osm['E_NO'] = emiss_per_cell(from_kTn_year_to_mol_day(E_NO_year, 30), 
                             total_vehicles, osm.veh.values) /cell_area
osm['E_NO2'] = emiss_per_cell(from_kTn_year_to_mol_day(E_NO2_year, 46), 
                             total_vehicles, osm.veh.values) / cell_area
osm['E_SO2'] = emiss_per_cell(from_kTn_year_to_mol_day(E_SO2_year, 64), 
                             total_vehicles, osm.veh.values) /cell_area
osm['E_PM'] = emiss_per_cell(from_kTn_year_to_mol_day(E_PM_year, 1), 
                             total_vehicles, osm.veh.values)
osm['E_PM'] = osm['E_PM'] / 3600. # to ug/s


# VOC
for voc, emi, mw in zip(voc_split.emi_spe, voc_split.total_emi, voc_split.mw):
    osm[voc] = emiss_per_cell(from_kTn_year_to_mol_day(emi, mw),
                              total_vehicles, osm.veh.values) / cell_area


times = ["{:02d}".format(t) + "Z" for t in range(0, 24)]

all_emiss = {}
all_emiss ={t:osm.copy()  for t in times}

def temporal_split(df, hrsplt_co, hrsplt_no):
    for sp in list(df.columns[5:]):
        if sp in ['E_NO', 'E_NO2', 'E_SO2', 'E_PM']:
            df[sp] = df[sp] * hrsplt_no
        else:
            df[sp] = df[sp] * hrsplt_co
    return df



for i in range(24):
    all_emiss[times[i]] = temporal_split(all_emiss[times[i]],
                                         hrsplt_co[i], 
                                         hrsplt_no[i])

    
all_emiss_df = pd.concat(all_emiss.values(), ignore_index=True)            
        
# Information based on
# Andrade et al. 2015, Vara-Vela et al. 2017 

# PM
f_fina_total = 0.670   # FRACAO FINA TOTAL
f_fina_so4 = 0.070   # FRACAO FINA DE SO4
f_fina_no3 = 0.016   # FRACAO FINA DE NO3
f_fina_org = 0.420   # FRACAO FINA DE ORG
f_fina_ec = 0.190   # FRACAO FINA DE EC
f_fina_pm25 = 0.304   # FRACAO FINA DE PM25
f_grossa_total = 0.330   # FRACAO GROSSA TOTAL
f_grossa_so4 = 0.000   # FRACAO GROSSA DE SO4
f_grossa_no3 = 0.000   # FRACAO GROSSA DE NO3
f_grossa_org = 0.000   # FRACAO GROSSA DE ORG
f_grossa_ec = 0.000   # FRACAO GROSSA DE EC
i_so4 = 0.136   # MODA AITKEN DE SO4
j_so4 = 0.864   # MODA ACUMULACAO DE SO4
i_no3 = 0.230   # MODA AITKEN DE NO3
j_no3 = 0.770   # MODA ACUMULACAO DE NO3
i_org = 0.190   # MODA AITKEN DE ORG
j_org = 0.810   # MODA ACUMULACAO DE ORG
i_ec = 0.940   # MODA AITKEN DE EC
j_ec = 0.060   # MODA ACUMULACAO DE EC
i_pm25 = 0.250   # MODA AITKEN DE PM25
j_pm25 = 0.750   # MODA ACUMULACAO DE PM25


def pm_split(df):
    df['E_PM25I'] = df['E_PM'].values * f_fina_total * f_fina_pm25 * i_pm25
    df['E_PM25J'] = df['E_PM'].values * f_fina_total * f_fina_pm25 * j_pm25
    df['E_SO4I'] = df['E_PM'].values * f_fina_total * f_fina_so4 * i_so4
    df['E_SO4J'] = df['E_PM'].values * f_fina_total * f_fina_so4 * j_so4
    df['E_NO3I'] = df['E_PM'].values * f_fina_total * f_fina_no3 * i_no3
    df['E_NO3J'] = df['E_PM'].values * f_fina_total * f_fina_no3 * j_no3
    df['E_ORGI'] = df['E_PM'].values * f_fina_total * f_fina_org * i_org
    df['E_ORGJ'] = df['E_PM'].values * f_fina_total * f_fina_org * j_org
    df['E_ECI'] = df['E_PM'].values * f_fina_total * f_fina_ec * i_ec
    df['E_ECJ'] = df['E_PM'].values * f_fina_total * f_fina_ec * j_ec
    df['E_NAAI'] = 0.0
    df['E_NAAJ'] = 0.0
    df['E_SO4C'] = df['E_PM'].values * f_grossa_so4
    df['E_NO3C'] = df['E_PM'].values * f_grossa_no3
    df['E_ORGC'] = df['E_PM'].values * f_grossa_org
    df['E_ECC'] = df['E_PM'].values * f_grossa_ec
    df['E_ORGI_A'] = 0.0
    df['E_ORGJ_A'] = 0.0
    df['E_ORGI_BB'] = 0.0
    df['E_ORGJ_BB'] = 0.0
    df['E_CO_A'] = 0.0
    df['E_CO_BB'] = 0.0
    df['E_PM_10'] = df['E_PM'] * f_grossa_total
    return df
    
    
all_emiss_df = pm_split(all_emiss_df)


def to_pychemis(df, emiss_sp, file_name):
    df['E_NH3'] = 0.0
    df['E_CSL'] = 0.0
    df['E_CO2'] = 0.0
    df['E_ORA2'] = 0.0
    col_names = ["x", "y"]
    col_to_pychemis = col_names + emiss_sp
    df[col_to_pychemis].to_csv(file_name, header=False,
                               sep="\t")


to_pychemis(all_emiss_df,emitted_species, 
            "simple_model_d01.txt")
