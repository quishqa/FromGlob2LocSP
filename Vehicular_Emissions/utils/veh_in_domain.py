#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:59:38 2020

@author: mgavidia
"""

import pandas as pd
import geopandas as gpd
import geobr as gb
import xarray as xr
import matplotlib.pyplot as plt
from shapely.geometry import Polygon



# States in domain
states_codes = ["SP", "RJ", "PR", "MG", 
                "ES", "SC", "MS"]

# Saving shapefile in a dict
states_in_domain = {code: gb.read_municipality(code_muni=code, year=2014)
                    for code in states_codes}


# Merging municipalities in domain
mun_gdf = pd.concat([state for code,state in states_in_domain.items()],
                    ignore_index=True)    


# Reading veicule number from
# https://www.gov.br/infraestrutura/pt-br/assuntos/transito/conteudo-denatran/frota-de-veiculos-2014
veic = pd.read_excel("./frota_por_municipio_e_tipo_out_2014.xlsx",
                     sheet_name="OUT_2014",
                     header=3)

# Clip shapefiles to WRF domain
wrf = xr.open_dataset("./geo_em.d01.nc")
xlat_c = wrf.XLAT_C.values[:, :, :]
xlon_c = wrf.XLONG_C.values[:, :, :]

xlat_min, xlat_max = xlat_c.min(), xlat_c.max()
xlon_min, xlon_max = xlon_c.min(), xlon_c.max()

wrf_dom = Polygon([[xlon_min, xlat_min], [xlon_max, xlat_min],
                   [xlon_max, xlat_max], [xlon_min, xlat_max]])

wrf_dom_gdf = gpd.GeoDataFrame([1], geometry=[wrf_dom], crs=mun_gdf.crs)

mun_clipped = gpd.clip(mun_gdf, wrf_dom_gdf)

# mun_clipped.plot()

# Merging municipalities with veic
muns_in_domain = mun_clipped.name_muni
veic_br = veic[["UF", "MUNICIPIO", "TOTAL"]].copy()

# The strategy is to combine muns_in_domains serie
# with veic_br, to do so, we need to remove the diacritics
# and turn all to lowers, we'll use unidecode modul

from unidecode import unidecode


muns_decode = (muns_in_domain
               .apply(lambda name: unidecode(name))
               .str.lower())

muns_to_merge = pd.DataFrame({"mun_name": muns_in_domain.values,
                              "mun_decode": muns_decode})
muns_to_merge.set_index("mun_decode", inplace=True)


veic_br["mun_decode"] = veic_br.MUNICIPIO.str.lower()
veic_br.set_index("mun_decode", inplace=True)

veic_in_domain = muns_to_merge.join(veic_br)

print("Muns total = ",len(veic_in_domain.index))
print("Muns with no data =",veic_in_domain.TOTAL.isna().sum() )
print("Vehicles in domain = ",veic_in_domain.TOTAL.sum())
