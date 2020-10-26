#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 17:06:58 2020

@author: mgavidia
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 03:04:19 2020

@author: mgavidia
"""

import osmnx as ox
import xarray as xr
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt

geo = xr.open_dataset("geo_em.d01.nc")

# Reading WRF cell corners
xlat_c = geo.XLAT_C.values[0, : ,:]
xlon_c = geo.XLONG_C.values[0, : ,:]

# Downloading OSM data
# We are going to downloading by selecting a box 
north = xlat_c.max()
south = xlat_c.min()
east = xlon_c.max()
west = xlon_c.min()

# Custom filter
# Based on https://github.com/gboeing/osmnx-examples/blob/master/notebooks/08-custom-filters-infrastructure.ipynb

street_main = ['primary', 'secondary', 'tertiary', 'motorway', 'trunk']
street_links = [st + "_link" for st in street_main]
all_streets = street_main + street_links

cf =( '[' +'"highway"' + "~" + 
     '"' + "|".join(all_streets) + '"' 
     + ']')

print("Starting download")
# Add a buffer zone of 0.2 to ensure WRF_Domain is inside
# traffic lines
SP = ox.graph_from_bbox(north, south, 
                       east, west,
                       network_type="drive",
                       custom_filter=cf)

print("Download complete, now saving")
# ox.save_graph_geopackage(SP, filepath='./data/geo_d02.gpkg') # Save as geopandas
ox.save_graphml(SP, filepath="./data/geo_d01_t.graphml") # Save as graph for OSMNx
# ox.save_graph_shapefile(SP, filepath='./data/geo_d01') # Save as shapefile
print("all done")
