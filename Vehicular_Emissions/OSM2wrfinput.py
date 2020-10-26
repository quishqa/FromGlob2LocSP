#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 03:04:19 2020

This script downloads openstreetmap street information into 
your WRF-Chem domain.

@author: mgavidia
"""

import osmnx as ox
import xarray as xr
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

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

# ox.save_graph_geopackage(SP, filepath='./data/geo_d02.gpkg') # Save as geopandas
ox.save_graphml(SP, filepath="./data/geo_d01.graphml") # Save as graph for OSMNx
# ox.save_graph_shapefile(SP, filepath='./data/geo_d01') # Save as shapefile

# Retrieving wrfinput cell corner coordinates
lat_c = xlat_c[:, 1]
lon_c = xlon_c[1, :]

# Creating a grid according to wrf cell corners
# Based on 
# https://gis.stackexchange.com/questions/269243/creating-polygon-grid-using-geopandas

poly = []

for j in range(len(lat_c) - 1):
    for i in range(len(lon_c) - 1):
        poly.append(Polygon([
            (lon_c[i], lat_c[j + 1]),
            (lon_c[i + 1], lat_c[j + 1]),
            (lon_c[i + 1], lat_c[j]),
            (lon_c[i], lat_c[j])
            ]))
grid_wrf = gpd.GeoDataFrame({'geometry':poly})
grid_wrf.to_file("grid_wrf_d01.shp")



# Extracting streets from SO graph file as geopandas
SP = ox.load_graphml("./data/geo_d02.graphml")
sp_gpd = ox.graph_to_gdfs(SP,  nodes=False, edges=True)


# Adding CRS to grid and adding an ID columns
grid = grid_wrf.set_crs(sp_gpd.crs)
grid['ID'] = range(0, len(grid))

# Subsetting requiered columns
gdf = sp_gpd[['highway', 'length', 'geometry']]

# Clipping to grid area
roads_clip = gpd.clip(gdf, grid)


# Plotting grid and roads
fig, ax = plt.subplots(figsize = (10, 10))
roads_clip.geometry.plot(linewidth=0.5, ax=ax, color="Black")
grid.geometry.plot(linewidth=0.5, edgecolor='k', color=None,ax=ax)
plt.savefig("a_plot_d01.png", dpi=300, bbox_inches="tight")


# Intercept roads_clip with grid
roads_int = gpd.overlay(roads_clip, grid, how='intersection')
roads_grid = roads_int.dissolve("ID") # Total_length (SUMlongKm)


# Calculating Main roads
roads_type = roads_int.copy()
roads_type = roads_type[((roads_type.highway != "primary_link") &
                         (roads_type.highway != "secondary_link") &
                         (roads_type.highway != "tertiary_link") &
                         (roads_type.highway != "trunk_link") &
                         (roads_type.highway != "motorway_link"))]
roads_grid_main = roads_type.dissolve("ID")

# Calculating UrbanGrau
# We use mercator projection to get the output in meters with precission
roads_grid['longKm'] = roads_grid.geometry.to_crs("EPSG:32733").length / 1000
roads_grid_main['mainKm'] = roads_grid_main.geometry.to_crs("EPSG:32733").length /1000

roads_grid['urban'] = (roads_grid_main.mainKm.values * 
                       (roads_grid.longKm.values - roads_grid_main.mainKm.values) /
                       (roads_grid_main.mainKm.values + 1))
roads_grid['mainKm'] = roads_grid_main.mainKm.values


# Creating input to LApat preprocessor
grid.set_index("ID", inplace=True)

grid_a4w = grid.join(roads_grid[["longKm", "mainKm", "urban"]])
grid_a4w = grid_a4w.fillna(0)
grid_a4w["centroid"] = grid_a4w.geometry.centroid # Same center as WRF domain

# Exporting to csv

s_df = grid_a4w.copy()
s_df["X"] = grid_a4w.centroid.geometry.x
s_df["Y"] = grid_a4w.centroid.geometry.y

s_df[["X", "Y", "longKm", "mainKm", "urban"]].to_csv("s3_test2.txt", sep=" ", header=False)




