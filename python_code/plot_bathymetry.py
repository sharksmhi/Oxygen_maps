#! /usr/bin/env python3

#===========================================================================
# script to make a nice figure for Bathymetry used in the DIVA tool to calculate oxygen area
# Author: Itzel Ruvalcaba
# data: 2025-11-28
#============================================================================
# load packages
import numpy as np
import xarray as xr
import netCDF4

# for plots and maps
import matplotlib
matplotlib.use("Agg")  # non-GUI backend (otherwise error when plotting)
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)

import cmocean

#=================
# paths
pathin = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data"

#=================
# file to open
mfile = "%s/bat_elevation_Baltic_Sea_masked.nc" %pathin

def set_up_cartopy_map(axis):
    """
    Create a map with cartopy in the given axis. 
    The map projections is PlateCarre to be able to plot gridlines and ticks.
    When creating the figures use ccrs.Mercator() to get a map with better ratio for our latitudes
    
    :param axis: matplotlib axis object
    """
    # Set map extent
    axis.set_extent(
        [
            9,
            29,
            53.5,
            61,
        ],
        crs=ccrs.PlateCarree()
    )

    # Coastlines
    axis.coastlines(resolution='10m', linewidth=0.5, color='gray')

    # Gridlines (meridians/parallels)
    gl = axis.gridlines(
        crs=ccrs.PlateCarree(), 
        draw_labels=True,
        linewidth=0.1,
        color='gray',
        alpha=0.5
    )
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 14}
    gl.ylabel_style = {'size': 14}

    return axis
#=================
#MAPS STARTS
mycmap=cmocean.cm.deep

#READ DATA
# First run (slow)
ds = xr.open_mfdataset(mfile, decode_times=False)
lon=ds.lon
lat=ds.lat
da=ds['bat'].compute()
# quick view
#da.plot.pcolormesh(ax=ax, x='lon', y='lat', transform=ccrs.PlateCarree(), cmap=mycmap)

# Mask land / zero values
da = da.where(da != 0) * -1

"""
# Om man vill spara som en nc fil och läsa igen, kanske snabbare
# 3. Wrap in a Dataset (keeps coordinates lon/lat)
ds_to_save = da.to_dataset(name="bat")

# 4. Save as NetCDF
ds_to_save.to_netcdf(f"bathymetry/bat_clean.nc", format="NETCDF4")

# Future runs (fast)
ds = xr.open_dataset(f"bathymetry/bat_clean.nc")
da = ds['bat']
lon = ds.lon
lat = ds.lat"""

# FIGURE
# make all text white (seashell)
plt.rcParams.update({
    "text.color": "seashell",       # general text color
    "axes.labelcolor": "seashell", # x/y axis labels
    "xtick.color": "seashell",        # x ticks
    "ytick.color": "seashell",        # y ticks
    "axes.titlecolor": "seashell"    # axes title
})

fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.Mercator()})
ax = set_up_cartopy_map(ax)

# Plot bathymetry
cs = da.plot.pcolormesh(
    ax=ax, x="lon", y="lat", 
    transform=ccrs.PlateCarree(), 
    vmin=50,
    vmax=200,
    cmap=mycmap, add_colorbar=False, rasterized=True)

# Define contour levels (depths in meters)
levels = [60, 100, 120]
# Add contour lines
contours = ax.contour(
    da.lon,
    da.lat,
    da,
    levels=levels,
    colors='black',
    linewidths=0.2,
    transform=ccrs.PlateCarree(),
    zorder=6
)

# Optional: label contours
ax.clabel(contours, inline=True, fontsize=10, fmt='%d')

cbar = fig.colorbar(cs, ax=ax, orientation="horizontal", pad=0.05, shrink=0.5)
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Bathymetry (m)", fontsize=18)

#plt.show()
plt.savefig('resultat/figures/bathymetry.png', transparent=True)
plt.savefig('resultat/figures/bathymetry.pdf', format="pdf", bbox_inches="tight")

plt.close()

