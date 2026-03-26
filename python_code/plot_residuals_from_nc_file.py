#! /usr/bin/env python3

#===========================================================================
# script to make a paper figure for residuals calculated by the the DIVA tool
# Author: Itzel Ruvalcaba
# data: 2025-17-28

# requires modeule cartopy!!!!
#============================================================================
# load packages
import numpy as np
import xarray as xr
import netCDF4
from pathlib import Path

# for plots and maps
import matplotlib
matplotlib.use("QtAgg")  # non-GUI backend (otherwise error when plotting)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as cm  # to get diff colormap with gray for 0
from matplotlib.lines import Line2D

#=================
# years of files to read
year = 1970

#choose seson of file to read
seasons = ['Spring', 'Autumn', 'Summer', 'Winter']

# paths to residual files
# *** Change to correct files
# /nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Baltic_Proper/20260114_1653/
results_dir = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results_lena_temp/Baltic_Proper/20260320_1659/")
levels = np.array([-150, -100, -50, 0, 50, 100, 150])
cmap = cm.get_cmap("viridis", len(levels) - 1)
norm = mcolors.BoundaryNorm(levels, cmap.N)

#=================
# file to open
for season in seasons:
    print('Opening file for season ', season)
    # change path to correct path
    mfile = results_dir / f"DIVArun/Oxygen_{year}_{season}_residuals.nc"

    #READ DATA
    # ds = xr.open_mfdataset(mfile, decode_times=True)
    ds = xr.open_dataset(mfile, engine='h5netcdf')
    lon=ds.obslon
    lat=ds.obslat
    depth=ds.obsdepth
    time=ds.obstime
    da=ds['Oxygen_residual'].compute()
    filtered_da = da.where((da < -50) | (da > 50))
    #!!!da=da.isel(....) #select below surface
    # da's dim = 1D
    # variables: 
    # time, depth, lat, lon (all in one single vector)
    
    # categorize depths into marker sizes:
    # create bin intervals
    depth_bins = np.arange(50, 250, 50)
    # keep 10 m above 100 m, 25 below and maybe 50 below 200m
    depth_labels = [ f"{depth_bins[i]}–{depth_bins[i+1]} m" for i in range(len(depth_bins) - 1)]
    # arrange depth into its categories
    depth_cat = np.digitize(depth, depth_bins)

    #=================================
    # PLOT RESIDUAL FIGURE FOR SEASONS
    #================================
    vmax = da.max() # drop surface too
    vmin = da.min()

    # reduce vmin and vmax so colors show better, removing too large values

    # Scale marker sizes (you can adjust the factor)
    sizes = depth_cat * 20 #np.ones(len(da.values)) * 10

    # Create legend elements: one per depth bin
    #bin_sizes = [(i+0.5)*10 for i in range(len(depth_bins)-1)]
    bin_sizes = np.arange(1, len(depth_bins)) * 10

    #Create legend element
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor='gray', markersize=np.sqrt(size),alpha=0.5)
        for label, size in zip(depth_labels, bin_sizes)
        ]

    len(legend_elements)
    #=================
    #MAPS STARTS
    fig, ax = plt.subplots(figsize=(15, 12), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)
    ax.set_extent([9, 32, 53, 66], ccrs.PlateCarree())
    # land in gray:
    ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor="none", zorder=0, alpha=0.5)
    # zorder ==> gives the llayer position (0 = background, 1 top layer, goes until 10). 
    # Coast
    ax.coastlines(resolution="10m", linewidth=1)
    # Grid
    gl = ax.gridlines(draw_labels=True, zorder=2) #to keep ticks inside map
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 18}#, "weight": "bold"}
    gl.ylabel_style = {"size": 18}#, "weight": "bold"}

    # Plot data

    # Because Matplotlib does not separate edge alpha from face alpha in scatter. We need two layers for the edges:
    cs = plt.scatter(
        x=lon,
        y=lat,
        c=filtered_da.values,
        s=sizes,
        transform=ccrs.PlateCarree(),
        marker='o',
        cmap=cmap,
        norm=norm,
        rasterized=True,
        edgecolors='none',
        alpha=1,
        zorder=3
    )
    # cs = plt.scatter(x=lon, y=lat, c = filtered_da.values, s = sizes, transform=ccrs.PlateCarree(), marker = 'o', cmap=mycmap, facecolor='none', edgecolors='black',  linewidth = 0.1, zorder=4, vmin=vmin, vmax=vmax)
    # rasterized=true ==>  good for large grids + PDFs

    cbar = fig.colorbar(
        cs,
        ax=ax,
        orientation="horizontal",
        pad=0.05,
        shrink=0.5,
        boundaries=levels,
        ticks=levels
    )
    cbar.ax.tick_params(labelsize=20)
    # what are the units: umol/l ?
    #cbar.set_label("Residuals before truncation (umol/L)", fontsize=18)
    # goal is to have:
    cbar.set_label("Difference between observations and interpolation at obs. point (umol/L)", fontsize=18)

    # add legend
    ax.legend(handles=legend_elements, title="Depth", loc="lower right")
    ax.set_title('%s for year %s' %(season,year), fontsize = 24, fontweight = 'bold')
    #plt.show()
    
    plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_{year}_{season}.png')
    plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_{year}_{season}.pdf', format="pdf", bbox_inches="tight")

    plt.close()

