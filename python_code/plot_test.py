import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from datetime import datetime
import pandas as pd
import json
# import cartopy
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import inset_locator
from mpl_toolkits.axes_grid1 import make_axes_locatable

def sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    # Create a Basemap instance
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')
    # Plot land borders from the bathymetry file

    # Plot data
    data = ds[parameter].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index]).squeeze()
    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)

    return m, pcm

def sub_plot_observations_basemap(ds, parameter, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None):
    
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    m, pcm = sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax)
    
    # Change vmin and vmax after the plot is created
    pcm.set_clim(vmin=vmin, vmax=vmax)
    # Create a custom discrete colormap
    #        -45     0        45    90 (2 ml/l)  135      180 (4 ml/l)
    colors = ['grey', 'brown', 'red',   'orange', 'yellow',  'pink',  'green']
    cmap = ListedColormap(colors)
    # Change the colormap after the plot is created
    pcm.set_cmap(cmap)

    # nedan är endast för att få med observationer
    # OBS att DIVAnd resultatet ligger under dimensionen 'Oxygen' i datasetet och obsevrationerna under 'Oxygen_data'
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']

    lon, lat = m(lon, lat)

    m.scatter(lon, lat, s=5, c=observations, cmap=cmap, edgecolors='k', linewidth=0.2, facecolor='none',
              vmin=vmin, vmax=vmax)
    
    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
    cbar.ax.tick_params(labelsize = 5)
    # cbar.set_label(label='µmol/l', fontsize = 10,  y=1.05)

    # Modify the colormap levels to control the step length
    levels = np.arange(vmin, vmax+1, 45)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    pcm.set_norm(norm)
    # Set the colorbar levels explicitly
    cbar.set_ticks(levels)

    # Add labels to the subplot
    axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m', fontsize = 10)
    axis.set_xlabel('Longitude', fontsize = 10)
    axis.set_ylabel('Latitude', fontsize = 10)

    return pcm

def sub_plot_area_at_threshold_basemap(ds, parameter, axis, year, vmin, vmax, threshold, unit='umol/l', bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))

    # Create a Basemap instance with Mercator projection
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')

    # Plot data
    try:
        data = ds[parameter].sel(time=ds['time'][time_index])
    except KeyError:
        print(f'No {parameter} in file')
        
    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    # Plot data using pcolormesh
    pcm = m.pcolormesh(lon, lat, data, cmap='ocean', vmin=vmin, vmax=vmax)
    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
    cbar.ax.tick_params(labelsize = 7)

    # Add labels to the subplot with adjusted text size
    axis.set_title(f'Area <= {threshold} {unit}', fontsize=10)
    axis.set_xlabel('Longitude', fontsize=8)
    axis.set_ylabel('Latitude', fontsize=8)


def sub_plot_error_area_at_threshold_basemap(ds, parameter, axis, year, vmin, vmax, threshold, unit='umol/l', bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))

    # Create a Basemap instance with Mercator projection
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')

    # Plot data
    try:
        data = ds[parameter].sel(time=ds['time'][time_index])
    except KeyError:
        print(f'No {parameter} in file')

    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    # Plot data using pcolormesh
    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)
    
    # Change vmin and vmax after the plot is created
    pcm.set_clim(vmin=0, vmax=1)
    # Create a custom discrete colormap
    #        0       0.3       0.5   
    colors = ['green', 'orange', 'red', 'red']
    cmap = ListedColormap(colors)
    # Change the colormap after the plot is created
    pcm.set_cmap(cmap)
    
    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax=cbax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)

    # Modify the colormap levels to control the step length
    levels = [0, 0.3, 0.5, 1]
    #norm = plt.Normalize(vmin=vmin, vmax=vmax)
    #pcm.set_norm(norm)
    # Set the colorbar levels explicitly
    cbar.set_ticks(levels)

    # Add labels to the subplot with adjusted text size
    axis.set_title(f'Error at <= {threshold} {unit} depth', fontsize=10)
    axis.set_xlabel('Longitude', fontsize=8)
    axis.set_ylabel('Latitude', fontsize=8)


def sub_plot_errorfields_basemap(ds, parameter, axis, year, show_depth, vmin, vmax):

    m, pcm = sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax)

    # Change vmin and vmax after the plot is created
    pcm.set_clim(vmin=0, vmax=1)
    # Create a custom discrete colormap
    #        0       0.3       0.5
    colors = ['green', 'orange', 'red', 'red']
    cmap = ListedColormap(colors)
    # Change the colormap after the plot is created
    pcm.set_cmap(cmap)

    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax=cbax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)

    # Modify the colormap levels to control the step length
    levels = [0, 0.3, 0.5, 1]
    # norm = plt.Normalize(vmin=vmin, vmax=vmax)
    # pcm.set_norm(norm)
    # Set the colorbar levels explicitly
    cbar.set_ticks(levels)

    axis.set_title(f'Errorfield at {show_depth} m results')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')

def plot(results_dir, netcdf_filename, year, season, ds):
    # in umol/l
    # 1 ml/l of O2 is approximately 43.570 µmol/kg
    # (assumes a molar volume of O2 of 22.392 l/mole and
    # a constant seawater potential density of 1025 kg/m3).
    # https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
    hypox = 90
    anox = 9
    unit = 'umol/l'
    # Select a time index from year and depth level from show_depth
    # choose depth
    show_depth = 70
    show_depth_1 = 100
    show_depth_2 = 125
    show_depth_3 = 150
    
    ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')

    print(f'producing maps for year {year} {season} at {show_depth} m')
    # Plot the data on a map
    plt.style.use('dark_background')

    # plot of areas at thresholds
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 4))
    # Adjust the spacing between subplots
    fig.tight_layout()

    sub_plot_error_area_at_threshold_basemap(ds, parameter='Relerr_per_grid_at_min_90_depth', axis=axs[0, 0], year=year,
                                             vmin=0, vmax=1, threshold=90)
    sub_plot_area_at_threshold_basemap(ds, parameter='Min_depth_90', axis=axs[1, 0], year=year, vmin=60, vmax=180,
                                       threshold=90)
    sub_plot_error_area_at_threshold_basemap(ds, parameter='Relerr_per_grid_at_min_0_depth', axis=axs[0, 1], year=year,
                                             vmin=0, vmax=1, threshold=0)
    sub_plot_area_at_threshold_basemap(ds, parameter='Min_depth_0', axis=axs[1, 1], year=year, vmin=60, vmax=180,
                                       threshold=0)
    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'Hypoxia and anoxia:  {year} {season}', fontsize=8)

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_{season}_areas_{netcdf_filename}.png', dpi=300,
                transparent=True)

    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 4, figsize=(10, 4))
    # Adjust the spacing between subplots
    fig.tight_layout() # (left, bottom, width, height)

    # Create a 1x4 grid of subplots
    #fig, axs = plt.subplots(1, 4, figsize=(18, 4))
    #ax_data = axs[0]
    #ax_error_field = axs[1]
    #ax_min_hypox = axs[2]
    #ax_min_anox = axs[3]
    vmin_o2 = -45
    vmax_o2 = 180+45
    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 0], year=year, show_depth=show_depth, vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 0], year=year, show_depth=show_depth,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 1], year=year, show_depth=show_depth_1,
                                        vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 1], year=year, show_depth=show_depth_1,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 2], year=year, show_depth=show_depth_2,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 2], year=year, show_depth=show_depth_2,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 3], year=year, show_depth=show_depth_3,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 3], year=year, show_depth=show_depth_3,
                                  vmin=0, vmax=0.5)

    # Adjust the spacing between subplots
    # fig.subplots_adjust(hspace=1, wspace = 0)  # You can adjust the value of hspace as needed


    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'Hypoxia and anoxia:  {year} {season}', fontsize=8)

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_{season}_{show_depth}m_{netcdf_filename}.png', dpi = 300, transparent=True)

## extract values that are within our limits, save to a new variable and nc-file. ####

def read_processed_nc(results_dir,file_list,year_list: json):
    for netcdf_filename in file_list:
        ds = xr.open_dataset(f"{results_dir}nc/processed/{netcdf_filename}")
       
        metadata_list = netcdf_filename.split('_')
        parameter = metadata_list[0]
        year_range = metadata_list[1]
        season = metadata_list[2]
        epsilon = metadata_list[3]
        corrlen = metadata_list[4]
        year_list = json.loads(year_list)
        for year in year_list:
            if str(year) not in year_range:
                   continue
            plot(results_dir, netcdf_filename, year, season, ds)

if __name__ == "__main__":
    # Result directory
    results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    # Open the JSON file
    with open(f"{results_dir}file_list.json", 'r') as file:
        # Load JSON data from the file
        file_list = json.load(file)

    year_list = json.dumps([1995, 1996])
    read_processed_nc(results_dir,file_list, year_list)