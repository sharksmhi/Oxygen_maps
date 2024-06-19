import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from datetime import datetime
import pandas as pd
import json
# import cartopy
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import inset_locator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['hatch.linewidth'] = 0.1

def create_custom_colormap(levels, colors):
    # Define color ranges and corresponding colors
    #ranges = np.linspace(0, 1, len(levels) + 1)
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(levels, cmap.N)

    return cmap, norm

def set_up_basemap(ds, axis):

    # Create a Basemap instance
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Rita kustlinjer
    m.drawcoastlines(linewidth=0.5, color='gray')

    # Rita kontinenter, men undvik att fylla sjöar
    #m.fillcontinents(color='white', lake_color=None, zorder=1)

    # Rita meridian- och parallelgränser med anpassade intervall
    m.drawmeridians(np.arange(9, 31, 4), labels=[True, False, False, True], linewidth=0.1, fontsize=3)
    m.drawparallels(np.arange(54, 67, 1), labels=[True, False, False, True], linewidth=0.1, fontsize=3)

    return m

def sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    # Create a Basemap instance
    m = set_up_basemap(ds, axis)
    # Plot land borders from the bathymetry file

    # Plot data
    data = ds[parameter].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index]).squeeze()
    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)

    return m, pcm

def sub_plot_only_observations(ds, axis, year, 
                                  show_depth = 0, vmin = 0, vmax = 180, observation_span=2, 
                                  colorbar = True, color = 'k'):
    
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    # Plocka ut observationer
    # OBS att DIVAnd resultatet ligger under dimensionen 'Oxygen' i datasetet och obsevrationerna under 'Oxygen_data'
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']

    m = set_up_basemap(ds, axis)
    lon, lat = m(lon, lat)

    if colorbar:
        levels = np.arange(vmin, vmax, 45)
        colors = ['black', 'brown', 'red', 'orange', 'yellow', 'green']
        # Create the custom colormap
        cmap, norm = create_custom_colormap(levels, colors)
        pcm = m.pcolormesh(lon, lat, observations, cmap='jet', vmin=vmin, vmax=vmax)

        # Change the colormap after the plot is created
        pcm.set_cmap(cmap)
        pcm.set_norm(norm)
        # Change vmin and vmax after the plot is created
        pcm.set_clim(vmin=vmin, vmax=vmax)

        m.scatter(lon, lat, s=5, c=observations, cmap=cmap, edgecolors='k', linewidth=0.2, facecolor='none',
              vmin=vmin, vmax=vmax)
        
        # Add a colorbar
        # Create an inset_axes for the colorbar
        cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
        cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
        cbar.ax.tick_params(labelsize = 4)
        # cbar.set_label(label='µmol/l', fontsize = 10,  y=1.05)

        # Modify the colormap levels to control the step length
        levels = np.arange(vmin, vmax+1, 45)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        pcm.set_norm(norm)
        # Set the colorbar levels explicitly
        cbar.set_ticks(levels)
        axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m', fontsize=8)
        # Add labels to the subplot
        #axis.set_xlabel('Longitude', fontsize=8)
        #axis.set_ylabel('Latitude', fontsize=8)
    else:
        m.scatter(lon, lat, s=5, edgecolors='k', linewidth=0.2, facecolor=color)

def sub_plot_observations_basemap(ds, parameter, axis, year, 
                                  show_depth = 0, vmin = 0, vmax = 180, observation_span=2, 
                                  colorbar = True, color = 'k'):
    
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    # Plocka ut observationer
    # OBS att DIVAnd resultatet ligger under dimensionen 'Oxygen' i datasetet och obsevrationerna under 'Oxygen_data'
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']

    m, pcm = sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax)
    lon, lat = m(lon, lat)

    if colorbar:
        levels = np.arange(vmin, vmax, 45)
        colors = ['black', 'brown', 'red', 'orange', 'yellow', 'green']
        # Create the custom colormap
        cmap, norm = create_custom_colormap(levels, colors)

        # Change the colormap after the plot is created
        pcm.set_cmap(cmap)
        pcm.set_norm(norm)
        # Change vmin and vmax after the plot is created
        pcm.set_clim(vmin=vmin, vmax=vmax)

        m.scatter(lon, lat, s=5, c=observations, cmap=cmap, edgecolors='k', linewidth=0.2, facecolor='none',
              vmin=vmin, vmax=vmax)
        
        # Add a colorbar
        # Create an inset_axes for the colorbar
        cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
        cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
        cbar.ax.tick_params(labelsize = 3)
        # cbar.set_label(label='µmol/l', fontsize = 10,  y=1.05)

        # Modify the colormap levels to control the step length
        levels = np.arange(vmin, vmax+1, 45)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        pcm.set_norm(norm)
        # Set the colorbar levels explicitly
        cbar.set_ticks(levels)
    else:
        m.scatter(lon, lat, s=5, edgecolors='k', linewidth=0.2, facecolor=color)

    # Add labels to the subplot
    axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m [µmol/l]', fontsize = 4)

    return pcm

def sub_plot_area_at_threshold_basemap(ds, parameter, axis, year, threshold, vmin = 0, vmax = 180, unit='umol/l', bath_file=None, colorbar=True, color = 'k', levels=[0.5, 1.5], hatches = []):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))

    # Create a Basemap instance with Mercator projection
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')
    m.drawmeridians(np.arange(9, 31, 4), labels=[True, False, False, True], linewidth=0.1, fontsize=3)
    m.drawparallels(np.arange(54, 67, 1), labels=[True, False, False, True], linewidth=0.1, fontsize=3)

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
    patches =[]
    if colorbar:
        pcm = m.pcolormesh(lon, lat, data, cmap='ocean', vmin=vmin, vmax=vmax)
        # Add a colorbar
        # Create an inset_axes for the colorbar
        cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
        cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
        cbar.ax.tick_params(labelsize = 6)
        #cbar.ax.patch.set_facecolor('white')
        axis.set_title(f'Area <= {threshold} {unit}', fontsize=6)
    else:
        # pcm = m.pcolormesh(lon, lat, data, color=color)
        pcm = m.contourf(lon, lat, data,  levels=levels, colors=[color], hatches = hatches)

def sub_plot_error_area_at_threshold_basemap(ds, parameter, axis, year, vmin, vmax, threshold, unit='umol/l', bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))

    # Create a Basemap instance with Mercator projection
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')
    m.drawmeridians(np.arange(9, 31, 4), labels=[True, False, False, True], linewidth=0.1, fontsize=3)
    m.drawparallels(np.arange(54, 67, 1), labels=[True, False, False, True], linewidth=0.1, fontsize=3)

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

    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax=cbax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)

    # Modify the colormap levels to control the step length
    levels = [0,0.3,0.5,1]
    # Create a custom discrete colormap
    #        0       0.3       0.5
    colors = ['green', 'orange', 'red']
    cmap, norm = create_custom_colormap(levels, colors)

    # Change the colormap after the plot is created
    pcm.set_cmap(cmap)
    pcm.set_norm(norm)

    cbar.set_ticks(levels)

    # Add labels to the subplot with adjusted text size
    axis.set_title(f'Error at <= {threshold} {unit} depth', fontsize=6)

def sub_plot_errorfields_basemap(ds, parameter, axis, year, show_depth, vmin, vmax):

    m, pcm = sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax)

    # Add a colorbar
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
    cbar = plt.colorbar(pcm, cax=cbax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)

    # Modify the colormap levels to control the step length
    levels = [0,0.3,0.5,1]
    # Create a custom discrete colormap
    #        0       0.3       0.5
    colors = ['green', 'orange', 'red']

    cmap, norm = create_custom_colormap(levels, colors)

    # Change the colormap after the plot is created
    pcm.set_cmap(cmap)
    pcm.set_norm(norm)

    # Change vmin and vmax after the plot is created
    pcm.set_clim(vmin=vmin, vmax=vmax)
    cbar.set_ticks(levels)
    axis.set_title(f'Errorfield at {show_depth} m', fontsize=6)

def plot(results_dir, netcdf_filename, year, season, ds, threshold_list):
    # 1 ml/l of O2 is approximately 43.570 µmol/kg
    # (assumes a molar volume of O2 of 22.392 l/mole and
    # a constant seawater potential density of 1025 kg/m3).
    # https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
    unit = 'umol/l'
    #Fullt sätt att fixa till den sträng som threshold_list verkar ha blivit när det var inne en sväng i julia.
    threshold_list = eval(threshold_list.replace("Any", ""))

    n_figs = len(threshold_list)
    print(n_figs)
    # plot of areas at thresholds
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, n_figs, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    fig.tight_layout()

    for index, threshold in enumerate(threshold_list):
        # Select a time index from year and depth level from show_depth
        ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')

        sub_plot_error_area_at_threshold_basemap(ds, parameter=f'Relerr_per_grid_at_min_{threshold}_depth', axis=axs[0, index], year=year,
                                                 vmin=0, vmax=1, threshold=threshold)
        sub_plot_area_at_threshold_basemap(ds, parameter=f'Min_depth_{threshold}', axis=axs[1, index], year=year, threshold=threshold, vmin=50, vmax=150,)

    # Add title and labels
    # Set the title for the whole figure
    if 0 in threshold_list:
        fig.suptitle(f'Hypoxia and anoxia:  {year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='bottom')
    else:
        fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=0.53, horizontalalignment='center',
                     verticalalignment='center')
    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_areas_{netcdf_filename}.png', dpi=300,
                transparent=False)
    plt.close()

    # plots of results at 4 different depths 10, 40, 50, 60
    fig, axs = plt.subplots(2, 4, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    fig.tight_layout()  # (left, bottom, width, height)
    if 0 in threshold_list:
        vmin_o2 = -45
        vmax_o2 = 180 + 45
    else:
        vmin_o2 = 90
        vmax_o2 = 360

    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 0], year=year, show_depth=10, vmin=vmin_o2,
                                  vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 0], year=year, show_depth=10,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 1], year=year, show_depth=40,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 1], year=year, show_depth=40,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 2], year=year, show_depth=50,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 2], year=year, show_depth=50,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 3], year=year, show_depth=60,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 3], year=year, show_depth=60,
                                 vmin=0, vmax=0.5)

    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_surf_{netcdf_filename}.png', dpi=300, transparent=False)
    plt.close()

    # plots of results at 4 different depths 60, 70, 80, 90
    fig, axs = plt.subplots(2, 4, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    fig.tight_layout() # (left, bottom, width, height)

    if 0 in threshold_list:
        vmin_o2 = -45
        vmax_o2 = 180 + 45
    else:
        vmin_o2 = 90
        vmax_o2 = 360

    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 0], year=year, show_depth=60, vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 0], year=year, show_depth=60,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 1], year=year, show_depth=70,
                                        vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 1], year=year, show_depth=70,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 2], year=year, show_depth=80,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 2], year=year, show_depth=80,
                                  vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 3], year=year, show_depth=90,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 3], year=year, show_depth=90,
                                  vmin=0, vmax=0.5)

    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_halo_{netcdf_filename}.png', dpi = 300, transparent=False)
    plt.close()

    # plots of results at 4 different depths 100, 110, 125, 150
    fig, axs = plt.subplots(2, 4, figsize=(10, 4.5)) #4
    # Adjust the spacing between subplots
    fig.tight_layout()  # (left, bottom, width, height)

    if 0 in threshold_list:
        vmin_o2 = -45
        vmax_o2 = 180 + 45
    else:
        vmin_o2 = 90
        vmax_o2 = 360
    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 0], year=year, show_depth=100, vmin=vmin_o2,
                                  vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 0], year=year, show_depth=100,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 1], year=year, show_depth=110,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 1], year=year, show_depth=110,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 2], year=year, show_depth=125,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 2], year=year, show_depth=125,
                                 vmin=0, vmax=0.5)

    sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, 3], year=year, show_depth=150,
                                  vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, 3], year=year, show_depth=150,
                                 vmin=0, vmax=0.5)

    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_deep_{netcdf_filename}.png', dpi=300, transparent=False)
    plt.close()


    # plots of results all observations and hypox area and with anox area overlayed
    fig, axs = plt.subplots(1, 1, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    fig.tight_layout()  # (left, bottom, width, height)

    # Vänder på threshold_list för att högst threshold skall hamna underst.
    threshold_list.reverse()
    color_list = ['lightgrey', 'grey', '#303030']
    hatches_list = [25 * '/', 25 * "\\",25*'|']

    for index, threshold in enumerate(threshold_list):
        #for threshold, index in threshold_list:
        sub_plot_area_at_threshold_basemap(ds, parameter=f'{threshold}_mask_firstlayer', axis=axs, year=year, colorbar=False, color = color_list[index],
                                       threshold=threshold)

    for index, threshold in enumerate(threshold_list):
        sub_plot_area_at_threshold_basemap(ds, parameter=f'Relerr_per_grid_at_min_{threshold}_depth', axis=axs, year=year, colorbar=False, color = 'none', hatches=hatches_list[index],threshold=threshold)
        print(hatches_list[index])

    sub_plot_only_observations(ds, axis=axs, year=year, colorbar=False, color = 'r',observation_span = 500)

    # Lägg till en "fejk" legend om syrefritt är med
    if 0 in threshold_list:
        fake_labels = [f'<{threshold_list[0]} µmol/l', f'<{threshold_list[1]} µmol/l', f'<{threshold_list[2]} µmol/l',
                       f'Error field <{threshold_list[0]} µmol/l', f'Error field <{threshold_list[1]} µmol/l',
                       f'Error field <{threshold_list[2]} µmol/l']
        fake_colors = [color_list[0], color_list[1], color_list[2],'none', 'none','none']
        fake_hatches = ['', '','', hatches_list[0], hatches_list[1],hatches_list[2]]
        # Skapa proxy-objekt för legenden
        patches = [mpatches.Patch(facecolor=color, hatch=hatch, label=label)
                   for color, hatch, label in zip(fake_colors, fake_hatches, fake_labels)]

        # Lägg till en "fejk" röd marker 'o'
        fake_marker = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markeredgecolor='k', markersize=2, label='Observations')
        # Lägg till legenden
        axs.legend(handles=patches + [fake_marker], loc='lower right', fontsize=7)
        # Add title and labels
        # Set the title for the whole figure
        fig.suptitle(f'Hypoxia and anoxia:  {year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                     verticalalignment='top')
    else: #om bottniska viken.
        fake_labels = [f'<{threshold_list[0]} µmol/l', f'<{threshold_list[1]} µmol/l', f'<{threshold_list[2]} µmol/l',
                       f'Error field <{threshold_list[0]} µmol/l', f'Error field <{threshold_list[1]} µmol/l', f'Error field <{threshold_list[2]} µmol/l']
        fake_colors = [color_list[0], color_list[1], color_list[2],'none', 'none','none']
        fake_hatches = ['', '','', hatches_list[0], hatches_list[1],hatches_list[2]]
        # Skapa proxy-objekt för legenden
        patches = [mpatches.Patch(facecolor=color, hatch=hatch, label=label)
                   for color, hatch, label in zip(fake_colors, fake_hatches, fake_labels)]

        # Lägg till en "fejk" röd marker 'o'
        fake_marker = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markeredgecolor='k',
                                 markersize=3, label='Observations')
        # Lägg till legenden
        axs.legend(handles=patches + [fake_marker], loc='lower right', fontsize=6)
        # Add title and labels
        # Set the title for the whole figure
        fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                     verticalalignment='top')
        if len(threshold_list) < 3:
            print("Bara två thresholds, ändra legenden!")
    #fake_marker =['none','none','none','none','o',]

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year}_overview_{netcdf_filename}.png', dpi=300, transparent=False)
    plt.close()

## extract values that are within our limits, save to a new variable and nc-file. ####

def read_processed_nc(results_dir,file_list,year_list: json):
    year_list = json.loads(year_list)

    for netcdf_filename in file_list:
        print(f'plto from {netcdf_filename}')
        ds = xr.open_dataset(f"{results_dir}/nc/processed/{netcdf_filename}", engine='h5netcdf')
        season = ds.attrs['season']
        epsilon = ds.attrs['epsilon']
        threshold_list = ds.attrs['threshold_list']
        corrlen = ds.attrs['horizontal correlation length m']
        start_year = ds.attrs['start year']
        end_year = ds.attrs['end year']

        print(repr(ds.attrs['threshold_list']))
        repr(ds.attrs['threshold_list'])

        for year in year_list:
            # str(year) testa om det inte funkar.
            if year not in range(int(start_year), int(end_year)+1):
                    pass
            print(f'plotting {year}')
            plot(results_dir, netcdf_filename, year, season, ds, threshold_list)

if __name__ == "__main__":
    print("running")
    # Result directory
    results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    results_dir = "C:/LenaV/code/DIVAnd/resultat/"
    results_dir = "C:/Work/DIVAnd/Oxygen_maps/resultat/Bothnian_Bay/"

    file_list = ["Oxygen_2000-2022_Summer_0.2_80000_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz_NBK.nc"]
    year_list = json.dumps([2001])
    ##ear_list = json.dumps([1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978,
    # 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997,
    # 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
    # 2017, 2018, 2019, 2020, 2021]);
    read_processed_nc(results_dir, file_list, year_list)