import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from datetime import datetime
import pandas as pd
import json
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import inset_locator, make_axes_locatable
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['hatch.linewidth'] = 0.1

def create_custom_colormap(levels, colors):
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(levels, cmap.N)

    return cmap, norm

def set_up_basemap(ds, axis):
    # Create a Basemap instance
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Rita kustlinjer
    m.drawcoastlines(linewidth=0.5, color='gray')

    # Rita meridian- och parallelgränser med anpassade intervall
    m.drawmeridians(np.arange(9, 31, 4), labels=[True, False, False, True], linewidth=0.1, fontsize=3)
    m.drawparallels(np.arange(54, 67, 1), labels=[True, False, False, True], linewidth=0.1, fontsize=3)

    return m

def sub_plot_parameter_basemap(ds, parameter, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None, time_index=0):
    
    depth_index = ds["depth"][:].values == show_depth

    # Create a Basemap instance
    m = set_up_basemap(ds, axis)

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
                                  colorbar = True, color = 'k', time_index=0, BG=False):

    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
    start_year = ds.attrs["start year"]
    end_year = ds.attrs["end year"]
    if not (int(start_year) <= int(year) <= int(end_year)):
        print(f"WARNING: {year} to plot not in range {start_dt.year}-{end_dt.year}")
        
    # Convert start/end years to pandas Timestamps
    start_dt = pd.Timestamp(f"{int(start_year)}-01-01")
    end_dt = pd.Timestamp(f"{int(end_year)}-12-31")
    # If we are plotting a background field we want observations from all years used
    # If it is not a background field we want observations from the year of the analysis
    # after the change to make one netcdf per analysed year, 
    if not BG:
        year = time_value.year
        start_dt = pd.Timestamp(f"{year}-01-01")
        end_dt = pd.Timestamp(f"{year}-12-31")
    # Plocka ut observationer
    # OBS att DIVAnd resultatet ligger under dimensionen 'Oxygen' i datasetet och obsevrationerna under 'Oxygen_data'
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = (
        (df.obsyear >= start_dt) &
        (df.obsyear <= end_dt) &
        (df.depth >= show_depth - observation_span) &
        (df.depth <= show_depth + observation_span)
    )

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

        m.scatter(lon, lat, s=5, c=observations, cmap=cmap, edgecolors='k', linewidth=0.2, alpha=0.5,
              vmin=vmin, vmax=vmax)
        
        # Add a colorbar
        # Create an inset_axes for the colorbar
        cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
        cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
        cbar.ax.tick_params(labelsize = 4)

        # Modify the colormap levels to control the step length
        levels = np.arange(vmin, vmax+1, 45)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        pcm.set_norm(norm)
        # Set the colorbar levels explicitly
        cbar.set_ticks(levels)
        axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m', fontsize=8)

    else:
        m.scatter(lon, lat, s=5, edgecolors='k', linewidth=0.2, facecolor=color)

def sub_plot_observations_basemap(ds, parameter, axis, year, 
                                  show_depth = 0, vmin = 0, vmax = 180, observation_span=2, 
                                  colorbar = True, color = 'k', time_index=0, BG = False):
    
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
    start_year = ds.attrs["start year"]
    end_year = ds.attrs["end year"]
    if not (int(start_year) <= int(year) <= int(end_year)):
        print(f"WARNING: {year} to plot not in range {start_dt.year}-{end_dt.year}")
        
    # Convert start/end years to pandas Timestamps
    start_dt = pd.Timestamp(f"{int(start_year)}-01-01")
    end_dt = pd.Timestamp(f"{int(end_year)}-12-31")
    # If we are plotting a background field we want observations from all years used
    # If it is not a background field we want observations from the year of the analysis
    # after the change to make one netcdf per analysed year, 
    if not BG:
        year = time_value.year
        start_dt = pd.Timestamp(f"{year}-01-01")
        end_dt = pd.Timestamp(f"{year}-12-31")    

    # Plocka ut observationer
    # OBS att DIVAnd resultatet ligger under dimensionen 'Oxygen' i datasetet och obsevrationerna under 'Oxygen_data'
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = (
        (df.obsyear >= start_dt) &
        (df.obsyear <= end_dt) &
        (df.depth >= show_depth - observation_span) &
        (df.depth <= show_depth + observation_span)
    )

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

        m.scatter(lon, lat, s=5, c=observations, cmap=cmap, edgecolors='k', linewidth=0.2, alpha=0.5,
              vmin=vmin, vmax=vmax)
        
        # Add a colorbar
        # Create an inset_axes for the colorbar
        cbax = inset_locator.inset_axes(axis, width="40%", height="3%", loc="lower right", bbox_to_anchor=(0, 0.15, 1, 1),
                                    bbox_transform=axis.transAxes)
        cbar = plt.colorbar(pcm, cax = cbax,  orientation = 'horizontal')
        cbar.ax.tick_params(labelsize = 3)

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

def sub_plot_area_at_threshold_basemap(ds, parameter, axis, year, threshold, vmin = 0, vmax = 180, unit='umol/l', bath_file=None, colorbar=True, color = 'k', levels=[0.5, 1.5], hatches = [], time_index=0):

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
        pcm = m.contourf(lon, lat, data,  levels=levels, colors=[color], hatches = hatches)

def sub_plot_error_area_at_threshold_basemap(ds, parameter, axis, year, vmin, vmax, threshold, unit='umol/l', bath_file=None, time_index=0):

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

    # Add a colorbar, Create an inset_axes for the colorbar
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

def plot(results_dir, netcdf_filename, year, season, ds, threshold_list, interval, time_index, BG):
    # 1 ml/l of O2 is approximately 43.570 µmol/kg
    # (assumes a molar volume of O2 of 22.392 l/mole and
    # a constant seawater potential density of 1025 kg/m3).
    # https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
    unit = 'umol/l'
    print(f"now plotting from: {netcdf_filename}")
    ############### FIRST plot, threshold results #################

    # Row 1: maps with error fields for each threshol
    # Row 2: maps of the depth of the onset of each threshold
    # Create a 2x2 grid of subplots
    n_figs = len(threshold_list)
    
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
        if season == "Winter":
            fig.suptitle(f"Hypoxia and anoxia: {'BG' if BG else f'{int(year)-1}-{int(year)}'} {season}", fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                         verticalalignment='bottom')
        else:
            fig.suptitle(f"Hypoxia and anoxia: {'BG' if BG else f'{int(year)}'} {season}", fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='bottom')
    else:
        if season == "Winter":
            fig.suptitle(f"{'BG' if BG else f'{int(year-1)}-{int(year)}'} {season}", fontsize=8, x=0.5, y=0.53, horizontalalignment='center',
                         verticalalignment='center')
        else:
            fig.suptitle(f"{'BG' if BG else f'{int(year)}'} {season}", fontsize=8, x=0.5, y=0.53, horizontalalignment='center',
                     verticalalignment='center')

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG_threshold_result{str(interval).replace(", ", "_")}_{season}.png', dpi=300,
                    transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/threshold_result{year}_{season}.png', dpi=300,
                    transparent=False)
        plt.close()

    ################### SECOND plot surface layer results ############################

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

    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    subplot_no = 0
    for show_depth in [10, 20, 30, 50]:
        try:
            sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, subplot_no], year=year, show_depth=show_depth, vmin=vmin_o2,
                                        vmax=vmax_o2, BG=BG)

            sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, subplot_no], year=year, show_depth=show_depth,
                                        vmin=0, vmax=0.5)
            subplot_no += 1
        except ValueError:
            continue

    # Add title and labels
    # Set the title for the whole figure
    if season == "Winter":
        fig.suptitle(f'{int(year)-1}-{int(year)} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                     verticalalignment='top')
    else:
        fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG_surf_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/surf_{year}_{season}.png', dpi=300, transparent=False)
        plt.close()

    ######################## THIRD plot halocline results #############################

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

    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    subplot_no = 0
    for show_depth in [60, 70, 80, 90]:
        try:
            sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, subplot_no], year=year, show_depth=show_depth, vmin=vmin_o2, vmax=vmax_o2, BG=BG)

            sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, subplot_no], year=year, show_depth=show_depth,
                                        vmin=0, vmax=0.5)
            subplot_no += 1
        except ValueError:
            continue

    # Add title and labels
    # Set the title for the whole figure
    if season == "Winter":
        fig.suptitle(f'{int(year)-1}-{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                     verticalalignment='top')
    else:
        fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG_halo_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/halo_{year}_{season}.png', dpi = 300, transparent=False)
        plt.close()

    ################# FOURTH plot results in the deep water ########################

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

    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    subplot_no = 0
    skipped_depth = []
    for show_depth in [100, 110, 125, 150]:
        try:
            sub_plot_observations_basemap(ds, parameter='Oxygen', axis=axs[0, subplot_no], year=year, show_depth=show_depth, vmin=vmin_o2,
                                    vmax=vmax_o2, BG=BG)

            sub_plot_errorfields_basemap(ds, parameter='Oxygen_relerr', axis=axs[1, subplot_no], year=year, show_depth=show_depth,
                                    vmin=0, vmax=0.5)
            subplot_no += 1

        except ValueError:
            print("Its not deep enough skip this depth...")
            skipped_depth.append(show_depth)
            continue

    if skipped_depth != [100, 110, 125, 150]:
        # Add title and labels
        # Set the title for the whole figure
        if season == "Winter":
            fig.suptitle(f'{int(year) - 1}-{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                         verticalalignment='top')
        else:
            fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center', verticalalignment='top')

        # Save the plot
        if BG:
            plt.savefig(f'{results_dir}/figures/BG_deep_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
            plt.close()
        else:
            plt.savefig(f'{results_dir}/figures/deep_{year}_{season}.png', dpi=300, transparent=False)
            plt.close()

    ################# FIFTH plot final result map  ###########################

    # plots of results all observations and hypox area and with anox area overlayed
    fig, axs = plt.subplots(1, 1, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    fig.tight_layout()  # (left, bottom, width, height)

    # Vänder på threshold_list för att högst threshold skall hamna underst.
    threshold_list.reverse()
    color_list = ['lightgrey', 'darkgrey', 'grey']
    hatches_list = [10 * '/', 10 * "\\",10*'|']

    for index, threshold in enumerate(threshold_list):
        #for threshold, index in threshold_list:
        sub_plot_area_at_threshold_basemap(ds, parameter=f'{threshold}_mask_firstlayer', axis=axs, year=year, colorbar=False, color = color_list[index],
                                       threshold=threshold)

    for index, threshold in enumerate(threshold_list):
        sub_plot_area_at_threshold_basemap(ds, parameter=f'Relerr_per_grid_at_min_{threshold}_depth', axis=axs, year=year, colorbar=False, color = 'none', hatches=[hatches_list[index]],threshold=threshold)

    sub_plot_only_observations(ds, axis=axs, year=year, colorbar=False, color = 'r',observation_span = 500, BG=BG)

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
        if season == "Winter":
            fig.suptitle(f'{int(year) - 1}-{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                         verticalalignment='top')
        else:
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
                                 markersize=4, label='Observations', markeredgewidth=0.5)
        # Lägg till legenden
        axs.legend(handles=patches + [fake_marker], loc='lower right', fontsize=6)
        # Add title and labels
        # Set the title for the whole figure
        if season == "Winter":
            fig.suptitle(f'{year - 1}-{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                         verticalalignment='top')
        else:
            fig.suptitle(f'{year} {season}', fontsize=8, x=0.5, y=1.0, horizontalalignment='center',
                     verticalalignment='top')
        if len(threshold_list) < 3:
            print("Bara två thresholds, ändra legenden!")

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG'
                    f'_result_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/final_result_{year}_{season}.png', dpi=300, transparent=False)
        plt.close()

## extract values that are within our limits, save to a new variable and nc-file. ####

def read_processed_nc(results_dir, file_list):
    
    for netcdf_filename in file_list:
        BG = False
        if 'Background' in netcdf_filename:
            BG = True
        
        ds = xr.open_dataset(f"{results_dir}/processed/{netcdf_filename}", engine='h5netcdf')
        
        ds_year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
        if not len(ds_year_list) == 1:
            print('Skipping file due to multiple years in files')
            continue
        ds_year = ds_year_list[0]    
        time_index = 0
        season = ds.attrs['season']
        epsilon = ds.attrs['epsilon']
        threshold_list = ds.attrs['threshold_list']
        # När attribute till nc sätts blir det av type Any, läses som en str från ds.attrs
        threshold_list = eval(threshold_list.replace("Any", ""))
        corrlen = ds.attrs['horizontal correlation length m']

        print(f"year in netcdf: {str(ds_year)}")
        interval = [int(ds_year) -1, int(ds_year) +1]   #Background year +/- 1
        plot(results_dir, netcdf_filename, ds_year, season, ds, threshold_list, interval, time_index, BG)

if __name__ == "__main__":
    print("running")
    # Result directory
    results_dir = "C:/LenaV/code/DIVAnd/resultat/"
    results_dir = "C:/Work/DIVAnd/Oxygen_maps/resultat/Baltic_Proper/20250613_0959/"

    file_list = ["Oxygen_2015-2015_Autumn_0.2_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz.nc",
                 "Oxygen_2015-2015_Spring_0.2_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz.nc",
                 "Oxygen_2015-2015_Summer_0.2_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz.nc",
                 "Oxygen_2015-2015_Winter_0.2_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked_varcorrlenz.nc",
                 "Background_Oxygen_10_year_Autumn_0.1_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked.nc",
                 "Background_Oxygen_10_year_Spring_0.1_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked.nc",
                 "Background_Oxygen_10_year_Summer_0.1_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked.nc",
                 "Background_Oxygen_10_year_Winter_0.1_80000.0_0.05_5.0_2.0_bat_elevation_Baltic_Sea_masked.nc"]

    year_list = json.dumps([2015])
    yearlist_background = json.dumps([[1960, 1969], [1970, 1979], [1980, 1989], [1990, 1999], [2000, 2009], [2010, 2019], [2020, 2024]])
    read_processed_nc(results_dir, file_list, year_list, yearlist_background)