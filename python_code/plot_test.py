import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import json
# import cartopy
from mpl_toolkits.basemap import Basemap

def sub_plot_observations_basemap(ds, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None):
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
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index]).squeeze()
    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)
    # Add a colorbar
    plt.colorbar(pcm, ax=axis, label='Oxygen umol/l', orientation='horizontal').ax.tick_params(labelsize=10)

    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']

    lon, lat = m(lon, lat)

    m.scatter(lon, lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
              vmin=vmin, vmax=vmax)

    # Add labels to the subplot
    axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')

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
    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)

    # Add labels to the subplot with adjusted text size
    axis.set_title(f'Area <= {threshold} {unit}', fontsize=12)
    axis.set_xlabel('Longitude', fontsize=10)
    axis.set_ylabel('Latitude', fontsize=10)


def sub_plot_errorfields_basemap(ds, axis, year, show_depth, vmin, vmax):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth

    # Create a Basemap instance with Mercator projection
    m = Basemap(projection='merc', llcrnrlat=ds['lat'].min(), urcrnrlat=ds['lat'].max(),
                llcrnrlon=ds['lon'].min(), urcrnrlon=ds['lon'].max(), resolution='l', ax=axis)

    # Plot land borders
    m.drawcoastlines(linewidth=0.5, color='gray')
    
    # plot data
    data = ds['Oxygen_relerr'].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index]).squeeze()

    lon = ds['lon'].values
    lat = ds['lat'].values

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    # Plot data using pcolormesh
    pcm = m.pcolormesh(lon, lat, data, cmap='jet', vmin=vmin, vmax=vmax)
    # Add labels to the 2nd subplot
    axis.set_title(f'Errorfield at {show_depth} m results')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')


"""
def sub_plot_observations_cartopy(ds, axis, year, show_depth, vmin, vmax, observation_span=2, bath_file=None):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    # Plot land borders
    axis.add_feature(cartopy.feature.COASTLINE, edgecolor='gray')
    axis.add_feature(cartopy.feature.BORDERS, linestyle=':', edgecolor='gray')

    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index])
    lon = ds['lon'].values
    lat = ds['lat'].values

    axis.pcolormesh(lon, lat, data, transform=cartopy.crs.PlateCarree(), cmap='jet', vmin=vmin, vmax=vmax)

    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']

    axis.scatter(x=lon, y=lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
                 vmin=vmin, vmax=vmax, transform=cartopy.crs.PlateCarree())

    # Add labels to the subplot
    axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')

    # Set the map projection to SWEREF 99
    axis.set_extent([10, 25, 55, 70], crs=cartopy.crs.PlateCarree())

    # Optionally, add gridlines
    axis.gridlines(draw_labels=True, linestyle='--', color='gray')
"""

def sub_plot_observations(ds, axis, year, show_depth, vmin, vmax, observation_span = 2):

    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
    # Plot land borders
    axis.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5, 0], colors="gray")
    
    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index])
    data.plot(ax=axis, x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth - observation_span) & (df.depth <= show_depth + observation_span))

    observations = df.loc[selection, 'Oxygen_data']
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']
    axis.scatter(x=lon, y=lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
                    vmin=vmin, vmax=vmax)
    # Add labels to the 1st subplot
    axis.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')

def sub_plot_errorfields(ds, axis, year, show_depth, vmin, vmax):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    depth_index = ds["depth"][:].values == show_depth
    # Plot land borders
    axis.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
    
    # plot data
    data = ds['Oxygen_relerr'].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index])
    data.plot(ax=axis, x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
    
    # Add labels to the 2nd subplot
    axis.set_title(f'Errorfield at {show_depth} m results')
    axis.set_xlabel('Longitude')
    axis.set_ylabel('Latitude')

def sub_plot_area_at_threshold(ds, parameter, axis, year, vmin, vmax, threshold, unit = 'umol/l'):
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))

    # Plot land borders
    axis.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
    # Plot data
    try:
        data = ds[parameter].sel(time=ds['time'][time_index])
        data.plot(ax=axis, x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
    except KeyError:
        print(f'no {parameter} in file')
    axis.set_title(f'area of oxygen <= {threshold} {unit}', fontsize=12)
    axis.set_xlabel('Longitude', fontsize=10)
    axis.set_ylabel('Latitude', fontsize=10)


def plot(results_dir, netcdf_filename, year, ds):
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

    depth_index = ds["depth"][:].values == show_depth
    depth_index_1 = ds["depth"][:].values == show_depth_1
    depth_index_2 = ds["depth"][:].values == show_depth_2
    depth_index_3 = ds["depth"][:].values == show_depth_3
    observation_span = 2

    # Extract year and month from the time value
    
    year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
    time_index = year_list.index(str(year))
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
    year_month = datetime.strftime(time_value, '%Y-%m')
    ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')

    print(f'producing maps for year {year} at {show_depth} m {year_month}')
    # Plot the data on a map
    plt.style.use('dark_background')
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 4, figsize=(10, 4))
    ax_data = axs[0, 0]
    ax_error_field = axs[0, 1]
    ax_min_hypox = axs[0, 2]
    ax_min_anox = axs[0, 3]
    ax_data_1 = axs[1, 0]
    ax_data_2 = axs[1, 1]
    ax_data_3 = axs[1, 2]
    ax_data_4 = axs[1, 3]

    # Create a 1x4 grid of subplots
    #fig, axs = plt.subplots(1, 4, figsize=(18, 4))
    #ax_data = axs[0]
    #ax_error_field = axs[1]
    #ax_min_hypox = axs[2]
    #ax_min_anox = axs[3]
    vmin_o2 = -180
    vmax_o2 = 180
    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, axis=ax_data, year=year, show_depth=show_depth, vmin=vmin_o2, vmax=vmax_o2)

    sub_plot_observations_basemap(ds, axis=ax_data_4, year=year, show_depth=show_depth, vmin=vmin_o2, vmax=vmax_o2)

    # 1,0 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, axis=ax_data_1, year=year, show_depth=show_depth_1, vmin=vmin_o2, vmax=vmax_o2)

    # 1,1 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscale
    sub_plot_observations_basemap(ds, axis=ax_data_2, year=year, show_depth=show_depth_2, vmin=vmin_o2, vmax=vmax_o2)

    # 1,2 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscala
    sub_plot_observations_basemap(ds, axis=ax_data_3, year=year, show_depth=show_depth_3, vmin=vmin_o2, vmax=vmax_o2)
    
    # 2222222 Plot the data on the 2nd subplot
    # plot relative error of results at the choosen depth
    sub_plot_errorfields_basemap(ds, axis=ax_error_field, year=year, show_depth=show_depth, vmin=0, vmax=0.5)
    
    # Plot 33333333333 the hypoxic min depth on the 3rd subplot
    # set common max, min limits for depth plots
    sub_plot_area_at_threshold_basemap(ds, parameter='Min_depth_hypoxia', axis=ax_min_hypox, year=year, vmin = 60, vmax = 150, threshold=hypox)

    # Plot the anoxic min depth  on the 4th subplot
    sub_plot_area_at_threshold_basemap(ds, parameter='Min_depth_anoxia', axis=ax_min_anox, year=year, vmin = 60, vmax = 150, threshold=anox)

    # Set a common title for each row of subplots
    """
        Modifications to the Groupers returned by get_shared_x_axes and get_shared_y_axes are deprecated. 
        In the future, these methods will return immutable views on the grouper structures. 
        Note that previously, calling e.g. join() would already fail to set up the correct structures for sharing axes; 
        use Axes.sharex or Axes.sharey instead.
    """
    ax_data.get_shared_x_axes().join(ax_data, ax_error_field)
    ax_data.get_shared_y_axes().join(ax_data, ax_error_field)

    # Adjust the spacing between subplots
    fig.tight_layout(rect=[0, 0, 1, 1]) # (left, bottom, width, height)

    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'maps of hypoxia and anoxia\n{year_month}')

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year_month}_{show_depth}_{netcdf_filename}.png', dpi = 300, transparent=True)

bath_file = xr.open_dataset("../Oxygen_maps/bathymetry/gebco_30sec_4.nc")
# Extract the required variables
b = bath_file["bat"]

## extract values that are within our limits, save to a new variable and nc-file. ####


# # choose year
start_year = 1995
stop_year = 1995
def read_processed_nc(results_dir,file_list,year_list):
    for netcdf_filename in file_list:

        print(netcdf_filename)
        ds = xr.open_dataset(f"{results_dir}nc/processed/{netcdf_filename}")
        """for year in year_list:
            print(year)
            plot(results_dir, netcdf_filename, year, ds)"""

        for year in range(start_year,stop_year+1):
            plot(results_dir, netcdf_filename, year, ds)

if __name__ == "__main__":
    # Result directory
    results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    # Open the JSON file
    with open(f"{results_dir}file_list.json", 'r') as file:
        # Load JSON data from the file
        file_list = json.load(file)

    year_list = json.dumps([1995])
    print(year_list)
    read_processed_nc(results_dir,file_list, year_list)