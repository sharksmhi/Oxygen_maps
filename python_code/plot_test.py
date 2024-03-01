import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import json

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

    # 1111111 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscala
    vmin_o2 = -180
    vmax_o2 = 180
    # Plot land borders
    ax_data.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index])
    data.plot(ax=ax_data, x='lon', y='lat', cmap='jet', vmin=vmin_o2, vmax=vmax_o2)
    # data = ds['Oxygen_data'].sel(obsdepth=ds['Oxygen_data'][40])
    # data.plot(ax=ax_data, x='lon', y='lat', cmap='jet')
    df = pd.DataFrame({'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'], 'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (df.depth >= show_depth-observation_span) & (df.depth <= show_depth+observation_span))
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']
    observations = df.loc[selection, 'Oxygen_data']
    # obs_val = df.loc[df.obsyear == datetime.strftime(time_value, '%Y'), 'Oxygen_data']
    # lat = xr.where(ds['obs_year']==datetime.strftime(time_value, '%Y'), ds['obslat'], np.nan)

    ax_data.scatter(x=lon, y=lat, s = 5, c = observations, cmap = 'jet', edgecolors = 'k', linewidth = 0.2, facecolor = 'none', vmin=vmin_o2, vmax=vmax_o2)
    # Add labels to the 1st subplot
    ax_data.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m')
    ax_data.set_xlabel('Longitude')
    ax_data.set_ylabel('Latitude')

    # 1,0 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscala

    vmin_o2 = -180
    vmax_o2 = 180
    # Plot land borders
    ax_data_1.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5, 0], colors="gray")
    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index_1])
    data.plot(ax=ax_data_1, x='lon', y='lat', cmap='jet', vmin=vmin_o2, vmax=vmax_o2)

    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
                df.depth >= show_depth_1 - observation_span) & (df.depth <= show_depth_1 + observation_span))
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']
    observations = df.loc[selection, 'Oxygen_data']

    ax_data_1.scatter(x=lon, y=lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
                    vmin=vmin_o2, vmax=vmax_o2)
    # Add labels to the 1st subplot
    ax_data_1.set_title(f'Oxygen at {show_depth_1} m\nobservation at +/- {observation_span} m')
    ax_data_1.set_xlabel('Longitude')
    ax_data_1.set_ylabel('Latitude')

    # 1,1 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscala

    vmin_o2 = -180
    vmax_o2 = 180
    # Plot land borders
    ax_data_2.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5, 0], colors="gray")
    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index_2])
    data.plot(ax=ax_data_2, x='lon', y='lat', cmap='jet', vmin=vmin_o2, vmax=vmax_o2)
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth_2 - observation_span) & (df.depth <= show_depth_2 + observation_span))
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']
    observations = df.loc[selection, 'Oxygen_data']


    ax_data_2.scatter(x=lon, y=lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
                    vmin=vmin_o2, vmax=vmax_o2)
    # Add labels to the 1st subplot
    ax_data_2.set_title(f'Oxygen at {show_depth_2} m\nobservation at +/- {observation_span} m')
    ax_data_2.set_xlabel('Longitude')
    ax_data_2.set_ylabel('Latitude')

    # 1,2 Plot the data on the 1st subplot
    # on the 1st and 2nd plot we show oxygen set min and max for colorscala

    vmin_o2 = -180
    vmax_o2 = 180
    # Plot land borders
    ax_data_3.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5, 0], colors="gray")
    # Plot data
    data = ds['Oxygen'].sel(time=ds['time'][time_index], depth=ds['depth'][depth_index_3])
    data.plot(ax=ax_data_3, x='lon', y='lat', cmap='jet', vmin=vmin_o2, vmax=vmax_o2)
    df = pd.DataFrame(
        {'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'],
         'depth': ds['obsdepth']})
    selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (
            df.depth >= show_depth_3 - observation_span) & (df.depth <= show_depth_3 + observation_span))
    lon = df.loc[selection, 'obslon']
    lat = df.loc[selection, 'obslat']
    observations = df.loc[selection, 'Oxygen_data']

    ax_data_3.scatter(x=lon, y=lat, s=5, c=observations, cmap='jet', edgecolors='k', linewidth=0.2, facecolor='none',
                    vmin=vmin_o2, vmax=vmax_o2)
    # Add labels to the 1st subplot
    ax_data_3.set_title(f'Oxygen at {show_depth_3} m\nobservation at +/- {observation_span} m')
    ax_data_3.set_xlabel('Longitude')
    ax_data_3.set_ylabel('Latitude')


    # 2222222 Plot the data on the 2nd subplot
    # plot relative error of results at the choosen depth
    data = ds['Oxygen_relerr'].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index])
    data.plot(ax=ax_error_field, x='lon', y='lat', cmap='jet', vmin=0, vmax=0.5)
    # Plot the relerrfield of the hypoxic depth layer
    # data = ds['Hypoxic_relerr_per_grid'].sel(time=ds['time'][time_index])
    # data.plot(ax=ax_error_field, x='lon', y='lat', cmap='jet', vmin=0, vmax=0.5)
    # Plot land borders
    ax_error_field.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")

    # Add labels to the 2nd subplot
    ax_error_field.set_title(f'Errorfield of {show_depth} m results')
    ax_error_field.set_xlabel('Longitude')
    ax_error_field.set_ylabel('Latitude')

    # Plot 33333333333 the hypoxic min depth on the 3rd subplot
    # set common max, min limits for depth plots
    vmin = 60
    vmax = 150
    # Plot land borders
    ax_min_hypox.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
    # Plot data
    try:
        data = ds['Min_depth_hypoxia'].sel(time=ds['time'][time_index])
        data.plot(ax=ax_min_hypox, x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
    except KeyError:
        print('no hypoxic area in file')
    ax_min_hypox.set_title(f'minimum depths where oxygen <= {hypox} {unit}')
    ax_min_hypox.set_xlabel('Longitude')
    ax_min_hypox.set_ylabel('Latitude')

    # Plot the anoxic min depth  on the 4th subplot
    # Plot land borders
    ax_min_anox.contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
    # Plot data
    try:
        data = ds['Min_depth_anoxia'].sel(time=ds['time'][time_index])
        data.plot(ax=ax_min_anox, x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
    except KeyError:
        print('no anoxic area in file')
    ax_min_anox.set_title(f'minimum depths where oxygen <= {anox} {unit}')
    ax_min_anox.set_xlabel('Longitude')
    ax_min_anox.set_ylabel('Latitude')

    # Set a common title for each row of subplots
    ax_data.get_shared_x_axes().join(ax_data, ax_error_field)
    ax_data.get_shared_y_axes().join(ax_data, ax_error_field)

    # Adjust the spacing between subplots
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    # Add title and labels
    # Set the title for the whole figure
    fig.suptitle(f'maps of hypoxia and anoxia\n{year_month}')

    # Save the plot
    plt.savefig(f'{results_dir}/figures/maps_{year_month}_{show_depth}_{netcdf_filename}.png', dpi = 300, transparent=True)

bath_file = xr.open_dataset("C:/Work/DIVAnd/Oxygen_maps/bathymetry/gebco_30sec_4.nc")
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