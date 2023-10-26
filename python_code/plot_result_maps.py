import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd

### open netcdf file ### 
location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
netcdf_filename = "Oxygen_1960-2021_Autumn_1_49000.0_0.05_5.0_2.0_with_backgroundfield_moredepths_gebco_30sec_4"
# choose depth
show_depth = 70
# choose year
show_year = 2009

ds = xr.open_dataset(f'{location}/resultat/nc/processed/{netcdf_filename}.nc')

bath_file = xr.open_dataset("./bathymetry/gebco_30sec_4.nc")
# Extract the required variables
b = bath_file["bat"]

## extract values that are within our limits, save to a new variable and nc-file. ####
# in ml/l
hypox = 2
anox = 0.1
unit = 'ml/l'
# in umol/l
# 1 ml/l of O2 is approximately 43.570 µmol/kg 
# (assumes a molar volume of O2 of 22.392 l/mole and 
# a constant seawater potential density of 1025 kg/m3).
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
hypox = 90
anox = 9
unit = 'umol/l'

# Select a specific time index and depth level
# time_index = 0  # Replace with the desired time index
# Extract year and month from the time value

year_list =[datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
time_index = year_list.index(str(show_year))
time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
year_month = datetime.strftime(time_value, '%Y-%m')
ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')


depth_index = ds["depth"][:].values == show_depth
observation_span = 2

print(f'producing maps for year {show_year} at {show_depth} m {year_month}')
# Plot the data on a map
plt.style.use('dark_background')
# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
ax_data = axs[0, 0]
ax_error_field = axs[0, 1]
ax_min_hypox = axs[1, 0]
ax_min_anox = axs[1, 1]
# Create a 1x4 grid of subplots
fig, axs = plt.subplots(1, 4, figsize=(18, 4))
ax_data = axs[0]
ax_error_field = axs[1]
ax_min_hypox = axs[2]
ax_min_anox = axs[3]

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

ax_data.scatter(x=lon, y=lat, s = 2, c = observations, cmap = 'jet', facecolor = 'none', vmin=vmin_o2, vmax=vmax_o2)
# Add labels to the 1st subplot
ax_data.set_title(f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m')
ax_data.set_xlabel('Longitude')
ax_data.set_ylabel('Latitude')

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
plt.savefig(f'{location}/resultat/figures/maps_{year_month}_{show_depth}_{netcdf_filename}.png', dpi = 300)

