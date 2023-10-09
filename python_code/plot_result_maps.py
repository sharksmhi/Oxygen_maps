import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd

### open netcdf file ### 
location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
netcdf_filename = "Oxygen_2009-2009_Autumn_2.5_7.0_gebco_30sec_4"
ds = xr.open_dataset(f'{location}/resultat/nc/processed/{netcdf_filename}.nc')

bath_file = xr.open_dataset("./bathymetry/gebco_30sec_4.nc")
print(bath_file)
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
time_index = 0  # Replace with the desired time index
# Extract year and month from the time value
print(repr(ds["time"][:].values))
print(repr(ds["time"][time_index].values))
time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
year_month = datetime.strftime(time_value, '%Y-%m')
ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')

print(repr(ds["depth"][:].values))
# choose depth
show_depth = 70
depth_index = ds["depth"][:].values == show_depth

# Plot the data on a map
# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
# 1111111 Plot the data on the 1st subplot
# Plot land borders
axs[0, 0].contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
# Plot data
data = ds['Oxygen'].sel(time=ds['time'][time_index], depth = ds['depth'][depth_index])
data.plot(ax=axs[0, 0], x='lon', y='lat', cmap='jet', vmin=0, vmax=180)
# data = ds['Oxygen_data'].sel(obsdepth=ds['Oxygen_data'][40])
# data.plot(ax=axs[0, 0], x='lon', y='lat', cmap='jet')
df = pd.DataFrame({'obsyear': ds['obsyear'], 'obslon': ds['obslon'], 'obslat': ds['obslat'], 'Oxygen_data': ds['Oxygen_data'], 'depth': ds['obsdepth']})
selection = ((df.obsyear == datetime.strftime(time_value, '%Y')) & (df.depth >= show_depth-5) & (df.depth <= show_depth+5))
lon = df.loc[selection, 'obslon']
lat = df.loc[selection, 'obslat']
observations = df.loc[selection, 'Oxygen_data']
# obs_val = df.loc[df.obsyear == datetime.strftime(time_value, '%Y'), 'Oxygen_data']
# lat = xr.where(ds['obs_year']==datetime.strftime(time_value, '%Y'), ds['obslat'], np.nan)
print(df.head())

axs[0,0].scatter(x=lon, y=lat, s = 2, c = observations, facecolor = 'none', vmin=0, vmax=180)
# Add labels to the 1st subplot
axs[0, 0].set_title(f'Oxygen at {show_depth} m\nobservation at +/- 5 m')
axs[0, 0].set_xlabel('Longitude')
axs[0, 0].set_ylabel('Latitude')

# 2222222 Plot the data on the 2nd subplot

# Plot the relerrfield of the hypoxuc depth layer
data = ds['Hypoxic_relerr_per_grid'].sel(time=ds['time'][time_index])
print(data.head())
data.plot(ax=axs[0, 1], x='lon', y='lat', cmap='jet', vmin=0, vmax=0.5)
# Plot land borders
axs[0, 1].contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")

# Add labels to the 2nd subplot
axs[0, 1].set_title(f'Errorfield of Hypoxic results')
axs[0, 1].set_xlabel('Longitude')
axs[0, 1].set_ylabel('Latitude')
# Set a common title for each row of subplots
fig.text(0.5, 0.94, f'{year_month} (Time index: {time_index})', ha='center')

# Select the last time index
last_time_index = 0
time_value = ds['time'][last_time_index].values.astype('datetime64[M]').item()
year_month = datetime.strftime(time_value, '%Y-%m')

# Plot the hypoxic min depth on the 3rd subplot
# set common max, min limits for depth plots
vmin = 60
vmax = 150
# Plot land borders
axs[1, 0].contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
# Plot data
data = ds['Min_depth_hypoxia'].sel(time=ds['time'][last_time_index])
data.plot(ax=axs[1, 0], x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
axs[1, 0].set_title(f'minimum depths where oxygen <= {hypox} {unit}')
axs[1, 0].set_xlabel('Longitude')
axs[1, 0].set_ylabel('Latitude')

# Plot the anoxic min depth  on the 4th subplot
# Plot land borders
axs[1, 1].contourf(bath_file["lon"], bath_file["lat"], -b, levels=[-1e5,0], colors="gray")
# Plot data
data = ds['Min_depth_anoxia'].sel(time=ds['time'][last_time_index])
data.plot(ax=axs[1, 1], x='lon', y='lat', cmap='jet', vmin=vmin, vmax=vmax)
axs[1, 1].set_title(f'minimum depths where oxygen <= {anox} {unit}')
axs[1, 1].set_xlabel('Longitude')
axs[1, 1].set_ylabel('Latitude')

# Set a common title for each row of subplots
axs[0, 0].get_shared_x_axes().join(axs[0, 0], axs[0, 1])
axs[0, 0].get_shared_y_axes().join(axs[0, 0], axs[0, 1])

# Set a common title for each row of subplots
fig.text(0.5, 0.47, f'{year_month} (Time index: {last_time_index})', ha='center')

# Adjust the spacing between subplots
fig.tight_layout(rect=[0, 0, 1, 0.95])

# Add title and labels
# Set the title for the whole figure
fig.suptitle(f'maps of hypoxia and anoxia')

# Save the plot
plt.savefig(f'{location}/resultat/figures/maps_{year_month}.png', dpi = 300)

