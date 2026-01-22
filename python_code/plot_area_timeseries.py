import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from datetime import datetime
import pandas as pd
from pathlib import Path

start_year = 1960
end_year = 2022

date = "20260114_1653"
freja_results_dir = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results/Baltic_Proper/{date}//")
on_freja = False
if Path(freja_results_dir).is_dir():
    results_dir = freja_results_dir
    on_freja = True

files = list(Path(results_dir).glob("area_data*.txt"))
df = pd.read_csv(files[0], sep="\t")
df["Area_90_km2_only"] = df["Area_90_km2"] - df["Area_0_km2"] 
df["Area_180_km2_only"] = df["Area_180_km2"] - df["Area_90_km2"] 

df_matlab = pd.read_csv('./data/area_volume_data_from_matlab.txt', sep = '\t')
df_matlab.rename(columns={"Anoxic_area_km2": "Area_0_km2", "Hypoxic_area_km2": "Area_90_km2"}, inplace=True)

mpl.rcParams['hatch.linewidth'] = 0.5 # sets the linewidth globally of all hatch lines
plt.rcParams['axes.axisbelow'] = True # makes the grid stay behind the bars
plt.rcParams["axes.formatter.limits"] = [-5, 4]

color_list = ['lightgrey', 'darkgrey', 'grey']

season_order = {'Winter': 0, 'Spring': 1, 'Summer': 2, 'Autumn': 3}
seasons = list(season_order.keys())
n_seasons = len(seasons)

# Timeseries one subplot for each season
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties

ax[0].set_ylim(0, 120000)

ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))

for season, subset in df.groupby('season'):
    subplot_row = season_order[season]
    # plot bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_180_km2"], width=1, color='lightgrey', edgecolor='black',linewidth=0.5, label='180 umol/l')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_180_km2'],bottom=subset['Area_180_km2'] - subset['Relerr_area_180_km2'], width=1, color='none', hatch='/////',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
    # plot hypoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_90_km2"], width = 1, color = 'darkgrey', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_90_km2'], bottom = subset['Area_90_km2']-subset['Relerr_area_90_km2'], width = 1, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
    # plot anoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset['Area_0_km2'], width = 1, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_0_km2'], bottom = subset['Area_0_km2']-subset['Relerr_area_0_km2'], width = 1, color = 'none', hatch='|||||', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")
    # set title, ylabel, grid on
    ax[subplot_row].set_title(season)
    ax[subplot_row].set_ylabel('area km$^2$') 
    ax[subplot_row].grid(which = 'major', linewidth=0.1, axis="y")
    ax[subplot_row].ticklabel_format(style='scientific', axis = 'y')

# turn on legend and set position
ax[0].legend(ncols=3, loc=(0.1,1.2))

# save figure
# TODO: fixa filnamnet!
temp_results_dir = "resultat"
fig.savefig(f'{temp_results_dir}/figures/{df.year.min()}-{df.year.max()}.png', dpi = 300)

# Timeseries group seasons together
fig, ax = plt.subplots(1, 1,  figsize=(15, 4), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties
ax.set_ylim(0, 120000)

ax.yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))

width = 0.2                 # bar width in year units
year_gap = width / 4        # gap between year groups
years = sorted(df['year'].unique())
# compute center offset for the group
group_offset = (n_seasons - 1) / 2 * width

# plot bars
# Loop through seasons
for i, season in enumerate(season_order.keys()):
    subset = df[df['season'] == season]

    # x positions for each bar: center the seasons around the year
    x_pos = subset['year'] + (i * width - group_offset)

    ax.bar(x=x_pos, height=subset[f"Area_180_km2"], width=width, color='lightgrey', edgecolor='black',linewidth=0.5, label='180 umol/l')
    ax.bar(x=x_pos, height=subset['Relerr_area_180_km2'],bottom=subset['Area_180_km2'] - subset['Relerr_area_180_km2'], width=width, color='none', hatch='/////',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
    # plot hypoxic bars
    ax.bar(x=x_pos, height=subset[f"Area_90_km2"], width = width, color = 'darkgrey', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
    ax.bar(x=x_pos, height=subset['Relerr_area_90_km2'], bottom = subset['Area_90_km2']-subset['Relerr_area_90_km2'], width = width, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
    # plot anoxic bars
    ax.bar(x=x_pos, height=subset['Area_0_km2'], width = width, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
    ax.bar(x=x_pos, height=subset['Relerr_area_0_km2'], bottom = subset['Area_0_km2']-subset['Relerr_area_0_km2'], width = width, color = 'none', hatch='|||||', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")

# set title, ylabel, grid on
ax.set_title("All seasons grouped")
ax.set_ylabel('area km$^2$') 
ax.ticklabel_format(style='scientific', axis = 'y')

# turn on legend and set position
ax.legend(ncols=3, loc=(0.1,1.2))

# save figure
# TODO: fixa filnamnet!
temp_results_dir = "resultat"
fig.savefig(f'{temp_results_dir}/figures/seasons_grouped_{df.year.min()}-{df.year.max()}.png', dpi = 300)


# Seasonal line plot
# List of variables to plot
variables = ['Area_0_km2', 'Area_90_km2', 'Area_180_km2']
variable_labels = ['Anoxic (0 µmol/l)', 'Hypoxic (90 µmol/l)', 'Normoxic (180 µmol/l)']

# Get all unique years
years = sorted(df['year'].unique())
n_years = len(years)
# Choose a colormap
cmap = cm.viridis  # or 'plasma', 'coolwarm', etc.
colors = cmap(np.linspace(0, 1, n_years))

# Normalize for colorbar
norm = Normalize(vmin=min(years), vmax=max(years))
# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(16,5), sharey=True)

for ax, var, label in zip(axes, variables, variable_labels):
    for i, year in enumerate(years):
        subset = df[df['year'] == year]
        # Make sure seasons are in correct order
        subset = subset.set_index('season').reindex(seasons)
        y = subset[var].values
        ax.plot(seasons, y, color=colors[i], linewidth=1)
    
    ax.set_title(label)
    ax.set_xlabel('Season')
    ax.grid(axis='y', linewidth=0.1)
    ax.ticklabel_format(style='scientific', axis='y')
    
axes[0].set_ylabel('Area (km²)')

# Add a single colorbar for all subplots
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label('Year')

temp_results_dir = "resultat"
fig.savefig(f'{temp_results_dir}/figures/seasonal_lineplot_{df.year.min()}-{df.year.max()}.png', dpi = 300)

# Create figure with 3 subplots
# List of variables to plot
variables = ['Area_0_km2', 'Area_90_km2_only', 'Area_180_km2_only']
variable_labels = ['Anoxic (0 µmol/l)', 'Hypoxic (90 µmol/l) only', 'Normoxic (180 µmol/l) only']
fig, axes = plt.subplots(1, 3, figsize=(16,5), sharey=True)

for ax, var, label in zip(axes, variables, variable_labels):
    for i, year in enumerate(years):
        subset = df[df['year'] == year]
        # Make sure seasons are in correct order
        subset = subset.set_index('season').reindex(seasons)
        y = subset[var].values
        ax.plot(seasons, y, color=colors[i], linewidth=1)
    
    ax.set_title(label)
    ax.set_xlabel('Season')
    ax.grid(axis='y', linewidth=0.1)
    ax.ticklabel_format(style='scientific', axis='y')
    
axes[0].set_ylabel('Area (km²)')

# Add a single colorbar for all subplots
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label('Year')

temp_results_dir = "resultat"
fig.savefig(f'{temp_results_dir}/figures/seasonal_lineplot_diff_{df.year.min()}-{df.year.max()}.png', dpi = 300)

# create figure object
fig, ax = plt.subplots(2, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties

ax[0].set_ylim(0, 100000)
ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))

df_merged = pd.merge(df.loc[df['season']=='Autumn'], df_matlab, on='year', how='outer', suffixes=['_DIVAnd', '_matlab'], sort=True)
# plot hypoxic bars
df_merged.plot.bar(ax=ax[0], x='year', y='Area_90_km2_matlab', width = 1, color = 'darkgrey', label='matlab')
df_merged.plot.bar(ax=ax[0], x='year', y='Area_90_km2_DIVAnd', width = 1, color = 'none', hatch='...', edgecolor = 'black', linewidth = 0.5, label='DIVAnd')

ax[0].set_title('Hypoxia')
ax[0].set_ylabel('area km$^2$') 
ax[0].grid(which = 'both')
# plot anoxic bars
df_merged.plot.bar(ax=ax[1], x='year', y='Area_0_km2_matlab', width = 1, color = 'grey', label='matlab')
df_merged.plot.bar(ax=ax[1], x='year', y='Area_0_km2_DIVAnd', width = 1, color = 'none', hatch='...', edgecolor = 'black', linewidth = 0.5, label='DIVAnd')
#df_merged.plot.bar(ax=ax[1], x='year', y='Anoxic_area_km2_DIVAnd', width = 1, color = 'r', label='DIVAnd')
#df_merged.plot.bar(ax=ax[1], x='year', y='Anoxic_area_km2_matlab', width = 1, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label='matlab')

# set title, ylabel, grid on
ax[1].set_title('Anoxia')
ax[1].set_ylabel('area km$^2$') 
ax[1].grid(which = 'both')

# turn on legend and set position
ax[0].legend(loc=(0.8,1.1))

# save figure
# TODO: fixa filnamnet
print(f'{results_dir}/figures/DIVAnd and matlab result.png')
temp_results_dir = "resultat"
fig.savefig(f'{temp_results_dir}/figures/DIVAnd and matlab result.png', dpi = 300)

"""
import matplotlib.pyplot as plt
# set the spacing between subplots
plt.subplot_tool()
plt.show()
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
plt.subplot_tool()
plt.show()
"""