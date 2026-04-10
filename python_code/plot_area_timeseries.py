import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from datetime import datetime
import pandas as pd
from pathlib import Path

def style_axis(ax, color="seashell"):
    ax.title.set_color(color)
    ax.xaxis.label.set_color(color)
    ax.yaxis.label.set_color(color)
    ax.tick_params(axis='both', colors=color)
    for spine in ax.spines.values():
        spine.set_edgecolor(color)

    leg = ax.get_legend()
    if leg:
        leg.get_frame().set_facecolor("none")
        for text in leg.get_texts():
            text.set_color(color)

start_year = 1960
end_year = 2025

results_dir = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results_lena_temp/Baltic_Proper/20260408_1932_high_res_1960_2025/")

#results_dir = Path(f"resultat/Baltic_Proper/{date}/")
print(results_dir)
files = list(Path(results_dir).glob("test_area_data*.txt"))
print(files)
df = pd.read_csv(files[0], sep="\t")
for basin in ["SEA-010_", "SEA-012_", ""]:
    df[f"Area_90_{basin}km2_only"] = df[f"Area_90_{basin}km2"] - df[f"Area_0_{basin}km2"] 
    df[f"Area_180_{basin}km2_only"] = df[f"Area_180_{basin}km2"] - df[f"Area_90_{basin}km2"] 

# List of variables to plot
basin = ""
variables = [f'Area_0_{basin}km2', f'Area_90_{basin}km2', f'Area_180_{basin}km2']
variable_labels = ['Anoxic (0 µmol/l)', 'Hypoxic (90 µmol/l)', 'Normoxic (180 µmol/l)']

df_matlab = pd.read_csv('./data/area_volume_data_from_matlab.txt', sep = '\t')
df_matlab.rename(columns={"Anoxic_area_km2": "Area_0_km2", "Hypoxic_area_km2": "Area_90_km2"}, inplace=True)

mpl.rcParams['hatch.linewidth'] = 0.5 # sets the linewidth globally of all hatch lines
plt.rcParams['axes.axisbelow'] = True # makes the grid stay behind the bars
plt.rcParams["axes.formatter.limits"] = [-5, 4]

color_list = ['lightgrey', 'darkgrey', 'grey']

season_order = {'Winter': 0, 'Spring': 1, 'Summer': 2, 'Autumn': 3}
seasons = list(season_order.keys())
n_seasons = len(seasons)

for season in seasons:
    fig, ax = plt.subplots(figsize=(10, 4))
    plt.subplots_adjust(top=0.8)

    subset = df.loc[df.season == season]
    # set axis properties
    if basin == "":
        ax.set_ylim(0, 120000)

    ax.yaxis.set_major_locator(ticker.MultipleLocator(20000))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10000))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # plot bars
    ax.bar(x=subset['year'], height=subset[f"Area_180_{basin}km2"], width=1, color='lightgrey', edgecolor='black',linewidth=0.5, label='180 umol/l')
    ax.bar(x=subset['year'], height=subset[f'Relerr_area_180_{basin}km2'],bottom=subset[f'Area_180_{basin}km2'] - subset[f'Relerr_area_180_{basin}km2'], width=1, color='none', hatch='/////',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
    # plot hypoxic bars
    ax.bar(x=subset['year'], height=subset[f"Area_90_{basin}km2"], width = 1, color = 'darkgrey', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
    ax.bar(x=subset['year'], height=subset[f'Relerr_area_90_{basin}km2'], bottom = subset[f'Area_90_{basin}km2']-subset[f'Relerr_area_90_{basin}km2'], width = 1, color = 'none', hatch=5* "\\", edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
    # plot anoxic bars
    ax.bar(x=subset['year'], height=subset[f'Area_0_{basin}km2'], width = 1, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
    ax.bar(x=subset['year'], height=subset[f'Relerr_area_0_{basin}km2'], bottom = subset[f'Area_0_{basin}km2']-subset[f'Relerr_area_0_{basin}km2'], width = 1, color = 'none', hatch='|||||', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")
    # set title, ylabel, grid on
    ax.set_title(f"{season}")
    ax.set_ylabel('area km$^2$') 
    ax.grid(which = 'major', linewidth=0.1, axis="y")
    ax.ticklabel_format(style='scientific', axis = 'y')

    # turn on legend and set position
    ax.legend(ncols=3, loc=(0.1,1.1))
    fig.suptitle(basin.strip("_"))
    style_axis(ax)
    # save figure
    # TODO: fixa filnamnet!
    fig.savefig(f'{results_dir}/figures/timeseries/{season}_{basin}{df.year.min()}-{df.year.max()}.png', dpi = 300, transparent=True)

# Timeseries one subplot for each season
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties

if basin == "":
    ax[0].set_ylim(0, 120000)

ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))

for season, subset in df.groupby('season'):
    subplot_row = season_order[season]
    # plot bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_180_{basin}km2"], width=1, color='lightgrey', edgecolor='black',linewidth=0.5, label='180 umol/l')
    ax[subplot_row].bar(x=subset['year'], height=subset[f'Relerr_area_180_{basin}km2'],bottom=subset[f'Area_180_{basin}km2'] - subset[f'Relerr_area_180_{basin}km2'], width=1, color='none', hatch='/////',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
    # plot hypoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_90_{basin}km2"], width = 1, color = 'darkgrey', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset[f'Relerr_area_90_{basin}km2'], bottom = subset[f'Area_90_{basin}km2']-subset[f'Relerr_area_90_{basin}km2'], width = 1, color = 'none', hatch=5* "\\", edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
    # plot anoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f'Area_0_{basin}km2'], width = 1, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset[f'Relerr_area_0_{basin}km2'], bottom = subset[f'Area_0_{basin}km2']-subset[f'Relerr_area_0_{basin}km2'], width = 1, color = 'none', hatch='|||||', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")
    # set title, ylabel, grid on
    ax[subplot_row].set_title(season)
    ax[subplot_row].set_ylabel('area km$^2$') 
    ax[subplot_row].grid(which = 'major', linewidth=0.1, axis="y")
    ax[subplot_row].ticklabel_format(style='scientific', axis = 'y')

# turn on legend and set position
ax[0].legend(ncols=3, loc=(0.1,1.2))
fig.suptitle(basin.strip("_"))
# save figure
# TODO: fixa filnamnet!
fig.savefig(f'{results_dir}/figures/timeseries/{basin}{df.year.min()}-{df.year.max()}.png', dpi = 300)

# Timeseries group seasons together
basins = ["", "SEA-007_", "SEA-009_", "SEA-010_", ]
fig, axs = plt.subplots(len(basins), 1,  figsize=(10, 4), layout="constrained", sharex=True)
#fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties
"""if basin == "":
   ax.set_ylim(0, 120000)
ax.yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))"""

width = 0.2                 # bar width in year units
year_gap = width / 4        # gap between year groups
years = sorted(df['year'].unique())
# compute center offset for the group
group_offset = (n_seasons - 1) / 2 * width

# plot bars
for basin_i, basin in enumerate(basins):
    # Loop through seasons
    ax = axs[basin_i]
    for i, season in enumerate(season_order.keys()):
        subset = df[df['season'] == season]

        # x positions for each bar: center the seasons around the year
        x_pos = subset['year'] + (i * width - group_offset)

        ax.bar(x=x_pos, height=subset[f"Area_180_{basin}km2"], width=width, color='lightgrey', edgecolor='black',linewidth=0.5, label='180 umol/l')
        ax.bar(x=x_pos, height=subset[f'Relerr_area_180_{basin}km2'],bottom=subset[f'Area_180_{basin}km2'] - subset[f'Relerr_area_180_{basin}km2'], width=width, color='none', hatch=10 * '/',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
        # plot hypoxic bars
        ax.bar(x=x_pos, height=subset[f"Area_90_{basin}km2"], width = width, color = 'darkgrey', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
        ax.bar(x=x_pos, height=subset[f'Relerr_area_90_{basin}km2'], bottom = subset[f'Area_90_{basin}km2']-subset[f'Relerr_area_90_{basin}km2'], width = width, color = 'none', hatch=10 * "\\", edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
        # plot anoxic bars
        ax.bar(x=x_pos, height=subset[f'Area_0_{basin}km2'], width = width, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
        ax.bar(x=x_pos, height=subset[f'Relerr_area_0_{basin}km2'], bottom = subset[f'Area_0_{basin}km2']-subset[f'Relerr_area_0_{basin}km2'], width = width, color = 'none', hatch=10*'|', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")
        ax.ticklabel_format(style='scientific', axis = 'y')
        ax.set_title(basin.strip("_"))
        ax.set_ylabel('area km$^2$')

# turn on legend and set position
#ax.legend(ncols=3, loc=(0.1,1.2))

fig.suptitle("all seasons grouped")
# save figure
# TODO: fixa filnamnet!
fig.savefig(f'{results_dir}/figures/timeseries/many_{basin}_seasons_grouped_{df.year.min()}-{df.year.max()}.png', dpi = 300)

# Seasonal line plot
# List of variables to plot
basin = "SEA-010_"
variables = [f'Area_0_{basin}km2', f'Area_90_{basin}km2', f'Area_180_{basin}km2']
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
if basin == "":
    axes[0].set_ylim(0, 120000)

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

fig.suptitle(basin)

# Add a single colorbar for all subplots
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label('Year')

fig.savefig(f'{results_dir}/figures/timeseries/{basin}_seasonal_lineplot_{df.year.min()}-{df.year.max()}.png', dpi = 300)

# Create figure with 3 subplots
# List of variables to plot
basin = ""
variables = [f'Area_0_{basin}km2', f'Area_90_{basin}km2_only', f'Area_180_{basin}km2_only']
variable_labels = ['Anoxic (0 µmol/l)', 'Hypoxic (90 µmol/l) only', 'Normoxic (180 µmol/l) only']
fig, axes = plt.subplots(1, 3, figsize=(16,5), sharey=True)
if basin == "":
    axes[0].set_ylim(0, 120000)

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
fig.suptitle(f"diff {basin}")
# Add a single colorbar for all subplots
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label('Year')

fig.savefig(f'{results_dir}/figures/timeseries/ALL_{basin}_seasonal_lineplot_diff_{df.year.min()}-{df.year.max()}.png', dpi = 300)
df_merged = pd.merge(df.loc[df['season']=='Autumn'], df_matlab, on='year', how='outer', suffixes=['_DIVAnd', '_matlab'], sort=True)

# create figure object
fig, ax = plt.subplots(2, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties
if basin == "":
    ax[0].set_ylim(0, 120000)

ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))
# plot hypoxic bars
df_merged.plot.line(ax=ax[0], x='year', y='Area_90_km2_matlab', color = 'deepskyblue', label='matlab')
df_merged.plot.line(ax=ax[0], x='year', y='Area_90_km2_DIVAnd', color = 'forestgreen', label='DIVAnd')

ax[0].set_title('Hypoxia')
ax[0].set_ylabel('area km$^2$') 
# plot anoxic bars
df_merged.plot.line(ax=ax[1], x='year', y='Area_0_km2_matlab', color = 'deepskyblue', label='matlab')
df_merged.plot.line(ax=ax[1], x='year', y='Area_0_km2_DIVAnd', color = 'forestgreen', label='DIVAnd')

# set title, ylabel, grid on
ax[1].set_title('Anoxia')
ax[1].set_ylabel('area km$^2$') 

# turn on legend and set position
ax[0].legend(loc=(0.8,1.1))

# save figure
# TODO: fixa filnamnet
print(f'{results_dir}/figures/timeseries/DIVAnd and matlab result.png')
fig.savefig(f'{results_dir}/figures/timeseries/DIVAnd and matlab result.png', dpi = 300)


# create figure object
# make all text white (seashell)
print(f"{plt.rcParams.keys()=}")
plt.rcParams.update({
    "text.color": "seashell",       # general text color
    "axes.labelcolor": "seashell", # x/y axis labels
    "axes.edgecolor": "seashell", 
    "xtick.color": "seashell",        # x ticks
    "ytick.color": "seashell",        # y ticks
    "axes.titlecolor": "seashell",   # axes title
    "axes.titlesize": 16,
    "legend.facecolor": "none"
})
fig, ax = plt.subplots(2, 1,  figsize=(10, 8), sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties
if basin == "":
    ax[0].set_ylim(0, 120000)

ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))
# plot hypoxic bars
df_merged.plot.line(ax=ax[0], x='year', y='Area_90_km2_matlab', ls="--", color = 'deepskyblue', label='Hypoxia Hansson et al 2011')
df_merged.plot.line(ax=ax[0], x='year', y='Area_90_km2_DIVAnd', color = 'deepskyblue', label='Hypoxia This study using DIVAnd')
df_merged.plot.line(ax=ax[0], x='year', y='Area_0_km2_matlab', ls="--", color = 'orangered', label='Anoxia Hansson et al 2011')
df_merged.plot.line(ax=ax[0], x='year', y='Area_0_km2_DIVAnd', color = 'orangered', label='Anoxia This study using DIVAnd')

# set title, ylabel, grid on
ax[0].set_title('Results')
ax[0].set_ylabel('area km$^2$') 

# turn on legend and set position
ax[0].legend(loc='upper left', ncol=2)

ax[1].set_ylim(-20000, 20000)
ax[1].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[1].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[1].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[1].xaxis.set_minor_locator(ticker.MultipleLocator(1))
# plot hypoxic bars
df_merged["diff_90"] = df_merged['Area_90_km2_matlab']-df_merged["Area_90_km2_DIVAnd"]
df_merged["diff_0"] = df_merged['Area_0_km2_matlab']-df_merged["Area_0_km2_DIVAnd"]
df_merged.plot.line(ax=ax[1], x='year', y='diff_90', color = 'deepskyblue', label='hypoxia' )
df_merged.plot.line(ax=ax[1], x='year', y='diff_0', color = 'orangered', label='anoxia')
ax[1].axhline(y=0, color = "seashell", lw=0.1)

# set title, ylabel, grid on
ax[1].set_title('Difference, Hansson et al 2011 - Results using DIVAnd')
ax[1].set_ylabel('delta area km$^2$') 

# turn on legend and set position
ax[1].legend(facecolor="none")

# save figure
# TODO: fixa filnamnet
print(f'{results_dir}/figures/timeseries/DIVAnd and matlab result.png')
fig.savefig(f'{results_dir}/figures/timeseries/DIVAnd and matlab result_together.png', dpi = 300, transparent = True)