import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from datetime import datetime
import pandas as pd

start_year = 1960
end_year = 2022

location = "//winfs-proj/proj/havgem/DIVA/syrekartor/" # or other location, like havgem path
df = pd.read_csv(f'{location}resultat/area_data_{start_year}_{end_year}.txt', sep = '\t')
df_matlab = pd.read_csv(f'{location}resultat/area_volume_data_from_matlab.txt', sep = '\t')
df_matlab.rename(columns={"Anoxic_area_km2": "Area_0_km2", "Hypoxic_area_km2": "Area_90_km2"}, inplace=True)

mpl.rcParams['hatch.linewidth'] = 0.5 # sets the linewidth globally of all hatch lines
plt.rcParams['axes.axisbelow'] = True # makes the grid stay behind the bars
plt.rcParams["axes.formatter.limits"] = [-5, 4]

subplot_order = {'Autumn': 0, 'Winter': 1, 'Spring': 2, 'Summer': 3}

# create figure object
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88, bottom=0.09, left=0.14, right=0.95, hspace=0.28)

# set axis properties

ax[0].set_ylim(0, 100000)

ax[0].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(1))

for season, subset in df.groupby('season'):
    subplot_row = subplot_order[season]
    # plot bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_180_km2"], width=1, color='b', edgecolor='black',linewidth=0.5, label='180 umol/l')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_180_km2'],bottom=subset['Area_180_km2'] - subset['Relerr_area_180_km2'], width=1, color='b', hatch='/////',edgecolor='black', linewidth=0.5, label="180 umol/l Relerror >0.5")
    # plot hypoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset[f"Area_90_km2"], width = 1, color = 'r', edgecolor = 'black', linewidth = 0.5, label = 'Hypoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_90_km2'], bottom = subset['Area_90_km2']-subset['Relerr_area_90_km2'], width = 1, color = 'r', hatch='/////', edgecolor = 'black', linewidth = 0.5, label="Hypoxic Relerror >0.5")
    # plot anoxic bars
    ax[subplot_row].bar(x=subset['year'], height=subset['Area_0_km2'], width = 1, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'Anoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Relerr_area_0_km2'], bottom = subset['Area_0_km2']-subset['Relerr_area_0_km2'], width = 1, color = 'grey', hatch='|||||', edgecolor = 'black', linewidth = 0.5, label="Anoxic Relerror >0.5")
    # set title, ylabel, grid on
    ax[subplot_row].set_title(season)
    ax[subplot_row].set_ylabel('area km$^2$') 
    ax[subplot_row].grid(which = 'both')
    ax[subplot_row].ticklabel_format(style='scientific', axis = 'y')

# turn on legend and set position
ax[0].legend(ncols=3, loc=(0.1,1.1))

# save figure
fig.savefig(f'{location}resultat/figures/{df.year.min()}-{df.year.max()}_barchart.png', dpi = 300)

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
df_merged.plot.bar(ax=ax[0], x='year', y='Area_90_km2_matlab', width = 1, color = 'red', label='matlab')
df_merged.plot.bar(ax=ax[0], x='year', y='Area_90_km2_DIVAnd', width = 1, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label='DIVAnd')

ax[0].set_title('Hypoxia')
ax[0].set_ylabel('area km$^2$') 
ax[0].grid(which = 'both')
# plot anoxic bars
df_merged.plot.bar(ax=ax[1], x='year', y='Area_0_km2_matlab', width = 1, color = 'red', label='matlab')
df_merged.plot.bar(ax=ax[1], x='year', y='Area_0_km2_DIVAnd', width = 1, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label='DIVAnd')
#df_merged.plot.bar(ax=ax[1], x='year', y='Anoxic_area_km2_DIVAnd', width = 1, color = 'r', label='DIVAnd')
#df_merged.plot.bar(ax=ax[1], x='year', y='Anoxic_area_km2_matlab', width = 1, color = 'none', hatch='/////', edgecolor = 'black', linewidth = 0.5, label='matlab')

# set title, ylabel, grid on
ax[1].set_title('Anoxia')
ax[1].set_ylabel('area km$^2$') 
ax[1].grid(which = 'both')

# turn on legend and set position
ax[0].legend(loc=(0.8,1.1))

# save figure
print(f'{location}resultat/figures/DIVAnd and matlab result.png')
fig.savefig(f'{location}resultat/figures/DIVAnd and matlab result.png', dpi = 300)

"""
import matplotlib.pyplot as plt
# set the spacing between subplots
plt.subplot_tool()
plt.show()
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
plt.subplot_tool()
plt.show()
"""