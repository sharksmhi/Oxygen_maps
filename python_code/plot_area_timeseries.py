import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from datetime import datetime
import pandas as pd

location = "//winfs-proj/proj/havgem/DIVA/syrekartor/" # or other location, like havgem path
df = pd.read_csv(f'{location}resultat/area_data.txt', sep = '\t')

mpl.rcParams['hatch.linewidth'] = 0.5 # sets the linewidth globally of all hatch lines
plt.rcParams['axes.axisbelow'] = True # makes the grid stay behind the bars

fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
fig.subplots_adjust(top=0.88,
                    bottom=0.09,
                    left=0.14,
                    right=0.95,
                    hspace=0.28)

subplot_row = 0
ax[subplot_row].set_ylim(0, 120000)
ax[subplot_row].yaxis.set_major_locator(ticker.MultipleLocator(20000))
ax[subplot_row].yaxis.set_minor_locator(ticker.MultipleLocator(10000))
ax[subplot_row].xaxis.set_major_locator(ticker.MultipleLocator(5))
ax[subplot_row].xaxis.set_minor_locator(ticker.MultipleLocator(1))
for season, subset in df.groupby('season'):
    # subset.plot.bar(ax=ax[subplot_row], x='year', y='Hypoxic_area_km2', width = 1, color = 'r')
    # subset.plot.bar(ax=ax[subplot_row], x='year', y='Hypoxic_relerr_area_km2', width = 1, color = 'grey')
    ax[subplot_row].bar(x=subset['year'], height=subset['Hypoxic_area_km2'], width = 1, color = 'r', edgecolor = 'black', linewidth = 0.5, label = 'hypoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Hypoxic_relerr_area_km2'], bottom = subset['Hypoxic_area_km2'], width = 1, color = 'r', hatch='////', edgecolor = 'black', linewidth = 0.5)
    ax[subplot_row].bar(x=subset['year'], height=subset['Anoxic_area_km2'], width = 1, color = 'grey', edgecolor = 'black', linewidth = 0.5, label = 'anoxic')
    ax[subplot_row].bar(x=subset['year'], height=subset['Anoxic_relerr_area_km2'], bottom = subset['Anoxic_area_km2'], width = 1, color = 'grey', hatch='////', edgecolor = 'black', linewidth = 0.5)
    ax[subplot_row].set_title(season)
    ax[subplot_row].set_ylabel('area km$^2$') 
    ax[subplot_row].grid(which = 'both')
    subplot_row += 1

ax[0].legend(loc=(0.8,1.1))
fig.savefig(f'{location}resultat/figures/timeseries_barchart.png', dpi = 300)



"""
import matplotlib.pyplot as plt
# set the spacing between subplots
plt.subplot_tool()
plt.show()
fig, ax = plt.subplots(4, 1,  figsize=(10, 8), sharey=True, sharex=True)
plt.subplot_tool()
plt.show()
"""