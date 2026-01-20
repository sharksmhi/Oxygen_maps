from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from datetime import datetime
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
from mpl_toolkits.axes_grid1 import inset_locator
import matplotlib as mpl
import matplotlib.patches as mpatches

mpl.rcParams['hatch.linewidth'] = 0.1

def add_created_stamp(ax, author="SMHI", fontsize=6):
    today = datetime.today().strftime("%Y-%m-%d")
    ax.text(
        0.02, 0.98,
        f"Created {today} by {author}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fontsize,
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"),
    )

def create_custom_colormap(levels, colors):
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(levels, cmap.N)

    return cmap, norm

def oxygen_cmap_norm(vmin, vmax, step=45):
    levels = np.arange(vmin, vmax + step, step)
    colors = ['black', 'brown', 'red', 'orange', 'yellow', 'green']
    cmap, _ = create_custom_colormap(levels, colors)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    return cmap, norm, levels

def create_colorbar(axis, plot_object, levels, title = None):
    # Create an inset_axes for the colorbar
    cbax = inset_locator.inset_axes(
        axis,
        width="40%",
        height="3%",
        loc="lower right",
        bbox_to_anchor=(0, 0.1, 1, 1),
        bbox_transform=axis.transAxes,
    )
    cbar = plt.colorbar(plot_object, cax=cbax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)
    cbar.set_ticks(levels)
    
    if title is not None:
        cbar.ax.set_title(title, fontsize=6, pad=2)

    return cbar

def set_up_cartopy_map(axis):
    """
    Create a map with cartopy in the given axis. 
    The map projections is PlateCarre to be able to plot gridlines and ticks.
    When creating the figures use ccrs.Mercator() to get a map with better ratio for our latitudes
    
    :param axis: matplotlib axis object
    """
    # Set map extent
    axis.set_extent(
        [
            9,
            29,
            53.5,
            61,
        ],
        crs=ccrs.PlateCarree()
    )

    # Coastlines
    axis.coastlines(resolution='10m', linewidth=0.5, color='gray')

    # Gridlines (meridians/parallels)
    gl = axis.gridlines(
        crs=ccrs.PlateCarree(), 
        draw_labels=True,
        linewidth=0.1,
        color='gray',
        alpha=0.5
    )
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}

    return axis

def plot_parameter(
    ds,
    parameter,
    axis,
    show_depth=None,
    vmin=None,
    vmax=None,
    levels=None,
    colors=None,
    time_index=0,
    cmap=None,
    norm=None,
    hatches=None,
    contourf=False
):
    """
    Plots a parameter from the xarray ds at the specified depth and time
    """
    axis = set_up_cartopy_map(axis)

    # select data to plot
    if "depth" in ds[parameter].dims and show_depth is not None:
        depth_index = ds["depth"].values == show_depth
        data = (
            ds[parameter]
            .isel(time=time_index)
            .sel(depth=ds["depth"][depth_index])
            .squeeze()
        )
    else:
        data = ds[parameter].isel(time=time_index).squeeze()

    lon = ds["lon"].values
    lat = ds["lat"].values

    lon2d, lat2d = np.meshgrid(lon, lat)

    # Build colormap & normalization if not provided
    if cmap is None or norm is None:
        if levels is not None and colors is not None:
            # Discrete colormap
            cmap, _ = create_custom_colormap(levels, colors)
            norm = plt.Normalize(vmin=min(levels), vmax=max(levels))
        elif vmin is not None and vmax is not None:
            # Continuous colormap
            cmap = plt.get_cmap('jet') if cmap is None else cmap
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
        else:
            # Default fallback
            cmap = plt.get_cmap('jet')
            norm = plt.Normalize(vmin=float(data.min()), vmax=float(data.max()))

    # add data to the axis and return axis and pcm to later change cmap
    if contourf:
        pcm = axis.contourf(
            lon,
            lat,
            data,
            levels=levels,
            colors=colors,
            hatches=hatches,
            transform=ccrs.PlateCarree(),
        )
    else:
        pcm = axis.pcolormesh(
            lon2d,
            lat2d,
            data,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree(),
        )
    

    return pcm

def overlay_observations(
    ds, axis, year,
    show_depth, observation_span,
    cmap, norm,
    time_index=0, BG=False
):
    """
    plots observations as scatter with colors from the observation values.
    """
    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()
    start_year = ds.attrs["start year"]
    end_year = ds.attrs["end year"]

    start_dt = pd.Timestamp(f"{int(start_year)}-01-01")
    end_dt = pd.Timestamp(f"{int(end_year)}-12-31")

    if not BG:
        year = time_value.year
        start_dt = pd.Timestamp(f"{year}-01-01")
        end_dt = pd.Timestamp(f"{year}-12-31")

    df = pd.DataFrame({
        'obsyear': ds['obsyear'],
        'obslon': ds['obslon'],
        'obslat': ds['obslat'],
        'value': ds['Oxygen_data'],
        'depth': ds['obsdepth'],
    })

    sel = (
        (df.obsyear >= start_dt) &
        (df.obsyear <= end_dt) &
        (df.depth >= show_depth - observation_span) &
        (df.depth <= show_depth + observation_span)
    )

    sc = axis.scatter(
        df.loc[sel, 'obslon'],
        df.loc[sel, 'obslat'],
        c=df.loc[sel, 'value'],
        s=6,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        linewidth=0.2,
        alpha=0.6,
        transform=ccrs.PlateCarree(),
    )

    return sc

def plot_parameter_with_observations(
    ds, parameter, axis, year,
    show_depth=0, observation_span=2,
    vmin=0, vmax=180,
    time_index=0, BG=False,
    add_colorbar=True
):
    """
    Combines a parameter at specified depth (using sub_plot_parameter)
    with observations using overlay_observations
    """
    cmap, norm, levels = oxygen_cmap_norm(vmin, vmax)

    pcm = plot_parameter(
        ds=ds,
        parameter=parameter,
        axis=axis,
        show_depth=show_depth,
        time_index=time_index,
        cmap=cmap,
        norm=norm,
    )

    sc = overlay_observations(
        ds, axis, year,
        show_depth, observation_span,
        cmap, norm,
        time_index, BG
    )

    if add_colorbar:
        create_colorbar(axis, pcm, levels, title="Oxygen umol/l")

    axis.set_title(
        f'Oxygen at {show_depth} m\nobservations ±{observation_span} m',
        fontsize=8
    )

    return pcm

def plot_only_observations(
    ds, axis, year,
    show_depth=0, observation_span=2,
    vmin=0, vmax=180,
    colorbar=True,
    color='k', time_index=0, BG=False
):
    """
    This plots observations from a selected depth 
    
    :param show_depth: integer, the depth that you want observations from, default 0 m (surface)
    :param observation_span: how large span around the show_depth you accept, default 2 m.
    :param BG: default False, set to True if the ds is a background fiels
    """

    time_value = ds['time'][time_index].values.astype('datetime64[M]').item()

    start_year = ds.attrs["start year"]
    end_year = ds.attrs["end year"]

    if not (int(start_year) <= int(year) <= int(end_year)):
        print(f"WARNING: {year} to plot not in range {start_year}-{end_year}")

    start_dt = pd.Timestamp(f"{int(start_year)}-01-01")
    end_dt = pd.Timestamp(f"{int(end_year)}-12-31")

    if not BG:
        year = time_value.year
        start_dt = pd.Timestamp(f"{year}-01-01")
        end_dt = pd.Timestamp(f"{year}-12-31")

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

    axis = set_up_cartopy_map(axis)

    if colorbar:
        levels = np.arange(vmin, vmax, 45)
        colors = ['black', 'brown', 'red', 'orange', 'yellow', 'green']
        cmap, norm = create_custom_colormap(levels, colors)

        sc = axis.scatter(
            lon,
            lat,
            s=5,
            c=observations,
            cmap=cmap,
            norm=norm,
            edgecolors='k',
            linewidth=0.2,
            alpha=0.5,
            transform=ccrs.PlateCarree()
        )

        create_colorbar(axis, sc, levels=np.arange(vmin, vmax+1, 45))

        axis.set_title(
            f'Oxygen at {show_depth} m\nobservation at +/- {observation_span} m',
            fontsize=8
        )

    else:
        sc = axis.scatter(
            lon,
            lat,
            s=5,
            edgecolors='k',
            linewidth=0.2,
            facecolor=color,
            transform=ccrs.PlateCarree()
        )
    return sc

def plot_errorfield(ds, parameter, axis, show_depth=None):

    # Modify the colormap levels to control the step length
    levels = [0,0.3,0.5,1]
    # Create a custom discrete colormap
    #        0       0.3       0.5
    colors = ['green', 'orange', 'red']
    cmap, norm = create_custom_colormap(levels, colors)

    pcm = plot_parameter(
        ds, parameter, axis, show_depth=show_depth, cmap=cmap, norm=norm
    )

    # Add a colorbar
    create_colorbar(axis=axis, plot_object=pcm, levels=levels)

    axis.set_title(f'Errorfield at {show_depth} m', fontsize=8)

def plot(results_dir, netcdf_filename, year, season, ds, threshold_list, interval, time_index, BG):
    # 1 ml/l of O2 is approximately 43.570 µmol/kg
    # (assumes a molar volume of O2 of 22.392 l/mole and
    # a constant seawater potential density of 1025 kg/m3).
    # https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html
    unit = 'umol/l'
    print(f"now plotting from: {netcdf_filename}")

    ############### FIRST plot, errors and min depths results #################
    # Row 1: maps with error fields for each threshol
    # Row 2: maps of the depth of the onset of each threshold
    # Create a 2x2 grid of subplots
    n_figs = len(threshold_list)
    
    fig, axs = plt.subplots(2, n_figs, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10, 4.5))
    # Adjust the spacing between subplots
    # fig.tight_layout()

    for index, threshold in enumerate(threshold_list):
        # Select a time index from year and depth level from show_depth

        # plot the relative error with a custom discrete colormap
        plot_errorfield(
            ds,
            parameter=f"Relerr_per_grid_at_min_{threshold}_depth",
            axis=axs[0, index],
        )
        axs[0, index].set_title(f'Error at <= {threshold} {unit} depth', fontsize=8)
        # second row plot the whole area each threshold with a colormap showing min_depth
        pcm = plot_parameter(ds, parameter=f"Min_depth_{threshold}",
            axis=axs[1, index], vmin=50,
            vmax=150,)
        axs[1, index].set_title(f'Area <= {threshold} {unit}', fontsize=8)
        create_colorbar(axs[1, index], pcm, levels=np.arange(50,150,20), title="depth [m]")

    # Add title and labels
    # Set the title for the whole figure
    if 0 in threshold_list:
        title_start = "Hypoxia and anoxia: "
    else:
        title_start = ""
    
    if season == "Winter":
        fig.suptitle(
            f"{title_start}{f'BG {int(year)-1}-{int(year)+1}' if BG else f'{int(year) - 1}-{int(year)}'} {season}",
            fontsize=10,
        )
    else:
        fig.suptitle(
            f"{title_start}{f'BG {int(year)-1}-{int(year)+1}' if BG else f'{int(year)}'} {season}",
            fontsize=10,
            )

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG_threshold_result{str(interval).replace(", ", "_")}_{season}.png', dpi=300,
                    transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/threshold_result{year}_{season}.png', dpi=300,
                    transparent=False)
        plt.close()

    # ################### SECOND to FOURTH plot result at different depthlayers ############################
    for name, depths in {"surf": [10, 20, 30, 50], "halo": [60, 70, 80, 90], "deep": [100, 110, 125, 150]}.items():

        # plots of results at 4 different depths 10, 40, 50, 60
        fig, axs = plt.subplots(2, 4, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10, 4.5), layout="constrained")

        if 0 in threshold_list:
            vmin_o2 = -45
            vmax_o2 = 180 + 45
        else:
            vmin_o2 = 90
            vmax_o2 = 360

        # on the 1st and 2nd plot we show oxygen set min and max for colorscale
        subplot_no = 0
        for show_depth in depths:
            plot_parameter_with_observations(
                    ds,
                    parameter="Oxygen",
                    axis=axs[0, subplot_no],
                    year=year,
                    show_depth=show_depth,
                    vmin=vmin_o2,
                    vmax=vmax_o2,
                    BG=BG,
                )
                
            # bottom row, plot error field at the selected depht
            plot_errorfield(
                ds,
                parameter="Oxygen_relerr",
                axis=axs[1, subplot_no],
                show_depth=show_depth
            )
            subplot_no += 1

        # Add title and labels
        # Set the title for the whole figure
        if season == "Winter":
            fig.suptitle(f'{int(year)-1}-{int(year)} {season}', fontsize=10,)
        else:
            fig.suptitle(f"{f'BG {int(year)-1}-{int(year)+1}' if BG else f'{int(year)}'} {season}", fontsize=10)

        # Save the plot
        if BG:
            plt.savefig(f'{results_dir}/figures/BG_{name}_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
            plt.close()
        else:
            plt.savefig(f'{results_dir}/figures/{name}_{year}_{season}.png', dpi=300, transparent=False)
            plt.close()

    ################# FIFTH plot final result map  ###########################

    # plots of results all observations and hypox area and with anox area overlayed
    fig, axs = plt.subplots(1, 1, subplot_kw={'projection': ccrs.Mercator()}, figsize=(10, 4.5), layout="constrained")

    # Vänder på threshold_list för att högst threshold skall hamna underst.
    threshold_list.reverse()
    color_list = ['lightgrey', 'darkgrey', 'grey']
    hatches_list = [10 * '/', 10 * "\\", 10*'|']

    for index, threshold in enumerate(threshold_list):
        pcm = plot_parameter(
            ds,
            parameter=f"{threshold}_mask_firstlayer",
            axis=axs,
            colors=[color_list[index]],
            hatches=[""], 
            levels=[0.5, 1.5], 
            contourf=True
        )
    
    for index, threshold in enumerate(threshold_list):
        # Mark areas with relative error >? with hatches
        pcm = plot_parameter(
            ds,
            parameter=f"Relerr_per_grid_at_min_{threshold}_depth",
            axis=axs,
            colors=["none"],
            hatches=[hatches_list[index]],
            levels=[0.5,1.5],
            contourf=True
        )

    # plot all points with observations in red
    sc = plot_only_observations(
        ds, axis=axs, year=year, colorbar=False, color="r", observation_span=500, BG=BG
    )

    # Skapa handle till legenden som motsvarar observationer i plotten
    sc.set_label("Observations")

    # Skapa colors and hatches objects for area handles
    area_colors = [color_list[0], color_list[1], color_list[2],'none', 'none','none']
    error_hatches = ['', '','', hatches_list[0], hatches_list[1],hatches_list[2]]

    # Skapa labels till patcherna 
    legend_labels_areas = [f'<{threshold_list[0]} µmol/l', f'<{threshold_list[1]} µmol/l', f'<{threshold_list[2]} µmol/l',
                       f'Error field <{threshold_list[0]} µmol/l', f'Error field <{threshold_list[1]} µmol/l',
                       f'Error field <{threshold_list[2]} µmol/l']
    
    # Skapa patch objekt till legenden
    patches = [mpatches.Patch(facecolor=color, hatch=hatch, label=label)
                   for color, hatch, label in zip(area_colors, error_hatches, legend_labels_areas)]

    # Lägg till en "fejk" legend om syrefritt är med
    if 0 not in threshold_list:
        # om bottniska viken.
        fake_labels = [f'<{threshold_list[0]} µmol/l', f'<{threshold_list[1]} µmol/l', f'<{threshold_list[2]} µmol/l',
                       f'Error field <{threshold_list[0]} µmol/l', f'Error field <{threshold_list[1]} µmol/l', f'Error field <{threshold_list[2]} µmol/l']
        fake_colors = [color_list[0], color_list[1], color_list[2],'none', 'none','none']
        fake_hatches = ['', '','', hatches_list[0], hatches_list[1],hatches_list[2]]
        # Skapa proxy-objekt för legenden
        patches = [mpatches.Patch(facecolor=color, hatch=hatch, label=label)
                   for color, hatch, label in zip(fake_colors, fake_hatches, fake_labels)]

        if len(threshold_list) < 3:
            print("Bara två thresholds, ändra legenden!")

    # Lägg till legenden
    axs.legend(handles=patches + [sc], loc='lower right', fontsize=8)
    add_created_stamp(axs)
    # Add title and labels
    # Set the title for the whole figure
    if season == "Winter":
        fig.suptitle(f'{int(year) - 1}-{year} {season}', fontsize=10)
    else:
        fig.suptitle(f'Hypoxia and anoxia: {f'BG {int(year)-1}-{int(year)+1}' if BG else f'{int(year)}'} {season}', fontsize=10)

    # Save the plot
    if BG:
        plt.savefig(f'{results_dir}/figures/BG'
                    f'_result_{str(interval).replace(", ", "_")}_{season}.png', dpi=300, transparent=False)
        plt.close()
    else:
        plt.savefig(f'{results_dir}/figures/final_result_{year}_{season}.png', dpi=300, transparent=False)
        plt.close()

## extract values that are within our limits, save to a new variable and nc-file. ####

def read_processed_nc(results_dir):
    
    processed_files = Path(results_dir) / "processed"
    for netcdf_filepath in list(processed_files.glob('*.nc')):#file_list:
        netcdf_filename = netcdf_filepath.name
        print(netcdf_filename)
        BG = False
        if 'Background' in netcdf_filename:
            BG = True
        
        ds = xr.open_dataset(f"{results_dir}/processed/{netcdf_filename}", engine='h5netcdf')
        ds['obsyear'] = ds['obstime'].values.astype('datetime64[Y]')
        
        ds_year_list = [datetime.strftime(timestr.astype('datetime64[M]').item(), '%Y') for timestr in ds["time"][:].values]
        if not len(ds_year_list) == 1:
            print('Skipping file due to multiple years in files')
            continue
        ds_year = ds_year_list[0]    
        time_index = 0
        season = ds.attrs['season']
        threshold_list = ds.attrs['threshold_list']
        # När attribute till nc sätts blir det av type Any, läses som en str från ds.attrs
        threshold_list = eval(threshold_list.replace("Any", ""))

        print(f"year in netcdf: {str(ds_year)}")
        interval = [int(ds_year) -1, int(ds_year) +1]   #Background year +/- 1
        plot(results_dir, netcdf_filename, ds_year, season, ds, threshold_list, interval, time_index, BG)

if __name__ == "__main__":
    print("running")
    # Result directory
    results_dir = "./results/Baltic_Proper/20260114_1653/"
    read_processed_nc(results_dir)