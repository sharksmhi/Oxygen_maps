import xarray as xr
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import matplotlib.colors as mcolors
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as cm 
import seaborn as sns

year = 2020
season = "Autumn"
results_dir = Path(f"/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/results_lena_temp/Baltic_Proper/20260320_1659/")

# Definiera djupintervall (label, lower, upper)
depth_ranges = [
    ("0–50 m", 0, 50),
    ("50–70 m", 50, 70),
    ("70–90 m", 70, 90),
    ("90–150 m", 90, 150),
    (">150 m", 150, float("inf")),
]

############### netcdf data ###################
pattern = f"Oxygen_{year}_{season}_0*.nc"
files = list((results_dir / "DIVArun").glob(pattern))
if not files:
    raise FileNotFoundError(f"No file found for {pattern}")
if len(files) > 1:
    print("Warning: multiple files found, using first one:", files[0])
# take the first match
observation_nc = files[0]

residuals_nc = results_dir / f"DIVArun/Oxygen_{year}_{season}_residuals.nc"
ds_obs = xr.open_dataset(observation_nc, engine='h5netcdf')
ds_res = xr.open_dataset(residuals_nc, engine='h5netcdf')

obs = ds_obs["Oxygen_data"].values
res = ds_res["Oxygen_residual"].values
lon = ds_obs["obslon"].values
lat = ds_obs["obslat"].values
depth = ds_obs["obsdepth"].values
time = ds_obs["obstime"].values

diva_result = obs - res

df = pd.DataFrame({
    "obs": ds_obs["Oxygen_data"].values,
    "res": ds_res["Oxygen_residual"].values,
    "lon": ds_obs["obslon"].values,
    "lat": ds_obs["obslat"].values,
    "depth": ds_obs["obsdepth"].values,
    "time": ds_obs["obstime"].values
})
df_clean = df.dropna(subset=["obs", "res"])
df_clean["diva_result"] = df_clean["obs"] - df_clean["res"]

############## histogram ##################
weights = np.ones_like(res) / len(res) * 100  # percent of total
filtered = res[(res < -50) | (res > 50)]
weights_filtered = np.ones_like(filtered) / len(res) * 100  # percent of total
percent_filtered = len(filtered) / len(res) *100

plt.hist(res, bins=20, weights=weights, edgecolor="grey")
plt.ylabel("Percent of total (%)")
plt.xlabel("Residuals, obs-interpolation (µmol/l)")
plt.suptitle(f"Residuals {year} {season}\n{np.round(percent_filtered,1)}% of |Residuals| > 50 µmol/l")
plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_histogram_percent_{year}_{season}.png')
plt.close()

plt.hist(filtered, bins=20, weights=weights_filtered, edgecolor="grey")
plt.suptitle(f"{year} {season} |Residuals| > 50\n{np.round(percent_filtered,1)}% of total")
plt.ylabel("Percent of total (%)")
plt.xlabel("Residuals, obs-interpolation (µmol/l)")
plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_filtered_histogram_percent_{year}_{season}.png')
plt.close()

data = [df_clean[(df_clean["depth"] >= dmin) & (df_clean["depth"] < dmax)]["res"] 
        for _, dmin, dmax in depth_ranges]
plt.hist(data, bins=20, stacked=True, label=[label for label, _, _ in depth_ranges], edgecolor='grey')
plt.ylabel("Count")
plt.xlabel("Residuals, obs-interpolation (µmol/l)")
plt.suptitle(f"{year}_{season}")
plt.legend()
plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_histogram_count_{year}_{season}.png')
plt.close()

for label, dmin, dmax in depth_ranges:
    subset = df_clean[(df_clean["depth"] >= dmin) & (df_clean["depth"] < dmax)]
    sns.kdeplot(subset["res"], label=label)

plt.xlabel("Residuals")
plt.ylabel("Density")
plt.legend()
plt.close()

bins = np.linspace(df_clean["res"].min(), df_clean["res"].max(), 21)
for label, dmin, dmax in depth_ranges:
    subset = df_clean[(df_clean["depth"] >= dmin) & (df_clean["depth"] < dmax)]
    weights = np.ones_like(subset["res"]) / len(subset) * 100
    
    plt.hist(subset["res"], bins=bins, weights=weights,
             histtype='step', linewidth=2, label=label)

plt.xlabel("Residuals, obs-interpolation (µmol/l)")
plt.ylabel("Percent (%)")
plt.legend()
plt.close()
threshold = 50

labels = []
percent_extreme = []

for label, dmin, dmax in depth_ranges:
    subset = df_clean[(df_clean["depth"] >= dmin) & (df_clean["depth"] < dmax)]
    extreme = subset[(subset["res"] < -threshold) | (subset["res"] > threshold)]
    
    pct = len(extreme) / len(subset) * 100
    labels.append(label)
    percent_extreme.append(pct)

plt.bar(labels, percent_extreme)
plt.ylabel("Extreme residuals (%)")
plt.title(f"Fraction of |residuals| > {threshold}\n{year}_{season}")
plt.savefig(f'{results_dir}/figures/residuals/fraction of |res| > {threshold}_{year}_{season}_.png')
plt.close()

df_clean["depth_group"] = pd.cut(
    df_clean["depth"],
    bins=[0, 50, 60, 100, 150, np.inf],
    labels=[d[0] for d in depth_ranges]
)

sns.boxplot(data=df_clean, x="depth_group", y="res")
plt.ylabel("Residuals")
plt.xlabel("Depth interval")
plt.title(f"{year}_{season}")
plt.savefig(f'{results_dir}/figures/residuals/boxplot_{year}_{season}_.png')
plt.close()

################### observation vs interpolation ##################
# Beräkna gemensamma axelgränser
xy_min = min(df_clean["obs"].min(), df_clean["diva_result"].min())
xy_max = max(df_clean["obs"].max(), df_clean["diva_result"].max())

# Skapa subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

# Loop över varje djupintervall och plotta
for i, (label, dmin, dmax) in enumerate(depth_ranges[1:]):
    # Filtrera data
    #mask = (depth > dmin) & (depth <= dmax)
    mask = (df_clean["depth"] >= dmin) & (df_clean["depth"] <= dmax) & (df_clean["res"].abs() > 50)
    obs_sel = df_clean["obs"][mask]
    interp_sel = df_clean["diva_result"][mask]
    depth_sel = df_clean["depth"][mask]

    # Gör scatterplot
    sc = axs[i].scatter(obs_sel, interp_sel, c=depth_sel, cmap='viridis')
    axs[i].set_title(f"{label}")
    axs[i].set_xlabel("Observation")
    axs[i].set_ylabel("Interpolation (DIVAnd result)")
    axs[i].set_xlim(xy_min, xy_max)
    axs[i].set_ylim(xy_min, xy_max)
    axs[i].set_aspect('equal')

    # Beräkna mått om data finns
    if len(obs_sel) > 0:
        rmse = np.sqrt(np.mean((obs_sel - interp_sel) ** 2))
        bias = np.mean(interp_sel - obs_sel)
        r2 = r2_score(obs_sel, interp_sel)

        textstr = f"RMSE: {rmse:.2f}\nBias: {bias:.2f}\nR²: {r2:.2f}"
        axs[i].text(0.05, 0.95, textstr, transform=axs[i].transAxes,
                    fontsize=10, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.8))

    # colorbar
    plt.colorbar(sc, ax=axs[i], label="Djup")

# Layout
plt.suptitle(f"{year}_{season}\nonly |residuals| > 50", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f'{results_dir}/figures/residuals/interp_vs_obs_{year}_{season}_only_worst_residuals.png')
plt.close()

#=================
#MAPS STARTS
for label, dmin, dmax in depth_ranges:
    subset = df_clean[(df_clean["depth"] >= dmin) & (df_clean["depth"] < dmax)]
    subset = subset.sort_values(by="res", key=np.abs, ascending=False)
    extreme = subset[(subset["res"] < -threshold) | (subset["res"] > threshold)]

    levels = np.array([-200, -150, -100, -50, 0, 50, 100, 150, 200])
    cmap = cm.get_cmap("jet", len(levels) - 1)
    norm = mcolors.BoundaryNorm(levels, cmap.N)

    abs_res = np.abs(extreme["res"])
    # Normalize to 0–1
    norm_sizes = (abs_res - abs_res.min()) / (abs_res.max() - abs_res.min())
    # Scale to a reasonable size range
    sizes = 10 + norm_sizes * 90   # sizes between 10 and 100

    fig, ax = plt.subplots(figsize=(15, 8), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)
    #fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Mercator()}, figsize=(10, 4.5))
    ax.set_extent([9, 32, 53, 62], ccrs.PlateCarree())
    # land in gray:
    ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor="none", zorder=0, alpha=0.5)
    # zorder ==> gives the llayer position (0 = background, 1 top layer, goes until 10). 
    # Coast
    ax.coastlines(resolution="10m", linewidth=1)
    # Grid
    gl = ax.gridlines(draw_labels=True, zorder=2) #to keep ticks inside map
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 18}#, "weight": "bold"}
    gl.ylabel_style = {"size": 18}#, "weight": "bold"}

    # Plot data
    cs = plt.scatter(
        x=extreme["lon"],
        y=extreme["lat"],
        c=extreme["res"],
        s=sizes,
        transform=ccrs.PlateCarree(),
        marker='o',
        cmap=cmap,
        norm=norm,
        rasterized=True,
        edgecolors='none',
        alpha=1,
        zorder=3
    )
    ax.text(
        0.02, 0.98,
        "Blue dots interpolation > observation\n"
        "Red dots observation > interpolation\n"
        "Larger size indicates larger |residual|\n"
        "map only shows |residual|>50",
        transform=ax.transAxes,   # <-- key difference
        fontsize=14,
        verticalalignment='top',
        bbox=dict(
            facecolor="white",
            edgecolor="black",
            boxstyle="round,pad=0.3",
            alpha=0.8
        ),
        zorder=5
    )
    cbar = fig.colorbar(
        cs,
        ax=ax,
        orientation="horizontal",
        pad=0.05,
        shrink=0.5,
        boundaries=levels,
        ticks=levels
    )
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label("Difference between observations and interpolation at obs. point (umol/L)", fontsize=18)
    plt.suptitle(f"Residuals {label} {year} {season}", fontsize=18)
    plt.savefig(f'{results_dir}/figures/residuals/DIVA_Oxygen_residual_{label}_{year}_{season}.png')
    
    plt.close()