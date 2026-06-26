import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

freja_input_dir = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/data/"
data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_260624.txt"
bath_fname = "bathymetry/emodnet_bathymetry_merged.nc"

out_flagged_fname = "measurements_probably_wrong_position.txt"
out_nominal_fname = "measurements_possible_nominal_position.txt"
out_clean_all_fname = "SHARK_SYKE_IOW_EMODNET_ICES_260624_cleaned_all.txt"

search_radius_m = 1000.0

df = pd.read_csv(
    freja_input_dir + data_fname,
    sep=r"\s+",
    header=None,
    usecols=[0, 1, 2, 3, 9, 10],
    names=["lon", "lat", "oxygen", "depth", "obstime", "obsid"],
    engine="python"
)

obs_lon = df["lon"].astype(float).values
obs_lat = df["lat"].astype(float).values

ds = xr.open_dataset(freja_input_dir + bath_fname)

lon = ds["lon"].values
lat = ds["lat"].values
bat = ds["bat"].values

# EMODnet elevation är negativ i havet, gör botten till positivt djup
bottom = -bat

bottom_da = xr.DataArray(
    bottom,
    coords={"lat": lat, "lon": lon},
    dims=("lat", "lon")
)

obs_points = xr.Dataset(
    {
        "lon": ("obs", obs_lon),
        "lat": ("obs", obs_lat),
    }
)

bottom_at_obs = bottom_da.sel(
    lon=obs_points["lon"],
    lat=obs_points["lat"],
    method="nearest"
).values

df["bottom_depth"] = bottom_at_obs
df["diff"] = df["depth"] - df["bottom_depth"]

# Klassning mot närmaste bottendjup
df["below_bottom"] = df["diff"] > 0
df["slightly_below_bottom"] = (df["diff"] > 0) & (df["diff"] <= 10)
df["clearly_below_bottom"] = df["diff"] > 10

# ------------------------------------------------------------
# Kolla om det finns djupare botten inom en radie runt varje obs
# ------------------------------------------------------------

max_bottom_depth_radius = np.full(len(df), np.nan)

below_indices = df.index[df["below_bottom"]].to_numpy()

print("Beräknar max bottendjup inom", search_radius_m, "m för observationer under botten...")

for idx in below_indices:
    olon = df.at[idx, "lon"]
    olat = df.at[idx, "lat"]

    dlat = search_radius_m / 111_000.0
    dlon = search_radius_m / (111_000.0 * np.cos(np.deg2rad(olat)))

    lat_sel = (lat >= olat - dlat) & (lat <= olat + dlat)
    lon_sel = (lon >= olon - dlon) & (lon <= olon + dlon)

    if not lat_sel.any() or not lon_sel.any():
        continue

    bottom_window = bottom[np.ix_(lat_sel, lon_sel)]

    lon_window = lon[lon_sel]
    lat_window = lat[lat_sel]

    lon2d, lat2d = np.meshgrid(lon_window, lat_window)

    dx = (lon2d - olon) * 111_000.0 * np.cos(np.deg2rad(olat))
    dy = (lat2d - olat) * 111_000.0
    dist = np.sqrt(dx**2 + dy**2)

    in_radius = dist <= search_radius_m

    vals = bottom_window[in_radius]
    vals = vals[np.isfinite(vals)]

    if vals.size > 0:
        max_bottom_depth_radius[idx] = np.nanmax(vals)

df["max_bottom_depth_1km"] = max_bottom_depth_radius
df["diff_to_max_bottom_1km"] = df["depth"] - df["max_bottom_depth_1km"]

df["possible_nominal_position"] = (
    df["below_bottom"]
    & df["max_bottom_depth_1km"].notna()
    & (df["depth"] <= df["max_bottom_depth_1km"])
)

bad = df[df["below_bottom"]].copy()

print("Antal observationer:", len(df))
print("Antal under botten:", len(bad))
print("Antal <= 10 m under botten:", df["slightly_below_bottom"].sum())
print("Antal > 10 m under botten:", df["clearly_below_bottom"].sum())
print("Antal möjliga nominella positioner:", df["possible_nominal_position"].sum())
print()
print("Diff-statistik, depth - bottom_depth:")
print(df["diff"].describe())

print(
    bad[[
        "lon", "lat", "depth", "bottom_depth", "diff",
        "max_bottom_depth_1km", "diff_to_max_bottom_1km",
        "possible_nominal_position", "obstime", "obsid"
    ]].head(30)
)

# Spara observationer som kan bero på nominell position
df_nominal_position = df[df["possible_nominal_position"]].copy()

df_nominal_position.to_csv(
    freja_input_dir + out_nominal_fname,
    sep="\t",
    index=False
)

print("Sparade möjliga nominella positioner till:", freja_input_dir + out_nominal_fname)

# Mättillfälle definieras av position + tid
df["lon_group"] = df["lon"].round(5)
df["lat_group"] = df["lat"].round(5)

group_cols = ["lon_group", "lat_group", "obstime"]

# Identifiera mättillfällen med flera djup som är tydligt under botten
# men exkludera observationer som kan förklaras av nominell position.
clearly_wrong_for_grouping = (
    df["clearly_below_bottom"]
    & ~df["possible_nominal_position"]
)

clearly_below_counts = (
    df[clearly_wrong_for_grouping]
    .groupby(group_cols)
    .size()
    .rename("n_depths_clearly_below_bottom")
    .reset_index()
)

wrong_position_events = clearly_below_counts[
    clearly_below_counts["n_depths_clearly_below_bottom"] > 1
]

df_wrong_position = df.merge(
    wrong_position_events[group_cols],
    on=group_cols,
    how="inner"
)

df_wrong_position.to_csv(
    freja_input_dir + out_flagged_fname,
    sep="\t",
    index=False
)

print()
print("Antal mättillfällen med trolig felposition:", len(wrong_position_events))
print("Antal rader i filen med trolig felposition:", len(df_wrong_position))
print("Sparade troliga felpositioner till:", freja_input_dir + out_flagged_fname)

# Lägg till nytt provtagningsdjup där diffen är <= 10 m under botten
df["sampling_depth_corrected"] = df["depth"]

mask_adjust_depth = df["slightly_below_bottom"]

df.loc[mask_adjust_depth, "sampling_depth_corrected"] = (
    df.loc[mask_adjust_depth, "bottom_depth"] + 1
)

print()
print("Antal observationer som fått nytt provtagningsdjup:", mask_adjust_depth.sum())

# Ta bort:
# 1. observationer där diff > 10 m och som inte kan förklaras av djupare område inom 1 km
# 2. hela mättillfällen där flera observationer har diff > 10 m och inte kan förklaras av nominell position
wrong_position_index = df.set_index(group_cols).index.isin(
    wrong_position_events.set_index(group_cols).index
)

mask_remove = (
    (df["clearly_below_bottom"] & ~df["possible_nominal_position"])
    | wrong_position_index
)

df_clean = df.loc[~mask_remove].copy()

df_clean = df_clean.drop(columns=["lon_group", "lat_group"])


df_clean.to_csv(
    freja_input_dir + out_clean_all_fname,
    sep="\t",
    index=False
)


# Skapa fil i samma format som används av DIVAnd/loadbigfile()
df_clean_original_format = pd.DataFrame({
    "obslon": df_clean["lon"],
    "obslat": df_clean["lat"],
    "obsval": df_clean["oxygen"],
    "obsdepth": df_clean["sampling_depth_corrected"].round(1),
    "obsdepth1": df_clean["sampling_depth_corrected"].round(1),
    "obsdepth2": df_clean["sampling_depth_corrected"].round(1),
    "obsdepth3": df_clean["sampling_depth_corrected"].round(1),
    "obsdepth4": df_clean["sampling_depth_corrected"].round(1),
    "obsdepth5": df_clean["sampling_depth_corrected"].round(1),
    "obstime": df_clean["obstime"],
    "obsid": df_clean["obsid"],
})


df_clean_original_format = df_clean_original_format.rename(
    columns={"sampling_depth_corrected": "depth"}
)

out_clean_original_no_header = "SHARK_SYKE_IOW_EMODNET_ICES_260624_cleaned.txt"
out_clean_original_with_header = "SHARK_SYKE_IOW_EMODNET_ICES_260624_cleaned_with_header.txt"

df_clean_original_format.to_csv(
    freja_input_dir + out_clean_original_no_header,
    sep="\t",
    index=False,
    header=False
)

df_clean_original_format.to_csv(
    freja_input_dir + out_clean_original_with_header,
    sep="\t",
    index=False,
    header=True
)

print("Sparade cleaned i originalformat utan header till:", freja_input_dir + out_clean_original_no_header)
print("Sparade cleaned i originalformat med header till:", freja_input_dir + out_clean_original_with_header)


print("Antal rader borttagna:", mask_remove.sum())
print("Antal rader kvar i rensad fil:", len(df_clean))
print("Andel borttagna rader i procent:", round(mask_remove.sum()/len(df_clean)*100)) 
print("Sparade rensad datafil till:", freja_input_dir + out_clean_original_no_header)

# ------------------------------------------------------------
# Plot 1: karta där botten och observationsdjup har samma färgskala
# ------------------------------------------------------------

bad_plot = bad.dropna(subset=["bottom_depth", "diff"]).copy()

vmin = 0
vmax = 200

plt.figure(figsize=(11, 9))

pcm = plt.pcolormesh(
    lon,
    lat,
    bottom,
    shading="auto",
    cmap="viridis",
    vmin=vmin,
    vmax=vmax
)

plt.scatter(
    bad_plot["lon"],
    bad_plot["lat"],
    c=bad_plot["depth"],
    s=10,
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    edgecolors="k",
    linewidths=0.2
)

plt.colorbar(pcm, label="Bottendjup / observationsdjup (m)")

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Observationer under botten: stationdjup och bottendjup med samma färgskala")
plt.tight_layout()

# ------------------------------------------------------------
# Plot 1b: karta över diff, bra för att se stora avvikelser
# ------------------------------------------------------------

plt.figure(figsize=(11, 9))

plt.pcolormesh(
    lon,
    lat,
    bottom,
    shading="auto",
    cmap="Greys",
    vmin=0,
    vmax=200
)

sc = plt.scatter(
    bad_plot["lon"],
    bad_plot["lat"],
    c=bad_plot["diff"],
    s=10,
    cmap="Reds",
    vmin=0,
    vmax=50
)

plt.colorbar(sc, label="Depth - bottom_depth (m)")

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Observationer under botten färgade efter diff")
plt.tight_layout()

# ------------------------------------------------------------
# Plot 2: zoomat histogram
# ------------------------------------------------------------

plt.figure(figsize=(10, 6))

plt.hist(
    df["diff"].dropna().clip(0, 50),
    bins=np.arange(0, 50.5, 0.5)
)

plt.axvline(0, color="k", linestyle="--", label="Botten")
plt.axvline(10, color="r", linestyle="--", label="10 m")

plt.xlabel("Depth - bottom_depth (m)")
plt.ylabel("Antal observationer")
plt.title("Fördelning av avstånd till botten (zoom)")
plt.legend()
plt.tight_layout()

# ------------------------------------------------------------
# Plot 3: större diff-intervall
# ------------------------------------------------------------

plt.figure(figsize=(10, 6))

max_diff = np.nanmax(df["diff"].values)

plt.hist(
    df["diff"].dropna(),
    bins=np.arange(0, max_diff + 1, 1)
)

plt.axvline(0, color="k", linestyle="--", label="Botten")
plt.axvline(10, color="r", linestyle="--", label="10 m")

plt.xlabel("Depth - bottom_depth (m)")
plt.ylabel("Antal observationer")
plt.title("Fördelning av avstånd till botten")
plt.legend()
plt.tight_layout()

plt.show()

# ------------------------------------------------------------
# Diagnostik
# ------------------------------------------------------------

test_lon = 15.3333
test_lat = 55.3833

nearest = bottom_da.sel(
    lon=test_lon,
    lat=test_lat,
    method="nearest"
)

print("Testpunkt:")
print("lon obs:", test_lon)
print("lat obs:", test_lat)
print("nearest lon:", float(nearest["lon"].values))
print("nearest lat:", float(nearest["lat"].values))
print("bottom depth:", float(nearest.values))

print("lon stigande:", np.all(np.diff(lon) > 0))
print("lat stigande:", np.all(np.diff(lat) > 0))
print("lat första/sista:", lat[0], lat[-1])
print("lon första/sista:", lon[0], lon[-1])

print("lon-dimension:", len(lon))
print("lat-dimension:", len(lat))

print("lon upplösning (grader):", np.unique(np.round(np.diff(lon), 6)))
print("lat upplösning (grader):", np.unique(np.round(np.diff(lat), 6)))

# Bornholmsbassängen, ungefärlig box
lon_min, lon_max = 13.5, 16.5
lat_min, lat_max = 54.5, 56.5

lon_mask = (lon >= lon_min) & (lon <= lon_max)
lat_mask = (lat >= lat_min) & (lat <= lat_max)

lon_sub = lon[lon_mask]
lat_sub = lat[lat_mask]
bottom_sub = bottom[np.ix_(lat_mask, lon_mask)]

plt.figure(figsize=(9, 7))

pcm = plt.pcolormesh(
    lon_sub,
    lat_sub,
    bottom_sub,
    shading="auto",
    cmap="viridis",
    vmin=0,
    vmax=120
)

cs = plt.contour(
    lon_sub,
    lat_sub,
    bottom_sub,
    levels=[40, 50, 60, 70, 80, 90, 100],
    colors="k",
    linewidths=0.5
)

plt.clabel(cs, fmt="%d m", fontsize=8)

plt.scatter(11.758, 56.958, s=40, color="red", label="Testpunkt")

plt.colorbar(pcm, label="Bottendjup (m)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Bathymetri på position")
plt.legend()
plt.tight_layout()

plt.show()