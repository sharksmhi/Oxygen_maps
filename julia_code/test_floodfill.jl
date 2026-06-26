# ### Oxygen maps
# ### 2023
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# calculates background fields for ten year periods, winter, spring, summer, autumn. Saves to netcdf
# Run seperate for each basin...

# #### Add necessary packages
#using Pkg
#Pkg.add("Plots")

using DIVAnd
using PyPlot
using NCDatasets
using Dates
using Statistics
using DelimitedFiles
using DataStructures
using Printf
using Missings
using JLD
using JSON
using DataFrames
using CSV

# ## Configuration
# * Define variabel and set horizontal, vertical and temporal resolutions.
### Oxygen ###
varname = "Oxygen"
savevar = "O2"
sdnp35 = "SDN:P35::EPC00002"
sdnp02 = "SDN:P02::DOXY"
doi = "10.6092/B5BB9EA1-4BDA-48A1-92B6-03645AC12FAE"
unit = "umol/l";

# ## Where to save the result. Create path if not there.
# File name based on the variable (but all spaces are replaced by _) _varlenz
# NC-files
location = "C:/LenaV/code/DIVAnd/Oxygen_maps/"
location = "C:/Work/DIVAnd/Oxygen_maps/"
location = ""

freja_location = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/"
freja_location = "freja/inprut_dur"
freja_location = "/nobackup/smhid20/proj/fouo/oxygen_indicator_2024/Oxygen_maps/"

on_freja = false
if isdir(freja_location)
    location = freja_location
    on_freja = true
end

outputdir = joinpath(location, "data/");
if !isdir(outputdir)
    @info("directory path for data does not exist $(outputdir)")
end

# ## Load data big files created by program "data_handling"
#data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_260325"
#data_fname = "SHARK_SYKE_IOW_EMODNET_ICES_250619"
#data_fname = "mat_file_1960_2024_reordered"
#@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(location, "data/$data_fname.txt"));

#Bottniska viken
#basin = "Gulf_of_Bothnia"

#Eg Östersjön o Kattegatt
#basin ="Baltic_Proper"

#basin ="test"

#Kattegatt
basin = "Kattegat"

# Läs in filens innehåll som en sträng
json_content = read("./settings.json", String)

# Parsar strängen som JSON
settings = JSON.parse(json_content)

# Exempel på att hämta specifika värden
lonr_range = settings[basin]["lonr"]
latr_range = settings[basin]["latr"]
dx = Float64.(settings[basin]["dx"])
dy = Float64.(settings[basin]["dy"])
# Skapa intervall (range) i Julia med angivet dx
lonr_min = lonr_range[1]
lonr_max = lonr_range[2]
latr_min = latr_range[1]
latr_max = latr_range[2]
lonr = lonr_min:dx:lonr_max  # 14:0.05:31
latr = latr_min:dy:latr_max

layers = settings[basin]["layers"]
depthr = Float64[layer[1] for layer in layers]
lenz_ = Float64[layer[2] for layer in layers]

lenf = Float64.(settings[basin]["lenf"])
threshold_list = settings[basin]["threshold_list"]

#BATHYMETRY
bathname = joinpath(location, "data", "bat_elevation_Baltic_Sea_masked.nc")
bath_file_name = split(bathname,"/")[end]
bathisglobal = true;
#@time bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

# ## MASK
@time xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

label = DIVAnd.floodfill(mmask);
@show size(label)

for k in axes(label, 3)
    labk = label[:, :, k]

    n_label1 = count(labk .== 1)
    labs_gt1 = filter(>(1), unique(labk))
    total_gt1 = count(labk .> 1)

    println("Depth=$(depthr[k]) m | label==1: $n_label1 | label>1: $total_gt1 celler i $(length(labs_gt1)) områden")
end


new_mask = (label .== 1);
not_used = (label .> 1);

# #### 4) mask for Skagerrak
grid_bx = [i for i in xmask, j in ymask];
grid_by = [j for i in xmask, j in ymask];
mask_edit = copy(new_mask);
#mask_edit = copy(mmask);                                                           #Maska bort:
sel_mask1 = (grid_by .>= 57.75) .& (grid_bx .<= 12.2);                              #Skagerrak
sel_mask2 = (grid_by .>= 57.4) .& (grid_by .< 57.75) .& (grid_bx .<= 10.4);         #Skagerrak
sel_mask3 = (grid_by .>= 57.) .& (grid_by .< 57.4) .& (grid_bx .<= 10.);            #Skagerrak

## If basin is....exclude other basins.
if basin == "Gulf_of_Bothnia"
    sel_mask4 = (grid_by .<= latr_min)        # keep all except north of min_lat (Åland Sea)
elseif basin == "test"
    sel_mask4 = (grid_by .>= latr_max)        # keep all except north of max_lat (Åland sea)
elseif basin == "Kattegat"
    sel_mask4 = (grid_bx .>= lonr_max);        # keep all except east of lonr_max (Arkona/Bornholm)
end
new_mask = mask_edit .* .!sel_mask1 .* .!sel_mask2 .* .!sel_mask3 .* .!sel_mask4;

NCDataset("new_mask.nc", "c") do ds
    defDim(ds, "lon", length(xmask))
    defDim(ds, "lat", length(ymask))
    defDim(ds, "depth", length(depthr))

    vlon = defVar(ds, "lon", Float64, ("lon",))
    vlat = defVar(ds, "lat", Float64, ("lat",))
    vdep = defVar(ds, "depth", Float64, ("depth",))
    vmask = defVar(ds, "new_mask", Int8, ("lon", "lat", "depth"))

    vlon[:] = collect(xmask)
    vlat[:] = collect(ymask)
    vdep[:] = collect(depthr)
    vmask[:, :, :] = Int8.(new_mask)
end

# ## Analysis parameters
# sz = (length(lonr), length(latr), length(depthr));
# lenx = fill(lenf,sz)   # 200 km
# leny = fill(lenf,sz)   # 200 km
# lenz =  [lenz_[k] for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]];
# len = (lenx, leny, lenz);

k = 10

figure(figsize=(10,8))
pcolormesh(
    xmask,
    ymask,
    permutedims(Float64.(label[:,:,k]))
)
colorbar()
title(" Depth=$(depthr[k]) m")
savefig("floodfill_labelmmask_$(depthr[k]).png")

figure(figsize=(10,8))
pcolormesh(
    xmask,
    ymask,
    permutedims(Float64.(new_mask[:,:,k]))
)
colorbar()
title("new_mask, depth=$(depthr[k]) m")
savefig("floodfill_new_mask_$(depthr[k]).png")

