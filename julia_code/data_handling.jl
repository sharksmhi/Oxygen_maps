#=
syrekartor_data_proc:
- Julia version: 
- Author: a001109
- Date: 2023-08-30

Program to handle O2 data before DIVAnd.
Emodnet, SHARKweb, ICES and SYKE data
-Load data
-Removes duplicates, with the help of DIVA-funktions.
-Handles both BTL and CTD-data. BTL-data is prioritiz
-Merge to one dataset
# ### Lena Viktorsson & Martin Hansson
=#
#using Pkg
#Pkg.add("PyCall")

# #### Add necessary packages
using DIVAnd
using PyPlot
using NCDatasets
using Dates
using Statistics
using DelimitedFiles
using DataStructures
using Printf
using Missings
using DataFrames
using CSV
using PyCall

function remove_matching_rows(data::DataFrame, bad_data::DataFrame, columns::Vector{Symbol})
    """
    Removes rows from `data` DataFrame where the specified `columns` match in the `bad_data` DataFrame.
    Args:
        data (DataFrame): Original dataset
        bad_data (DataFrame): Dataset describing rows to be removed
        columns (Vector{Symbol}): List of column names to use for comparison
    Returns:
        DataFrame: Filtered DataFrame with rows removed
    """
    # Skapa en Set av rader från `bad_data` för de angivna kolumnerna
    rows_to_remove = Set([NamedTuple(bad_data[i, columns]) for i in 1:nrow(bad_data)])

    # Filtrera bort rader i `data` som matchar `rows_to_remove`
    filtered_rows = [!(NamedTuple(data[i, columns]) in rows_to_remove) for i in 1:nrow(data)]
    filtered_data = data[filtered_rows, :]

    return filtered_data
end

# ## Configuration
### Oxygen ###
varname = "Oxygen"
savevar = "O2"
sdnp35 = "SDN:P35::EPC00002"
sdnp02 = "SDN:P02::DOXY"
doi = "10.6092/B5BB9EA1-4BDA-48A1-92B6-03645AC12FAE"
unit = "umol/l";

# ## Where to save the result. Create path if not there.
# File name based on the variable (but all spaces are replaced by _) _varlenz
# data-files
#location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
location = "C:/Work/DIVAnd/Oxygen_maps/"
outputdir = joinpath(location,"data/");
if !isdir(outputdir)
    mkpath(outputdir)
end
# Figures
figdir = joinpath(location,"resultat/figures/");

if !isdir(figdir)
    mkpath(figdir)
    #mkpath(joinpath(figdir, "test"))
end

# ## Load data big files
# emodnet BTL + CTD
@show("Loading SHARK BTL/CTD...")
fname_shark = joinpath(location, "data/all_baltic/sharkweb_btlctd_02_241107.txt")
@time obsval_shark,obslon_shark,obslat_shark,obsdepth_shark,obstime_shark,obsid_shark = loadbigfile(fname_shark);

@show("Loading IOW BTL/CTD...")
fname_iow = joinpath(location, "data/all_baltic/IOW_1959_2024.txt")
@time obsval_iow,obslon_iow,obslat_iow,obsdepth_iow,obstime_iow,obsid_iow = loadbigfile(fname_iow);
obsid_iow = string.("IOW-",obsid_iow)

@show("Loading EMODNET BTL...")
datafile_emod_btl = joinpath(location, "data/EMODNET_2024/BTL_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt")
@time obsval_emod_btl,obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl,obsid_emod_btl = ODVspreadsheet.load(Float64,[datafile_emod_btl],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
# obsid_emod_btl = LOCAL_CDI_ID + EMOD_CODE
obsid_emod_btl = string.("emod_btl-",obsid_emod_btl)

@show("Loading EMODNET CTD...")
datafile_emod_ctd = joinpath(location, "data/EMODNET_2024/CTD_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt")
@time obsval_emod_ctd,obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd,obsid_emod_ctd = ODVspreadsheet.load(Float64,[datafile_emod_ctd],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
obsid_emod_ctd = string.("emod_ctd-",obsid_emod_ctd)
@show(obsid_emod_ctd[1])

@show("Loading ICES...")
datafile_ices_btlctd = joinpath(location, "data/all_baltic/ICES_btl_lowres_ctd_02_241107.txt")
@time obsval_ices,obslon_ices,obslat_ices,obsdepth_ices,obstime_ices,obsid_ices = loadbigfile(datafile_ices_btlctd);
obsid_ices = string.("ICES-", obsid_ices)
@show(obsid_ices[1])

@show("Loading SYKE data...")
datafile_syke_btlctd = joinpath(location, "data/all_baltic/syke_data_no_header_241107.txt")
@time obsval_syke,obslon_syke,obslat_syke,obsdepth_syke,obstime_syke,obsid_syke = loadbigfile(datafile_syke_btlctd);
@show(obsid_syke[1])

@show("Loading BAD-data file...data to be removed")
datafile_bad_data = joinpath(location, "data/all_baltic/bad_data.txt")
bad_data = CSV.read(datafile_bad_data, DataFrame)

@show(length(obsval_shark));
@show(length(obsval_iow));
@show(length(obsval_emod_btl));
@show(length(obsval_emod_ctd));
@show(length(obsval_ices));
@show(length(obsval_syke));

# ## EMODNET: Remove low res CTD when BTL is available.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.01
#Time separation: 60 minute.
time_sep = 60
#obsval difference: 90 µmol/l correspond to ~2 ml/l
obsval_diff = 90

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl), obsval_emod_btl,
    (obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd),obsval_emod_ctd,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_emod_ctd) * 100; digits=2);
@info("Number of possible duplicates emodnet BTL/CTD: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon_emod = [obslon_emod_btl; obslon_emod_ctd[newpoints]];
obslat_emod = [obslat_emod_btl; obslat_emod_ctd[newpoints]];
obsdepth_emod = [obsdepth_emod_btl; obsdepth_emod_ctd[newpoints]];
obstime_emod = [obstime_emod_btl; obstime_emod_ctd[newpoints]];
obsval_emod = [obsval_emod_btl; obsval_emod_ctd[newpoints]];
obsid_emod = [obsid_emod_btl; obsid_emod_ctd[newpoints]];

# ## Remove SYKE_data when SHARK_data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.1
#Time separation: 10 minute.
time_sep = 10
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 100

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_shark,obslat_shark,obsdepth_shark,obstime_shark), obsval_shark,
    (obslon_syke,obslat_syke,obsdepth_syke,obstime_syke),obsval_syke,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_syke) * 100; digits=2);
@info("Number of possible duplicates in SHARKweb/SYKE: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_syke = isempty.(dupl);
@info("Number of new points: $(sum(newpoints_syke)))")

obslon_shark_syke = [obslon_shark; obslon_syke[newpoints_syke]];
obslat_shark_syke = [obslat_shark; obslat_syke[newpoints_syke]];
obsdepth_shark_syke = [obsdepth_shark; obsdepth_syke[newpoints_syke]];
obstime_shark_syke = [obstime_shark; obstime_syke[newpoints_syke]];
obsval_shark_syke = [obsval_shark; obsval_syke[newpoints_syke]];
obsid_shark_syke = [obsid_shark; obsid_syke[newpoints_syke]];




# ## Remove IOW_data when SHARK_SYKE_data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.1
#Time separation: 10 minute.
time_sep = 10
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 100

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_shark_syke,obslat_shark_syke,obsdepth_shark_syke,obstime_shark_syke), obsval_shark_syke,
    (obslon_iow,obslat_iow,obsdepth_iow,obstime_iow),obsval_iow,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_iow) * 100; digits=2);
@info("Number of possible duplicates in SHARKweb_SYKE/IOW: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_iow = isempty.(dupl);
@info("Number of new points: $(sum(newpoints_iow)))")

obslon_shark_syke_iow = [obslon_shark_syke; obslon_iow[newpoints_iow]];
obslat_shark_syke_iow = [obslat_shark_syke; obslat_iow[newpoints_iow]];
obsdepth_shark_syke_iow = [obsdepth_shark_syke; obsdepth_iow[newpoints_iow]];
obstime_shark_syke_iow = [obstime_shark_syke; obstime_iow[newpoints_iow]];
obsval_shark_syke_iow = [obsval_shark_syke; obsval_iow[newpoints_iow]];
obsid_shark_syke_iow = [obsid_shark_syke; obsid_iow[newpoints_iow]];

# ## Remove EMOD_data when SHARK_SYKE_IOW data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.1
#Time separation: 10 minute.
time_sep = 10
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 100

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_shark_syke_iow,obslat_shark_syke_iow,obsdepth_shark_syke_iow,obstime_shark_syke_iow), obsval_shark_syke_iow,
    (obslon_emod,obslat_emod,obsdepth_emod,obstime_emod),obsval_emod,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_emod) * 100; digits=2);
@info("Number of possible duplicates in SHARK_SYKE_IOW/EMOD: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_emod = isempty.(dupl);
@info("Number of new points: $(sum(newpoints_emod)))")

obslon_shark_syke_iow_emod = [obslon_shark_syke_iow; obslon_emod[newpoints_emod]];
obslat_shark_syke_iow_emod = [obslat_shark_syke_iow; obslat_emod[newpoints_emod]];
obsdepth_shark_syke_iow_emod = [obsdepth_shark_syke_iow; obsdepth_emod[newpoints_emod]];
obstime_shark_syke_iow_emod = [obstime_shark_syke_iow; obstime_emod[newpoints_emod]];
obsval_shark_syke_iow_emod = [obsval_shark_syke_iow; obsval_emod[newpoints_emod]];
obsid_shark_syke_iow_emod = [obsid_shark_syke_iow; obsid_emod[newpoints_emod]];

# ## Remove ICES data when SHARK_SYKE_IOW_EMOD is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.05
# Vertical separation: 0.01 m depth
depth_dist= 1
#Time separation: 10 minute.
time_sep = 60
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 100

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_shark_syke_iow_emod,obslat_shark_syke_iow_emod,obsdepth_shark_syke_iow_emod,obstime_shark_syke_iow_emod), obsval_shark_syke_iow_emod,
    (obslon_ices,obslat_ices,obsdepth_ices,obstime_ices),obsval_ices,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_ices) * 100; digits=2);
@info("Number of possible duplicates in SHARKweb_SYKE_IOW_EMOD/ICES: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_ICES = isempty.(dupl);
@info("Number of new points: $(sum(newpoints_ICES)))")

obslon = [obslon_shark_syke_iow_emod; obslon_ices[newpoints_ICES]];
obslat = [obslat_shark_syke_iow_emod; obslat_ices[newpoints_ICES]];
obsdepth = [obsdepth_shark_syke_iow_emod; obsdepth_ices[newpoints_ICES]];
obstime = [obstime_shark_syke_iow_emod; obstime_ices[newpoints_ICES]];
obsval = [obsval_shark_syke_iow_emod; obsval_ices[newpoints_ICES]];
obsid = [obsid_shark_syke_iow_emod; obsid_ices[newpoints_ICES]];

# # ## Remove SYKE data when EMODnet_SHARK_ICES data is available.
# # Remove true duplicates, hence when exactly the same data is found in both datasets.
# # ## Criteria (can be adapted according to the application):
# # Horizontal distance: 0.01 degree (about 1km)
# xy_dist = 0.05
# # Vertical separation: 0.01 m depth
# depth_dist= 1
# #Time separation: 10 minute.
# time_sep = 60
# #obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
# obsval_diff = 100
#
# @time dupl = DIVAnd.Quadtrees.checkduplicates(
#     (obslon_emodsharkices,obslat_emodsharkices,obsdepth_emodsharkices,obstime_emodsharkices), obsval_emodsharkices,
#     (obslon_syke_btlctd,obslat_syke_btlctd,obsdepth_syke_btlctd,obstime_syke_btlctd),obsval_syke_btlctd,
#     (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);
#
# # ## Find the indices of the possible duplicates:
# index = findall(.!isempty.(dupl));
# ndupl = length(index);
# pcdupl = round(ndupl / length(obslon_ices_btlctd) * 100; digits=2);
# @info("Number of possible duplicates in emodnetSHARKwebICES/SYKE: $ndupl")
# @info("Percentage of duplicates: $pcdupl%")
# # ## If you decide to combine the 2 (or more) datasets:
# newpoints_SYKE = isempty.(dupl);
# @info("Number of new points: $(sum(newpoints)))")
#
# obslon = [obslon_emodsharkices; obslon_syke_btlctd[newpoints_SYKE]];
# obslat = [obslat_emodsharkices; obslat_syke_btlctd[newpoints_SYKE]];
# obsdepth = [obsdepth_emodsharkices; obsdepth_syke_btlctd[newpoints_SYKE]];
# obstime = [obstime_emodsharkices; obstime_syke_btlctd[newpoints_SYKE]];
# obsval = [obsval_emodsharkices; obsval_syke_btlctd[newpoints_SYKE]];
# obsid = [obsid_emodsharkices; obsid_syke_btlctd[newpoints_SYKE]];


# ## Do a range check and remove outliers µomol/l
min_value = -500.0
max_value = 600.0

# Perform the range check
in_range = (obsval .>= min_value) .& (obsval .<= max_value)

# Extract values within the range (valid data)
valid_obsval = obsval[in_range]

# Extract values outside the range (invalid data)
out_of_range = .!in_range
invalid_obsval = obsval[out_of_range]

# Store the invalid values in a separate file
invalid_values_file = joinpath(location, "data/all_baltic/invalid_obsval.txt")
open(invalid_values_file, "w") do file
    for value in invalid_obsval
        write(file, "$value\n")
    end
end

# Save the valid data back to obsval_shark
obsval = valid_obsval
obslon = obslon[in_range];
obslat = obslat[in_range];
obsdepth = obsdepth[in_range];
obstime = obstime[in_range];
obsid = obsid[in_range];

@show("sets all neg o2 to 0.44661 µmol/l")
obsval[obsval .<= 0] .= -1
obsval[obsval .<= 0] .= 0.44661

# Output some statistics
@show sum(in_range)
@show sum(out_of_range)

# ## Hantera negativa värden i obsdepth_shark och gör dem positiva???
obsdepth_shark = abs.(obsdepth_shark)

# ## Create a plot showing the additional data points:
figure("Additional-Data")
aspectratio = 1/cos(mean(54:0.05:61) * pi/180)
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
ylim(53, 66);
xlim(8, 31);
#contourf(bx, by, permutedims(Float64.(mask_edit[:,:,1]),[2,1]),
#    levels=[-1e5,0],cmap="binary");
plot(obslon_ices[newpoints_ICES], obslat_ices[newpoints_ICES], "ro", mfc="none",
   markersize=.2, label="Additional data\nfrom ICES BTL lowres CTD")
plot(obslon_emod[newpoints_emod], obslat_emod[newpoints_emod], "go",
   markersize=.2, label="Additional data\nfrom EMODnet")
plot(obslon_iow[newpoints_iow], obslat_iow[newpoints_iow], "co", mfc="none",
   markersize=.2, label="Additional data\nfrom IOW")
plot(obslon_syke[newpoints_syke], obslat_syke[newpoints_syke], "yo", mfc="none",
   markersize=.2, label="Additional data\nfrom SYKE")
plot(obslon_shark, obslat_shark, "bo", markersize=.2, label="SHARK")

legend(loc=4, fontsize=4)
gca().set_aspect(aspectratio)
figname = "additional_data_SHARK_SYKE_IOW_EMOD_ICES.png"
PyPlot.savefig(joinpath(figdir,"$(figname)"), dpi=300);
PyPlot.close_figs()

# Skapar en dataframe för alla data
data  = DataFrame(obslon=obslon,obslat=obslat,obsval=obsval,obsdepth=obsdepth,obsdepth1=obsdepth,obsdepth2=obsdepth,obsdepth3=obsdepth,obsdepth4=obsdepth,obsdepth5=obsdepth,obstime=obstime,obsid=obsid)

@show describe(data)
@show describe(bad_data)
@show("Removing bad data...from det bad-data file..")

columns = [:obslon, :obslat, :obsdepth, :obstime]

# Anropa funktionen för att ta bort bad data
filtered_data = remove_matching_rows(data, bad_data, columns)

# Visa resultatet
println("Originaldata, antal rader:")
println(nrow(data))
println("\nData att ta bort, antal rader:")
println(nrow(bad_data))
println("\nFiltrerat data, antal rader:")
println(nrow(filtered_data))

# Skriver data till fil
filename = "SHARK_SYKE_IOW_EMODNET_ICES_250619"
CSV.write(joinpath(outputdir, "$(filename).txt"), filtered_data, delim="\t", writeheader=false)
CSV.write(joinpath(outputdir, "$(filename)_with_header.txt"), filtered_data, delim="\t", writeheader=true)
#DIVAnd.saveobs(joinpath(outputdir, "$(filename).nc"),varname, obsval, (obslon,obslat,obsdepth,obstime),obsid)
@show "Important: Make new background files and remove your weighting file before you do a new DIVAnd analysis!"
# Formaterar det till önskat format för Johannas/fredriks QC....
#formatted_times = [Dates.format(time, "yyyymmddHHMM") for time in obstime]
#obsval = obsval.* (22.414./1000) #Tillbaka till ml/l för QC?
#QC = ones(Int, length(formatted_times))
#NN = fill(NaN, length(formatted_times))
##                   tid,                lat,    lat_qc, lon,          lon_qc,    djup,       djup_qc,      tryck,     tryck_qc, temp, temp_qc, salt, salt_qc,chl,  chl_qc, turb, turb_qc, syre,       syre_qc, PAr, PAR_qc, ljudhast, ljudhast_qc, ID (Skall eg inte va med)
#               40014	8002	88002	8003	88003	8178	88178	8168	88168	8179	88179	8181	88181	8063	88063	8174	88174	8191	88191	8175	88175	8187	88187
#df_QC  = DataFrame(formatted_times=formatted_times, obslat=obslat, QC1=QC, obslon=obslon, QC2=QC, obsdepth=obsdepth, QC3=QC, obsdepth_tryck=obsdepth, QC4=QC, NN5=NN, QC5=QC, NN6=NN, QC6=QC, NN7=NN, QC7=QC, NN8=NN, QC8=QC, obsval=obsval, QC9=QC, NN10=NN, QC10=QC, NN11=NN, QC11=QC, obsid=obsid)
#filename = "EMODNET_SHARK_ICES_SYKE_241107_QC"
#CSV.write(joinpath(outputdir, "$(filename).txt"), df_QC, delim="\t", writeheader=false)