#=
syrekartor_data_proc:
- Julia version: 
- Author: a001109
- Date: 2023-08-30

Program to handle O2 data before DIVAnd.
Emodnet, SHARKweb anb ICES data
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
@show("Loading SHARKWEB...")
fname_shark = joinpath(location, "data/all_baltic/sharkweb_btlctd_02_240603.txt")
@time obsval_shark,obslon_shark,obslat_shark,obsdepth_shark,obstime_shark,obsid_shark = loadbigfile(fname_shark);

@show("Loading EMODNET BTL...")
datafile_emod_btl = joinpath(location, "data/all_baltic/BTLdata_from_ALL_emodnet-chem_v2023_2022.txt")
@time obsval_emod_btl,obslon_emod_btl,obslat_emod_btl,obsdepth_emod_btl,obstime_emod_btl,obsid_emod_btl = ODVspreadsheet.load(Float64,[datafile_emod_btl],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );

@show("Loading EMODNET CTD...")
datafile_emod_ctd = joinpath(location, "data/all_baltic/CTDdata_from_ALL_emodnet-chem_v2023_2022.txt")
@time obsval_emod_ctd,obslon_emod_ctd,obslat_emod_ctd,obsdepth_emod_ctd,obstime_emod_ctd,obsid_emod_ctd = ODVspreadsheet.load(Float64,[datafile_emod_ctd],
                           ["Water body dissolved oxygen concentration"]; nametype = :localname );
@show("Loading ICES...")
datafile_ices_btlctd = joinpath(location, "data/all_baltic/ICES_btl_lowres_ctd_02_NEW.txt")
@time obsval_ices_btlctd,obslon_ices_btlctd,obslat_ices_btlctd,obsdepth_ices_btlctd,obstime_ices_btlctd,obsid_ices_btlctd = loadbigfile(datafile_ices_btlctd);

@show(length(obsval_emod_btl));
@show(length(obsval_emod_ctd));
@show(length(obsval_shark));
@show(length(obsval_ices_btlctd));

# ## Remove low res CTD when BTL is available.
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

# ## Remove SHARK_data when EMODnet_data is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.01
# Vertical separation: 0.01 m depth
depth_dist= 0.01
#Time separation: 10 minute.
time_sep = 10
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 1

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emod,obslat_emod,obsdepth_emod,obstime_emod), obsval_emod,
    (obslon_shark,obslat_shark,obsdepth_shark,obstime_shark),obsval_shark,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_shark) * 100; digits=2);
@info("Number of possible duplicates in emodnet/SHARKweb: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints_shark = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon_emodshark = [obslon_emod; obslon_shark[newpoints_shark]];
obslat_emodshark = [obslat_emod; obslat_shark[newpoints_shark]];
obsdepth_emodshark = [obsdepth_emod; obsdepth_shark[newpoints_shark]];
obstime_emodshark = [obstime_emod; obstime_shark[newpoints_shark]];
obsval_emodshark = [obsval_emod; obsval_shark[newpoints_shark]];
obsid_emodshark = [obsid_emod; obsid_shark[newpoints_shark]];

# ## Remove ICES data when EMODnet_SHARK is available.
# Remove true duplicates, hence when exactly the same data is found in both datasets.
# ## Criteria (can be adapted according to the application):
# Horizontal distance: 0.01 degree (about 1km)
xy_dist = 0.05
# Vertical separation: 0.01 m depth
depth_dist= 1
#Time separation: 10 minute.
time_sep = 60
#obsval difference: 1 µmol/l correspond to ~pyttelite ml/l
obsval_diff = 1

@time dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon_emodshark,obslat_emodshark,obsdepth_emodshark,obstime_emodshark), obsval_emodshark,
    (obslon_ices_btlctd,obslat_ices_btlctd,obsdepth_ices_btlctd,obstime_ices_btlctd),obsval_ices_btlctd,
    (xy_dist,xy_dist,depth_dist,time_sep/(24*60)),obsval_diff);

# ## Find the indices of the possible duplicates:
index = findall(.!isempty.(dupl));
ndupl = length(index);
pcdupl = round(ndupl / length(obslon_ices_btlctd) * 100; digits=2);
@info("Number of possible duplicates in emodnetSHARKweb/ICES: $ndupl")
@info("Percentage of duplicates: $pcdupl%")
# ## If you decide to combine the 2 (or more) datasets:
newpoints = isempty.(dupl);
@info("Number of new points: $(sum(newpoints)))")

obslon = [obslon_emodshark; obslon_ices_btlctd[newpoints]];
obslat = [obslat_emodshark; obslat_ices_btlctd[newpoints]];
obsdepth = [obsdepth_emodshark; obsdepth_ices_btlctd[newpoints]];
obstime = [obstime_emodshark; obstime_ices_btlctd[newpoints]];
obsval = [obsval_emodshark; obsval_ices_btlctd[newpoints]];
obsid = [obsid_emodshark; obsid_ices_btlctd[newpoints]];

# ## Create a plot showing the additional data points:
figure("Additional-Data")
aspectratio = 1/cos(mean(54:0.05:61) * pi/180)
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
ylim(53, 66);
xlim(8, 31);
#contourf(bx, by, permutedims(Float64.(mask_edit[:,:,1]),[2,1]),
#    levels=[-1e5,0],cmap="binary");
plot(obslon_emod, obslat_emod, "bo", markersize=.2, label="Emodnet")
plot(obslon_shark[newpoints_shark], obslat_shark[newpoints_shark], "go",
   markersize=.2, label="Additional data\nfrom SHARKweb")
plot(obslon_ices_btlctd[newpoints], obslat_ices_btlctd[newpoints], "ro", mfc="none",
   markersize=.2, label="Additional data\nfrom ICES BTL lowres CTD")

legend(loc=3, fontsize=4)
gca().set_aspect(aspectratio)
figname = "additional_data.png"
@show joinpath(figdir,"$(figname)")
PyPlot.savefig(joinpath(figdir,"$(figname)"), dpi=300);
PyPlot.close_figs()

df  = DataFrame(obslon=obslon,obslat=obslat,obsval=obsval,obsdepth=obsdepth,obsdepth1=obsdepth,obsdepth2=obsdepth,obsdepth3=obsdepth,obsdepth4=obsdepth,obsdepth5=obsdepth,obstime=obstime,obsid=obsid)
filename = "EMODNET_SHARK_ICES_240620"
CSV.write(joinpath(outputdir, "$(filename).txt"), df, delim="\t", writeheader=false)
#DIVAnd.saveobs(joinpath(outputdir, "$(filename).nc"),varname, obsval, (obslon,obslat,obsdepth,obstime),obsid)


