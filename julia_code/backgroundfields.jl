# ### Oxygen maps
# ### 2023
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# calculates background fields for ten year periods, winter, spring, summer, autumn. Saves to netcdf

# #### Add necessary packages
#using Pkg
#Pkg.add("JLD")

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
location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
outputdir = joinpath(location, "resultat/nc/$(savevar)/");
if !isdir(outputdir)
    mkpath(outputdir)
end
# Figures
figdir = "./resultat/figures/$(savevar)/";
if !isdir(figdir)
    mkpath(figdir)
    mkpath(joinpath(figdir, "test"))
end

# ## Load data big files created by program "syrekartor_data_proc"
#fname = "SHARK_EMODNET.txt"
data_fname = "EMODNET_SHARK_ICES.txt"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(location, "data/$data_fname"));

dx, dy = 0.05, 0.05         #~5km?
lonr = 9.:dx:31.
latr = 53.5:dy:61.0
timerange = [Date(1960,1,1),Date(2021,12,31)];
depthr = [0., 10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];

#BACKGROUND
# year and month-list for background analysis
yearlist = [1960:1969,1970:1979,1980:1989,1990:1999,2000:2009,2010:2021];
month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];  # Seasonal climatology
TSbackground = DIVAnd.TimeSelectorYearListMonthList(yearlist,month_list);
seasons=["Winter","Spring","Summer","Autumn"]
months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

"BATHYMETRY"
bathname = joinpath(location, "bathymetry/gebco_30sec_4.nc")
#bathname = joinpath(location, "bathymetry/bat_elevation_Baltic_Sea_masked.nc")
bath_file_name = split(bathname,"/")[end]
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

# ## MASK
xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

label = DIVAnd.floodfill(mmask);
new_mask = (label .== 1);

# #### 4) mask for Skagerrak
grid_bx = [i for i in xmask, j in ymask];
grid_by = [j for i in xmask, j in ymask];
mask_edit = copy(new_mask);
#mask_edit = copy(mmask);
sel_mask1 = (grid_by .>= 57.75) .& (grid_bx .<= 12.2);
sel_mask2 = (grid_by .>= 57.4) .& (grid_by .< 57.75) .& (grid_bx .<= 10.4);
sel_mask3 = (grid_by .>= 57.) .& (grid_by .< 57.4) .& (grid_bx .<= 10.);
sel_mask4 = (grid_by .>= 60.2) .& (grid_bx .>= 15.) .& (grid_bx .<= 25.);;
new_mask = mask_edit .* .!sel_mask1 .* .!sel_mask2 .* .!sel_mask3 .* .!sel_mask4;

# ## Analysis parameters
sz = (length(lonr), length(latr), length(depthr));
lenx = fill(50_000.,sz)   # 200 km
leny = fill(50_000.,sz)   # 200 km
#lenz = [min(max(10.,depthr[k]/150),300.) for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]]
lenz = fill(10,sz);      # 25 m
len = (lenx, leny, lenz);
epsilon = 0.1;

w_depth = 5.
w_days = 2.

rdiag_jldfile = joinpath(location, "data/$(data_fname)_weighted_$(w_depth)_$(w_days).jld")

if isfile(rdiag_jldfile)
    @load rdiag_jldfile rdiag
    @info "Loading saved rdiag file with w_depth = $(w_depth) and w_days = $(w_days)"
else
    @info "Calculating rdiag with w_depth = $(w_depth) and w_days = $(w_days)!"
    @time rdiag = 1 ./ DIVAnd.weight_RtimesOne((obslon,obslat,obsdepth,float.(Dates.dayofyear.(obstime))),(0.10,0.10,w_depth,w_days));
    @save rdiag_jldfile rdiag
end

@show maximum(rdiag),mean(rdiag)
epsilon_weighted = epsilon * rdiag;


# To include December & Nov from previous year in the analyse
obstime_shifted = copy(obstime)
obstime_shifted[Dates.month.(obstime) .== 12 .& Dates.month.(obstime) .== 11] .+= Dates.Year(1)

filenamebackground = joinpath(outputdir, "$(replace(varname,' '=>'_'))_background_weighted_0.05_field_$(bath_file_name)")

dbinfo = @time diva3d((lonr,latr,depthr,TSbackground),
           (obslon,obslat,obsdepth,obstime_shifted), obsval,
           len, epsilon_weighted,
           filenamebackground,varname,
           bathname=bathname,
           mask = new_mask,
           fitcorrlen = false,
           niter_e = 1,
           solver = :direct,
           MEMTOFIT = 120,
       );


