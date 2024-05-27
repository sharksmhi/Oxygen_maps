# ### Oxygen maps
# ### 2023-2024
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# read dataset
# set up metadata info
# run DIVAnd and save to netcdf.

# #### Add necessary packages
#using Pkg
#Pkg.add("PyCall")

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
using PyCall

args= ARGS
input_dir = args[1]
results_dir = args[2]
data_fname = args[3]
year_list = JSON.parse(args[4])
month_list = JSON.parse(args[5])
seasons = JSON.parse(args[6])
lenf = JSON.parse(args[7])
epsilon = JSON.parse(args[8])
dx = JSON.parse(args[9])
bath_file_name = args[10]
w_depth = JSON.parse(args[11])
w_days = JSON.parse(args[12])

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

# Figures
figdir = "./resultat/figures/$(savevar)/";
if !isdir(figdir)
    mkpath(figdir)
    mkpath(joinpath(figdir, "test"))
end

# ## Load data big files created by program "syrekartor_data_proc"
#data_fname = "EMODNET_SHARK_ICES"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(input_dir, "$(data_fname)"))

# Sätt horisontell uppplösning i grader
# 0.05 motsvarra ca 5km
#dx, dy = 0.125, 0.125  #Karin dx, dy = 0.1, 0.1
dy = dx #0.05, 0.05  #Karin dx, dy = 0.1, 0.1
lonr = 9.:dx:31.
latr = 53.5:dy:61.

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

# ## Extract the bathymetry
# It is used to delimit the domain where the interpolation is performed.
# Modify bathname according to the resolution required.
# Bathymetry needs to be namned "bat"
#bath_file_name = "bat_elevation_Baltic_Sea_masked"
#bath_file_name = "gebco_30sec_4"
bathname = joinpath(input_dir, "$(bath_file_name).nc")
@show(bathname)
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

#Load background field
#filenamebackground = joinpath(input_dir, "$(replace(varname,' '=>'_'))_background_weighted_0.05_field_$(bath_file_name).nc")
filenamebackground = joinpath(input_dir, "$(replace(varname,' '=>'_'))_background_weighted_$(dx)_field_$(bath_file_name).nc")

# year and month-list for background analysis
year_list_background = [1960:1969,1970:1979,1980:1989,1990:1999,2000:2009,2010:2021];
TSbackground = DIVAnd.TimeSelectorYearListMonthList(year_list_background,month_list);

# ## MASK
# #### 1) Create a mask from bathymetry and your selected depth vector
# First set the depth resolotion, this will be the depths used in DIVArun
# Then set the horizontal correlation length (should be twice the resolution)
depthr = [0.,  10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];
lenz_ =  [20., 20., 20., 10., 10., 10., 20., 20., 10., 10., 10., 10,  10., 10., 10., 10., 10.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50., 50., 50.,100., 100.,100.];
xmask,ymask,mmask = load_mask(bathname,true,lonr,latr,depthr);

# #### 2) apply floodfill.
# Karin: Varför depthr? borde gå med 0?
# It gives you an array which has the same size has the mask, and contains values equal to 0, 1, 2...
# 1 corresponds to the largest area, label 2 the 2nd largest and so on... fråga karin

# #### 3) Create a new mask
# which means it will have the value true only in the main sea area, not on land or in the small isolated pixels.
# I also attach a plot with the 2 masks (before and after the floodfill). You can see that your problem
# with the pixel in the middle of the island has disappeared. Note also that the northwesternmost part of the
# domain has also been masked, as this area was the 2nd largest sea area.
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

# Ökas slen så minskas avståndet från kusten där len är mindre.
sz = (length(lonr),length(latr),length(depthr));
lenx = fill(lenf,sz)
leny = fill(lenf,sz)
#Horisontell korrelationslängd
lenz =  [lenz_[k] for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]];
#lenz = fill(10,sz);      # 10m 25 m

len = (lenx, leny, lenz);
lx = lenf
ly = lenf

# To include December from previous year in the analyse
obstime_shifted = copy(obstime)
# Get all dates with dec and nov and add one year to these specific dates.
# Winter period will be jan-feb + nov-dec from previous year.
obstime_shifted[Dates.month.(obstime) .== 12 .| Dates.month.(obstime) .== 11] .+= Dates.Year(1)

# Settings for DIVAnd-------------------------------------------------------------------------------
error_thresholds = [("L1", 0.3), ("L2", 0.5)];

# Modify data weight
# Compute the new weights that takes into account close points.
# If the dataset is large, this can take a few minutes.
# The maximal and mean values provide an indication of the spatial proximity between the data.
# If you apply this technique, you need to adapt epsilon2:

#Make a test with different values for depth and days.
# 60m och 30 dagar - funkar ej
# 5m och 7 dagar - funkar ej
# 20m och 7 dagar - funkar ej
# 10m och 7 dagar - funkar ej
# 2.5m och 7 dagar - funkar bäst just nu
# 5m och 2 dagar - funkar bäst just nu
# 2.5m och 2 dagar
# 2.5m och 30 dagar

#w_depth = 5.
#w_days = 2.

@show data_fname
@show split(data_fname,".")[1]
rdiag_jldfile = joinpath(input_dir, "$(split(data_fname,".")[1])_weighted_$(w_depth)_$(w_days).jld")

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

# One metadata set up per season
metadata=Array{DataStructures.OrderedDict{String,Any}}(undef,4) ;

#A list of created files
file_list = []

for monthlist_index in 1:length(month_list)
    season = seasons[monthlist_index]

    @info("Creating metadata dicitonary for the season $(season)")
    metadata_season = OrderedDict(
        # set attributes for DVIA run from our settings
        "season" => season,
        "epsilon" => "$epsilon",
        "horizontal correlation length m" => "$lx",
        "start year" => string(year_list[1]),
        "end year" => string(year_list[end]),
        # Name of the project (SeaDataCloud, SeaDataNet, EMODNET-chemistry, ...)
        "project" => "EMODNET-chemistry",
        # URN code for the institution EDMO registry,
        # e.g. SDN:EDMO::1579
        "institution_urn" => "SDN:EDMO::545",
        # Production group
        "production" => "SMHI",
        # Name and emails from authorsseasons
        "Author_e-mail" => ["Martin Hansson <martin.hansson@smhi.se>"],
        # Source of the observation
        "source" => "observational data from SeaDataNet/EMODnet Chemistry Data Network",
        # Additional comment
        "comment" => "Every year of the time dimension corresponds to a 1-year centred average for one season.",

        # SeaDataNet Vocabulary P35 URN
        # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p35
        # example: SDN:P35::WATERTEMP
        "parameter_keyword_urn" => "$sdnp35", # Water body phosphate

        # List of SeaDataNet Parameter Discovery Vocabulary P02 URNs
        # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=p02
        # example: ["SDN:P02::TEMP"]
        "search_keywords_urn" => ["$sdnp02"], # Water body phosphate

        # List of SeaDataNet Vocabulary C19 area URNs
        # SeaVoX salt and fresh water body gazetteer (C19)
        # http://seadatanet.maris2.nl/v_bodc_vocab_v2/search.asp?lib=C19
        # example: ["SDN:C19::3_1"]
        "area_keywords_urn" => ["SDN:C19::2"],
        "product_version" => "2.0",
        "product_code" => "SMHI-Baltic Sea-$(replace(varname,' '=>'_'))-v2023-ANA",

        # bathymetry source acknowledgement
        # see, e.g.
        # * EMODnet Bathymetry Consortium (2016): EMODnet Digital Bathymetry (DTM).
        # https://doi.org/10.12770/c7b53704-999d-4721-b1a3-04ec60c87238
        #
        # taken from
        # http://www.emodnet-bathymetry.eu/data-products/acknowledgement-in-publications
        #
        # * The GEBCO Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003
        #
        # taken from
        # https://www.bodc.ac.uk/projects/data_management/international/gebco/gebco_digital_atlas/copyright_and_attribution/

        "bathymetry_source" => "The GEBCO 30sec Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003",

        # NetCDF CF standard name
        # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
        # example "standard_name" = "sea_water_temperature",
        "netcdf_standard_name" => "$(replace(varname,' '=>'_'))",
        "netcdf_long_name" => "$varname",
        "netcdf_units" => "$unit",

        # Abstract for the product
        #"abstract" => "...",

        # This option provides a place to acknowledge various types of support for the
        # project that produced the data
        "acknowledgement" => "Aggregated data products are generated by SMHI with the support of SWAM",
        "documentation" => "https://doi.org/10.6092/A8CFB472-10DB-4225-9737-5A60DA9AF523",

        # Digital Object Identifier of the data product
        "doi" => "$doi",

        "DIVAnd_source" => "https://github.com/gher-ulg/DIVAnd.jl",
        "DIVAnd_version" => "2.7.5",
        "DIVA_code_doi" => "10.5281/zenodo.4715361",
        "DIVA_references" => "Barth, A.; Beckers, J.-M.; Troupin, C.; Alvera-Azcárate, A. & Vandenbulcke, L.
    divand-1.0: n-dimensional variational data analysis for ocean observations
    Geoscientific Model Development, 2014, 7, 225-241. DOI: :10.5194/gmd-7-225-2014"); 

    metadata[monthlist_index]=metadata_season

    @info("starting DIVAnd computations for $(seasons[monthlist_index])")
    @info(Dates.now())

    # Time selection for the analyse. This was already defined together with yearlist, month_list, seasons
    TS = DIVAnd.TimeSelectorYearListMonthList(year_list,month_list[monthlist_index:monthlist_index])
    #@show TS;

    # File name based on the variable (but all spaces are replaced by _)
    nc_filename = "$(replace(varname,' '=>'_'))_$(minimum(year_list))-$(maximum(year_list))_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_$(bath_file_name)_varcorrlenz.nc"
    nc_filepath = joinpath("$(results_dir)nc/O2", nc_filename)

    #Append the created files to file_list
    push!(file_list, nc_filename)

    if isfile(nc_filepath)
       rm(nc_filepath) # delete the previous analysis
    end
    @info("Will write results in $nc_filepath")
    # create attributes for the netcdf file (need an internet connexion)
    ncglobalattrib,ncvarattrib = SDNMetadata(metadata_season,nc_filepath,varname,lonr,latr)

    @time dbinfo = diva3d((lonr,latr,depthr,TS),
              (obslon,obslat,obsdepth,obstime_shifted),
              obsval,
              len,
              epsilon_weighted,
              nc_filepath, varname,
              bathname = bathname,
              bathisglobal = bathisglobal,
              ncvarattrib = ncvarattrib,
              ncglobalattrib = ncglobalattrib,
              timeorigin = timeorigin,
              #transform = Anam.loglin(epsilon),
              #transform = Anam.loglin(-5.),
              mask = new_mask,
              solver = :direct,
              niter_e = 1,
              background = DIVAnd.backgroundfile(filenamebackground,varname,TSbackground),
              error_thresholds = error_thresholds,
              surfextend = true,
              alphabc = 0,
              MEMTOFIT = 250
       );

    # Save the observation metadata in the NetCDF file
    #DIVAnd.saveobs(nc_filepath,(obslon,obslat,obsdepth,obstime),obsid,used = dbinfo[:used])
    DIVAnd.saveobs(nc_filepath,"Oxygen_data", obsval, (obslon,obslat,obsdepth,obstime),obsid, used = dbinfo[:used])

end

# Serialize the list of filenames to JSON and print it
json_string = JSON.json(file_list)
# Write the JSON string to a file
open(joinpath(results_dir, "file_list.json"), "w") do file
    write(file, json_string)
end