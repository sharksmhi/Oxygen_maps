# ### Oxygen maps
# ### 2023
# ### Lena Viktorsson & Martin Hansson
# ### Based on copy from Karin Wesslander and DIVA workshop notebooks

# read dataset
# set up metadata info
# run DIVAnd and save to netcdf.

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

args= ARGS
input_dir = args[1]
result_dir = args[2]

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
#fname = "SHARK_EMODNET.txt"
data_fname = "EMODNET_SHARK_ICES"
@time obsval,obslon,obslat,obsdepth,obstime,obsid = loadbigfile(joinpath(input_dir, "$(data_fname).txt"))

#Load background field
filenamebackground = joinpath(input_dir, "$(replace(varname,' '=>'_'))_background_weighted_0.05_field.nc")
# year and month-list for background analysis
yearlist = [1960:1969,1970:1979,1980:1989,1990:1999,2000:2009,2010:2021];
month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];  # Seasonal climatology
TSbackground = DIVAnd.TimeSelectorYearListMonthList(yearlist,month_list);

# Sätt horisontell uppplösning
#dx, dy = 0.125, 0.125  #Karin dx, dy = 0.1, 0.1
dx, dy = 0.05, 0.05  #Karin dx, dy = 0.1, 0.1
lonr = 9.:dx:31.
latr = 53.5:dy:61.

#New lists for analysis
yearlist = [1993];
#yearlist = [1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021];
#yearlist = [1960,1965,1970,1975,1980,1994,1995,1998,2000,2003,2004,2005,2010,2014,2015,2018,2020,2021];
month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
seasons=["Winter","Spring","Summer","Autumn"]
months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];

#month_list = [[6,7,8], [8,9,10]];
#seasons=["Summer","Autumn"]
#months=["(June-Aug)","(Aug-Oct)"];

#month_list = [[3,4,5]];
#seasons=["Spring"]
#months=["(Mar-May)"];

# Time origin for the NetCDF file
timeorigin = DateTime(1900,1,1,0,0,0);
aspect_ratio = 1/cos(mean(latr) * pi/180);

# #### Set negative oxygen to 0
# Sätter alla syrevärden som är noll eller mindre än noll till ngt väldigt nära noll, dock ej noll.
# Om syre är noll eller negativt får vi inte ut några värden till analysen. Lite knepigt att den inte klarar noll dock.
# Kanske fråga Charles när vi har tid. 
# %% <0.2 ml/l är från SMHIs ctd för syrefritt vilket motsvarar 8.93 µmol/l.

#if varname == "Oxygen"
#    sel_q = (obsval .<= -250.);
#    obsval[sel_q] .= -250.;
#    sel_q = (obsval .> 600.);
#    obsval[sel_q] .= 600.;
#end

# ## Extract the bathymetry
# It is used to delimit the domain where the interpolation is performed.
# Modify bathname according to the resolution required.

bath_file_name = "gebco_30sec_4"
bathname = joinpath(input_dir, "$(bath_file_name).nc")
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,lonr,latr);

# ## MASK
# #### 1) Create a mask from bathymetry and your selected depth vector

# First set the depth resolotion, this will be the depths used in DIVArun
# Then set the horizontal correlation length (should be twice the resolution)
depthr = [0.,  10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70.,  75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];
# lenz_ =  [10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,  5.,  10., 10., 10., 10., 10.,  10.,  10.,  10.,  10.,   5.,  10.,  10.]; #används ej

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

aspectratio = 1/cos(mean(54:0.05:61) * pi/180)
fig=figure(1)
ax = subplot(1,1,1)
ax.tick_params("both",labelsize=6)
#pyplot.plot(xmask, ymask, new_mask[:,:,1]');
pcolor(xmask, ymask, new_mask[:,:,1]');
gca().set_aspect(aspectratio)
gcf()

# #### Plot

# ## Analysis parameters
# ## Horizontal correlation length
# From https://github.com/gher-uliege/DIVAnd.jl/issues/121 
#
# Here is a piece of code (from Emodnet Chemistry), which could be used to create a correlation length lenfilled
# decreasing to a fifth of the normal value lenf when approaching a cost.
# The distance at which the change happens is defined by slen (here in meters, so metrics pmn are in /m). sz is size(mask)
#
# The mask is 2D.
#
# Ökas slen så minskas avståndet från kusten där len är mindre.

# lenf should be the normal value for horizontal correlation length
# Sätt korrelationslängd (m)
# Vi bör köra med lite längre lenf troligen 80_000km då vi har ca 40nm mellan våra station i eg.Östersjön
lenf = 80_000.#49_000. #78_000, 50_000 25_000

sz = (length(lonr),length(latr),length(depthr));
lenx = fill(lenf,sz)   # 100 km
leny = fill(lenf,sz)   # 100 km
lenz = fill(10,sz);      # 25 m
len = (lenx, leny, lenz);

lx = lenf #78_000. #50_000.
ly = lenf #78_000.

# ### Plotting function
# Define a plotting function that will be applied for each time index and depth level.     
# All the figures will be saved in a selected directory.
# This function is within the time loop. Try to move outside and understand the how the input argumetns are sent to plotres from DIVArun

# If you do not want to generate plots but print the time index at every time slice
# you can use the function `plotres_timeindex`.
# function plotres_timeindex(timeindex,sel,fit,erri)
#    @show timeindex
#end

# # Run DIVA3d
# ### Background field is calculated in notebook: 
# filename_background = "Background_file_6y_20230203_$(savevar).nc"
# TS_back = DIVAnd.TimeSelectorYearListMonthList(yearlist_back, month_list)
# background = DIVAnd.backgroundfile(filename_background,varname,TS_back)

# To include December from previous year in the analyse
obstime_shifted = copy(obstime)
# Get all dates with dec and nov and add one year to these specific dates.
# Winter period will be jan-feb + nov-dec from previous year.
obstime_shifted[Dates.month.(obstime) .== 12 .| Dates.month.(obstime) .== 11] .+= Dates.Year(1)

# Settings for DIVAnd-------------------------------------------------------------------------------
error_thresholds = [("L1", 0.3), ("L2", 0.5)];
epsilon = 0.2; #1., 0.1, 10

# low epsilon means higher noise in data and result is more smoothed
# high epsilon means lower noise in data and result is less smoothed and each observation is seen more

# Modify data weight
# ⌛⌛⌛
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

w_depth = 5.
w_days = 2.

rdiag_jldfile = joinpath(input_dir, "$(data_fname)_weighted_$(w_depth)_$(w_days).jld")

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
#epsilon_weighted = epsilon
#---------------------------------------------------------------------------------------------------
# sätter min, max för färgskalan i plottarna
vmin_ = 0
vmax_ = 400
# %%
# One metadata set up per season
metadata=Array{DataStructures.OrderedDict{String,Any}}(undef,4) ;



for monthlist_index in 1:length(month_list)
    season = seasons[monthlist_index]
    @show monthlist_index
    @show seasons
    @info("Creating metadata dicitonary for the season $(season)")
    metadata_season = OrderedDict(
        # set attributes for DVIA run from our settings
        "season" => season,
        "epsilon" => "$epsilon",
        "horizontal correlation length m" => "$lx",
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
        "acknowledgement" => "Aggregated data products are generated by EMODnet Chemistry under the support of DG MARE Call for Tenders EASME/EMFF/2016/006-lot4, EASME/2019/OP/0003-lot4.",
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

    function plotres(timeindex,selected,fit,erri)
        year=yearlist[timeindex]
        @show timeindex
        tmp = copy(fit)
        L2 = 0.5
        L1 = 0.3

        #Välj djup som skall plottas
        selected_depth = 70.

        depth_index = findall(depthr .== selected_depth)
        @info("making figures for $(year) in $(season) including months $(month_list[monthlist_index]) $(selected_depth) m")
        #tmp[erri .> rel_error] .= NaN;
        figure(figsize = (10,8))
        title("test")

        #Plot1
        subplot(2,2,1)
        title("$(year),$(season),$(selected_depth)m")

        # select the data in a depthrange of +/- 5 m from selected depth
        sel_obs = selected .& (obsdepth .> selected_depth-5) .& (obsdepth .< selected_depth+5)

        # om man vill sätta max och min för färgskalen efter valda data.
        # nu sätter vi vmin, vmax utanför loopen så det blir samma färgskala på alla.
        #vmin = vmin_ #minimum(obsval[selsurface])
        #vmax = vmax #maximum(obsval[selsurface])
        # plot the data
        xlim(minimum(lonr),maximum(lonr))
        ylim(minimum(latr),maximum(latr))
        # Land
        contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
        # Djupnivå
        contour(bx,by,permutedims(b,[2,1]), levels = [selected_depth], linewidths=[1], colors = [[.15,.25,.35]])
        gca().set_aspect(aspect_ratio)
        # Observationer
        scatter(obslon[sel_obs],obslat[sel_obs],6,obsval[sel_obs];
                vmin = vmin_, vmax = vmax_, cmap="jet_r")  #"PuBuGn","viridis"
        colorbar(extend="max")

        # plot2 the analysis
        subplot(2,2,2)
        @show lenz[depth_index[1]]
        title("lx=$(lx)[m],ly=$(ly)[m], lz=$(lenz[1,1,depth_index[1]])[m],e=$(epsilon)")

        #title("L=$(L)[m], lenx=$(lx)[m], leny=$(ly)[m]")
        @show depth_index[1]
        pcolor(lonr,latr,permutedims(tmp[:,:,depth_index[1]],[2,1]);
            vmin = vmin_, vmax = vmax_, cmap="jet_r")
        colorbar(extend="max")
        #Land
        contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
        gca().set_aspect(aspect_ratio)

        #Plot 3
        subplot(2,2,3)
        title("Error threshold L2>$(L2)")
        #Resultat
        tmp_L2=copy(tmp)
        tmp_L2[erri .< L2] .= NaN
        pcolor(lonr,latr,permutedims(tmp_L2[:,:,depth_index[1]],[2,1]);
            vmin = vmin_, vmax = vmax_, cmap="jet_r")
        colorbar(extend="max")
        #Land
        contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
        gca().set_aspect(aspect_ratio)

        # plot the analysis
        #Plot 4
        subplot(2,2,4)
        title("Error threshold L1>$(L1)")
        #Resultat
        tmp_L1=copy(tmp)
        tmp_L1[erri .< L1] .= NaN
        pcolor(lonr,latr,permutedims(tmp_L1[:,:,depth_index[1]],[2,1]);
            vmin = vmin_, vmax = vmax_, cmap="jet_r")
        colorbar(extend="max")
        #Land
        contourf(bx,by,permutedims(b,[2,1]), levels = [-1e5,0],colors = [[.5,.5,.5]])
        gca().set_aspect(aspect_ratio)

        # save the figure
        figname = savevar * @sprintf("_%04d_", year) * season * @sprintf("_%d",selected_depth) * 
        @sprintf("_%.2f",epsilon) * @sprintf("_%05d",lx) * @sprintf("_%d",lenz[1,1,depth_index[1]]) *
        @sprintf("_%03d.png",timeindex)
        PyPlot.savefig(joinpath(figdir, figname), dpi=300, bbox_inches="tight");
        PyPlot.close_figs()
    end

    # Time selection for the analyse. This was already defined together with yearlist, month_list, seasons
    TS = DIVAnd.TimeSelectorYearListMonthList(yearlist,month_list[monthlist_index:monthlist_index])
    @show TS;

    # File name based on the variable (but all spaces are replaced by _)
    nc_filename = "$(replace(varname,' '=>'_'))_$(minimum(yearlist))-$(maximum(yearlist))_$(season)_$(epsilon)_$(lx)_$(dx)_$(w_depth)_$(w_days)_with_backgroundfield_moredepths_$(bath_file_name).nc"
    nc_filepath = joinpath(result_dir, nc_filename)

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
              plotres = plotres,
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
