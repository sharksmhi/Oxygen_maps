"""
julia julia_code\oxygen_analysis.jl, år, säsong, DIVAsettings
python python_code\calculata_areas.py, år, säsong, DIVASettings
"""
import subprocess
import json
from python_code import calculate_areas as calculate_areas
from python_code import plot_result
from pathlib import Path
from sys import stdout

def run_julia_function(args):
    try:
        # Run Julia script using subprocess and pass function name and arguments
        # process = subprocess.Popen(['julia', 'julia_code/my_julia.jl', function_name, arg1, arg2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # ropa på vår oxygen analysis funktion i oxygen analysis juila scriptet. 
        # Lägg till fler argument efter outputdir med komma mellan
        subprocess.run(args, check=True)

    except FileNotFoundError:
        print("Julia executable not found. Make sure Julia is installed and added to the system PATH.")
        
if __name__ == "__main__":

    # Specify function name and arguments
    # function_name = "my_julia_function"
    # argument1 = "value1"

    # Data input directory
    input_dir = "data"
    # Result directory
    #results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    # results_dir = "C:/LenaV/code/DIVAND/resultat/"
    #results_dir = "C:/Work/DIVAnd/Oxygen_maps/resultat/"
    # Input data filename
    #data_fname = "EMODNET_SHARK_ICES_240620.txt"
    data_fname = "EMODNET_SHARK_ICES_SYKE_240913.txt"

   #Definiera basins
    basin_dict = [
        {
            "Basin": "Baltic Proper",
            "lonr": '14.:dx:31',
            #"lonr": '9.:dx:31',
            "latr": '53.5:dy:61',
            "lenf": json.dumps(80000),
            "epsilon": json.dumps(0.2),
            "dx": json.dumps(0.05),
            "threshold_list" : json.dumps([0, 90, 180]),
            "depthr": json.dumps([0., 10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110.,
              115., 120., 125., 130., 135., 140., 145., 150., 175., 200., 250., 300.]),
            "lenz_" : json.dumps([20., 20., 20., 10., 10., 10., 20., 20., 10., 10., 10., 10., 10., 10., 10., 10., 10., 50., 50., 50., 50.,
             50., 50., 50., 50., 50., 50., 50., 50., 100., 100., 100.])
        },
        {
            "Basin": "Bothnian Bay",
            #lon = json.dumps([16.5, 25.5]),
            "lonr": '16.5:dx:27',
            "latr": '60:dy:66.',
            "lenf": json.dumps(80000),
            "epsilon": json.dumps(0.2),
            "dx": json.dumps(0.05),
            #                            3ml/l, 4ml/l 5ml/l
            "threshold_list": json.dumps([135, 180, 225]),
            "depthr": json.dumps(
                [0., 10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110.,
                 115., 120., 125., 130., 135., 140., 145., 150., 175., 200.]),
            "lenz_": json.dumps(
                [20., 20., 20., 10., 10., 10., 20., 20., 10., 10., 10., 10., 10., 10., 10., 10., 10., 50., 50., 50.,
                 50.,
                 50., 50., 50., 50., 50., 50., 50., 50., 100.])
        },
        {
            "Basin": "Kattegat",
            #"lonr": '9.:dx:13.1',
            "lonr": '9.:dx:15',
            "latr": '53.8:dy:57.75',
            "lenf": json.dumps(40000.0),
            "epsilon": json.dumps(0.2),
            "dx": json.dumps(0.05),
            "threshold_list": json.dumps([0, 90, 180]),
            "depthr": json.dumps(
                [0., 5., 10., 12.5, 15., 17.5, 20., 22.5,  25., 27.5, 30., 32.5, 35., 40., 50., 55., 60., 65., 70., 75., 80.]),
            "lenz_": json.dumps(
                [10., 10., 10., 5., 5.,   5.,  5.,   5.,    5.,  5.,   5.,  5.,  10., 20., 10., 10., 10., 10., 10., 10., 10.])
        }
    ]
    # Initialize variables
    basin = []
    lonr = []
    latr = []
    lenf = []
    epsilon = []
    dx = []
    threshold_list = []
    depthr = []
    lenz_ = []

    # Här bestämmer vi vilka bassäng som skall analyseras
    # Hämta info från vald bassäng
    for item in basin_dict:
        if item["Basin"] == "Kattegat":
        #if item["Basin"] == "Bothnian Bay":
        #if item["Basin"] == "Baltic Proper":
            basin = item["Basin"]
            lonr = item["lonr"]
            latr = item["latr"]
            lenf = item.get("lenf")
            epsilon = item.get("epsilon")
            dx = item.get("dx")
            threshold_list = item.get("threshold_list", [])
            depthr = item.get("depthr", [])
            lenz_ = item.get("lenz_", [])

    # Visa extracted information
    print("Show the setup for the choosen basin: ", basin)
    print("lonr:", lonr)
    print("latr:", latr)
    print("lenf:", lenf)
    print("epsilon:", epsilon)
    print("dx:", dx)
    print("threshold list:", threshold_list)
    print("depthr list:", depthr)
    print("lenz_ list:", lenz_)

    #Samlar ihop resultaten:
    results_dir = Path(f"C:/Work/DIVAnd/Oxygen_maps/resultat/{basin.replace(' ', '_')}/")
    results_dir.mkdir(parents=True, exist_ok=True)
    Path(results_dir, "figures/").mkdir(parents=True, exist_ok=True)
    Path(results_dir, "nc/").mkdir(parents=True, exist_ok=True)
    Path(results_dir, "nc/O2/").mkdir(parents=True, exist_ok=True)
    Path(results_dir, "nc/processed/").mkdir(parents=True, exist_ok=True)
    #print(results_dir)
    #results_dir = str(results_dir) + '\\'
    #exit()
    # Years, month and seasons to be analysed
    #year_list = json.dumps([1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978,
    #1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997,
    #1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
    #2017, 2018, 2019, 2020, 2021, 2022])
    year_list = json.dumps([2000, 2001])
    #year_list = json.dumps([2000])

    seasons_dict = {
                "Winter": [11,12,1,2],
                "Spring": [3,4,5],
                "Summer": [6,7,8],
                "Autumn": [8,9,10]
               }
    seasons = []
    month_list = []
    for season, months in seasons_dict.items():
        seasons.append(season)
        month_list.append(months)
    seasons = json.dumps(seasons)
    month_list = json.dumps(month_list)

    # Correlation length
    # Vi bör köra med lite längre lenf troligen 80_000km då vi har ca 40nm mellan våra station i eg.Östersjön
    #lenf = json.dumps(80000)    #Km
    # Resolution
    #dx = json.dumps(0.05)       #deg
    # Signal to noise ratio
    # low epsilon means higher noise in data and result is more smoothed
    # high epsilon means lower noise in data and result is less smoothed and each observation is seen more
    #epsilon = json.dumps(0.2)
    #Bathymetry file
    bath_file_name = "bat_elevation_Baltic_Sea_masked"

    print("Bathymetry file: ", bath_file_name)
    #Thresholds to analyse in µmol/l oxygen (0, 2, 4 ml/l)
    #threshold_list = [0, 90, 180]   #µmol/l
    # Modify data weight
    w_depth = json.dumps(5.)
    w_days = json.dumps(2.)
    #Set True if you want to save area_data to file (for time-series bar plots)
    save_area_data=True

    print("running DIVAnd in Julia...")
    args = ['julia', 'julia_code/oxygen_analysis.jl', input_dir, results_dir, data_fname, year_list, month_list, seasons, lenf, epsilon, dx, bath_file_name, w_depth, w_days, depthr, lenz_, lonr, latr, basin, threshold_list]
    # Call the function and save a json-file with a file_list containing the results. That we can send to the calculate_areas function.
    run_julia_function(args)

    file_list = []

    for season in json.loads(seasons):

        file_list.append(f"Oxygen_{min(json.loads(year_list))}-{max(json.loads(year_list))}_{season}_{json.loads(epsilon)}_{lenf}_{json.loads(dx)}_{w_depth}_{w_days}_{bath_file_name}_varcorrlenz.nc")
    
    # file_list = [
    #     "Oxygen_background_weighted_0.03_field_gebco_30sec_4.nc",
    #     "Oxygen_background_weighted_0.03_field_bat_elevation_Baltic_Sea_masked.nc",
    #     "Oxygen_background_weighted_0.05_field_bat_elevation_Baltic_Sea_masked.nc",
    #     "Oxygen_background_weighted_0.05_field.nc"
    # ]
    #Calculate areas from DIVA-results and save in a new nc-file. Results in file_list
    print("calculating areas...")

    calculate_areas.calculate_areas(results_dir, file_list, json.loads(threshold_list), save_area_data)

    #Read and plot areas in file_list
    print("plotting...")
    plot_result.read_processed_nc(results_dir,file_list,year_list)

### extract values that are within our limits, save to a new variable and nc-file. ####
# 1 ml/l of O2 is approximately 43.570 µmol/kg
# (assumes a molar volume of O2 of 22.392 l/mole and
# a constant seawater potential density of 1025 kg/m3).
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html

