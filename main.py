"""
julia julia_code\oxygen_analysis.jl, år, säsong, DIVAsettings
python python_code\calculata_areas.py, år, säsong, DIVASettings
"""
import subprocess
import json
from python_code import calculate_areas as calculate_areas
from python_code import plot_test
from sys import stdout

"""# year: list, seasons: list, DIVAsettings: dict
# location = "//winfs-proj/proj/havgem/DIVA/syrekartor/"
# outputdir = joinpath(location, "resultat/nc/$(savevar)/");
# data_fname = "EMODNET_SHARK_ICES.txt"
# year_list_background = [1960:1969,1970:1979,1980:1989,1990:1999,2000:2009,2010:2021];
# month_list_background = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
# filenamebackground
# year_list = []
# month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
# seasons=["Winter","Spring","Summer","Autumn"]
# months=["(Nov-Feb)","(Mar-May)","(June-Aug)","(Aug-Oct)"];
# depthr = [0.,  10., 20., 25., 30., 35., 40., 50., 55., 60., 65., 70.,  75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145.,150.,175.,200.,250.,300.];
"""

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

    #Data input directory
    input_dir = "data"
    #Result directory
    results_dir = "//winfs-proj/proj/havgem/DIVA/syrekartor/resultat/"
    #Input data filename
    data_fname = "EMODNET_SHARK_ICES.txt"
    #Years, month and seasons to be analysed
    year_list = json.dumps([1996])
    #month_list = [ [11,12,1,2], [3,4,5], [6,7,8], [8,9,10]];
    month_list = json.dumps([[11, 12, 1, 2]])
    seasons = json.dumps(["Winter"])
    #seasons=["Winter","Spring","Summer","Autumn"]
    #Correlation length
    lenf = 80_000   #Km
    #Thresholds to analyse
    threshold_list = [0, 90, 180]

    args = ['julia', 'julia_code/oxygen_analysis.jl', input_dir, results_dir, data_fname, year_list, month_list, seasons]
    # Call the function and save a json-file with a file_list containing the results. That we can send to the calculate_areas function.
    run_julia_function(args)

    # Open the JSON file
    with open(f"{results_dir}file_list.json", 'r') as file:
        # Load JSON data from the file
        file_list = json.load(file)

    #Calculate areas from DIVA-results and save in a new nc-file. Results in file_list
    calculate_areas.calculate_areas(results_dir, file_list, threshold_list)

    #Read and plot areas in file_list
    plot_test.read_processed_nc(results_dir,file_list,year_list)

### extract values that are within our limits, save to a new variable and nc-file. ####
# 1 ml/l of O2 is approximately 43.570 µmol/kg
# (assumes a molar volume of O2 of 22.392 l/mole and
# a constant seawater potential density of 1025 kg/m3).
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html

