import scipy.io
import numpy as np

data_folder = "havgem/DIVA/syrekartor/data/"

# Load all variables from MATLAB file
# mat_vars = scipy.io.loadmat(data_folder+"o2_1960_2022.mat")

# Load information about variables in MATLAB file
mat_info = scipy.io.whosmat(data_folder+"o2_1960_2022.mat")

# Print information about each variable
# Combine variables into a single array, skipping variables with different length
data_list = []
header = []
first_len = None
for var in mat_info:
    print(f"Variable name: {var[0]}")
    print(f"Variable shape: {var[1]}")
    print(f"Variable data type: {var[2]}")

    # Load variable data from MATLAB file
    data = scipy.io.loadmat(data_folder+"o2_1960_2022.mat", variable_names=var[0])

    # Remove any dimensions of size 1 using np.squeeze()
    x = np.squeeze(data[var[0]])

    # Convert array elements to strings and remove brackets and quotes
    x_str = np.array([str(elem).replace('[','').replace(']','').replace('\'','').replace('\"','') for elem in x])

    # Check if all arrays have the same length
    if first_len is None:
        first_len = len(x_str)
    elif len(x_str) != first_len:
        print(f"Skipping variable {var[0]} due to different length ({len(x_str)} instead of {first_len})")
        continue

    data_list.append(x_str)
    header.append(var[0])

    # Save variable data to text file
    # np.savetxt(f"{data_folder}mat_file_as_txt/{var[0]}.txt", x_str, fmt='%s')
    
if len(data_list) == 0:
    print("No variables with same length found. Stopping.")
else:
    # Stack data arrays
    data = np.column_stack(data_list)

    # Reorder columns as required by DIVAnd bigfile reader
    data_reordered = np.zeros((data.shape[0], 11), dtype=np.object)
    data_reordered[:, 0] = data[:, 2]  # lon
    data_reordered[:, 1] = data[:, 1]  # lat
    data_reordered[:, 2] = data[:, 4]  # value /O2
    data_reordered[:, 3] = data[:, 3]  # depth
    # columns 4-8 can contain any data, but cannot be left empty,
    # fill columns 4-6 with information from mat-file, leave colummns 7-8 with only zeros
    data_reordered[:, 4] = data[:, 5]  # Yr, can contain any data, but cannot be left empty
    data_reordered[:, 5] = data[:, 6]  # Month, can contain any data, but cannot be left empty
    data_reordered[:, 6] = data[:, 0]  # Stno, can contain any data, but cannot be left empty
    data_reordered[:, 9] = data[:, 7]  # date
    data_reordered[:, 10] = data[:, 8] # id /pltf

# Save the reordered data to a file with headers. Headers have to be removed before reading the file in DIVAnd
header_reordered = ['lon', 'lat', 'oxygen', 'depth', 'year', 'month', 'station_no', 'extra1', 'extra2', 'date', 'pltf']
header_str = '\t'.join(var_name for var_name in header_reordered)
np.savetxt(f"{data_folder}mat_file_as_txt/output_file_reordered.txt", data_reordered, header=header_str, comments='', delimiter='\t', fmt='%s')

# order of columns in the outputfiler have to be:
# 0 lon
# 1 lat
# 2 value
# 3 depth
# index 4 to 8 can contain any data, but cannot be left empty
# 9 date
# 10 id

# order in .mat file
# 0 St_No
# 1 lat
# 2 lon
# 3 depth
# 4 value
# 5 Yr
# 6 Month
# 7 date
# 8 id


