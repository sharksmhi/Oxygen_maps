# ##Program to convert a sharkwebexport to a "big file" readable for DIVAnd
# ##Export with decoding: UTF-8
# ##

import pandas as pd
import numpy as np

# Define file path and data types for each column
file_path = "C:/Work/DIVAnd/Oxygen_maps/data/all_baltic/"
file_name = "sharkweb_data_1960_2024_utf8.txt"

# Load data from CSV file into a pandas DataFrame
#df = pd.read_csv(file_path+file_name, dtype=dtype_dict)
df = pd.read_table(file_path+file_name, delimiter="\t")

# Rename the '%Platform' column to 'ID'
#df.rename(columns={"Sampling platform (code)": "ID"}, inplace=True)

# Replace invalid hours with '0000'
df.loc[df["Sampling time (start)"].isin(['-9', '2400', '',np.nan]), "Sampling time (start)"] = '00:00'

# Get day from "Sampling date"
df['Day'] = df["Sampling date"].apply(lambda x: x.split("-")[2])

# Combine the date and hour columns into a new column 'date_time'
df['date_time'] = df["Sampling date"].astype(str) + "T" + df["Sampling time (start)"]
#df['date_time'] = pd.to_datetime(df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2)).dt.strftime("%Y-%m-%dT%H:%M")
# df['date'] = df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2) #+ ":00"
# # format date column to YYYY-MM-DDTHH:MM
# df['date_time'] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%dT%H:%M")

# Make a new combined O2 variable.
# Base new variable on O2 bottle and add O2 CTD when O2 bottle is missing.
def get_doxy(row):
    if row["Hydrogen sulphide H2S (umol/l)"] and not np.isnan(row["Hydrogen sulphide H2S (umol/l)"]):
        return row["Hydrogen sulphide H2S (umol/l)"] * -0.04488
    if row['Dissolved oxygen O2 bottle (ml/l)'] and not np.isnan(row['Dissolved oxygen O2 bottle (ml/l)']):
        return row['Dissolved oxygen O2 bottle (ml/l)']
    if row['Dissolved oxygen O2 CTD (ml/l)'] and not np.isnan(row['Dissolved oxygen O2 CTD (ml/l)']):
        return row['Dissolved oxygen O2 CTD (ml/l)']
    return np.nan

df["DOXY_ml"] = df.apply(get_doxy, axis=1)

# Convert O2 ml/l to µmol/l-> 1ml/l = 44.661 µmol/l
#df['DOXY_umol'] = df['DOXY_ml'].multiply(44.661)
df['DOXY_umol'] = df['DOXY_ml'].apply(lambda x: x*44.661)
print(df.dtypes)
# Filter out rows with missing data
#df_filtered = df.loc[df["DOXY_umol"] != np.nan]
df_filtered = df.dropna(subset=["DOXY_umol"])

# Define the order of columns for the output files
#column_list = ["Longdeg", "Latdeg", "DOXY(umol/l)", "Pressure(Dbars)", "TEMP", "PSAL", "H2SX(umol/l)", "St_No", "Year", "date_time", "ID", "Mnth", "Dy", "Hr"]
column_list = ["Sample longitude (DD)", "Sample latitude (DD)", "DOXY_umol", "Sampling depth (m)", "Temperature bottle (C)", "Salinity bottle (o/oo psu)", "Hydrogen sulphide H2S (umol/l)", "Visit event identifier", "Year", "date_time","Sampling platform (code)", "Month", "Day", "Sampling time (start)"]

# Write the filtered data to two output files, one with headers and one without
# df_filtered[column_list].to_csv(file_path+"bot.txt", index=False, sep='\t')
df_filtered[column_list].to_csv(file_path+"sharkweb_btlctd_02_240603.txt", index=False, header=False, sep='\t')



