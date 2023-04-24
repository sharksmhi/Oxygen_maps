import pandas as pd

# Define file path and data types for each column
file_path = "havgem/DIVA/syrekartor/data/"
file_name = "bot.csv"
dtype_dict = {
    "%Platform": str,
    "St_No": str,
    "Latdeg": str,
    "Longdeg": str,
    "Year": str,
    "Mnth": str,
    "Dy": str,
    "Hr": str,
    "Sound": int,
    "Pressure(Dbars)": float,
    "TEMP": float,
    "PSAL": float,
    "DOXY(umol/l)": float,
    "H2SX(umol/l)": float
}

# Load data from CSV file into a pandas DataFrame
df = pd.read_csv(file_path+file_name, dtype=dtype_dict)

# Rename the '%Platform' column to 'ID'
df.rename(columns={"%Platform": "ID"}, inplace=True)

# Replace invalid hours with '0000'
df.loc[df["Hr"].isin(['-9', '2400']), "Hr"] = '0000'

# Combine the date and hour columns into a new column 'date_time'
df['date_time'] = pd.to_datetime(df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2)).dt.strftime("%Y-%m-%dT%H:%M")
# df['date'] = df["Year"] + "-" + df["Mnth"] + "-" + df["Dy"] + " " + df["Hr"].str.slice(0,2) #+ ":00"
# # format date column to YYYY-MM-DDTHH:MM
# df['date_time'] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%dT%H:%M")

# Define the order of columns for the output files
column_list = ["Longdeg", "Latdeg", "DOXY(umol/l)", "Pressure(Dbars)", "TEMP", "PSAL", "H2SX(umol/l)", "St_No", "Year", "date_time", "ID", "Mnth", "Dy", "Hr"]

# Filter out rows with missing data
df_filtered = df.loc[df["DOXY(umol/l)"] != -9]

# Write the filtered data to two output files, one with headers and one without
# df_filtered[column_list].to_csv(file_path+"bot.txt", index=False, sep='\t')
# df_filtered[column_list].to_csv(file_path+"bot_no_header.txt", index=False, header=False, sep='\t')



