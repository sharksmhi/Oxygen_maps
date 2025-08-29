
import polars as pl
from pathlib import Path
import pandas as pd
import numpy as np

def _format_date(df):
    df['date'] = pd.to_datetime(df['date'], format="mixed", utc=True)
    # standardize as string with desired format, e.g. "YYYY-MM-DD HH:MM"
    # df['date'] = df['date'].dt.strftime('%Y-%m-%dT%H:%M')
    # df['date'] = df['date'].astype(str).str.strip()
    # df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_shark_data(file_path):
    """Loads shark tab-delimited dataset."""
    cols = ['depth', 'Dissolved oxygen O2 bottle (ml/l)', 'Q-flag Dissolved oxygen O2 bottle',	'Dissolved oxygen O2 CTD (ml/l)', 'Q-flag Dissolved oxygen O2 CTD',	'Hydrogen sulphide H2S (umol/l)', 'Q-flag Hydrogen sulphide H2S', 'DOXY_ml', 'value', 'NEG_DOXY_umol']

    df = pd.read_csv(file_path, sep="\t", names=cols)
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    df = _format_date(df)
    return df

def load_ices_btl_data(file_path):
    """Loads ICES BTL dataset as downloaded from ices.dk."""
    """
    ICES: remove all Swedish and Russian data. 
            Swedish data because we have more reliable data in shark.
            Russian because it is generally unreliable quality.
    """

    # Read file
    df = pd.read_csv(file_path, sep=",", skiprows=21)
    # Select columns to keep
    df = df[
        ['Longitude [degrees_east]',
            'Latitude [degrees_north]',
            'Oxygen (DOXYZZXX_UMLL) [ml/l]',
            'QV:ODV:Oxygen (DOXYZZXX_UMLL)',
            'Depth (ADEPZZ01_ULAA) [m]',
            'yyyy-mm-ddThh:mm:ss.sss',
            'Cruise',
            'Station'
            ]
        ]
    str_cols = ["Cruise", "Station", "yyyy-mm-ddThh:mm:ss.sss", "QV:ODV:Oxygen (DOXYZZXX_UMLL)"]
    df[str_cols] = df[str_cols].astype(str)
    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Oxygen (DOXYZZXX_UMLL) [ml/l]": "value",
            "Depth (ADEPZZ01_ULAA) [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    # Create id colulmn
    df['id'] = "ices_btl-" + df['Cruise'] + "-" + df['Station'] + '-' + df['date']

    # Remove any rows with NaN Oxygen
    df = df.dropna(subset=["value"])
    # Change unit from ml/l to µmol/l
    df['value'] = df['value'].apply(lambda x: x * 44.661)

    # Remove any data from russia (90), which is considered always bad and Sweden (77) which at ICES contain bad and unflagged data
    # Instead we add our own Sweden data from SHARK, much better and accurate"""
    remove_countries = ['90', '77']
    # Remove bad or questionable data; QV= 4 or 8 in ICES-data. 
    df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"] = df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"].str.strip()
    removed_ICES_data = df[df["Cruise"].str.startswith(tuple(remove_countries)) | df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"].isin(['4', '8'])]

    # Only keep rows that should not be removed
    df = df[~df.index.isin(removed_ICES_data.index)]

    # Save the removed rows to a new CSV-fil
    removed_ICES_data.to_csv("removed_ices_btl_data.txt", sep="\t", index=False)

    df = _format_date(df)

    return df

def load_ices_ctd_data(file_path):
    """Loads ICES CTD dataset with specific format."""
    df = pd.read_csv(file_path, sep=",", skiprows=9)
    df = df[
        ['Longitude [degrees_east]',
        'Latitude [degrees_north]',
        'Oxygen (DOXYZZXX_UMLL) [ml/l]',
        'QV:ODV:Oxygen (DOXYZZXX_UMLL)',
        'Depth (ADEPZZ01_ULAA) [m]',
        'yyyy-mm-ddThh:mm:ss.sss',
        'Cruise',
        'Station']
        ]
    str_cols = ["Cruise", "Station", "yyyy-mm-ddThh:mm:ss.sss", "QV:ODV:Oxygen (DOXYZZXX_UMLL)"]
    df[str_cols] = df[str_cols].astype(str)
    df['id'] = "ices_ctd-" + df['Cruise'].astype(str) + "-" + df['Station'].astype(str) + '-' + df['yyyy-mm-ddThh:mm:ss.sss'].astype(str)
    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Oxygen (DOXYZZXX_UMLL) [ml/l]": "value",
            "Depth (ADEPZZ01_ULAA) [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    #Remove any rows with NaN Oxygen
    df = df.dropna(subset=["value"])

    """Remove bad or questionable data; QV= 4 or 8 in ICES-data. 
    Remove any data from russia (90), which is considered always bad and Sweden (77) which at ICES contain bad and unflagged data
    Instead we add our own Sweden data from SHARK, much better and accurate"""
    remove = ['90', '77']

    # Filter our row that should be removed
    df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"] = df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"].astype(str).str.strip()
    removed_ICES_data = df[df["Cruise"].str.startswith(tuple(remove)) | df["QV:ODV:Oxygen (DOXYZZXX_UMLL)"].isin(['4', '8'])]
    
    # Only keep rows that should not be removed
    df = df[~df.index.isin(removed_ICES_data.index)]

    # Save the removed rows to a new CSV-fil
    removed_ICES_data.to_csv("removed_ices_ctd_data.txt", sep="\t", index=False)

    # Change unit from ml/l to µmol/l
    df['value'] = df['value'].apply(lambda x: x * 44.661)

    df = _format_date(df)

    return df


def ICES_to_big_file(file_path):
    # Convert a ICES-file to a "big file" readable for DIVAnd
    # ICES data from ices.dk 
    # BTL and CTD low res dataset

    file_name = file_path.name

    dtype_dict = {
        "Cruise": str,
        "Station": str,
        "Latitude [degrees_north]": str,
        "Longitude [degrees_east]": str,
        "mon/day/yr":str,
        "hh:mm":str,
        "Depth (ADEPZZ01_ULAA) [m]": float,
        "Temperature (TEMPPR01_UPAA) [degC]": float,
        "Salinity (PSALPR01_UUUU) [dmnless]": float,
        "Oxygen (DOXYZZXX_UMLL) [ml/l]": float,
        "Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]": float
    }

    # Load data from CSV file into a pandas DataFrame
    #df = pd.read_csv(file_path+file_name, dtype=dtype_dict)
    #df = pd.read_table(file_path+file_name, delimiter="\t",low_memory=False)
    df = pd.read_table(file_path+file_name, delimiter="\t",dtype=dtype_dict)

    # Rename the '%Platform' column to 'ID'
    #df.rename(columns={"Sampling platform (code)": "ID"}, inplace=True)

    # Konvertera datumkolumnen till datetime
    df['datum'] = pd.to_datetime(df['mon/day/yr'], format='%m/%d/%Y')
    # Extrahera månaden, dag, år
    df['Month'] = df['datum'].dt.month.astype(str)
    df['Day'] = df['datum'].dt.day.astype(str)
    df['Year'] = df['datum'].dt.year.astype(str)

    df["Visit event identifier"] = df["Year"] + '-' + df["Cruise"] + '-' + df["Station"]

    # Replace invalid hours with '00:00'
    df.loc[df["hh:mm"].isin(['-9', '24', '',np.nan]), "hh:mm"] = '00:00'
    #df.loc[df["hh:mm"].isin(['-9', '60', '',np.nan]), "hh:mm"] = '00:00'

    #df["Time"] = df["hh:mm"][0:1].str.zfill(2) + ":" + df["hh:mm"][3:4].str.zfill(2)

    # Combine the date and hour columns into a new column 'date_time'
    df['date_time'] = df["Year"] + "-" + df["Month"].str.zfill(2) + "-" + df["Day"].str.zfill(2) + "T" + df["hh:mm"]

    # Make a new combined O2 variable.
    # Base new variable on O2 bottle and add O2 CTD when O2 bottle is missing.
    def get_doxy(row):
        if row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] and not np.isnan(row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"]):
            return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * -0.04488
            #return row["Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]"] * 0 + 0.01
        if row['Oxygen (DOXYZZXX_UMLL) [ml/l]'] and not np.isnan(row['Oxygen (DOXYZZXX_UMLL) [ml/l]']):
            return row['Oxygen (DOXYZZXX_UMLL) [ml/l]']
        return np.nan

    df["DOXY_ml"] = df.apply(get_doxy, axis=1)

    # Convert O2 ml/l to µmol/l-> 1ml/l = 44.661 µmol/l
    df['DOXY_umol'] = df['DOXY_ml'].apply(lambda x: x*44.661)
    print(df.dtypes)
    # Filter out rows with missing data
    #df_filtered = df.loc[df["DOXY_umol"] != np.nan]
    df_filtered = df.dropna(subset=["DOXY_umol"])

    # Ta bort alla data som är ryska (90) och svenska data 77
    remove_countries = ['90', '77']

    # Filtrera ut rader som matchar något av mönstren
    removed_ICES_data = df_filtered[df_filtered["Cruise"].str.startswith(tuple(remove_countries))]

    # Behåll endast rader som inte matchar mönstren
    df_filtered = df_filtered[~df_filtered["Cruise"].str.startswith(tuple(remove_countries))]

    # Spara de borttagna raderna till en ny CSV-fil (valfritt)
    removed_ICES_data.to_csv(file_path+"removed_data.csv", index=False)

    # Define the order of columns for the output files
    column_list = ["Longitude [degrees_east]", "Latitude [degrees_north]", "DOXY_umol", "Depth (ADEPZZ01_ULAA) [m]", "Temperature (TEMPPR01_UPAA) [degC]", "Salinity (PSALPR01_UUUU) [dmnless]", "Hydrogen Sulphide (H2SXZZXX_UPOX) [umol/l]", "Cruise", "Year", "date_time", "Visit event identifier", "Month", "Day", "hh:mm"]

    # Write the filtered data to two output files, one with headers and one without
    df_filtered[column_list].to_csv(file_path+"ICES_btl_lowres_ctd_02_250328.txt", index=False, header=False, sep='\t')

    return df

def load_emodnet_btl_data(file_path):
    """Loads EMODNET dataset with specific format."""
    df = pd.read_csv(file_path, sep="\t", skiprows=26)
    df = df[
        ['Longitude [degrees_east]',
            'Latitude [degrees_north]',
            'Water body dissolved oxygen concentration [umol/l]',
            'Depth [m]',
            'yyyy-mm-ddThh:mm:ss.sss',
            'LOCAL_CDI_ID',
            'EDMO_code'
            ]
        ]
    str_cols = ["EDMO_code", "LOCAL_CDI_ID", "yyyy-mm-ddThh:mm:ss.sss"]
    df[str_cols] = df[str_cols].astype(str)
    df['id'] = "emod_btl-" + df['EDMO_code'] + "-" + df['LOCAL_CDI_ID']
    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [umol/l]": "value",
            "Depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    # Remove any rows with NaN Oxygen
    df = df.dropna(subset=["value"])

    df = _format_date(df)

    return df

def load_emodnet_ctd_data(file_path):
    """Loads EMODNET CTD dataset with specific format."""
    df = pd.read_csv(file_path, sep="\t", skiprows=26)
    df = df[
        ['Longitude [degrees_east]',
        'Latitude [degrees_north]',
        'Water body dissolved oxygen concentration [umol/l]',
        'Depth [m]',
        'yyyy-mm-ddThh:mm:ss.sss',
        'LOCAL_CDI_ID',
        'EDMO_code']
        ]

    # Fyll på metadata (alla kolumner utom de som varierar per rad) för varje batch
    metadata_columns = ['Longitude [degrees_east]', 'Latitude [degrees_north]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID', 'EDMO_code']
    df[metadata_columns] = df[metadata_columns].ffill()
    str_cols = ["EDMO_code", "LOCAL_CDI_ID", "yyyy-mm-ddThh:mm:ss.sss"]
    df[str_cols] = df[str_cols].astype(str)
    df['id'] = "emod_ctd-" + df['EDMO_code'] + "-" + df['LOCAL_CDI_ID']

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [umol/l]": "value",
            "Depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    df = _format_date(df)

    return df


def load_big_file_format(file_path):
    """Loads generic tab-delimited dataset."""
    cols = ['lon', 'lat', 'value', 'depth', 'drop1', 'drop2', 'drop3', 'drop4', 'drop5', 'date', 'id']
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.columns = cols + df.columns[len(cols):].tolist()
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    df = _format_date(df)
    
    return df
