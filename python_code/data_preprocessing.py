import pandas as pd
import polars as pl
import numpy as np
from duplicate_check import find_duplicate_candidates, remove_duplicates
import os
from pathlib import Path

# order of columns in the outputfile have to be:
# ['lon', 'lat', 'value', 'depth', 'date', 'id']
# 0 lom
# 1 lat
# 2 value
# 3 depth
# index 4 to 8 can contain any data, but cannot be left empty
# 9 date 
# 10 id
# 
# date column have to be a string that har the format "%Y-%m-%dT%H:%M"

def range_check(df, column, min_val, max_val):
    """Removes values outside the specified range."""
    return df[(df[column] >= min_val) & (df[column] <= max_val)]

def set_neg_depth_to(df, value):
    "Sets all negative depth to value "
    df["depth"] = df["depth"].astype(float)
    df.loc[df["depth"] < 0, "depth"] = value
    return df

def set_neg_oxyg_to(df, near_zero):
    "Sets all negative oxygen to value"
    df["value"] = df["value"].astype(float)
    df.loc[df["value"] < 0, "value"] = near_zero
    return df

def remove_shallow_low_ox_values(df, column, shallow_ox_threshold = 180, shallow_depth = 20):
    """Removes values shallower than given depth and with values below shallow_ox_threshold"""
    return df[~((df[column] <= shallow_ox_threshold) & (df["depth"] <= shallow_depth))]

# Function to downscale the precision of timestamps to the least precise one
def downscale_to_least_precise(ts1, ts2):
    # Compare the precision (e.g., 'ns' for nanosecond precision)
    if ts1.microsecond == 0 and ts2.microsecond == 0:
        ts1 = ts1.normalize()  # Remove time (downscale to day precision)
        ts2 = ts2.normalize()  # Remove time (downscale to day precision)
    elif ts1.minute == 0 and ts2.minute == 0:
        ts1 = ts1.replace(second=0, microsecond=0)  # Remove seconds and microseconds
        ts2 = ts2.replace(second=0, microsecond=0)  # Remove seconds and microseconds
    return ts1, ts2

def remove_matching_rows(df, bad_data, column, log_file):
    """Removes rows in df that match rows in bad_data based on id column and logs removed rows."""
    removed_rows = df[df['id'].isin(bad_data['id'])]
    cleaned_df = df[~df['id'].isin(bad_data['id'])]

    if not removed_rows.empty:
        with open(log_file, "a") as f:
            removed_rows.to_csv(f, sep="\t", index=False, header=f.tell()==0)
    
    return cleaned_df

def format_output(df):
    """Ensures the output dataset follows the specified 11-column format."""
    df['date'] = pd.to_datetime(df['date']).dt.strftime('%Y-%m-%dT%H:%M')
    filler_cols = pd.DataFrame(np.nan, index=df.index, columns=['dummy_4', 'dummy_5', 'dummy_6', 'dummy_7', 'dummy_8'])
    df = pd.concat([df[['lon', 'lat', 'value', 'depth']], filler_cols, df[['date', 'id']]], axis=1)
    return df

def load_emodnet_btl_data(file_path):
    """Loads EMODNET dataset with specific format."""
    df = pd.read_csv(file_path, sep="\t", skiprows=26)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Water body dissolved oxygen concentration [umol/l]', 'Depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID', 'EDMO_code']]
    df['id'] = "emod_btl-" + df['EDMO_code'].astype(str) + "-" + df['LOCAL_CDI_ID'].astype(str)
    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [umol/l]": "value",
            "Depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    df['date'] = df['date'].str[:16]  # Trim milliseconds
    # Convert the 'date' column to datetime using the specified format
    df['date'] = pd.to_datetime(df['date'], format="mixed")

    return df

def load_emodnet_ctd_data(file_path):
    """Loads EMODNET CTD dataset with specific format."""
    df = pd.read_csv(file_path, sep="\t", skiprows=26)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Water body dissolved oxygen concentration [umol/l]', 'Depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID', 'EDMO_code']]

    # Fyll på metadata (alla kolumner utom de som varierar per rad) för varje batch
    metadata_columns = ['Longitude [degrees_east]', 'Latitude [degrees_north]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID', 'EDMO_code']
    df[metadata_columns] = df[metadata_columns].ffill()
    df['id'] = "emod_ctd-" + df['EDMO_code'].astype(str) + "-" + df['LOCAL_CDI_ID'].astype(str)

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [umol/l]": "value",
            "Depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_ices_btl_data(file_path):
    """Loads ICES BTL dataset with specific format."""
    df = pd.read_csv(file_path, sep=",", skiprows=21)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Oxygen (DOXYZZXX_UMLL) [ml/l]', 'QV:ODV:Oxygen (DOXYZZXX_UMLL)', 'Depth (ADEPZZ01_ULAA) [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'Cruise', 'Station']]

    df['id'] = "ices_btl-" + df['Cruise'].astype(str) + "-" + df['Station'].astype(str) + '-' + df['yyyy-mm-ddThh:mm:ss.sss'].astype(str)

    #Remove any rows with NaN Oxygen
    df = df.dropna(subset=["Oxygen (DOXYZZXX_UMLL) [ml/l]"])

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
    removed_ICES_data.to_csv("removed_ices_btl_data.txt", sep="\t", index=False)

    # Change unit from ml/l to µmol/l
    df['Oxygen(DOXYZZXX_UMLL)[µmol/l]'] = df['Oxygen (DOXYZZXX_UMLL) [ml/l]'].apply(lambda x: x * 44.661)

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Oxygen(DOXYZZXX_UMLL)[µmol/l]": "value",
            "Depth (ADEPZZ01_ULAA) [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )

    df['date'] = pd.to_datetime(df['date'], format="mixed")
    df['date'] = df['date'].astype(str).str.strip()
    df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_ices_ctd_data(file_path):
    """Loads ICES CTD dataset with specific format."""
    df = pd.read_csv(file_path, sep=",", skiprows=9)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Oxygen (DOXYZZXX_UMLL) [ml/l]', 'QV:ODV:Oxygen (DOXYZZXX_UMLL)', 'Depth (ADEPZZ01_ULAA) [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'Cruise', 'Station']]

    df['id'] = "ices_ctd-" + df['Cruise'].astype(str) + "-" + df['Station'].astype(str) + '-' + df['yyyy-mm-ddThh:mm:ss.sss'].astype(str)

    #Remove any rows with NaN Oxygen
    df = df.dropna(subset=["Oxygen (DOXYZZXX_UMLL) [ml/l]"])

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
    df['Oxygen(DOXYZZXX_UMLL)[µmol/l]'] = df['Oxygen (DOXYZZXX_UMLL) [ml/l]'].apply(lambda x: x * 44.661)

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Oxygen(DOXYZZXX_UMLL)[µmol/l]": "value",
            "Depth (ADEPZZ01_ULAA) [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )

    df['date'] = pd.to_datetime(df['date'], format="mixed")
    df['date'] = df['date'].astype(str).str.strip()
    df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_iow_data(file_path):
    """Loads IOW CTD/BTL dataset with specific format.
        DOXY6TID [ml/l] = BTL (inkl neg o2)
        DOXY4IDD [ml/l] = CTD
    """
    df = pd.read_csv(file_path, sep="\t", skiprows=20)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'DOXY4IDD [ml/l]', 'DOXY6TID [ml/l]', 'depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'Cruise', 'Station','Type']]

    df['id'] = "iow-" + df['yyyy-mm-ddThh:mm:ss.sss'].astype(str) + "-" + df['Station'].astype(str)

    def get_doxy(row):
        if row['DOXY6TID [ml/l]'] and not np.isnan(row['DOXY6TID [ml/l]']):
            return row['DOXY6TID [ml/l]']
        if row['DOXY4IDD [ml/l]'] and not np.isnan(row['DOXY4IDD [ml/l]']):
            return row['DOXY4IDD [ml/l]']
        return np.nan

    df['Water body dissolved oxygen concentration [ml/l]'] = df.apply(get_doxy, axis=1)
    # Change unit from ml/l to µmol/l
    df['Water body dissolved oxygen concentration [µmol/l]'] = df['Water body dissolved oxygen concentration [ml/l]'].apply(lambda x: x * 44.661)

    df.rename(columns=
        {
            "Longitude [degrees_east]": "lon",
            "Latitude [degrees_north]": "lat",
            "Water body dissolved oxygen concentration [µmol/l]": "value",
            "depth [m]": "depth",
            "yyyy-mm-ddThh:mm:ss.sss": "date",
        }, inplace=True
    )
    df['date'] = df['date'].str[:16]  # Trim milliseconds

    return df

def load_shark_data(file_path):
    """Loads shark tab-delimited dataset."""
    cols = ['depth', 'Dissolved oxygen O2 bottle (ml/l)', 'Q-flag Dissolved oxygen O2 bottle',	'Dissolved oxygen O2 CTD (ml/l)', 'Q-flag Dissolved oxygen O2 CTD',	'Hydrogen sulphide H2S (umol/l)', 'Q-flag Hydrogen sulphide H2S', 'DOXY_ml', 'value', 'NEG_DOXY_umol']

    df = pd.read_csv(file_path, sep="\t", names=cols)
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    return df

def load_data(file_path):
    """Loads generic tab-delimited dataset."""
    cols = ['lon', 'lat', 'value', 'depth', 'drop1', 'drop2', 'drop3', 'drop4', 'drop5', 'date', 'id']
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.columns = cols + df.columns[len(cols):].tolist()
    df = df[['lon', 'lat', 'value', 'depth', 'date', 'id']]
    
    return df

def process_dataset(file_name, is_emodnet_btl=False, is_emodnet_ctd=False, is_iow=False, is_ices_btl=False, is_ices_ctd=False):
    if is_emodnet_btl:
        df = load_emodnet_btl_data(file_name)
    elif is_emodnet_ctd:
        df = load_emodnet_ctd_data(file_name)
    elif is_iow:
        df = load_iow_data(file_name)
    elif is_ices_btl:
        df = load_ices_btl_data(file_name)
    elif is_ices_ctd:
        df = load_ices_ctd_data(file_name)
    else: 
        df = load_data(file_name)
    
    df['date'] = pd.to_datetime(df['date'], format="mixed")
    df = range_check(df, 'value', -500, 600)
    """No the below mess up the Kattegatt analysis /MH"""
    #df = remove_shallow_low_ox_values(df, 'value')

    """Lägg till check på att depth inte är neg. Sätt neg till noll. """
    df = set_neg_depth_to(df, 0)

    """Lägg till check på att alla negativa värden sätts till ngt lågt...
    obsval[obsval .<= 0] .= 0.44661"""
    df = set_neg_oxyg_to(df, 0.44661)

    if not df.loc[~df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)].empty:
        print("Warning")
        print("\tFound nan values in one or more of the required columns:")
        print(f"\t{df.loc[~df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)]}")
        print("\tRemoving these from the dataset")
        df.dropna(subset=['lon', 'lat', 'value', 'depth','date'], inplace=True)

    return df

def main():
    work_dir = "./data"
    data_dir = Path("data/all_baltic")
    output_dir = os.path.join(data_dir, "processed/")
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "removed_rows_log.txt")
    
    # Load and process EMODNET and ICES BTL & CTD datasets
    # print("Load EMODnet btl & ctd....")
    # emodnet_btl_df = process_dataset(os.path.join(data_dir, "BTL_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_emodnet_btl=True)
    # emodnet_ctd_df = process_dataset(os.path.join(data_dir, "CTD_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_emodnet_ctd=True)

    # print("Load IOW data....")
    # iow_df = process_dataset(os.path.join(data_dir, "IOW_data_from_Spreadsheet_Collection_2025-03-25T12-45-01.txt"), is_iow=True)
    print("Load ICES BTL data....")
    data_file = data_dir / "all_ices_data_utf-8.txt"
    print(data_file.resolve())  # shows the absolute path
    # ices_btl_df = process_dataset(data_file, is_ices_btl=True)
    ices_df = process_dataset(os.path.join(data_dir, "ICES_btl_lowres_ctd_02_241107.txt"))
    # print("Load ICES CTD data....")
    # ices_ctd_df = process_dataset(os.path.join(data_dir, "1c9c37d5-2236-42a8-8755-f553c1f4fa0f.csv"), is_ices_ctd=True)
    # print("Load SYKE data...")
    # syke_df = process_dataset(os.path.join(data_dir, "syke_data_no_header_241107.txt"))
    ## SHARK problem, data has no lat, lom, move shark preprocessing from excel_to_tact_sharkweb to this code
    print("Load SHARK data...")
    shark_df = process_dataset(os.path.join(data_dir, "sharkweb_btlctd_02_241107.txt"))

    _, ices_duplicates = find_duplicate_candidates(shark_df, ices_df, suffix="ices")

    with pl.Config(tbl_cols=-1):
        print(ices_duplicates.head(20))
        print(ices_duplicates.height)
    
    # remove duplicates in ices df
    print("remove shark duplicates from ices")
    cleaned_ices = remove_duplicates(shark_df, ices_df)
    
    return
    print("\nCheck_duplicates for emodnet btl and ctd")
    emodnet_df, emod_net_duplicates = find_duplicate_candidates(emodnet_btl_df, emodnet_ctd_df, suffix='emodnet_ctd')
    emodnet_df.to_csv(os.path.join(output_dir, "emodnet.txt"), sep="\t", index=False)
    emod_net_duplicates.to_csv(os.path.join(output_dir, "emodnet_duplicates.txt"), sep="\t", index=False)

    print("\nCheck_duplicates for ices btl and ctd")
    ices_df, ices_net_duplicates = find_duplicate_candidates(ices_btl_df, ices_ctd_df, suffix='ices_ctd')
    ices_df.to_csv(os.path.join(output_dir, "ices.txt"), sep="\t", index=False)
    ices_net_duplicates.to_csv(os.path.join(output_dir, "ices_duplicates.txt"), sep="\t", index=False)

    print("\nCheck_duplicates for emodnet btl/ctd and ices btl/ctd")
    merged_emodnet_ices, emodnet_ices_duplicates = find_duplicate_candidates(emodnet_df, ices_df, suffix='ices')
    merged_emodnet_ices.to_csv(os.path.join(output_dir, "merged_emodnet_ices.txt"), sep="\t", index=False)
    emodnet_ices_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_btl_ices_btl_duplicatecheck_ices_removed.txt"), sep="\t", index=False)

    print("\nCheck_duplicates for syke and emodnet/ices")
    # merged_emodnet_ices = pd.concat([emodnet_df, ices_df], ignore_index=True)
    merged_emodnet_ices_syke, emodnet_ices_syke_duplicates = find_duplicate_candidates(syke_df, merged_emodnet_ices, suffix='emodnet_ices')
    merged_emodnet_ices_syke.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke.txt"), sep="\t", index=False)
    emodnet_ices_syke_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_duplicatecheck_emodnetisec_removed.txt"), sep="\t", index=False)

    print("\nCheck_duplicates for shark vs emodnet/ices/syke")
    # merged_emodnet_ices_syke = pd.concat([merged_emodnet_ices, syke_df], ignore_index=True)
    merged_emodnet_ices_syke_shark, emodnet_ices_syke_shark_duplicates = find_duplicate_candidates(shark_df, merged_emodnet_ices_syke, time_threshold=4, vertical_threshold=2, oxygen_value_threshold=50, name1='shark', name2='emodnet_ices_syke')
    merged_emodnet_ices_syke_shark.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke_shark.txt"), sep="\t", index=False)
    emodnet_ices_syke_shark_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_shark_duplicatecheck_emodneticessyke_removed.txt"), sep="\t", index=False)

    print("\nCheck_duplicates for iow vs emodnet/ices/syke/shark")
    merged_emodnet_ices_syke_shark_iow, emodnet_ices_syke_shark_iow_duplicates = find_duplicate_candidates(iow_df, merged_emodnet_ices_syke_shark, time_threshold=4, vertical_threshold=2, oxygen_value_threshold=50, name1='iow', name2='emodnet_ices_syke_shark')
    merged_emodnet_ices_syke_shark_iow.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke_shark_iow.txt"), sep="\t", index=False)
    emodnet_ices_syke_shark_iow_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_shark_iow_duplicatecheck_emodneticessykeshark_removed.txt"), sep="\t",index=False)

    formatted_df = format_output(merged_emodnet_ices_syke_shark_iow)
    formatted_df.to_csv(os.path.join(output_dir, "duplicate_checked.txt"), sep="\t", index=False)
    
    # formatted_df = load_data(os.path.join(output_dir, "merged_cleaned.txt"))
    # Load bad data
    bad_data_df = pd.read_csv(os.path.join(data_dir, "bad_data.txt"), sep = '\t')
    bad_data_df['date'] = pd.to_datetime(bad_data_df['obstime'], format="mixed")
    bad_data_df['id'] = bad_data_df['obsid']
    bad_data_df['lon'] = bad_data_df['obslon']
    bad_data_df['lat'] = bad_data_df['obslat']
    bad_data_df['value'] = bad_data_df['obsvalue']
    bad_data_df['depth'] = bad_data_df['obsdepth']

    # Remove bad data after merging and log removed rows
    
    cleaned_df = remove_matching_rows(formatted_df, bad_data_df, ['lon', 'lat', 'depth', 'date', 'value', 'id'], log_file)
    # Format output
    df = format_output(cleaned_df)

    # Save cleaned dataset
    df.to_csv(os.path.join(output_dir, "cleaned_df_using_id_col.txt"), sep="\t", index=False)

    cleaned_df_by_dupl_check, bad_data_removed = find_duplicate_candidates(formatted_df, bad_data_df,
                     suffix ='bad_data')
    
    # Format output
    formatted_df = format_output(cleaned_df_by_dupl_check)

    # Save cleaned dataset
    formatted_df.to_csv(os.path.join(output_dir, "cleaned_df_using_dupl_250619.txt"), sep="\t", index=False)
    bad_data_removed.to_csv(os.path.join(output_dir, "bad_data_removed_using_dupl.txt"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
