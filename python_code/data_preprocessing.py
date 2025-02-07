import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.neighbors import BallTree
from datetime import datetime, timedelta


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

def check_duplicates(df1, df2, distance_threshold, time_threshold, name1 = "df1", name2 = "df2"):
    # Combine the lat/lon distance and time difference in the same space (Euclidean distance)
    # Assuming normalized time (in hours) and normalized lat/lon in meters
    combined_threshold = np.sqrt(distance_threshold**2 + time_threshold**2)
    # Convert the 'date' column to datetime (if it's not already in datetime format)
    df1['date'] = pd.to_datetime(df1['date'])
    df2['date'] = pd.to_datetime(df2['date'])

    # Convert the 'date' column to a numerical format (seconds since epoch)
    df1['date_num'] = (df1['date'] - pd.Timestamp("1900-01-01T00:00")) // pd.Timedelta('1h')
    df2['date_num'] = (df2['date'] - pd.Timestamp("1900-01-01T00:00")) // pd.Timedelta('1h')

    # Convert lat, lon to meters (using rough approximations)
    # 1 degree of latitude ≈ 111,000 meters
    df1['lat_meters'] = df1['lat'] * 111000  # Convert latitude to meters
    df2['lat_meters'] = df2['lat'] * 111000  # Convert latitude to meters

    # 1 degree of longitude ≈ 111,000 meters * cos(latitude)
    mean_latitude = df1['lat'].mean()  # Use the mean latitude for an approximation
    longitude_conversion_factor = 111000 * np.cos(np.radians(mean_latitude))

    df1['lon_meters'] = df1['lon'] * longitude_conversion_factor  # Convert longitude to meters
    df2['lon_meters'] = df2['lon'] * longitude_conversion_factor  # Convert longitude to meters

    # Now we have lat, lon, time, and depth all in comparable units (meters or hours)
    # We will normalize the columns
    dataset1 = df1[['lon_meters', 'lat_meters', 'date_num', 'depth']].values
    dataset2 = df2[['lon_meters', 'lat_meters', 'date_num', 'depth']].values

    # Create a BallTree for dataset1
    tree1 = BallTree(dataset1)

    # Query dataset2 points against the BallTree (find nearest neighbors)
    distances, indices = tree1.query(dataset2, k=1)  # k=1 to find the nearest point in dataset1

    # Set a distance threshold to define duplicates (e.g., points considered "duplicates" if they are very close)
    distance_threshold = 0.1  # Adjust this based on your data's precision

    # Find duplicates based on distance threshold
    duplicates = distances.flatten() < combined_threshold

    # Get indices of duplicate points in df2
    duplicate_indices_df2 = np.where(duplicates)[0]

    # Get indices of non-duplicate points in df2
    non_duplicate_indices_df2 = np.where(~duplicates)[0]

    # Create df2_filtered (removing duplicates from df2)
    df2_filtered = df2.iloc[non_duplicate_indices_df2]

    # Create df3 (duplicates from df1 and df2)
    # We'll create a new column 'source' to indicate which dataset the point came from
    duplicates_df1 = df1.iloc[indices[duplicate_indices_df2].flatten()]  # Use indices from BallTree for df1
    duplicates_df1.loc[:, 'source'] = name1

    duplicates_df2 = df2.iloc[duplicate_indices_df2]
    duplicates_df2.loc[:, 'source'] = name2

    print(f"\nDuplicates found {name2} (unique lon, lat, date): {duplicates_df2.drop_duplicates(subset=['lon', 'lat', 'date']).shape[0]}")
    print(f"Duplicates found {name1} (unique lon, lat, date): {duplicates_df1.drop_duplicates(subset=['lon', 'lat', 'date']).shape[0]}")

    # Combine the duplicates from df1 and df2 into df3
    duplicates = pd.concat([duplicates_df1, duplicates_df2])

    # Reset index of df3
    duplicates.reset_index(drop=True, inplace=True)

    merged = pd.concat([df1, df2_filtered], ignore_index=True)

    return merged, duplicates


# def _depr_check_duplicates(df1, df2, columns, tolerance, time_tolerance):
#     """Checks for duplicates within 0.01 degrees distance and a 10-minute timespan."""
#     df2_filtered = df2.copy()
#     df2_filtered['keep'] = True
#     df2_filtered[['lon', 'lat']] = df2_filtered[['lon', 'lat']].replace([np.inf, -np.inf], np.nan) # Replace inf with NaN
#     # Drop rows where any of the specified columns have NaN values
#     df2_filtered.dropna(subset=['lon', 'lat'], inplace=True)
    
#     coords1 = np.radians(df1[['lon', 'lat']].values)
#     coords2 = np.radians(df2_filtered[['lon', 'lat']].values)
#     tree = KDTree(coords1)
    
#     indices = tree.query_ball_point(coords2, np.radians(tolerance))
    
#     time_threshold = pd.Timedelta(minutes=time_tolerance)
    
#     for i, nearby in enumerate(indices):
#         if nearby:
#             # time_diffs = [abs((datetime.strptime(df1.iloc[idx]['date'], "%Y-%m-%dT%H:%M") - datetime.strptime(df2.iloc[i]['date'], "%Y-%m-%dT%H:%M")).total_seconds()) for idx in nearby]
#             for idx in nearby:
#                 try:
#                     abs(df1.iloc[idx]['date'] - df2_filtered.iloc[i]['date']) 
#                 except TypeError:
#                     print(f"df1.iloc[idx]['date'] {df1.iloc[idx]['date']}")
#                     print(f"df2.iloc[i]['date'] {df2_filtered.iloc[i]['date']}")
#                     exit()
#             time_diffs = [abs(df1.iloc[idx]['date'] - df2_filtered.iloc[i]['date']) for idx in nearby]
#             if any(diff <= time_threshold for diff in time_diffs):  # 600 seconds = 10 minutes
#                 df2_filtered.at[i, 'keep'] = False
    
#     return df2_filtered[df2_filtered['keep']].drop(columns=['keep']), df2_filtered[~df2_filtered['keep']].drop(columns=['keep'])

def remove_matching_rows(df, bad_data, columns, log_file):
    """Removes rows in df that match rows in bad_data based on specific columns and logs removed rows."""
    merged = df.merge(bad_data, on=columns, how='left', indicator=True)
    removed_rows = merged.query('_merge == "both"').drop('_merge', axis=1)
    cleaned_df = merged.query('_merge == "left_only"').drop('_merge', axis=1)
    
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
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Water body dissolved oxygen concentration [umol/l]', 'Depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID']]
    df.columns = ['lon', 'lat', 'value', 'depth', 'date', 'id']
    df['date'] = df['date'].str[:16]  # Trim milliseconds
    # Convert the 'date' column to datetime using the specified format
    df['date'] = pd.to_datetime(df['date'], format="mixed")

    return df

def load_odv_data(file_path):
    """Loads EMODNET CTD dataset with specific format."""
    df = pd.read_csv(file_path, sep="\t", skiprows=26)
    df = df[['Longitude [degrees_east]', 'Latitude [degrees_north]', 'Water body dissolved oxygen concentration [umol/l]', 'Depth [m]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID']]

    # Fyll på metadata (alla kolumner utom de som varierar per rad) för varje batch
    metadata_columns = ['Longitude [degrees_east]', 'Latitude [degrees_north]', 'yyyy-mm-ddThh:mm:ss.sss', 'LOCAL_CDI_ID']
    df[metadata_columns] = df[metadata_columns].ffill()

    df.columns = ['lon', 'lat', 'value', 'depth', 'date', 'id']
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

def process_dataset(file_name, is_emodnet_btl=False, is_odv=False):
    if is_emodnet_btl:
        df = load_emodnet_btl_data(file_name)
    elif is_odv:
        df = load_odv_data(file_name)
    else: 
        df = load_data(file_name)
    
    df['date'] = pd.to_datetime(df['date'], format="mixed")
    df = range_check(df, 'value', -500, 600)
    if not df.loc[df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)].empty:
        print("Warning")
        print("\tFound nan values in one or more of the required columns:")
        print(f"\t{df.loc[df.index.isin(df.dropna(subset=['lon', 'lat', 'value', 'depth','date']).index)]}")
        print("\tRemoving these from the dataset")
        df.dropna(subset=['lon', 'lat', 'value', 'depth','date'], inplace=True)

    return df

def main():
    data_dir = "Oxygen_maps/data/all_baltic/"
    output_dir = os.path.join(data_dir, "processed/")
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "removed_rows_log.txt")
    
    # Load and process EMODNET and ICES BTL & CTD datasets
    print("Load EMODnet....")
    emodnet_btl_df = process_dataset(os.path.join(data_dir, "BTL_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_emodnet_btl=True)
    emodnet_ctd_df = process_dataset(os.path.join(data_dir, "CTD_data_from_EMODnet_Eutrophication_European_2024_unrestricted.txt"), is_odv=True)
    # print("check_duplicates for emodnet btl and ctd")
    # emodnet_df, emod_net_duplicates = check_duplicates(emodnet_btl_df, emodnet_ctd_df, time_threshold=4, distance_threshold=1000)
    # emodnet_df.to_csv(os.path.join(output_dir, "emodnet.txt"), sep="\t", index=False)
    # emod_net_duplicates.to_csv(os.path.join(output_dir, "emodnet_duplicates.txt"), sep="\t", index=False)

    print("Load ICES....")
    ices_df = process_dataset(os.path.join(data_dir, "ICES_btl_lowres_ctd_02_241107.txt"))
    print("\nCheck_duplicates for emodnet btl and ices low res")
    merged_emodnet_ices, emodnet_ices_duplicates = check_duplicates(emodnet_btl_df, ices_df, time_threshold=4, distance_threshold=1000, name1='emodnet_btl', name2='ices_low_res')
    merged_emodnet_ices.to_csv(os.path.join(output_dir, "merged_emodnet_ices.txt"), sep="\t", index=False)
    emodnet_ices_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_duplicates.txt"), sep="\t", index=False)
    print("merge emdodnet and ices")
    # merged_emodnet_ices = pd.concat([emodnet_df, ices_df], ignore_index=True)
    
    syke_df = process_dataset(os.path.join(data_dir, "syke_data_no_header_241107.txt"))
    merged_emodnet_ices_syke, emodnet_ices_syke_duplicates = check_duplicates(syke_df, merged_emodnet_ices,  time_threshold=4, distance_threshold=1000, name1='syke', name2='emodnet_ices')
    merged_emodnet_ices_syke.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke.txt"), sep="\t", index=False)
    emodnet_ices_syke_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_duplicates.txt"), sep="\t", index=False)
    print("merge emdodnet and ices")
    # merged_emodnet_ices_syke = pd.concat([merged_emodnet_ices, syke_df], ignore_index=True)
    
    ## SHARK problem, data has no lat, lom, move shark preprocessing from excel_to_tact_sharkweb to this code
    shark_df = process_dataset(os.path.join(data_dir, "sharkweb_btlctd_02_241107.txt"))
    merged_emodnet_ices_syke_shark, emodnet_ices_syke_shark_duplicates = check_duplicates(shark_df, merged_emodnet_ices_syke, time_threshold=4, distance_threshold=1000, name1='shark', name2='emodnet_ices_syke')
    merged_emodnet_ices_syke_shark.to_csv(os.path.join(output_dir, "merged_emodnet_ices_syke_shark.txt"), sep="\t", index=False)
    emodnet_ices_syke_shark_duplicates.drop_duplicates(subset=['lon', 'lat', 'date', 'source']).to_csv(os.path.join(output_dir, "emodnet_ices_syke_shark_duplicates.txt"), sep="\t", index=False)
    
    formatted_df = format_output(merged_emodnet_ices_syke_shark)
    formatted_df.to_csv(os.path.join(output_dir, "merged_cleaned.txt"), sep="\t", index=False)
    
    exit()
    # Load bad data
    bad_data_df = load_data(os.path.join(data_dir, "bad_data.txt"))
    
    # Remove bad data after merging and log removed rows
    cleaned_df = remove_matching_rows(merged_emodnet_ices_syke_shark, bad_data_df, ['lon', 'lat', 'depth', 'date', 'value', 'id'], log_file)
    
    # Format output
    formatted_df = format_output(cleaned_df)

    # Save cleaned dataset
    formatted_df.to_csv(os.path.join(output_dir, "merged_cleaned.txt"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
