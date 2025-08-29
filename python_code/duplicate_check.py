import polars as pl
import datetime as dt

def add_coarse_bins(df: pl.DataFrame) -> pl.DataFrame:
    # Set a common startyear (default is 1970, but windows cant use that for years <1970)
    epoch = dt.datetime(1970, 1, 1)#, tzinfo=dt.timezone.utc)
    # Step 1: cast date column properly
    # df = df.with_columns([
    #     pl.col("date").str.to_datetime(time_zone="utc", format="%Y-%m-%dT%H:%M").alias("date")
    # ])
    # print(df.select(["date"]))

    return (
        df
        .with_columns([
            # timebin 25 hours to catch midnight errors
            (((pl.col("date").cast(pl.Datetime("us")) - pl.lit(epoch)).dt.total_seconds() / 3600 ) // 25).cast(pl.Int64).alias("date_bin"), 
            (pl.col("lat").round(2)).alias("lat_bin"),          # spatial bin ~1 km
            (pl.col("lon").round(2)).alias("lon_bin"),          # spatial bin ~1 km
            (pl.col("depth").round(1)).alias("depth_bin")
        ])
    )

def find_duplicate_candidates(best_data, not_so_good_data, suffix="not_so_good_data"):
    # print(best_data["date"])
    # print(not_so_good_data["date"])
    # not_so_good_data['date'] = not_so_good_data['date'].dt.tz_convert(None)  # make naive datetime
    # best_data['date'] = best_data['date'].dt.tz_convert(None)  # make naive datetime
    best_data = add_coarse_bins(pl.from_pandas(best_data))
    not_so_good_data = add_coarse_bins(pl.from_pandas(not_so_good_data))

    candidates = (
        best_data.join(
            not_so_good_data,
            on=["date_bin", "lat_bin", "lon_bin", "depth_bin"],
            how="inner",
            suffix=f"_{suffix}"
        )
    )
    # print(candidates.select(["date", f"date_{suffix}"]))
    candidates = candidates.with_columns([
        (abs((pl.col("date") - pl.col(f"date_{suffix}")).dt.total_minutes())).alias("dt_minutes"),
        ((pl.col("lat") - pl.col(f"lat_{suffix}"))*111_000).alias("dy_m"),    # latdiff i meter
        ((pl.col("lon") - pl.col(f"lon_{suffix}"))*70_000).alias("dx_m"),     # approx för Sverige
        ((pl.col("depth") - pl.col(f"depth_{suffix}"))).alias("dz_m"),        # diff i djupled 
        (abs((pl.col("value") - pl.col(f"value_{suffix}"))).alias("dO2_umol")),
    ])
    # total horisontell distans (pytagoras i plan approx)
    candidates = candidates.with_columns(
        ((pl.col("dx_m")**2 + pl.col("dy_m")**2).sqrt()).alias("dist_horiz_m")
    )

    # slutgiltig filtrering
    duplicates = candidates.filter(
        (pl.col("dt_minutes") <= 90) &      # ≤ 90 min
        (pl.col("dist_horiz_m") <= 1000) &  # ≤ 1000 m
        (pl.col("dz_m") <= 1) &             # ≤ 1 m
        (pl.col("dO2_umol") <= 0.446)       # 0.446 umol/l ca 0.1 ml/l
    )

    return candidates, duplicates

def remove_duplicates(best_data, not_so_good_data, dt = 90, dz = 1, do2 = 0.446, d_horiz = 1000):
    """
    Choose thresholds for 
    dt time in minutes
    dz depth in meters
    do2 oxygen in umol/l
    d_horiz horizontal distance between positions
    """
    
    candidates, _ = find_duplicate_candidates(best_data, not_so_good_data)
    # slutgiltig filtrering
    duplicates = candidates.filter(
        (pl.col("dt_minutes") <= dt) &    # ≤ 90 min
        (pl.col("dist_horiz_m") <= d_horiz) &    # ≤ 1000 m
        (pl.col("dz_m") <= dz) &   # ≤ 1 m
        (pl.col("dO2_umol") >= do2)
    )

    # all ids in df2 that are duplicates
    duplicate_ids = duplicates.select([
        pl.col("id_not_so_good_data").alias("id"),
        pl.col("depth_not_so_good_data").alias("depth")
    ])
    print(f"remove {duplicate_ids.height} rows of {pl.from_pandas(not_so_good_data).height} rows in original")
    print(not_so_good_data.columns)

    cleaned_data = pl.from_pandas(not_so_good_data).join(duplicate_ids, on=["id", "depth"], how="anti")
    removed_data = pl.from_pandas(not_so_good_data).join(duplicate_ids, on=["id", "depth"], how="semi")

    return cleaned_data.to_pandas(), removed_data.to_pandas()