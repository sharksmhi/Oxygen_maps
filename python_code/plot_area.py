import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json

def area_bar_plot(results_dir, year_list):
    print("plottar area bars av resultat...")
    year_list = json.loads(year_list)
    # Läs in tab-separerad datafil
    df = pd.read_csv(f'{results_dir}/area_data_{year_list[0]}_{year_list[-1]}.txt', sep="\t")

    # Sortera datan så att säsongerna alltid kommer i ordning
    season_order = ["Winter", "Spring", "Summer", "Autumn"]
    df["season"] = pd.Categorical(df["season"], categories=season_order, ordered=True)
    df = df.sort_values(["year", "season"])

    # Skapa unika år och säsonger
    years = df["year"].unique()
    n_years = len(years)
    season_labels = season_order

    # Skapa positioner för staplar
    bar_width = 0.2

    # Skapa färg för varje säsong
    colors = ["skyblue", "lightgreen", "salmon", "orange"]

    plt.figure(figsize=(18, 6))

    group_width = 0.8  # hur bred en hel årsgrupp är
    season_step = group_width / len(season_labels)  # avstånd mellan säsonger
    x = np.arange(n_years) * 1.5  # större mellanrum mellan åren

    for i, season in enumerate(season_labels):
        df_season = df[df["season"] == season]
        x_pos = x + (i - 1.5) * bar_width  # position för alla staplar i denna säsong

        # Plotta först Area_90_km2 (så Area_0_km2 hamnar ovanpå)
        plt.bar(x_pos,
                df_season["Area_90_km2"],
                width=bar_width,
                color=colors[i],
                alpha=0.5,
                label=f"{season} - Area 90 km²")

        # Plotta Area_0_km2 exakt ovanpå samma x_pos
        plt.bar(x_pos,
                df_season["Area_0_km2"],
                width=bar_width,
                color=colors[i],
                alpha=1.0,
                label=f"{season} - Area 0 km²")

    # X-tick: År
    plt.xticks(x, years, rotation=45)
    plt.xlabel("År")
    plt.ylabel("Area (km²)")
    plt.title("Area per säsong och år")
    plt.legend(ncol=4, bbox_to_anchor=(0.5, -0.15), loc="upper center")
    plt.tight_layout()
    plt.show()

    ###.......###
    print("plottar area bars av background fields...")

    # Läs in tab-separerad datafil
    df = pd.read_csv(f'{results_dir}/Background_area_data.txt', sep="\t")

    # Sortera datan så att säsongerna alltid kommer i ordning
    season_order = ["Winter", "Spring", "Summer", "Autumn"]
    df["season"] = pd.Categorical(df["season"], categories=season_order, ordered=True)
    df = df.sort_values(["year", "season"])

    # Skapa unika år och säsonger
    years = df["year"].unique()
    n_years = len(years)
    season_labels = season_order

    # Skapa positioner för staplar
    bar_width = 0.2

    # Skapa färg för varje säsong
    colors = ["skyblue", "lightgreen", "salmon", "orange"]

    plt.figure(figsize=(18, 6))

    group_width = 0.8  # hur bred en hel årsgrupp är
    season_step = group_width / len(season_labels)  # avstånd mellan säsonger
    x = np.arange(n_years) * 1.5  # större mellanrum mellan åren

    for i, season in enumerate(season_labels):
        df_season = df[df["season"] == season]
        x_pos = x + (i - 1.5) * bar_width  # position för alla staplar i denna säsong

        # Plotta först Area_90_km2 (så Area_0_km2 hamnar ovanpå)
        plt.bar(x_pos,
                df_season["Area_90_km2"],
                width=bar_width,
                color=colors[i],
                alpha=0.5,
                label=f"{season} - Area 90 km²")

        # Plotta Area_0_km2 exakt ovanpå samma x_pos
        plt.bar(x_pos,
                df_season["Area_0_km2"],
                width=bar_width,
                color=colors[i],
                alpha=1.0,
                label=f"{season} - Area 0 km²")

    # X-tick: År
    plt.xticks(x, years, rotation=45)
    plt.xlabel("År")
    plt.ylabel("Area (km²)")
    plt.title("Bakgrundsarea per säsong och mittår i varje period")
    plt.legend(ncol=4, bbox_to_anchor=(0.5, -0.15), loc="upper center")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    print("Plotting barplot of that thing...")
    # Result directory
    results_dir = "C:/Work/DIVAnd/Oxygen_maps/resultat/Baltic_Proper/20250617_0908/"
    year_list = json.dumps([2015])

    area_bar_plot(results_dir, year_list)