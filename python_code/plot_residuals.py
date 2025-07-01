import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score

# Läs in tabbseparerad fil
filnamn = r"C:\Work\DIVAnd\Oxygen_maps\results\Baltic_Proper\20250625_1144\DIVArun\Oxygen_residual.txt"
df = pd.read_csv(filnamn, sep="\t")

#obsval	diva	residual	lat	long	depth	time	id
# Extrahera variabler
obs = df["obsval"]
interp = df["diva"]
residual = df["residual"]
depth = df["obsdepth"]

# Beräkna gemensamma axelgränser
xy_min = min(obs.min(), interp.min())
xy_max = max(obs.max(), interp.max())

# Definiera djupintervall (label, lower, upper)
depth_ranges = [
    ("0–50 m", 0, 50),
    ("50–100 m", 50, 100),
    ("100–150 m", 100, 150),
    (">150 m", 150, float("inf")),
]

# Skapa subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

# Loop över varje djupintervall och plotta
for i, (label, dmin, dmax) in enumerate(depth_ranges):
    # Filtrera data
    mask = (depth > dmin) & (depth <= dmax)
    obs_sel = obs[mask]
    interp_sel = interp[mask]
    depth_sel = depth[mask]

    # Gör scatterplot
    sc = axs[i].scatter(obs_sel, interp_sel, c=depth_sel, cmap='viridis')
    axs[i].set_title(f"{label}")
    axs[i].set_xlabel("Observation")
    axs[i].set_ylabel("Interpolerat")
    axs[i].set_xlim(xy_min, xy_max)
    axs[i].set_ylim(xy_min, xy_max)
    axs[i].set_aspect('equal')

    # Beräkna mått om data finns
    if len(obs_sel) > 0:
        rmse = np.sqrt(np.mean((obs_sel - interp_sel) ** 2))
        bias = np.mean(interp_sel - obs_sel)
        r2 = r2_score(obs_sel, interp_sel)

        textstr = f"RMSE: {rmse:.2f}\nBias: {bias:.2f}\nR²: {r2:.2f}"
        axs[i].text(0.05, 0.95, textstr, transform=axs[i].transAxes,
                    fontsize=10, verticalalignment='top',
                    bbox=dict(facecolor='white', alpha=0.8))

    # Färgbar
    plt.colorbar(sc, ax=axs[i], label="Djup")

# Layout
plt.suptitle("Observation vs Interpolerat värde per djupintervall", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
