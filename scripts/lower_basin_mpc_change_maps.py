import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.backends.backend_pdf import PdfPages

# -----------------------------
# USER SETTINGS
# -----------------------------
shapefile_path = r"E:\git_projects\biebrza-data\kopec_and_slawik_2020\MPC_1966_2015.shp"

# Base MPC columns
col_1997 = "mpc_1997"
col_2006 = "mpc_2006"
col_2015 = "mpc_2015"

# Output mode: "pdf" or "images"
output_mode = "images"      # change to "images" if desired

# For PDF mode
output_pdf = "output/kopec_and_slawik_2020_analysis/biebrza_mpc_changes_and_histograms.pdf"

# For image mode (one file per change period & plot type)
map_image_template = "output/kopec_and_slawik_2020_analysis/biebrza_mpc_change_map_{label}.png"
hist_image_template = "output/kopec_and_slawik_2020_analysis/biebrza_mpc_change_hist_{label}.png"

# Plot appearance
map_cmap_name = "RdBu_r"  # diverging colormap; red=increase, blue=decrease
figsize_map = (8, 8)
figsize_hist = (8, 5)
dpi = 300

# Histogram settings
hist_min = -100
hist_max = 100
bin_width = 10
bins = np.arange(hist_min, hist_max + bin_width, bin_width)

# -----------------------------
# LOAD DATA
# -----------------------------
gdf = gpd.read_file(shapefile_path)

needed_cols = [col_1997, col_2006, col_2015]
missing = [c for c in needed_cols if c not in gdf.columns]
if missing:
    raise ValueError(f"Missing MPC columns in shapefile: {missing}")

# -----------------------------
# COMPUTE CHANGE COLUMNS
# -----------------------------
changes = {
    "1997_2015": (col_1997, col_2015),
    "1997_2006": (col_1997, col_2006),
    "2006_2015": (col_2006, col_2015),
}

change_cols = []

for label, (col_start, col_end) in changes.items():
    change_col = f"mpc_change_{label}"
    gdf[change_col] = gdf[col_end] - gdf[col_start]
    change_cols.append(change_col)

# -----------------------------
# OUTER BOUNDARY (NO GRID LINES)
# -----------------------------
outer_geom = gdf.unary_union
outer_boundary_gs = gpd.GeoSeries([outer_geom.boundary], crs=gdf.crs)

# -----------------------------
# GLOBAL COLOR SCALE FOR CHANGES
# -----------------------------
min_change = float(gdf[change_cols].min().min())
max_change = float(gdf[change_cols].max().max())
abs_max = max(abs(min_change), abs(max_change))

try:
    norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)
except Exception:
    norm = Normalize(vmin=-abs_max, vmax=abs_max)

cmap = get_cmap(map_cmap_name)

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

# -----------------------------
# FUNCTIONS TO MAKE FIGURES
# -----------------------------
def make_map_figure(change_col: str, label: str):
    """Spatial map of MPC change."""
    fig, ax = plt.subplots(figsize=figsize_map, dpi=dpi)

    gdf.plot(
        column=change_col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        edgecolor="none",   # no internal square borders
        linewidth=0
    )

    # Outer park boundary only
    outer_boundary_gs.plot(
        ax=ax,
        edgecolor="black",
        linewidth=0.6,
        facecolor="none"
    )

    start, end = label.split("_")
    ax.set_title(
        f"Biebrza National Park — Change in Tree & Shrub MPC {start}–{end}",
        fontsize=14
    )
    ax.set_axis_off()

    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.03)
    cbar.set_label("Change in Tree & Shrub Mean Percentage Coverage [%-points]")

    plt.tight_layout()
    return fig, ax


def make_hist_figure(change_col: str, label: str):
    """Histogram of MPC change distribution."""
    data = gdf[change_col].dropna()

    fig, ax = plt.subplots(figsize=figsize_hist, dpi=dpi)

    weights = np.ones_like(data) * 100.0 / len(data)

    ax.hist(
        data,
        bins=bins,
        range=(hist_min, hist_max),
        weights=weights,
        edgecolor="black"
    )

    # Vertical line at 0 (no change)
    ax.axvline(0, linestyle="--")

    start, end = label.split("_")
    ax.set_title(
        f"Distribution of Change in Tree & Shrub MPC {start}–{end}",
        fontsize=14
    )
    ax.set_xlabel("Change in MPC [%-points]")
    ax.set_ylabel("Number of squares")

    ax.set_xlim(hist_min, hist_max)

    # Tick marks every 20 units for readability (optional)
    ax.set_xticks(np.arange(hist_min, hist_max + 20, 20))

    plt.tight_layout()
    return fig, ax

# -----------------------------
# OUTPUT: PDF OR IMAGES
# -----------------------------
if output_mode.lower() == "pdf":
    with PdfPages(output_pdf) as pdf:
        for change_col, label in zip(change_cols, changes.keys()):
            # Map page
            fig_map, ax_map = make_map_figure(change_col, label)
            pdf.savefig(fig_map, dpi=dpi, bbox_inches="tight")
            plt.close(fig_map)

            # Histogram page
            fig_hist, ax_hist = make_hist_figure(change_col, label)
            pdf.savefig(fig_hist, dpi=dpi, bbox_inches="tight")
            plt.close(fig_hist)

    print(f"Created multi-page PDF with maps and histograms: {output_pdf}")

elif output_mode.lower() == "images":
    for change_col, label in zip(change_cols, changes.keys()):
        # Map image
        fig_map, ax_map = make_map_figure(change_col, label)
        map_name = map_image_template.format(label=label)
        fig_map.savefig(map_name, dpi=dpi, bbox_inches="tight")
        plt.close(fig_map)
        print(f"Saved change map: {map_name}")

        # Histogram image
        fig_hist, ax_hist = make_hist_figure(change_col, label)
        hist_name = hist_image_template.format(label=label)
        fig_hist.savefig(hist_name, dpi=dpi, bbox_inches="tight")
        plt.close(fig_hist)
        print(f"Saved change histogram: {hist_name}")

else:
    raise ValueError('output_mode must be either "pdf" or "images".')
