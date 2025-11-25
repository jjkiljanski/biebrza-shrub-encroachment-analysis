import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.backends.backend_pdf import PdfPages

# -----------------------------
# USER SETTINGS
# -----------------------------
shapefile_path = r"E:\git_projects\biebrza-data\kopec_and_slawik_2020\MPC_1966_2015.shp"

# Columns with Tree & Shrub MPC
mpc_cols = ["mpc_1966", "mpc_1980", "mpc_1997",
            "mpc_2006", "mpc_2010", "mpc_2015"]

# Output mode: "pdf" or "images"
output_mode = "images"      # change to "images" if desired

# For PDF mode
output_pdf = "output/kopec_and_slawik_2020_analysis/biebrza_mpc_all_years.pdf"

# For image mode (one file per year)
image_template = "output/kopec_and_slawik_2020_analysis/biebrza_mpc_{year}.png"

# Plot appearance
cmap_name = "YlGn"
figsize = (8, 8)
dpi = 300

# -----------------------------
# LOAD DATA
# -----------------------------
gdf = gpd.read_file(shapefile_path)

missing = [c for c in mpc_cols if c not in gdf.columns]
if missing:
    raise ValueError(f"Missing columns in shapefile: {missing}")

# -----------------------------
# OUTER BOUNDARY (NO GRID LINES)
# -----------------------------
# This creates a single geometry for the whole park,
# then we'll plot just its boundary as an outline.
outer_geom = gdf.unary_union
# Wrap in GeoSeries so we can plot it easily
outer_boundary_gs = gpd.GeoSeries([outer_geom.boundary], crs=gdf.crs)

# -----------------------------
# GLOBAL COLOR SCALE
# -----------------------------
vmin = float(gdf[mpc_cols].min().min())
vmax = float(gdf[mpc_cols].max().max())

norm = Normalize(vmin=vmin, vmax=vmax)
cmap = get_cmap(cmap_name)

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

# -----------------------------
# FUNCTION TO MAKE ONE FIGURE
# -----------------------------
def make_figure_for_column(col: str):
    year = col.split("_")[-1]

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot the squares filled by MPC, WITHOUT their borders
    gdf.plot(
        column=col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        edgecolor="none",   # no square borders
        linewidth=0
    )

    # Overlay only the outer boundary
    outer_boundary_gs.plot(
        ax=ax,
        edgecolor="black",
        linewidth=0.6,
        facecolor="none"
    )

    ax.set_title(f"Biebrza National Park — Tree & Shrub MPC {year}", fontsize=14)
    ax.set_axis_off()

    # Colorbar
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.03)
    cbar.set_label("Tree & Shrub Mean Percentage Coverage [%]")

    plt.tight_layout()
    return fig, ax, year

# -----------------------------
# OUTPUT: PDF OR IMAGES
# -----------------------------
if output_mode.lower() == "pdf":
    with PdfPages(output_pdf) as pdf:
        for col in mpc_cols:
            fig, ax, year = make_figure_for_column(col)
            pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
    print(f"Created multi-page PDF: {output_pdf}")

elif output_mode.lower() == "images":
    for col in mpc_cols:
        fig, ax, year = make_figure_for_column(col)
        out_name = image_template.format(year=year)
        fig.savefig(out_name, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved image: {out_name}")

else:
    raise ValueError('output_mode must be either "pdf" or "images".')