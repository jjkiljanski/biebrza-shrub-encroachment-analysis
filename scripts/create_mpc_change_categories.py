import os
import numpy as np
import geopandas as gpd

# =============================
# USER SETTINGS
# =============================

# Input shapefile from Kopeć & Sławik (2020)
shapefile_path = r"E:\git_projects\biebrza-data\kopec_and_slawik_2020\MPC_1966_2015.shp"

# Output directory for per-category shapefiles
output_dir_categories = r"E:\git_projects\biebrza-data\kopec_and_slawik_2020\mpc_change_categories"

# Column names for Mean Percentage Coverage (MPC) of trees + shrubs
col_1997 = "mpc_1997"
col_2006 = "mpc_2006"   # not used for final categories, but kept for completeness
col_2015 = "mpc_2015"

# Thresholds for translating MPC into cover states
wetland_threshold = 10.0   # MPC < 10% -> wetland (open marsh / meadow)
tree_threshold    = 60.0   # MPC >= 60% -> trees (closed canopy)

# Minimum %-point change to treat as a "real" change
delta_threshold   = 20.0   # |mpc_2015 - mpc_1997| < 10% -> stable


# =============================
# PREPARATION
# =============================

os.makedirs(output_dir_categories, exist_ok=True)

print("Reading shapefile...")
gdf = gpd.read_file(shapefile_path)

# Check required columns
required_cols = [col_1997, col_2006, col_2015]
missing = [c for c in required_cols if c not in gdf.columns]
if missing:
    raise ValueError(f"Missing MPC columns in shapefile: {missing}")

# Compute MPC change 1997–2015 (in %-points)
change_1997_2015_col = "mpc_change_1997_2015"
gdf[change_1997_2015_col] = gdf[col_2015] - gdf[col_1997]


# =============================
# HELPER FUNCTIONS
# =============================

def cover_state(mpc: float,
                wet_th: float = wetland_threshold,
                tree_th: float = tree_threshold) -> str:
    """
    Map MPC value to a discrete state: 'wetland', 'shrubs', 'trees', or 'no_data'.
    """
    if mpc is None or np.isnan(mpc):
        return "no_data"
    if mpc < wet_th:
        return "wetland"
    elif mpc < tree_th:
        return "shrubs"
    else:
        return "trees"


def classify_trajectory(row) -> str:
    """
    Classify 1997 -> 2015 trajectory based on discrete states + MPC change magnitude.
    """
    s0 = row["state_1997"]
    s1 = row["state_2015"]
    dMPC = row[change_1997_2015_col]

    if s0 == "no_data" or s1 == "no_data":
        return "no_data"

    # If states are the same OR the change is very small, call it stable
    if s0 == s1 or abs(dMPC) < delta_threshold:
        return f"stable_{s1}"

    # Otherwise treat as a transition: s0 -> s1
    return f"{s0}_to_{s1}"


# "Nice" categories you primarily care about; others will be grouped as "other"
allowed_categories = {
    "stable_wetland",
    "wetland_to_shrubs",
    "wetland_to_trees",
    "stable_shrubs",
    "shrubs_to_trees",
    "stable_trees",
}


def simplify_trajectory(traj: str) -> str:
    """
    Map the full trajectory label to a simpler set of key categories.
    """
    if traj == "no_data":
        return "no_data"
    if traj in allowed_categories:
        return traj
    return "other"


def sanitize_name(name: str) -> str:
    """
    Make a category string safe for use in a filename.
    """
    name = name.replace("->", "_to_")
    name = name.replace(" ", "_")
    return name


# =============================
# ASSIGN STATES & TRAJECTORIES
# =============================

print("Assigning cover states for each year...")

# States for 1997, 2006, 2015 (we mainly use 1997 & 2015 for trajectories)
gdf["state_1997"] = gdf[col_1997].apply(cover_state)
gdf["state_2006"] = gdf[col_2006].apply(cover_state)
gdf["state_2015"] = gdf[col_2015].apply(cover_state)

print("Classifying trajectories 1997 -> 2015...")
gdf["traj_1997_2015"] = gdf.apply(classify_trajectory, axis=1)
gdf["traj_simple_97_15"] = gdf["traj_1997_2015"].apply(simplify_trajectory)


# =============================
# SUMMARY STATISTICS (PRINT)
# =============================

print("\n=== Full trajectory classes (traj_1997_2015) ===")
vc_full = gdf["traj_1997_2015"].value_counts().sort_index()
for cat, cnt in vc_full.items():
    print(f"{cat:20s} : {cnt}")

print("\n=== Simplified trajectory classes (traj_simple_97_15) ===")
vc_simple = gdf["traj_simple_97_15"].value_counts().sort_index()
for cat, cnt in vc_simple.items():
    print(f"{cat:20s} : {cnt}")


# =============================
# EXPORT ONE SHAPEFILE PER CATEGORY
# =============================

category_col = "traj_simple_97_15"  # use simplified categories for outputs

print(f"\nWriting per-category shapefiles to: {output_dir_categories}\n")

unique_categories = sorted(gdf[category_col].unique())

for category in unique_categories:
    # Optionally skip no_data / other
    if category in ("no_data", "other"):
        print(f"Skipping category '{category}'")
        continue

    subset = gdf[gdf[category_col] == category].copy()
    if subset.empty:
        print(f"Category '{category}' has 0 features, skipping.")
        continue

    safe_cat = sanitize_name(category)
    out_path = os.path.join(
        output_dir_categories,
        f"biebrza_MPC_{safe_cat}.shp"
    )

    subset.to_file(out_path)
    print(f"Written shapefile for category '{category}' -> {out_path}")

print("\nDone.")