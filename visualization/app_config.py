"""
Configuration and constants for the Biebrza encroachment Streamlit app.

Edit this file to:
- change class names / order / colors,
- add more rasters (e.g. a 2005–2023 window, delta maps),
- tweak feature list.
"""

# Class index -> name mapping (dominant_class band)
CLASS_INDEX_TO_NAME = {
    0: "wetland_to_woody",
    1: "shrubs_to_trees",
    2: "stable_wetland",
    3: "stable_trees",
    4: "stable_shrubs",
}

# Canonical ordered list of class names
CLASS_NAMES = [
    "wetland_to_woody",
    "shrubs_to_trees",
    "stable_wetland",
    "stable_trees",
    "stable_shrubs",
]

# Colors per class (for maps and plots) - tweak freely
CLASS_COLORS = {
    "wetland_to_woody": "#d73027",  # red
    "shrubs_to_trees": "#fc8d59",   # orange
    "stable_wetland": "#1a9850",    # green
    "stable_trees": "#4575b4",      # blue
    "stable_shrubs": "#e6ab02",     # yellow-ish
}

# Raster configuration for easy extension
# Assumption: rasters are in EPSG:4326 (lat/lon). If not, reproject beforehand.
RASTER_CONFIG = {
    "1997–2015": {
        "path": "data/mosaic_1997_2015_biebrza_web.tif",
        "type": "mosaic",
        "bands": {
            "dominant_class": 1,   # 1-based indices for rasterio.read
            "p_encroachment": 2,
            "uncertainty": 3,
        },
    },
    # Example future entries:
    # "2005–2023": {
    #     "path": "data/mosaic_2005_2023.tif",
    #     "type": "mosaic",
    #     "bands": {
    #         "dominant_class": 1,
    #         "p_encroachment": 2,
    #         "uncertainty": 3,
    #     },
    # },
    # "Delta encroachment": {
    #     "path": "data/delta_p_encroachment.tif",
    #     "type": "delta",
    #     "bands": {
    #         "delta_p_encroachment": 1,
    #     },
    # },
}

FEATURES = ["NDMI", "NBR", "NDVI", "NIR", "SWIR1", "SWIR2"]