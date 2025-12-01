"""
Merges all GeoTIFF tiles in a folder into a single mosaic and then clips
the mosaic to the study-area boundary defined in a GeoJSON file.
"""

import os
from glob import glob

import rasterio
from rasterio.merge import merge
from rasterio.io import MemoryFile
from rasterio.mask import mask
import fiona

# --------------------------------------------------------------
# USER SETTINGS
# --------------------------------------------------------------

# Path to folder containing *all* your prediction tiles
# (e.g. 300–500 TIFF files)
TILES_DIR = r"E:\git_projects\biebrza-data\model_predictions\tiles"

# Path to GeoJSON boundary used for clipping
BOUNDARY_GEOJSON = r"E:\git_projects\biebrza-shrub-encroachment-analysis\data\biebrza_boundary_clean.geojson"

# Output (clipped) mosaic file
OUTPUT_TIF = r"E:\git_projects\biebrza-shrub-encroachment-analysis\data\mosaic_1997_2015.tif"

# --------------------------------------------------------------


def main():
    # Find all .tif tiles in folder
    tile_paths = sorted(glob(os.path.join(TILES_DIR, "*.tif")))
    print(f"Found {len(tile_paths)} tiles.")

    if len(tile_paths) == 0:
        print("No TIFF files found. Check TILES_DIR.")
        return

    # Open all tiles as rasterio objects
    src_files = [rasterio.open(p) for p in tile_paths]

    # Merge them
    print("Merging tiles (this may take a few minutes)...")
    mosaic, out_transform = merge(src_files)

    # Use profile from the first tile
    out_profile = src_files[0].profile.copy()
    for src in src_files:
        src.close()

    # Update profile for mosaic
    out_profile.update(
        height=mosaic.shape[1],
        width=mosaic.shape[2],
        transform=out_transform,
        compress="lzw",
    )

    # Read clipping geometry from GeoJSON
    print(f"Reading clip geometry from:\n  {BOUNDARY_GEOJSON}")
    with fiona.open(BOUNDARY_GEOJSON, "r") as boundary_src:
        # Collect all geometries in case there are multiple features
        geoms = [feature["geometry"] for feature in boundary_src]

        # NOTE: This assumes the GeoJSON is in the same CRS as the raster tiles.
        # If not, reproject the vector layer beforehand or add reprojection here.

    # Clip the mosaic using the boundary geometry
    print("Clipping mosaic to boundary geometry...")

    # Use an in-memory dataset so we don't have to write the full mosaic to disk
    with MemoryFile() as memfile:
        with memfile.open(**out_profile) as dataset:
            dataset.write(mosaic)

            # mask() returns the clipped array and updated transform
            clipped, clipped_transform = mask(
                dataset,
                geoms,
                crop=True
            )

    # Update profile for clipped mosaic
    out_profile.update(
        height=clipped.shape[1],
        width=clipped.shape[2],
        transform=clipped_transform,
    )

    # Save clipped mosaic
    print(f"Writing clipped mosaic to:\n  {OUTPUT_TIF}")
    with rasterio.open(OUTPUT_TIF, "w", **out_profile) as dest:
        dest.write(clipped)

    print("Done! Clipped mosaic written successfully.")


if __name__ == "__main__":
    main()
