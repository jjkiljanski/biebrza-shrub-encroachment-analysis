"""
Prepare lighter rasters for the Streamlit web app:

1. Clip the original mosaic to the Biebrza boundary.
2. Downsample the clipped raster to a smaller size for web display.

Inputs
------
RASTER_IN  : full mosaic (e.g. for the whole tile set or region)
BOUNDARY   : Biebrza NP boundary as a GeoJSON
Outputs
-------
RASTER_CLIPPED : mosaic clipped to boundary
RASTER_WEB     : further downsampled version of RASTER_CLIPPED for web

You can adjust:
- MAX_WEB_PIXELS: how many pixels the web version should roughly have
                  (H * W <= MAX_WEB_PIXELS)
- COMPRESS / BIGTIFF options as you like.
"""

import os

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.enums import Resampling
from affine import Affine


# ---------------------------------------------------------------------
# Paths (adjust if needed)
# ---------------------------------------------------------------------
RASTER_IN = "data/mosaic_1997_2015.tif"
BOUNDARY = r"data/biebrza_boundary_clean.geojson"
RASTER_CLIPPED = "data/mosaic_1997_2015_biebrza_clipped.tif"
RASTER_WEB = "data/mosaic_1997_2015_biebrza_web.tif"

# Rough target for web raster pixel count (H * W).
# 4 million pixels -> e.g. 2000x2000; tweak as you like.
MAX_WEB_PIXELS = 4_000_000


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def clip_raster_to_boundary(raster_in, boundary_geojson, raster_out):
    """
    Clip raster_in to boundary_geojson and save as raster_out.

    Keeps all bands, nodata, CRS. Uses rasterio.mask.mask with crop=True.
    """
    print(f"Clipping raster:\n  {raster_in}\nusing boundary:\n  {boundary_geojson}")

    # Load boundary
    gdf = gpd.read_file(boundary_geojson)

    with rasterio.open(raster_in) as src:
        # Reproject boundary to match raster CRS, if needed
        if gdf.crs is None:
            raise ValueError("Boundary GeoJSON has no CRS defined.")
        if gdf.crs != src.crs:
            gdf = gdf.to_crs(src.crs)

        geom = gdf.unary_union  # merge all geometries
        geoms = [geom]

        # Clip
        out_image, out_transform = mask(
            src,
            geoms,
            crop=True,
            nodata=src.nodata,
        )

        out_meta = src.meta.copy()
        out_meta.update(
            {
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )

    # Write clipped raster
    os.makedirs(os.path.dirname(raster_out), exist_ok=True)
    with rasterio.open(raster_out, "w", **out_meta) as dst:
        dst.write(out_image)

    print(f"Clipped raster written to:\n  {raster_out}")


def create_web_raster(raster_in, raster_out, max_pixels=MAX_WEB_PIXELS):
    """
    Downsample a raster for web visualization.

    - Uses nearest-neighbour resampling (good for class labels in band 1).
    - Keeps all bands, nodata, CRS.
    - Scales resolution so that H * W ≲ max_pixels.
    """
    print(f"Creating web raster from:\n  {raster_in}")
    with rasterio.open(raster_in) as src:
        height, width = src.height, src.width
        bands = src.count

        total_pixels = height * width
        print(f"Original clipped size: {width} x {height} ({total_pixels:,} pixels)")

        if total_pixels <= max_pixels:
            print("Raster is already small enough; copying without resampling.")
            data = src.read()
            out_transform = src.transform
        else:
            # Compute scale factor so that H*W ~ max_pixels
            scale = (total_pixels / max_pixels) ** 0.5
            out_height = max(1, int(height / scale))
            out_width = max(1, int(width / scale))
            print(
                f"Downsampling to approx: {out_width} x {out_height} "
                f"({out_width * out_height:,} pixels)"
            )

            data = src.read(
                out_shape=(bands, out_height, out_width),
                resampling=Resampling.nearest,  # preserves class labels
            )

            # Compute new transform
            scale_x = src.width / out_width
            scale_y = src.height / out_height
            out_transform = src.transform * Affine.scale(scale_x, scale_y)

        # Prepare metadata for writing
        out_meta = src.meta.copy()
        out_meta.update(
            {
                "height": data.shape[1],
                "width": data.shape[2],
                "transform": out_transform,
                # Optionally add compression to shrink file size:
                "compress": "lzw",
                "tiled": True,
            }
        )

    os.makedirs(os.path.dirname(raster_out), exist_ok=True)
    with rasterio.open(raster_out, "w", **out_meta) as dst:
        dst.write(data)

    print(f"Web raster written to:\n  {raster_out}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
if __name__ == "__main__":
    # 1) Clip original mosaic to Biebrza boundary
    clip_raster_to_boundary(RASTER_IN, BOUNDARY, RASTER_CLIPPED)

    # 2) Create downsampled web version from the clipped raster
    create_web_raster(RASTER_CLIPPED, RASTER_WEB, max_pixels=MAX_WEB_PIXELS)

    print("Done.")
