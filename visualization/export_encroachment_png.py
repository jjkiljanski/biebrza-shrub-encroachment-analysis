"""
Export a zoomed + shifted PNG of encroachment probability.

- Zoom factor = 3×
- Pan: shift window towards upper-left (north-west)
- Same colormap as the Streamlit app (plasma)
"""

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from matplotlib import colormaps
from PIL import Image

# -----------------------------
# CONFIG
# -----------------------------
RASTER_PATH = "data/mosaic_1997_2015.tif"
P_ENCROACHMENT_BAND = 2
CMAP_NAME = "plasma"
OUTPUT = "encroachment_zoomed.png"

# Zoom: higher = zoom more
ZOOM_FACTOR = 2.0

# Pan offsets in FRACTIONS of the raster size (0.1 = 10% of width/height)
# Negative dx moves left; negative dy moves up
DX = -0.20     # move viewport 10% left
DY = 0.17     # move viewport 10% up

# Optional: final PNG size
MAX_SIZE = 1200


# -----------------------------
# MAIN
# -----------------------------
def main():

    with rasterio.open(RASTER_PATH) as src:
        nodata = src.nodata

        # Original full bounds
        left, bottom, right, top = src.bounds
        width = right - left
        height = top - bottom

        print("Full bounds:")
        print(" left, bottom, right, top =", src.bounds)

        # --------------------------------------------------
        # 1) Compute zoomed window bounds
        # --------------------------------------------------
        new_width = width / ZOOM_FACTOR
        new_height = height / ZOOM_FACTOR

        # Start in center of raster
        center_x = (left + right) / 2
        center_y = (top + bottom) / 2

        # Apply panning (in geographic units)
        center_x += DX * width
        center_y += DY * height

        # Compute box
        win_left   = center_x - new_width / 2
        win_right  = center_x + new_width / 2
        win_bottom = center_y - new_height / 2
        win_top    = center_y + new_height / 2

        print("Zoom window:")
        print(" left, bottom, right, top =",
              (win_left, win_bottom, win_right, win_top))

        # --------------------------------------------------
        # 2) Read that bounding box
        # --------------------------------------------------
        window = from_bounds(win_left, win_bottom, win_right, win_top, src.transform)

        # Adjust output size to something manageable
        # We target MAX_SIZE for the longer dimension
        h = int((window.height))
        w = int((window.width))

        scale = min(MAX_SIZE / w, MAX_SIZE / h, 1.0)
        out_w = int(w * scale)
        out_h = int(h * scale)

        # Read downsampled
        p_enc = src.read(
            P_ENCROACHMENT_BAND,
            window=window,
            out_shape=(out_h, out_w),
            resampling=Resampling.average,
        ).astype(float)

    # --------------------------------------------------
    # 3) Build RGBA with colormap
    # --------------------------------------------------
    valid = np.isfinite(p_enc)
    if nodata is not None:
        valid &= p_enc != nodata

    norm = np.clip(p_enc, 0, 1)
    cmap = colormaps.get_cmap(CMAP_NAME)

    h, w = norm.shape
    rgba = np.zeros((h, w, 4), dtype=np.uint8)

    # Apply colormap to valid pixels only
    colors = cmap(norm[valid])

    r = np.zeros_like(norm)
    g = np.zeros_like(norm)
    b = np.zeros_like(norm)
    a = np.zeros_like(norm)

    r[valid], g[valid], b[valid], a[valid] = colors.T

    rgba[..., 0] = (255 * r).astype(np.uint8)
    rgba[..., 1] = (255 * g).astype(np.uint8)
    rgba[..., 2] = (255 * b).astype(np.uint8)
    rgba[..., 3] = 0
    rgba[valid, 3] = 255

    img = Image.fromarray(rgba, mode="RGBA")
    img.save(OUTPUT)

    print(f"Saved zoomed PNG: {OUTPUT}")


if __name__ == "__main__":
    main()
