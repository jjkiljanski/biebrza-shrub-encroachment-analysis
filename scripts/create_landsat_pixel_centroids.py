import rasterio
from rasterio.transform import xy
import geopandas as gpd
from shapely.geometry import Point
import numpy as np

tiff_path = r"E:\git_projects\biebrza-data\landsat_photos\biebrza_1997.tif"
output_path = r"E:\git_projects\biebrza-data\shapefiles\landast_pixel_centroids.geojson"

# Open raster
with rasterio.open(tiff_path) as src:
    transform = src.transform
    width = src.width
    height = src.height
    crs = src.crs
    band1 = src.read(1)  # read first band (optional, for value attribute)

# Prepare lists for geometries and attributes
geoms = []
values = []

for row in range(height):
    for col in range(width):
        # Skip NoData pixels if you want
        if band1[row, col] == src.nodata:
            continue
        
        # Get center coordinate of pixel
        x, y = xy(transform, row, col, offset='center')
        geoms.append(Point(x, y))
        values.append(band1[row, col])

# Build GeoDataFrame
gdf = gpd.GeoDataFrame(
    {"value": values},
    geometry=geoms,
    crs=crs
)

# Save as GeoJSON (or .shp, .gpkg, etc.)
gdf.to_file(output_path, driver="GeoJSON")
print(f"Saved {len(gdf)} pixel centroids to {output_path}")