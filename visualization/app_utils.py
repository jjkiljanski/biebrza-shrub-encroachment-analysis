"""
Utility functions for the Biebrza encroachment Streamlit app.

Includes:
- cached data loading (rasters, CSVs),
- helpers to convert rasters to RGBA images for folium overlays,
- plotting functions for representative examples and class trajectories,
- a small legend renderer for the UI.
"""

import numpy as np
import pandas as pd
import streamlit as st

import rasterio
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

import plotly.graph_objects as go

from .app_config import (
    CLASS_INDEX_TO_NAME,
    CLASS_NAMES,
    CLASS_COLORS,
    FEATURES,
)


# --------------------------------------------------------------------------
# Basic helpers
# --------------------------------------------------------------------------


def hex_to_rgb(hex_color: str):
    """Convert '#RRGGBB' -> (R, G, B) in 0-255."""
    rgb = mcolors.to_rgb(hex_color)  # returns (0-1)
    return tuple(int(255 * v) for v in rgb)


# --------------------------------------------------------------------------
# Data loaders (cached)
# --------------------------------------------------------------------------


@st.cache_data(show_spinner=False)
def load_mosaic(path):
    """
    Load a GeoTIFF mosaic into memory.

    Returns
    -------
    data : np.ndarray
        Array of shape (bands, height, width).
    bounds : rasterio.coords.BoundingBox
        Bounds in the raster CRS (assumed EPSG:4326 for folium).
    crs : str
        CRS as WKT or EPSG string.
    transform : affine.Affine
        Geotransform.
    nodata : float or int or None
        NoData value from the raster.
    """
    with rasterio.open(path) as src:
        data = src.read()  # shape (bands, H, W)
        bounds = src.bounds
        crs = src.crs
        transform = src.transform
        nodata = src.nodata
    return data, bounds, str(crs), transform, nodata


@st.cache_data(show_spinner=False)
def load_representative_examples(path="data/representative_examples.csv"):
    """Load representative pixel examples CSV."""
    return pd.read_csv(path)


@st.cache_data(show_spinner=False)
def load_class_trajectories(path="data/class_trajectories.csv"):
    """Load class-level trajectories CSV."""
    return pd.read_csv(path)


# --------------------------------------------------------------------------
# Raster -> RGBA helpers
# --------------------------------------------------------------------------


def create_class_rgba_image(class_band, nodata=None):
    """
    Convert a dominant_class band (H, W) into an RGBA image (H, W, 4)
    using CLASS_COLORS.

    Pixels == nodata (if nodata is not None) get alpha=0 (transparent).
    """
    class_band = class_band.astype(int)
    height, width = class_band.shape
    rgba = np.zeros((height, width, 4), dtype=np.uint8)

    # Default: fully transparent
    rgba[..., 3] = 0

    # Determine valid mask
    if nodata is not None:
        valid = class_band != nodata
    else:
        valid = np.ones_like(class_band, dtype=bool)

    for idx, cls_name in CLASS_INDEX_TO_NAME.items():
        color = CLASS_COLORS.get(cls_name, "#000000")
        r, g, b = hex_to_rgb(color)
        mask = (class_band == idx) & valid
        rgba[mask, 0] = r
        rgba[mask, 1] = g
        rgba[mask, 2] = b
        rgba[mask, 3] = 200  # semi-opaque

    # Any other integer codes (if present) will stay transparent
    return rgba


def create_continuous_rgba_image(
    band_data,
    nodata=None,
    cmap_name="viridis",
    vmin=None,
    vmax=None,
    alpha=200,
    threshold=None,
):
    """
    Convert a continuous band (H, W) into an RGBA image (H, W, 4) using a colormap.

    Parameters
    ----------
    band_data : np.ndarray
        2D array (H, W)
    nodata : float/int/None
        nodata value (optional)
    cmap_name : str
        Name of matplotlib colormap.
    vmin, vmax : float or None
        Min/max for normalization. If None, use data min/max (excluding nodata).
    alpha : int
        Base alpha (0-255) for valid pixels.
    threshold : float or None
        If set, pixels < threshold will be fully transparent (alpha=0).
    """
    arr = band_data.astype(float)
    height, width = arr.shape
    rgba = np.zeros((height, width, 4), dtype=np.uint8)

    # Valid mask: finite and not nodata
    valid = np.isfinite(arr)
    if nodata is not None:
        valid &= arr != nodata

    if not np.any(valid):
        # Nothing valid: return fully transparent
        return rgba

    if vmin is None:
        vmin = np.nanmin(arr[valid])
    if vmax is None:
        vmax = np.nanmax(arr[valid])
    if vmax <= vmin:
        vmax = vmin + 1e-6

    # Normalize to 0-1
    norm = (arr - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)

    cmap = cm.get_cmap(cmap_name)

    # Apply colormap only on valid pixels
    r, g, b, a = np.zeros_like(arr), np.zeros_like(arr), np.zeros_like(arr), np.zeros_like(arr)
    r[valid], g[valid], b[valid], a[valid] = cmap(norm[valid]).T

    rgba[..., 0] = (255 * r).astype(np.uint8)
    rgba[..., 1] = (255 * g).astype(np.uint8)
    rgba[..., 2] = (255 * b).astype(np.uint8)
    rgba[..., 3] = alpha

    # Apply threshold mask: make < threshold transparent
    if threshold is not None:
        mask_threshold = (arr < threshold) & valid
        rgba[mask_threshold, 3] = 0

    return rgba


# --------------------------------------------------------------------------
# Time labels & plotting helpers
# --------------------------------------------------------------------------


def get_time_labels_from_df(df: pd.DataFrame):
    """
    Infer ordered time labels from the representative_examples DataFrame.

    Assumption: all features share the same set of time labels, and column
    names are like 'NDMI_1997_1998', 'NDMI_1999_2000', ...

    We take NDMI_* as a template.
    """
    ndmi_cols = [c for c in df.columns if c.startswith("NDMI_")]
    time_labels = [c.split("_", 1)[1] for c in ndmi_cols]
    # Sort by starting year (before underscore)
    time_labels = sorted(time_labels, key=lambda s: int(s.split("_")[0]))
    return time_labels

def show_continuous_legend(cmap_name, vmin=0.0, vmax=1.0, label="Value"):
    """
    Display a simple horizontal colorbar legend for continuous rasters.

    Intended for encroachment probability and uncertainty (0–1),
    but vmin/vmax are configurable.
    """
    gradient = np.linspace(vmin, vmax, 256).reshape(1, -1)

    fig, ax = plt.subplots(figsize=(4, 0.5))
    im = ax.imshow(gradient, aspect="auto", cmap=cmap_name)
    ax.set_axis_off()

    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.8, pad=0.3)
    cbar.set_label(label)
    # Optional nicer ticks for 0–1
    cbar.set_ticks([vmin, (vmin + vmax) / 2.0, vmax])
    cbar.set_ticklabels([f"{vmin:.1f}", f"{(vmin + vmax) / 2.0:.1f}", f"{vmax:.1f}"])

    st.pyplot(fig, use_container_width=False)


def plot_example_timeseries(row: pd.Series, time_labels, features=None):
    """Plot all feature trajectories over time for a single row using Plotly."""
    if features is None:
        features = FEATURES

    fig = go.Figure()

    for feat in features:
        ys = [row[f"{feat}_{t}"] for t in time_labels]
        fig.add_trace(
            go.Scatter(
                x=time_labels,
                y=ys,
                mode="lines+markers",
                name=feat,
            )
        )

    fig.update_layout(
        xaxis_title="Time",
        yaxis_title="Normalized value",
        margin=dict(l=40, r=20, t=40, b=40),
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_example_probabilities(row: pd.Series):
    """Plot class probability bar chart for a single row."""
    probs = [row.get(f"prob_{cls}", np.nan) for cls in CLASS_NAMES]
    x = CLASS_NAMES
    pred_class = row.get("pred_class_str", None)

    colors = []
    for cls in CLASS_NAMES:
        if cls == pred_class:
            colors.append("#d73027")  # highlight predicted class in red
        else:
            colors.append("#a6bddb")  # muted blue

    fig = go.Figure(
        go.Bar(
            x=x,
            y=probs,
            marker_color=colors,
        )
    )
    fig.update_layout(
        xaxis_title="Class",
        yaxis_title="Predicted probability",
        yaxis=dict(range=[0, 1]),
        margin=dict(l=40, r=20, t=40, b=40),
    )
    st.plotly_chart(fig, use_container_width=True)


def plot_class_trajectories(df_traj: pd.DataFrame, feature: str, selected_classes):
    """Plot mean ± std trajectories per class using Plotly."""
    df_f = df_traj[df_traj["feature"] == feature].copy()
    df_f = df_f[df_f["class_str"].isin(selected_classes)]

    if df_f.empty:
        st.info("No data for this feature / class selection.")
        return

    # Determine global sorted time order
    time_labels = sorted(df_f["time_label"].unique(), key=lambda s: int(s.split("_")[0]))

    fig = go.Figure()

    for cls in selected_classes:
        sub = df_f[df_f["class_str"] == cls].copy()
        sub = sub.set_index("time_label").reindex(time_labels)
        mean = sub["mean"].values
        std = sub["std"].values

        color = CLASS_COLORS.get(cls, "#666666")

        # Main mean line
        fig.add_trace(
            go.Scatter(
                x=time_labels,
                y=mean,
                mode="lines+markers",
                name=cls,
                line=dict(color=color),
            )
        )

        # Shaded std band: upper then lower with fill='tonexty'
        upper = mean + std
        lower = mean - std

        # Convert base color to rgba string with alpha=0.2 for the fill
        r, g, b, a = mcolors.to_rgba(color, alpha=0.2)
        fill_rgba = f"rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, {a})"

        fig.add_trace(
            go.Scatter(
                x=time_labels,
                y=upper,
                mode="lines",
                line=dict(width=0),
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=time_labels,
                y=lower,
                mode="lines",
                line=dict(width=0),
                fill="tonexty",
                fillcolor=fill_rgba,  # <- FIX HERE
                showlegend=False,
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        xaxis_title="Time",
        yaxis_title="Normalized value",
        margin=dict(l=40, r=20, t=40, b=40),
        hovermode="x unified",
    )
    st.plotly_chart(fig, use_container_width=True)

def show_class_legend():
    """Render a simple class legend in Streamlit (colored squares + labels)."""
    st.markdown("**Class legend**")
    for cls in CLASS_NAMES:
        color = CLASS_COLORS.get(cls, "#000000")
        st.markdown(
            f"<span style='display:inline-block;width:16px;height:16px;"
            f"background-color:{color};margin-right:8px;border-radius:2px;'></span> {cls}",
            unsafe_allow_html=True,
        )
