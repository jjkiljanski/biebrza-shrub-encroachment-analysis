# app.py
"""
Streamlit app entry point for exploring shrub/woody encroachment Conv1D
model outputs in Biebrza National Park.

Run with:
    streamlit run app.py
"""

import streamlit as st
import folium
from folium.raster_layers import ImageOverlay
from streamlit_folium import st_folium

from visualization.app_config import (
    CLASS_NAMES,
    RASTER_CONFIG,
)
from visualization.app_utils import (
    load_mosaic,
    load_representative_examples,
    load_class_trajectories,
    create_class_rgba_image,
    create_continuous_rgba_image,
    get_time_labels_from_df,
    plot_example_timeseries,
    plot_example_probabilities,
    plot_class_trajectories,
    show_class_legend,
    show_continuous_legend
)


# -----------------------------------------------------------------------------
# PAGE CONFIG
# -----------------------------------------------------------------------------

st.set_page_config(
    page_title="Biebrza Encroachment Explorer",
    layout="wide",
)

st.title("Shrub Encroachment Explorer – Biebrza National Park")

# Tabs: add a landing / About tab as requested
tab_about, tab_map, tab_examples, tab_trajectories = st.tabs(
    ["🏠 About", "🗺️ Map", "🔬 Representative examples", "📈 Class-level trajectories"]
)

# -----------------------------------------------------------------------------
# TAB 0: ABOUT / LANDING PAGE
# -----------------------------------------------------------------------------
with tab_about:
    st.subheader("About this project")

    # 🔁 Replace this demo text with your own description later
    st.markdown(
        """

# 🌿 Monitoring Shrub Encroachment in Biebrza National Park Using Landsat Time Series and Deep Learning

Shrub and woody vegetation encroachment is one of the most pressing ecological threats to **Biebrza National Park**, one of Europe’s largest and most pristine peatland complexes. The Biebrza wetlands support a unique mosaic of marshes, sedge meadows, fens, and open wet grasslands shaped by centuries of extensive, low-intensity agriculture. Over recent decades, socio-economic changes have led to the abandonment of traditional mowing and grazing practices. As these open wetlands lose regular disturbance, **shrubs and young trees spread rapidly**, altering hydrology, reducing biodiversity, and transforming habitat structure for species of conservation concern.

Understanding *where* and *how fast* encroachment occurs is essential for conservation management, restoration planning, and evaluating ecological resilience.

---

## 📘 Background: MPC Trajectories From Kopeć & Sławik (2020)

A key reference point for this project is the work of **Kopeć & Sławik (2020)**, who developed a robust framework for analyzing long-term vegetation change using the **Mean Percentage Coverage (MPC)** concept. They derived shrub/tree cover estimates for **125 × 125 m grid cells** in the **lower Biebrza basin** using multi-temporal Landsat imagery. These MPC trajectories enabled a classification of land-cover change processes such as *wetland stability*, *shrub encroachment*, and *succession to forest*.

In this project, their published MPC dataset is used as **weak supervision** to train a machine-learning model at a finer spatial resolution. Each MPC square provides a coarse trajectory label, which is then transferred to all 30-m Landsat pixels within that square.

**Full citation (APA):**

Kopeć, D., & Sławik, Ł. (2020). How to effectively use long-term remotely sensed data to analyze the process of tree and shrub encroachment into open protected wetlands. Applied Geography, 125, 102345. https://doi.org/10.1016/j.apgeog.2020.102345

---

## 🛠️ Technical Overview

### **Satellite Data**
- Landsat 5/7/8 surface reflectance (1997–2015)
- Biannual composites (10 time steps)
- Per-pixel spectral & index features (6 channels):
  - **NDMI**, **NBR**, **NDVI**, **NIR**, **SWIR1**, **SWIR2**

### **Model Input**
Each pixel is represented by a **(T=10, C=6)** time series.

### **Label Space (5 classes)**
After merging related categories, the model predicts:
- **wetland_to_woody**
- **shrubs_to_trees**
- **stable_wetland**
- **stable_trees**
- **stable_shrubs**

### **Tech Stack**
- **Google Earth Engine (GEE)** — for generating composites & sampling pixels  
- **Python (Colab)** — for all Machine Learning workflows  
- **PyTorch** — model development & training  
- **rasterio / numpy / pandas** — data preparation  
- **Streamlit** — interactive visualization platform

---

## 🧹 Data Cleaning & Preparation

Because MPC labels apply to large 125 m squares while training is done per 30 m pixel, the dataset is inherently **noisy**. Several cleaning steps were applied:

1. **Correcting label inconsistencies**  
   Merging closely related MPC transitions into a unified class (`wetland_to_woody`).

2. **Stratified split by MPC square**  
   Ensures that all pixels from the same square fall into the same train/val/test partition.

3. **Normalization using training statistics only**  
   Mean/std computed per band per time step.

4. **Removing ambiguous or low-information pixels**  
   Pixels that were extremely ambiguous (flat spectra, no temporal pattern, or highly contradictory labels) were dropped to improve robustness.

This cleaning substantially improved prediction quality, particularly for transition classes.

---

## 🤖 Model Training and Selection

Two model architectures were investigated:

### **1. 1D Convolutional Neural Network (Conv1D)**  
- Learns short- and medium-term temporal patterns.
- Fast to train and efficient for inference.

### **2. LSTM Recurrent Neural Network**  
- Better theoretical capacity for long-time dependencies.

Both models were trained on the cleaned dataset using class weighting to compensate for class imbalance.

### **Results**
Final test accuracies:
- **Conv1D:** 0.741  
- **LSTM:** 0.755

Despite the LSTM scoring slightly higher overall, the **Conv1D model demonstrated clearly better precision and recall for the most important target: the `wetland_to_woody` encroachment class**. Since conservation relevance prioritizes detecting actual encroachment over marginal gains in total accuracy, the **Conv1D classifier was chosen as the final model**.

---

## 🌍 Generating Spatial Predictions

The final Conv1D model is applied in GEE to the entire Biebrza landscape.  
To optimize performance, the area is processed in **252 × 252 pixel tiles**, each batch containing:

- the full 10-step time series for 6 features  
- inference on ~63,500 pixels per tile  
- careful memory management to ensure smooth operation

This approach balances:
- computational feasibility,
- reliable batching,
- and consistent tile-based normalization,

while enabling high-resolution predictions across large wetland areas.

### **Prediction Outputs**
Each raster mosaic includes:

1. **Dominant class** (argmax over 5 classes)  
2. **Probability of encroachment** (P(wetland_to_woody))  
3. **Prediction uncertainty** (1 – max probability)

These maps provide actionable insights into spatial patterns of past woody encroachment and hotspots of ongoing transition.

---

## 🎯 Purpose of the Streamlit App

The app allows users to:

- **Explore spatial predictions** over Biebrza using interactive maps  
- **Inspect representative pixel-level time series** (both correct and misclassified cases)  
- **Understand class-level temporal signatures** (NDMI/NBR/NDVI/NIR/SWIR1/SWIR2 mean trajectories)  
- Build ecological intuition about how different land-cover trajectories express themselves in spectral time-series data.

This is intended as both a communication tool and a scientific companion to the modeling workflow.

"""
    )

# -----------------------------------------------------------------------------
# TAB 1: MAP
# -----------------------------------------------------------------------------
with tab_map:
    st.subheader("Prediction mosaic map")

    col_controls, col_map = st.columns([1, 2])

    with col_controls:
        st.markdown("### Map settings")

        period = st.selectbox(
            "Mosaic period",
            options=list(RASTER_CONFIG.keys()),
            index=0,
        )
        raster_conf = RASTER_CONFIG[period]
        raster_path = raster_conf["path"]

        map_type = st.radio(
            "Map type",
            options=["Dominant class", "Encroachment probability", "Uncertainty"],
        )

        basemap_choice = st.selectbox(
            "Basemap",
            options=["OpenStreetMap", "ESRI WorldImagery", "CartoDB Positron"],
            index=2 if map_type == "Encroachment probability" else 0,
        )

        # Control how strongly the raster stands out against the basemap
        default_opacity = 0.9 if map_type == "Encroachment probability" else 0.7
        overlay_opacity = st.slider(
            "Overlay opacity",
            min_value=0.0,
            max_value=1.0,
            value=default_opacity,
            step=0.05,
            help="Increase to make the raster stand out more; decrease to see more of the basemap.",
        )

        threshold = None
        if map_type == "Encroachment probability":
            threshold = st.slider(
                "Probability threshold (τ)",
                min_value=0.0,
                max_value=1.0,
                value=0.5,
                step=0.05,
            )

        st.markdown(
            """
**Notes**

- Dominant class: argmax over 5 classes.
- Encroachment prob.: probability of *wetland_to_woody*.
- Uncertainty: `1 - max_class_probability`.
"""
        )

    with col_map:
        # Load mosaic from cache
        try:
            data, bounds, crs, transform, nodata = load_mosaic(raster_path)
        except Exception as e:
            st.error(f"Could not load raster '{raster_path}': {e}")
            st.stop()

        # Select band & colormap based on map type
        band_key = None
        cmap_name = "viridis"
        if map_type == "Dominant class":
            band_key = "dominant_class"
        elif map_type == "Encroachment probability":
            band_key = "p_encroachment"
            cmap_name = "plasma"  # dark purple -> yellow
        elif map_type == "Uncertainty":
            band_key = "uncertainty"
            cmap_name = "Greys"

        band_index = raster_conf["bands"][band_key]  # 1-based
        band_data = data[band_index - 1, :, :]  # 2D

        # Create RGBA overlay
        if map_type == "Dominant class":
            rgba = create_class_rgba_image(band_data, nodata=nodata)
        else:
            # Continuous; apply threshold only for encroachment.
            # Use full alpha here and control visibility via overlay_opacity.
            thr = threshold if map_type == "Encroachment probability" else None
            rgba = create_continuous_rgba_image(
                band_data,
                nodata=nodata,
                cmap_name=cmap_name,
                alpha=255,   # fully opaque in the image
                threshold=thr,
            )

        # Folium map setup
        # Assumption: CRS is EPSG:4326 (lat/lon). If not, reproject raster beforehand.
        center_lat = (bounds.top + bounds.bottom) / 2.0
        center_lon = (bounds.left + bounds.right) / 2.0

        if basemap_choice == "OpenStreetMap":
            m = folium.Map(
                location=[center_lat, center_lon],
                zoom_start=11,
                tiles="OpenStreetMap",
            )
        elif basemap_choice == "ESRI WorldImagery":
            m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles=None)
            folium.TileLayer(
                tiles=(
                    "https://server.arcgisonline.com/ArcGIS/rest/services/"
                    "World_Imagery/MapServer/tile/{z}/{y}/{x}"
                ),
                attr="Esri WorldImagery",
                name="Esri World Imagery",
            ).add_to(m)
        else:  # CartoDB Positron
            m = folium.Map(
                location=[center_lat, center_lon],
                zoom_start=11,
                tiles="CartoDB positron",
            )

        img_bounds = [[bounds.bottom, bounds.left], [bounds.top, bounds.right]]

        ImageOverlay(
            image=rgba,
            bounds=img_bounds,
            opacity=overlay_opacity,
            name=map_type,
            interactive=False,
            cross_origin=False,
        ).add_to(m)

        folium.LayerControl().add_to(m)

        st_folium(m, width="100%", height=600)

        # ---- LEGENDS (this is the bit that was easy to break) ----
        if map_type == "Dominant class":
            show_class_legend()
        else:
            if map_type == "Encroachment probability":
                label = "Probability of wetland_to_woody encroachment"
            else:
                label = "Uncertainty (1 - max class probability)"

            st.markdown("### Legend")
            show_continuous_legend(
                cmap_name=cmap_name,
                vmin=0.0,
                vmax=1.0,
                label=label,
            )

# -----------------------------------------------------------------------------
# TAB 2: REPRESENTATIVE EXAMPLES
# -----------------------------------------------------------------------------
with tab_examples:
    st.subheader("Representative pixel trajectories")

    df_examples = load_representative_examples()
    time_labels = get_time_labels_from_df(df_examples)

    col_controls, col_plots = st.columns([1, 2])

    with col_controls:
        st.markdown("### Filters")

        # True class selection (keep order from CLASS_NAMES, filter to those present)
        classes_in_data = [c for c in CLASS_NAMES if c in df_examples["class_str"].unique()]
        class_choice = st.selectbox(
            "True class (class_str)",
            options=classes_in_data,
        )

        example_types = sorted(df_examples["example_type"].dropna().unique())
        example_type_choice = st.selectbox(
            "Example type",
            options=example_types,
        )

        # Filter dataframe
        df_filtered = df_examples[
            (df_examples["class_str"] == class_choice)
            & (df_examples["example_type"] == example_type_choice)
        ]

        if df_filtered.empty:
            st.warning("No examples found for this combination of class and example type.")
            example_row = None
        else:
            example_ids = sorted(df_filtered["example_id"].astype(str).unique())
            example_id_choice = st.selectbox(
                "Example ID",
                options=example_ids,
            )
            example_row = df_filtered[df_filtered["example_id"] == example_id_choice].iloc[0]

        st.markdown("---")
        st.markdown("### Metadata")

        if example_row is not None:
            st.write(f"**Example ID:** `{example_row['example_id']}`")
            st.write(f"**True class (class_str):** `{example_row['class_str']}`")
            st.write(f"**Predicted class (pred_class_str):** `{example_row['pred_class_str']}`")
            st.write(f"**Example type:** `{example_row['example_type']}`")

            # Optional fields, if present
            for optional_col in ["numark", "pixel_id", "orig_row_idx"]:
                if optional_col in example_row.index:
                    st.write(f"**{optional_col}:** `{example_row[optional_col]}`")

    with col_plots:
        if example_row is None:
            st.info("Select a combination with available examples to see plots.")
        else:
            st.markdown("### Time-series of features")
            plot_example_timeseries(example_row, time_labels)

            st.markdown("### Class probability vector")
            plot_example_probabilities(example_row)

# -----------------------------------------------------------------------------
# TAB 3: CLASS-LEVEL TRAJECTORIES
# -----------------------------------------------------------------------------
with tab_trajectories:
    st.subheader("Class-level trajectories (mean ± std)")

    df_traj = load_class_trajectories()

    col_controls, col_plot = st.columns([1, 2])

    with col_controls:
        st.markdown("### Settings")

        features_in_data = sorted(df_traj["feature"].dropna().unique())
        # Keep canonical ordering where possible
        features_ordered = [f for f in features_in_data if f in features_in_data]
        # In practice FEATURES from config covers them, but we keep it simple here
        feature_choice = st.selectbox(
            "Feature",
            options=features_in_data,
        )

        classes_in_traj = [c for c in CLASS_NAMES if c in df_traj["class_str"].unique()]
        selected_classes = st.multiselect(
            "Classes to show",
            options=classes_in_traj,
            default=classes_in_traj,
        )

        st.markdown(
            """
Each line shows the **mean normalized value** per class over time.
The shaded band is ±1 standard deviation.
"""
        )

    with col_plot:
        if not selected_classes:
            st.info("Select at least one class to display trajectories.")
        else:
            plot_class_trajectories(df_traj, feature_choice, selected_classes)
