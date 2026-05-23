############################################################
# Nuclear GFP intensity quantification from 2-channel ND2 files
#
# Channel default:
#   Channel 1 = GFP
#   Channel 2 = DAPI
#
# Python uses zero-based channel indexing:
#   GFP_CHANNEL = 0 means channel 1
#   DAPI_CHANNEL = 1 means channel 2
############################################################

############################
# USER SETTINGS
############################

SAVE_EXCEL = True
SAVE_PLOTS_PDF = True
SKIP_FILES_ALREADY_ANALYZED = False
REMOVE_FINAL_NUMBER_FOR_GROUPING = True

RUN_STATS = True
CONTROL_MATCH_TEXT = "DMSO"
STATS_TEST = "mannwhitney"  # "mannwhitney" or "ttest"
STATS_LEVEL = "file"        # "file" or "cell"

############################
# CHANNEL SETTINGS
############################

GFP_CHANNEL = 0
DAPI_CHANNEL = 1

############################
# NUCLEUS SEGMENTATION
############################

MIN_NUC_AREA = 1000
NUC_INT_PCT = 1

############################
# NON-NUCLEAR GFP SETTINGS
############################

EXCLUDE_NUCLEAR_BORDER_PIXELS = 3
NON_NUCLEAR_EXCLUDE_HIGH_PCT = 99.5

############################
# IMPORTS
############################

import os
import re
import time
import warnings
import logging

import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tkinter import Tk, filedialog
from scipy.stats import mannwhitneyu, ttest_ind

from skimage import exposure, measure, morphology
from skimage.segmentation import clear_border
from skimage.morphology import disk

from stardist.models import StarDist2D
from csbdeep.utils import normalize
from nd2reader import ND2Reader

############################
# CLEANUP
############################

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.getLogger("tensorflow").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

############################
# SELECT FOLDER
############################

def select_folder():
    root = Tk()
    root.withdraw()
    root.lift()
    root.attributes("-topmost", True)
    root.focus_force()

    folder = filedialog.askdirectory(
        parent=root,
        title="Select folder containing ND2 files"
    )

    root.destroy()

    if folder == "":
        raise ValueError("No folder selected.")

    return folder


print("\nSelect folder containing ND2 files...\n")
input_folder = select_folder()
print(f"Selected folder:\n{input_folder}\n")

############################
# OUTPUT FOLDERS
############################

analysis_folder = os.path.join(input_folder, "Analysis")
overlay_folder = os.path.join(analysis_folder, "Overlays")
plot_folder = os.path.join(analysis_folder, "Plots")

os.makedirs(analysis_folder, exist_ok=True)
os.makedirs(overlay_folder, exist_ok=True)
os.makedirs(plot_folder, exist_ok=True)

############################
# FIND FILES
############################

nd2_files = sorted([
    f for f in os.listdir(input_folder)
    if f.lower().endswith(".nd2")
])

if len(nd2_files) == 0:
    raise ValueError("No ND2 files found.")

print(f"Found {len(nd2_files)} ND2 files.\n")

############################
# LOAD STARDIST
############################

print("Loading StarDist model...")
model = StarDist2D.from_pretrained("2D_versatile_fluo")
print("StarDist loaded.\n")

############################
# STORAGE
############################

summary_records = []

############################
# HELPER FUNCTIONS
############################

def get_group_name(filename):
    base = os.path.splitext(filename)[0]

    if REMOVE_FINAL_NUMBER_FOR_GROUPING:
        base = re.sub(r"([_-])\d+$", "", base)

    return base


def p_to_stars(p):
    if pd.isna(p):
        return "NA"
    if p <= 0.0001:
        return "****"
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return "ns"


def load_nd2_channel_max_projection(path, channel_index):
    max_proj = None

    with ND2Reader(path) as reader:

        axes = reader.sizes

        if "c" in axes:
            reader.iter_axes = "z"
            reader.default_coords["c"] = channel_index
        else:
            reader.iter_axes = "z"

        for frame in reader:
            img = frame.astype(np.float32)

            if max_proj is None:
                max_proj = img
            else:
                max_proj = np.maximum(max_proj, img)

    if max_proj is None:
        raise ValueError(f"Could not read image: {path}")

    return max_proj


def get_stats_input_table(summary_df, per_replicate_summary):
    if STATS_LEVEL == "cell":
        return summary_df.copy()

    if STATS_LEVEL == "file":
        stats_df = per_replicate_summary.copy()

        stats_df = stats_df.rename(
            columns={
                "Mean_Nuclear_GFP": "Nuclear_GFP_Mean",
                "Mean_NonNuclear_GFP": "NonNuclear_GFP_Mean",
                "Mean_Nuclear_to_NonNuclear_GFP_Ratio": "Nuclear_to_NonNuclear_GFP_Ratio"
            }
        )

        return stats_df

    raise ValueError("STATS_LEVEL must be 'file' or 'cell'")


def run_stats_against_control(df, metric):
    control_groups = [
        g for g in df["Group"].unique()
        if CONTROL_MATCH_TEXT.lower() in g.lower()
    ]

    if len(control_groups) == 0:
        raise ValueError(
            f"No control group found containing: {CONTROL_MATCH_TEXT}"
        )

    if len(control_groups) > 1:
        print(
            f"Warning: multiple control groups found: {control_groups}. "
            f"Using first one: {control_groups[0]}"
        )

    control_group = control_groups[0]

    control_values = df.loc[
        df["Group"] == control_group,
        metric
    ].dropna()

    results = []

    for group in sorted(df["Group"].unique()):

        test_values = df.loc[
            df["Group"] == group,
            metric
        ].dropna()

        if group == control_group:
            p_value = np.nan
            statistic = np.nan
            stars = "Control"

        elif len(control_values) < 2 or len(test_values) < 2:
            p_value = np.nan
            statistic = np.nan
            stars = "NA"

        else:
            if STATS_TEST == "mannwhitney":
                stat_result = mannwhitneyu(
                    control_values,
                    test_values,
                    alternative="two-sided"
                )

                statistic = stat_result.statistic
                p_value = stat_result.pvalue

            elif STATS_TEST == "ttest":
                stat_result = ttest_ind(
                    control_values,
                    test_values,
                    equal_var=False
                )

                statistic = stat_result.statistic
                p_value = stat_result.pvalue

            else:
                raise ValueError(
                    "STATS_TEST must be 'mannwhitney' or 'ttest'"
                )

            stars = p_to_stars(p_value)

        results.append({
            "Stats_Level": STATS_LEVEL,
            "Metric": metric,
            "Control_Group": control_group,
            "Compared_Group": group,
            "Test": STATS_TEST,
            "N_Control": len(control_values),
            "N_Group": len(test_values),
            "Control_Mean": control_values.mean(),
            "Group_Mean": test_values.mean(),
            "Statistic": statistic,
            "P_Value": p_value,
            "Stars": stars
        })

    return pd.DataFrame(results)

############################
# PROCESS FILE
############################

def process_file(fname, idx, total):
    print(f"[{idx}/{total}] Processing {fname}")

    start = time.perf_counter()

    path = os.path.join(input_folder, fname)
    base_name = os.path.splitext(fname)[0]
    group_name = get_group_name(fname)

    overlay_path = os.path.join(
        overlay_folder,
        f"{base_name}_nucleus_overlay.tif"
    )

    if SKIP_FILES_ALREADY_ANALYZED and os.path.exists(overlay_path):
        print("  Skipping existing overlay.\n")
        return

    ############################
    # LOAD CHANNELS
    ############################

    gfp = load_nd2_channel_max_projection(
        path,
        GFP_CHANNEL
    )

    dapi = load_nd2_channel_max_projection(
        path,
        DAPI_CHANNEL
    )

    gfp_8 = exposure.rescale_intensity(
        gfp,
        in_range="image",
        out_range=(0, 255)
    ).astype(np.uint8)

    dapi_8 = exposure.rescale_intensity(
        dapi,
        in_range="image",
        out_range=(0, 255)
    ).astype(np.uint8)

    ############################
    # SEGMENT NUCLEI ON DAPI
    ############################

    dapi_norm = normalize(
        dapi,
        1,
        99.8,
        axis=(0, 1)
    )

    labels, _ = model.predict_instances(dapi_norm)

    labels = clear_border(labels)

    labels = morphology.remove_small_objects(
        labels,
        min_size=MIN_NUC_AREA
    )

    nuc_props = measure.regionprops(
        labels,
        intensity_image=dapi_8
    )

    if len(nuc_props) == 0:
        print("  No nuclei detected.\n")
        return

    dapi_means = np.array([
        p.mean_intensity
        for p in nuc_props
    ])

    cutoff = np.percentile(
        dapi_means,
        NUC_INT_PCT
    )

    nuclei = [
        p for p in nuc_props
        if p.mean_intensity >= cutoff
    ]

    print(f"  Nuclei detected: {len(nuclei)}")

    ############################
    # NON-NUCLEAR GFP MASK
    ############################

    all_nuclei_mask = labels > 0

    non_nuclear_mask = ~morphology.binary_dilation(
        all_nuclei_mask,
        disk(EXCLUDE_NUCLEAR_BORDER_PIXELS)
    )

    non_nuclear_pixels = gfp[non_nuclear_mask]

    if len(non_nuclear_pixels) == 0:
        non_nuclear_gfp_mean = np.nan
    else:
        high_cutoff = np.percentile(
            non_nuclear_pixels,
            NON_NUCLEAR_EXCLUDE_HIGH_PCT
        )

        non_nuclear_pixels_filtered = non_nuclear_pixels[
            non_nuclear_pixels <= high_cutoff
        ]

        non_nuclear_gfp_mean = np.mean(
            non_nuclear_pixels_filtered
        )

    ############################
    # OVERLAY IMAGE
    ############################

    overlay = np.zeros(
        (*gfp_8.shape, 3),
        dtype=np.uint8
    )

    overlay[..., 1] = gfp_8
    overlay[..., 2] = dapi_8

    ############################
    # MEASURE EACH NUCLEUS
    ############################

    for p in nuclei:

        nucleus_mask = labels == p.label

        nuclear_gfp_pixels = gfp[nucleus_mask]

        nuclear_gfp_mean = np.mean(nuclear_gfp_pixels)

        nuclear_gfp_median = np.median(nuclear_gfp_pixels)

        nuclear_gfp_integrated_density = np.sum(
            nuclear_gfp_pixels
        )

        nuclear_area = np.sum(nucleus_mask)

        if (
            pd.isna(non_nuclear_gfp_mean)
            or non_nuclear_gfp_mean == 0
        ):
            nuclear_to_non_nuclear_ratio = np.nan
        else:
            nuclear_to_non_nuclear_ratio = (
                nuclear_gfp_mean / non_nuclear_gfp_mean
            )

        summary_records.append({
            "File": fname,
            "Group": group_name,
            "Cell_ID": p.label,
            "Nucleus_Area": float(nuclear_area),
            "DAPI_Mean_Intensity": float(p.mean_intensity),
            "Nuclear_GFP_Mean": float(nuclear_gfp_mean),
            "Nuclear_GFP_Median": float(nuclear_gfp_median),
            "Nuclear_GFP_Integrated_Density": float(nuclear_gfp_integrated_density),
            "NonNuclear_GFP_Mean": float(non_nuclear_gfp_mean),
            "Nuclear_to_NonNuclear_GFP_Ratio": float(nuclear_to_non_nuclear_ratio)
        })

        border = (
            morphology.binary_dilation(
                nucleus_mask,
                disk(1)
            )
            ^ nucleus_mask
        )

        overlay[border] = [255, 0, 0]

    overlay_bgr = cv2.cvtColor(
        overlay,
        cv2.COLOR_RGB2BGR
    )

    cv2.imwrite(
        overlay_path,
        overlay_bgr
    )

    elapsed = time.perf_counter() - start

    print(f"  Done in {elapsed:.2f}s\n")

############################
# RUN ALL FILES
############################

for i, fname in enumerate(nd2_files, start=1):
    process_file(fname, i, len(nd2_files))

############################
# SAVE TABLES
############################

summary_df = pd.DataFrame(summary_records)

if len(summary_df) == 0:
    raise ValueError("No nuclei were quantified.")

per_cell_path = os.path.join(
    analysis_folder,
    "nuclear_gfp_per_cell_summary.csv"
)

summary_df.to_csv(
    per_cell_path,
    index=False
)

per_replicate_summary = (
    summary_df
    .groupby(["Group", "File"], as_index=False)
    .agg(
        N_Cells=("Cell_ID", "count"),
        Mean_Nuclear_GFP=("Nuclear_GFP_Mean", "mean"),
        SEM_Nuclear_GFP=("Nuclear_GFP_Mean", "sem"),
        Mean_NonNuclear_GFP=("NonNuclear_GFP_Mean", "mean"),
        SEM_NonNuclear_GFP=("NonNuclear_GFP_Mean", "sem"),
        Mean_Nuclear_to_NonNuclear_GFP_Ratio=("Nuclear_to_NonNuclear_GFP_Ratio", "mean"),
        SEM_Nuclear_to_NonNuclear_GFP_Ratio=("Nuclear_to_NonNuclear_GFP_Ratio", "sem"),
        Mean_Nuclear_GFP_Integrated_Density=("Nuclear_GFP_Integrated_Density", "mean"),
        SEM_Nuclear_GFP_Integrated_Density=("Nuclear_GFP_Integrated_Density", "sem"),
        Mean_Nucleus_Area=("Nucleus_Area", "mean"),
        SEM_Nucleus_Area=("Nucleus_Area", "sem")
    )
)

per_replicate_path = os.path.join(
    analysis_folder,
    "nuclear_gfp_per_replicate_summary.csv"
)

per_replicate_summary.to_csv(
    per_replicate_path,
    index=False
)

group_summary = (
    per_replicate_summary
    .groupby("Group", as_index=False)
    .agg(
        N_Files=("File", "nunique"),
        Total_Cells=("N_Cells", "sum"),
        Mean_Nuclear_GFP=("Mean_Nuclear_GFP", "mean"),
        SEM_Nuclear_GFP=("Mean_Nuclear_GFP", "sem"),
        Mean_NonNuclear_GFP=("Mean_NonNuclear_GFP", "mean"),
        SEM_NonNuclear_GFP=("Mean_NonNuclear_GFP", "sem"),
        Mean_Nuclear_to_NonNuclear_GFP_Ratio=("Mean_Nuclear_to_NonNuclear_GFP_Ratio", "mean"),
        SEM_Nuclear_to_NonNuclear_GFP_Ratio=("Mean_Nuclear_to_NonNuclear_GFP_Ratio", "sem"),
        Mean_Nuclear_GFP_Integrated_Density=("Mean_Nuclear_GFP_Integrated_Density", "mean"),
        SEM_Nuclear_GFP_Integrated_Density=("Mean_Nuclear_GFP_Integrated_Density", "sem"),
        Mean_Nucleus_Area=("Mean_Nucleus_Area", "mean"),
        SEM_Nucleus_Area=("Mean_Nucleus_Area", "sem")
    )
)

group_path = os.path.join(
    analysis_folder,
    "nuclear_gfp_group_summary.csv"
)

group_summary.to_csv(
    group_path,
    index=False
)

############################
# STATISTICS
############################

stats_df = pd.DataFrame()

if RUN_STATS:

    stats_input_df = get_stats_input_table(
        summary_df,
        per_replicate_summary
    )

    stats_metrics = [
        "Nuclear_GFP_Mean",
        "NonNuclear_GFP_Mean",
        "Nuclear_to_NonNuclear_GFP_Ratio"
    ]

    stats_tables = []

    for metric in stats_metrics:
        stats_tables.append(
            run_stats_against_control(
                stats_input_df,
                metric
            )
        )

    stats_df = pd.concat(
        stats_tables,
        ignore_index=True
    )

    print("\n==============================")
    print("STATISTICS AGAINST CONTROL")
    print("==============================")
    print(f"Control match text: {CONTROL_MATCH_TEXT}")
    print(f"Test used: {STATS_TEST}")
    print(f"Stats level: {STATS_LEVEL}\n")
    print(stats_df.to_string(index=False))

############################
# EXCEL OUTPUT
############################

excel_path = os.path.join(
    analysis_folder,
    "nuclear_gfp_analysis_with_statistics.xlsx"
)

if SAVE_EXCEL:

    with pd.ExcelWriter(
        excel_path,
        engine="openpyxl"
    ) as writer:

        summary_df.to_excel(
            writer,
            sheet_name="Per_Cell",
            index=False
        )

        per_replicate_summary.to_excel(
            writer,
            sheet_name="Per_Replicate",
            index=False
        )

        group_summary.to_excel(
            writer,
            sheet_name="Grouped",
            index=False
        )

        if RUN_STATS:
            stats_input_df.to_excel(
                writer,
                sheet_name="Stats_Input",
                index=False
            )

            stats_df.to_excel(
                writer,
                sheet_name="Stats_vs_Control",
                index=False
            )

    print(f"\nExcel file saved:\n{excel_path}")

############################
# PLOTS
############################

def get_stars_for_metric(metric, group):
    if not RUN_STATS or len(stats_df) == 0:
        return ""

    match = stats_df.loc[
        (stats_df["Metric"] == metric)
        & (stats_df["Compared_Group"] == group),
        "Stars"
    ]

    if len(match) == 0:
        return ""

    stars = str(match.iloc[0])

    if stars in ["Control", "NA"]:
        return ""

    return stars


def grouped_plot(value_col, ylabel, filename):

    metric_column = value_col

    file_means = (
        per_replicate_summary
        .rename(
            columns={
                "Mean_Nuclear_GFP": "Nuclear_GFP_Mean",
                "Mean_NonNuclear_GFP": "NonNuclear_GFP_Mean",
                "Mean_Nuclear_to_NonNuclear_GFP_Ratio": "Nuclear_to_NonNuclear_GFP_Ratio"
            }
        )
    )

    group_means = (
        file_means
        .groupby("Group", as_index=False)
        .agg(
            Mean=(metric_column, "mean"),
            SEM=(metric_column, "sem")
        )
    )

    groups = list(group_means["Group"])
    x = np.arange(len(groups))

    fig, ax = plt.subplots(
        figsize=(max(6, len(groups) * 1.2), 5)
    )

    ax.bar(
        x,
        group_means["Mean"],
        yerr=group_means["SEM"],
        capsize=5,
        edgecolor="black"
    )

    ymax = 0

    for i, group in enumerate(groups):

        vals = file_means.loc[
            file_means["Group"] == group,
            metric_column
        ].values

        if len(vals) > 0:
            ymax = max(ymax, np.nanmax(vals))

        jitter = np.random.normal(
            0,
            0.05,
            len(vals)
        )

        ax.scatter(
            np.repeat(i, len(vals)) + jitter,
            vals,
            s=40,
            edgecolor="black",
            zorder=3
        )

    y_range = max(
        ymax,
        group_means["Mean"].max()
    ) * 0.15

    if y_range == 0:
        y_range = 1

    for i, group in enumerate(groups):

        stars = get_stars_for_metric(
            value_col,
            group
        )

        if stars == "":
            continue

        mean_value = group_means.loc[
            group_means["Group"] == group,
            "Mean"
        ].iloc[0]

        sem_value = group_means.loc[
            group_means["Group"] == group,
            "SEM"
        ].iloc[0]

        if pd.isna(sem_value):
            sem_value = 0

        y_pos = mean_value + sem_value + y_range * 0.25

        ax.text(
            i,
            y_pos,
            stars,
            ha="center",
            va="bottom",
            fontsize=14
        )

    ax.set_ylim(
        0,
        max(
            ymax,
            group_means["Mean"].max()
        ) + y_range
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        groups,
        rotation=45,
        ha="right"
    )

    ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    outpath = os.path.join(
        plot_folder,
        filename
    )

    plt.savefig(outpath)
    plt.close()

    print(f"Plot saved:\n{outpath}")


if SAVE_PLOTS_PDF:

    grouped_plot(
        "Nuclear_GFP_Mean",
        "Nuclear GFP mean intensity",
        "nuclear_gfp_mean_intensity.pdf"
    )

    grouped_plot(
        "Nuclear_to_NonNuclear_GFP_Ratio",
        "Nuclear / non-nuclear GFP ratio",
        "nuclear_to_non_nuclear_gfp_ratio.pdf"
    )

############################
# FINAL
############################

print("\n✅ All done!")
print(f"\nPer-cell CSV:\n{per_cell_path}")
print(f"\nPer-replicate CSV:\n{per_replicate_path}")
print(f"\nGrouped CSV:\n{group_path}")

if SAVE_EXCEL:
    print(f"\nExcel with statistics:\n{excel_path}")

print(f"\nResults saved in:\n{analysis_folder}")