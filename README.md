# Nuclear GFP Quantification from ND2 Files

This repository contains a Python script to quantify nuclear GFP intensity from multi-channel `.nd2` microscopy images.

The script segments nuclei using the DAPI channel and measures GFP intensity inside each nucleus.

It was designed for fast, reproducible analysis of fluorescence microscopy datasets with multiple replicates and experimental conditions.

---

# What the script does

For each `.nd2` image, the script:

1. Opens a folder selection window for the user to choose the input folder.
2. Creates an `Analysis` folder automatically.
3. Loads GFP and DAPI channels from `.nd2` files.
4. Creates maximum-intensity projections across z-slices.
5. Segments nuclei using StarDist on the DAPI channel.
6. Measures GFP intensity inside each nucleus.
7. Measures non-nuclear GFP intensity outside nuclei.
8. Calculates nuclear enrichment ratios.
9. Groups replicate files automatically based on file names.
10. Performs optional statistics against a control condition.
11. Saves overlay images, tables, Excel files, and PDF plots.

---

# Channel setup

Default channel configuration:

```python
GFP_CHANNEL = 0
DAPI_CHANNEL = 1
```

Python uses zero-based indexing:

| Biological channel | Python index |
|---|---|
| Channel 1 | 0 |
| Channel 2 | 1 |
| Channel 3 | 2 |

If your acquisition order is different, simply change these values at the top of the script.

---

# Input

Place all `.nd2` files inside a single folder.

Example:

```text
Experiment/
├── DMSO_001.nd2
├── DMSO_002.nd2
├── DrugA_001.nd2
├── DrugA_002.nd2
├── DrugB_001.nd2
└── DrugB_002.nd2
```

When the script starts, a Windows folder selection window will open.

Select the folder containing the `.nd2` files.

---

# Automatic grouping

Replicates are grouped automatically from file names.

Example:

```text
DMSO_001.nd2
DMSO_002.nd2
DMSO_003.nd2
```

becomes:

```text
DMSO
```

This allows automatic grouping of biological or technical replicates.

---

# Output

The script creates:

```text
Experiment/
└── Analysis/
    ├── Overlays/
    ├── Plots/
    ├── nuclear_gfp_per_cell_summary.csv
    ├── nuclear_gfp_per_replicate_summary.csv
    ├── nuclear_gfp_group_summary.csv
    └── nuclear_gfp_analysis_with_statistics.xlsx
```

---

# Overlay images

Overlay images are saved in:

```text
Analysis/Overlays/
```

The overlays show:

- GFP signal
- DAPI signal
- Segmented nuclei outlines

This allows fast visual quality control of segmentation.

---

# Output tables

## nuclear_gfp_per_cell_summary.csv

One row per detected nucleus.

Main columns:

| Column | Description |
|---|---|
| File | Original ND2 file |
| Group | Experimental condition |
| Cell_ID | Nucleus ID |
| Nucleus_Area | Nuclear area in pixels |
| DAPI_Mean_Intensity | Mean DAPI intensity |
| Nuclear_GFP_Mean | Mean GFP intensity inside nucleus |
| Nuclear_GFP_Median | Median GFP intensity inside nucleus |
| Nuclear_GFP_Integrated_Density | Sum of GFP signal inside nucleus |
| NonNuclear_GFP_Mean | Mean GFP intensity outside nuclei |
| Nuclear_to_NonNuclear_GFP_Ratio | Nuclear enrichment ratio |

---

## nuclear_gfp_per_replicate_summary.csv

One row per ND2 file.

This is useful when each file is considered one biological replicate.

Statistics are performed using this table by default.

---

## nuclear_gfp_group_summary.csv

One row per experimental group.

Contains grouped averages and SEM values across replicates.

---

# Metrics

## Nuclear GFP mean intensity

```text
Mean GFP intensity inside the segmented nucleus
```

Useful for measuring nuclear accumulation.

---

## Non-nuclear GFP mean intensity

```text
Mean GFP intensity outside nuclei
```

This acts as a non-nuclear/background reference.

The script excludes nuclear borders and extreme bright outlier pixels to avoid artifacts.

---

## Nuclear-to-non-nuclear GFP ratio

```text
Nuclear GFP mean / Non-nuclear GFP mean
```

This measures nuclear enrichment of GFP signal.

Values:

| Ratio | Interpretation |
|---|---|
| >1 | GFP enriched in nucleus |
| ~1 | Similar nuclear and non-nuclear signal |
| <1 | GFP depleted from nucleus |

---

## Nuclear GFP integrated density

```text
Sum of all GFP pixel intensities inside the nucleus
```

This reflects total nuclear GFP signal.

Unlike mean intensity, this metric is affected by nuclear size.

---

# Statistics

Statistics are optional and controlled at the top of the script.

Current default:

```python
RUN_STATS = True
CONTROL_MATCH_TEXT = "DMSO"
STATS_TEST = "mannwhitney"
STATS_LEVEL = "file"
```

This means:

- Any group containing `"DMSO"` is used as control.
- Statistics are performed using one value per ND2 file.
- Mann–Whitney U test is used.

---

# Statistics level

## File-level statistics (recommended)

```python
STATS_LEVEL = "file"
```

Each ND2 file is treated as one replicate.

This is usually biologically correct because cells within the same image are not fully independent measurements.

---

## Cell-level statistics

```python
STATS_LEVEL = "cell"
```

Each nucleus is treated as one datapoint.

This increases sample size but may inflate significance if many cells come from the same image.

---

# Significance stars

Plots and Excel tables include significance stars:

```text
ns      p > 0.05
*       p ≤ 0.05
**      p ≤ 0.01
***     p ≤ 0.001
****    p ≤ 0.0001
```

---

# PDF plots

If enabled:

```python
SAVE_PLOTS_PDF = True
```

PDF plots are saved in:

```text
Analysis/Plots/
```

Plots include:

- Group means
- SEM error bars
- Individual replicate points
- Statistical significance stars

---

# Main user settings

General settings:

```python
SAVE_EXCEL = True
SAVE_PLOTS_PDF = True
SKIP_FILES_ALREADY_ANALYZED = False
REMOVE_FINAL_NUMBER_FOR_GROUPING = True

RUN_STATS = True
CONTROL_MATCH_TEXT = "DMSO"

STATS_TEST = "mannwhitney"
STATS_LEVEL = "file"
```

Channel settings:

```python
GFP_CHANNEL = 0
DAPI_CHANNEL = 1
```

Nucleus segmentation settings:

```python
MIN_NUC_AREA = 1000
NUC_INT_PCT = 1
```

Non-nuclear GFP settings:

```python
EXCLUDE_NUCLEAR_BORDER_PIXELS = 3
NON_NUCLEAR_EXCLUDE_HIGH_PCT = 99.5
```

---

# Installation

Create a Python environment and install required packages:

```bash
pip install numpy pandas matplotlib scipy scikit-image nd2reader stardist csbdeep tensorflow openpyxl opencv-python
```

Depending on the operating system, TensorFlow and StarDist may require compatible Python versions.

---

# How to run

Run the script from Python or PyCharm.

A folder selection window will open.

Select the folder containing the `.nd2` files.

Results will automatically be saved in:

```text
SelectedFolder/Analysis/
```

---

# Notes

- The script uses maximum-intensity projection across z-slices.
- All measurements are in pixel units unless externally calibrated.
- The script assumes that DAPI labels nuclei clearly.
- StarDist is used for nucleus segmentation.
- Detection settings may require adjustment depending on microscope settings and image quality.
- The script was designed for fluorescence microscopy datasets with multiple replicates and experimental conditions.

---
