# SNAP

SNAP (Spot and Nuclei Analysis Pipeline) is a MATLAB toolkit for quantitative analysis of fluorescence microscopy images. It is designed for robust, interactive spot localization with optional nuclei segmentation and per-nucleus spot assignment.

SNAP is built for practical experimental workflows:
- interactive parameter tuning with deterministic preview updates
- high-throughput batch processing
- classifier training with channel-specific feature sets
- modular extension points for community contributions
- microscopy-appropriate localization methods (Gaussian-family fits, radial symmetry, fit-quality filtering)

## Table of Contents

1. [What Is Included](#what-is-included)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Application Guide](#application-guide)
5. [Input Conventions](#input-conventions)
6. [Outputs](#outputs)
7. [Extension and Modularity](#extension-and-modularity)
8. [Contributing](#contributing)
9. [Troubleshooting](#troubleshooting)

## What Is Included

| Application | Purpose | Typical Use |
|---|---|---|
| `SNAP` | Main interactive analysis UI | Tune preprocessing, maxima detection, fitting, filtering, classification, and export |
| `SNAP_prepare` | Convert multi-channel microscopy libraries into SNAP-ready folder structure | Build clean datasets for batch runs |
| `SNAP_batch` | Run parameterized analysis across many samples | Large studies, reproducible processing |
| `SNAP_classify` | Manual spot labeling and classifier training from exported candidates | Curate high-quality training labels visually |
| `SNAP_train` | Channel-wise SVM training from image+CSV labels or scripted inputs | Scalable, per-channel classifier training |

Supporting components:
- `+snap_helpers/`: core algorithms and UI helpers
- `+snap_modules/`: shared modular runtime engines
- `external_plugins/`: optional external plugin entry points
- `+snap_contrib/` and `examples/`: non-core contribution examples

## Installation

### Requirements

- MATLAB R2020a or newer (recommended)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox (for SVM workflows)
- Bio-Formats Toolbox (optional, for `SNAP_prepare`)

Bio-Formats installation (if needed):
- MATLAB Home -> Add-Ons -> Get Add-Ons -> search `Bio-Formats`

### Download

Option 1 (recommended):
```bash
git clone https://github.com/marcorojas-cessa/SNAP.git
```

Option 2:
- Download ZIP from GitHub
- Extract to a local folder

### Add SNAP to MATLAB path

In MATLAB:
```matlab
addpath('/absolute/path/to/SNAP');
savepath; % optional
```

`/absolute/path/to/...` means the full path from your filesystem root, for example:
- macOS/Linux: `/Users/you/projects/SNAP`
- Windows: `C:\Users\you\projects\SNAP`

### Launch check

```matlab
SNAP
```

If the UI opens without errors, installation is complete.

## Quick Start

### Workflow A: Interactive analysis on one sample

1. Run `SNAP`.
2. Load nuclei and signal channel images.
3. Configure per-channel processing and fitting settings.
4. Click `Update Previews` to process and inspect results.
5. Use `Abort Processing` if a run is too heavy.
6. Export data and parameter files.

### Workflow B: Prepare libraries and batch process

1. Run `SNAP_prepare` and read your microscopy library.
2. Map channel identities (DIC, nuclei, fluorescence channels).
3. Export selected images into SNAP folder structure.
4. Run `SNAP_batch` with your parameter file.

Example CLI:
```matlab
SNAP_batch('/path/to/prepared_dataset', '/path/to/parameters.mat');
```

### Workflow C: Train and apply classifiers

Two training routes:
- `SNAP_classify`: visual/manual labeling workflow
- `SNAP_train`: per-channel training from directories and labels

After training, load classifiers in `SNAP` or `SNAP_batch` for automated spot/noise classification.

## Application Guide

### `SNAP` (main UI)

`SNAP` is the full interactive environment for analysis and parameter tuning.

Key behavior:
- Preview updates are manual (`Update Previews`), not automatic.
- Long operations are abort-aware (`Abort Processing`).
- Channels are processed through shared stage logic for consistent behavior.
- Main `SNAP` UI supports up to five fluorescence channels, plus nuclei and DIC inputs.

Main stage groups available in the UI:
- signal preprocessing (deconvolution, denoising, background correction)
- maxima detection (regional, extended/H-maxima, LoG)
- localization fitting (1D, 2D+1D, 3D, distorted 3D, radial symmetry)
- fit filtering (quality and parameter bounds)
- optional classification and nuclei-based filtering

### `SNAP_prepare`

Use `SNAP_prepare` when source data comes from complex multi-channel formats (for example CZI, ND2, LIF).

What it does:
- reads supported formats (Bio-Formats)
- shows image inventory
- maps channel identities
- exports SNAP-compatible per-sample folders

### `SNAP_batch`

`SNAP_batch` runs reproducible analysis on many samples using one parameter file.

Modes:
- GUI mode: `SNAP_batch`
- command-line mode: `SNAP_batch(inputDir, paramFile, ...)`

Example:
```matlab
SNAP_batch('/path/to/data', '/path/to/params.mat', 'OutputDir', '/path/to/results');
```

### `SNAP_classify`

`SNAP_classify` is for manual curation and training from exported candidates.

Typical flow:
1. Load SNAP export MAT data
2. Select classification features
3. Label spots as real/noise
4. Train and export classifier

Useful shortcuts include real/noise labeling and spot navigation directly in the UI.

### `SNAP_train`

`SNAP_train` is the scalable training interface for per-channel classifiers.

GUI mode:
```matlab
SNAP_train
```

Capabilities:
- reads parameter file to infer active channels and fitting context
- trains one classifier per selected channel
- per-channel training and validation directories
- per-channel match distance
- optional FIJI coordinate conversion per channel
- manual hyperparameters or validation sweep optimization
- channel-specific base features and custom expressions
- `Load Expression Pack...` support for applying saved feature sets

Programmatic mode is also supported:
```matlab
SNAP_train(exportFiles, labelFiles, '/path/to/classifier.mat', 'MatchDistance', 2);
```

## Input Conventions

### Batch folder layout

`SNAP_batch` expects sample subfolders, for example:

```text
input_directory/
  Sample_001/
    nuclei.tif
    channel1.tif
    channel2.tif
    dic.tif
  Sample_002/
    ...
```

`SNAP_prepare` can create this structure automatically.

### `SNAP_train` image/label pairing

For directory-based discovery, `SNAP_train` matches by same base filename:
- image: `.tif` or `.tiff`
- label: `.csv`

Example pair:
- `cell_01_ch1.tif`
- `cell_01_ch1.csv`

### CSV label columns

Accepted coordinate headers are case-insensitive and may use:
- `x`, `y`, optional `z` or `slice`
- or maxima/fitted variants (`maxima_x`, `fitted_y`, etc.)

Label columns can be `label`, `class`, or `is_real`.

### FIJI coordinate conversion

If labels were generated in FIJI/ImageJ indexing, enable `Convert FIJI Coords` for that channel in `SNAP_train`.

Applied conversion:
- swap `x/y` to MATLAB row/column convention
- add `+1` index offset

## Outputs

SNAP exports analysis-ready data in CSV and MAT formats.

Common outputs include:
- nuclei measurements
- per-channel signal localization measurements
- nuclei-signal composition summaries
- parameter/configuration MAT files
- classifier MAT files

Signal MAT exports include fit-window context fields used for advanced downstream analysis, including:
- `rawDataWindow`
- `fitWindowDimensions`
- `fitWindowOrigin`
- `localMaximaInWindow`
- `originalMaximaCoords`

## Extension and Modularity

SNAP is compartmentalized so core behavior and contributed methods can evolve without becoming tightly coupled.

### Core runtime paths

- signal engine: `+snap_modules/+signal/`
- prepare engine: `+snap_modules/+prepare/`
- built-in plugins: `+snap_modules/+plugins/`

### External extension paths

- signal plugins: `external_plugins/signal/`
- prepare reader plugins: `external_plugins/prepare/readers/`
- prepare exporter plugins: `external_plugins/prepare/exporters/`

This naming separation is intentional:
- built-ins live under `+snap_modules/+plugins/...`
- community or lab-specific add-ons live under `external_plugins/...`

### SVM expression packs

SNAP includes an example contribution for expression-pack generation:
- `+snap_contrib/+svm/`
- `examples/svm_feature_pack/`

Build and save a pack:
```matlab
pack = snap_contrib.svm.buildExpressionPack('ParameterFile', '/absolute/path/to/params.mat');
snap_contrib.svm.saveExpressionPack(pack, '/absolute/path/to/snap_svm_expression_pack.mat');
```

Then in `SNAP_train`, use `Load Expression Pack...`.

Compatibility guardrails in `SNAP_train`:
- pack contents are validated against channel context
- incompatible features/expressions are auto-pruned with log messages
- all-NaN features are auto-dropped during feature-matrix build
- if needed, training falls back to automatic non-position base features

## Contributing

Contributions are welcome for both core improvements and optional modules.

Recommended process:
1. Open an issue (`bug report`, `feature request`, or `contribution proposal`).
2. Fork and create a focused branch.
3. Keep contributions modular:
   - core changes in core paths only when necessary
   - optional methods in `external_plugins/`, `+snap_contrib/`, and `examples/`
4. Include documentation and usage examples with your PR.
5. Submit a PR using the repository template.

Repository governance aids are already included:
- CODEOWNERS
- PR template
- issue templates
- CI and contribution guard workflows

## Troubleshooting

### `SNAP_prepare` cannot read my file
- Install Bio-Formats Toolbox in MATLAB Add-Ons.

### Training fails with zero valid samples
- Review selected features/custom expressions.
- In `SNAP_train`, check log lines for auto-pruned incompatible features.
- Verify label CSV columns and coordinate convention.

### Batch run finds no channels or nuclei
- Confirm folder naming and channel mapping.
- Verify parameter file and channel count are aligned with dataset.

### Processing is slow on large stacks
- Use preview tuning on representative images first.
- Reduce expensive options while searching for stable parameters.
- Use batch mode after parameters are fixed.
