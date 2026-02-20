# SNAP - Spot & Nuclei Analysis Pipeline

<p align="center">
  <strong>A comprehensive MATLAB toolkit for fluorescent signal localization and nuclei segmentation in microscopy images</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#workflow">Workflow</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#citation">Citation</a>
</p>

---

## Overview

**SNAP** (Spot & Nuclei Analysis Pipeline) is a powerful, interactive MATLAB application designed for quantitative analysis of fluorescence microscopy images. It provides a complete end-to-end workflow for:

- **Nuclei segmentation** using advanced adaptive thresholding algorithms
- **Sub-pixel localization** of fluorescent signals via Gaussian fitting
- **Machine learning classification** to distinguish real spots from noise
- **High-throughput batch processing** for large datasets

Current UI behavior is intentionally manual and deterministic:
- Previews are updated only when you click **Update Previews**
- Processing can be interrupted with **Abort Processing**, including long maxima fitting runs

SNAP is particularly suited for applications in:
- Single-molecule localization microscopy
- FISH (Fluorescence In Situ Hybridization) analysis
- Protein localization studies
- Cell cycle analysis
- Any application requiring precise spot detection within cellular compartments

---

## Features

### 🔬 **Multi-Channel Analysis**
- Support for up to 5 fluorescence channels plus nuclei and DIC
- Independent parameter tuning per channel
- Cross-channel signal composition analysis

### 🎯 **Advanced Spot Detection**
- **Simple Regional Maxima** - Fast neighborhood-based detection
- **Extended Maxima (H-maxima)** - Noise-robust peak finding
- **Laplacian of Gaussian (LoG)** - Scale-space blob detection with anisotropic 3D support

### 📐 **Sub-Pixel Gaussian Fitting**
- **1D (X,Y,Z)** - Independent axis fitting for speed
- **2D (XY) + 1D (Z)** - Hybrid approach for anisotropic PSFs
- **3D Gaussian** - Full volumetric fitting
- **Distorted 3D Gaussian** - With correlation terms (ρ_xy, ρ_xz, ρ_yz)
- **Radial Symmetry** - Ultra-fast centroid refinement

### 🧬 **Intelligent Nuclei Segmentation**
- Multiple thresholding strategies: Absolute, Statistical (Mean/Median ± offset)
- **Auto Local Threshold** algorithms matching ImageJ implementation:
  - Bernsen, Mean, Median, MidGrey
  - Niblack, Phansalkar, Sauvola, Otsu
- Morphological filtering: Size, Circularity, Solidity
- Edge nuclei exclusion

### 🤖 **Machine Learning Classification**
- SVM-based spot classification (Real vs. Noise)
- Interactive training interface with visual spot labeling
- Custom feature expressions for advanced classification
- Z-score normalization with stored parameters for reproducibility

### 📊 **Comprehensive Export**
- Nuclei morphological measurements (CSV, MAT)
- Signal localization data with fit parameters
- Nuclei-signal composition analysis
- Annotated visualizations

### ⚡ **Performance Optimized**
- Efficient batch processing with progress tracking
- Preview/runtime caching for responsive interaction
- Abort-aware maxima fitting for safer long runs

---

## Installation

### Requirements

- **MATLAB R2020a or later** (App Designer UI components)
- **Image Processing Toolbox**
- **Statistics and Machine Learning Toolbox** (for SVM classification)
- **Bio-Formats Toolbox** (optional, for SNAP_prepare - install from MATLAB Add-Ons)

### Setup

1. Clone or download this repository:
   ```bash
   git clone https://github.com/yourusername/SNAP.git
   ```

2. Add SNAP to your MATLAB path:
   ```matlab
   addpath('/path/to/SNAP');
   ```

3. Launch SNAP:
   ```matlab
   SNAP
   ```

---

## Quick Start

### Interactive GUI Mode

```matlab
% Launch the main SNAP interface
SNAP

% Steps:
% 1. Load nuclei image (Browse → Nuclei)
% 2. Load fluorescence channel(s) (Browse → Channel 1, 2, ...)
% 3. Configure processing parameters in each tab
% 4. Click "Update Previews" to run processing and refresh all previews
% 5. Use "Abort Processing" any time to stop long runs
% 6. Export data when satisfied
```

### Batch Processing

```matlab
% GUI mode - configure and run batch jobs
SNAP_batch

% Command-line mode
SNAP_batch('path/to/data/', 'parameters.mat')
SNAP_batch('path/to/data/', 'parameters.mat', 'OutputDir', 'results/')
```

### Prepare Multi-Channel Libraries

```matlab
% Convert proprietary formats (CZI, ND2, LIF, etc.) to SNAP-compatible folders
SNAP_prepare
```

### Train a Classifier

```matlab
% Open the classification training interface
SNAP_classify
```

---

## Important Behavior

- **No auto preview updates**: SNAP does not auto-refresh previews while parameters change.
- **Abort is active during fitting**: abort requests are checked during chunked fitting and inside per-maxima fitting loops.
- **No autosave on close**: export **Parameters** to save a reproducible configuration.

---

## Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           SNAP PROCESSING PIPELINE                          │
└─────────────────────────────────────────────────────────────────────────────┘

     ┌──────────────┐
     │ IMAGE INPUT  │
     │ ─────────────│
     │ • Nuclei     │
     │ • DIC        │
     │ • Channels   │
     └──────┬───────┘
            │
            ▼
┌───────────────────────────────────────────────────────────────────────────┐
│                              PREPROCESSING                                 │
├───────────────────────────────────┬───────────────────────────────────────┤
│         NUCLEI CHANNEL            │         SIGNAL CHANNELS               │
│  ┌─────────────────────────────┐  │  ┌─────────────────────────────────┐  │
│  │ 1. Deconvolution (optional) │  │  │ 1. Deconvolution (optional)     │  │
│  │    • Lucy-Richardson        │  │  │    • Lucy-Richardson            │  │
│  │    • Wiener                 │  │  │    • Wiener                     │  │
│  │    • Blind                  │  │  │    • Blind                      │  │
│  ├─────────────────────────────┤  │  ├─────────────────────────────────┤  │
│  │ 2. Smoothing/Denoising      │  │  │ 2. Smoothing/Denoising          │  │
│  │    • Gaussian               │  │  │    • Gaussian                   │  │
│  │    • Median                 │  │  │    • Median                     │  │
│  │    • Non-Local Means        │  │  │    • Non-Local Means            │  │
│  │    • Wavelet                │  │  │    • Wavelet                    │  │
│  ├─────────────────────────────┤  │  ├─────────────────────────────────┤  │
│  │ 3. Background Correction    │  │  │ 3. Background Correction        │  │
│  │    • Rolling Ball           │  │  │    • Rolling Ball               │  │
│  │    • Top-Hat                │  │  │    • Top-Hat                    │  │
│  │    • Gaussian               │  │  │    • Gaussian                   │  │
│  └─────────────────────────────┘  │  └─────────────────────────────────┘  │
└───────────────────────────────────┴───────────────────────────────────────┘
            │                                       │
            ▼                                       ▼
┌───────────────────────────┐       ┌───────────────────────────────────────┐
│   NUCLEI SEGMENTATION     │       │          SIGNAL DETECTION             │
├───────────────────────────┤       ├───────────────────────────────────────┤
│ • Absolute threshold      │       │ LOCAL MAXIMA DETECTION                │
│ • Mean/Median ± offset    │       │ • Simple Regional                     │
│ • Auto Local Threshold:   │       │ • Extended Maxima (H-max)             │
│   - Bernsen               │       │ • Laplacian of Gaussian               │
│   - Niblack               │       ├───────────────────────────────────────┤
│   - Phansalkar            │       │ GAUSSIAN FITTING                      │
│   - Sauvola               │       │ • 1D (X,Y,Z)                          │
│   - Otsu                  │       │ • 2D (XY) + 1D (Z)                    │
│ ──────────────────────────│       │ • 3D Gaussian                         │
│ FILTERING                 │       │ • Distorted 3D (with correlations)    │
│ • Size (pixels/µm²/µm³)   │       │ • Radial Symmetry                     │
│ • Circularity/Sphericity  │       ├───────────────────────────────────────┤
│ • Solidity                │       │ FIT FILTERING                         │
│ • Edge exclusion          │       │ • R² threshold                        │
└───────────────────────────┘       │ • Sigma sum range                     │
            │                       │ • Amplitude range                     │
            │                       │ • Intensity range                     │
            │                       ├───────────────────────────────────────┤
            │                       │ CLASSIFICATION (optional)             │
            │                       │ • SVM: Real vs. Noise                 │
            └───────────┬───────────┴───────────────────────────────────────┘
                        │
                        ▼
            ┌───────────────────────┐
            │   SIGNAL COMPOSITION  │
            │   (Spots per Nucleus) │
            └───────────┬───────────┘
                        │
                        ▼
            ┌───────────────────────┐
            │        EXPORT         │
            ├───────────────────────┤
            │ • Nuclei data (CSV)   │
            │ • Signal data (CSV)   │
            │ • Composition (CSV)   │
            │ • Parameters (MAT)    │
            │ • Visualizations      │
            └───────────────────────┘
```

---

## Modular Architecture

SNAP is compartmentalized so core UIs remain stable while algorithmic methods can evolve through pluggable modules.

### Design Intent

The architecture is built around scientific reproducibility:
- Keep preprocessing, detection, fitting, and classification behavior explicit and auditable
- Preserve quantitative output schemas unless intentionally versioned
- Keep channel behavior homogeneous by running channels through shared pipeline stages
- Allow extension without modifying `SNAP.m`, `SNAP_batch.m`, `SNAP_prepare.m`, or `SNAP_train.m`

### Shared Runtime Engines

- **Signal pipeline engine**: `+snap_modules/+signal/runPipeline.m`
  - Used by `SNAP_batch` and `SNAP_train`
  - Executes one selected plugin per stage (highest-priority enabled plugin wins)
  - Stage order:
    - `signal_processing`
    - `maxima_detection`
    - `gaussian_fitting`
    - `fit_filtering`
    - `classification`
    - `nuclei_filtering`

- **Prepare engines**:
  - `+snap_modules/+prepare/readLibrary.m`
  - `+snap_modules/+prepare/exportMappedChannels.m`
  - Used by `SNAP_prepare` for reader/exporter selection

### Built-In Module Sources

- Signal built-ins: `+snap_modules/+plugins/+signal/*.m`
- Prepare built-ins: `+snap_modules/+plugins/+prepare/*.m`

Built-ins wrap the established SNAP helper methods to preserve current behavior while exposing clear extension boundaries.

### External Extension Points

- Signal plugins:
  - `external_plugins/signal/`
  - `SNAP_EXTERNAL_SIGNAL_PLUGIN_PATH` (path-separated directories)
- SNAP_prepare reader plugins:
  - `external_plugins/prepare/readers/`
  - `SNAP_EXTERNAL_PREPARE_PROVIDER_PATH`
- SNAP_prepare exporter plugins:
  - `external_plugins/prepare/exporters/`
  - `SNAP_EXTERNAL_PREPARE_EXPORTER_PATH`

Naming is intentionally distinct:
- Built-in runtime plugins live only under `+snap_modules/+plugins/...`
- Community/external plugins live only under `external_plugins/...`

Included templates:
- `external_plugins/signal/template_signal_plugin.m`
- `external_plugins/prepare/readers/template_prepare_reader.m`
- `external_plugins/prepare/exporters/template_prepare_exporter.m`

### Collaboration Guide (In-Repo)

This README is the canonical contributor guide for modular development.

#### Engineering Rules

1. Keep app entrypoints thin (`SNAP.m`, `SNAP_batch.m`, `SNAP_prepare.m`, `SNAP_train.m`, `SNAP_classify.m`).
2. Add/replace algorithm behavior through plugins, not by hard-coding app logic.
3. Maintain backward compatibility for parameter files and exported analysis formats.
4. Keep progress logging informative but non-blocking.
5. Avoid hidden global state in plugins.

#### Signal Plugin Contract

Signal plugin function returns a struct with:
- `id` (unique string)
- `stage` (one of: `signal_processing`, `maxima_detection`, `gaussian_fitting`, `fit_filtering`, `classification`, `nuclei_filtering`)
- `displayName` (string)
- `priority` (numeric; highest enabled plugin wins within a stage)
- `version` (string)
- `supportsFcn` (optional function handle)
- `isEnabledFcn` (optional function handle)
- `run` (function handle with signature `state = run(state, context)`)

Common `state` fields:
- `rawImage`
- `processedImage`
- `maximaCoords`
- `fitResults`
- `aborted`

Common `context` fields:
- `channelIdx`
- `params`
- `handles`
- `mode`
- `progressCallback`

#### SNAP_prepare Reader Contract

Reader provider function returns a struct with:
- `id`, `displayName`, `priority`, `version`
- `canReadFcn(filePath) -> logical`
- `readFcn(filePath, progressCb) -> struct`

Reader output must include:
- `imageData` (cell)
- `tableData` (cell)
- `numChannels` (scalar numeric)

#### SNAP_prepare Exporter Contract

Exporter provider function returns a struct with:
- `id`, `displayName`, `priority`, `version`
- `canExportFcn(imageInfo, mapping, outputPath) -> logical`
- `exportFcn(imageInfo, mapping, outputPath, progressCb)`

#### Quality Checklist for Plugin PRs

1. Plugin loads without warnings.
2. Default built-in behavior is preserved in `SNAP_batch`, `SNAP_prepare`, and `SNAP_train`.
3. Output schemas remain stable unless intentionally versioned and documented.
4. Progress callbacks never crash or block the pipeline.
5. Current SNAP `.mat` parameter files and classifiers remain loadable/usable across SNAP apps.

### Example Contribution: SVM Feature Pack

An example, non-core contribution is included to demonstrate high-value SVM feature engineering and modular extensibility:

- `+snap_contrib/+svm/buildExpressionPack.m`
- `+snap_contrib/+svm/saveExpressionPack.m`
- `+snap_contrib/+svm/augmentFitResultsWithModelStats.m`
- `examples/svm_feature_pack/README.md`
- `examples/svm_feature_pack/create_example_expression_pack.m`
- `external_plugins/signal/template_signal_modelstats_plugin.m`

This pack provides a literature-informed custom expression set and optional model-selection statistics (AIC/BIC and residual descriptors) for advanced programmatic workflows.

#### How to use the custom pack

1. Export fitted channel data from SNAP or SNAP_batch (MAT exports now include `signals(i).rawDataWindow` and fit-window metadata).
2. Build and save a per-channel pack:
   ```matlab
   pack = snap_contrib.svm.buildExpressionPack('ParameterFile', 'svm_parameters.mat');
   snap_contrib.svm.saveExpressionPack(pack, 'snap_svm_expression_pack.mat');
   ```
3. In `SNAP_train`, use **Select Features...** per channel and apply the corresponding feature/expression set from the saved pack.
4. If you train programmatically, pass the channel-specific set into `SNAP_train(..., 'SelectedFeatures', ..., 'CustomExpressions', ...)`.

#### How to keep contributions out of SNAP core on GitHub

1. Keep core runtime in `+snap_helpers`, `+snap_modules`, and top-level app files (`SNAP.m`, `SNAP_batch.m`, `SNAP_prepare.m`, `SNAP_classify.m`, `SNAP_train.m`).
2. Put optional community extensions under `+snap_contrib`, `examples`, and `external_plugins`.
3. Ensure core code does not hard-depend on `+snap_contrib`; contributions should be opt-in only.
4. Require each contribution PR to include:
   - A short README in its folder
   - Example invocation
   - Clear statement of required inputs/outputs and any added fields
5. Use a `contribution` label in GitHub PRs so these remain visibly separate from core maintenance PRs.
6. Any PR that touches core paths must carry the `core` label (enforced by workflow).

#### GitHub mediation settings to enable in repo UI

1. Protect `main` (Settings → Branches):
   - Require pull request before merging
   - Require approvals
   - Require conversation resolution
   - Disable direct pushes to `main`
2. Mark these checks as required on `main`:
   - `MATLAB Smoke Checks`
   - `guard` (from `Contribution Scope Guard`)
3. Keep CODEOWNERS enabled so core/contribution paths are auto-routed for review.
4. Run `Sync Repository Labels` once (Actions tab) so template/labeler labels exist.

---

## Documentation

### Project Structure

```
SNAP/
├── SNAP.m                    # Main interactive GUI
├── SNAP_batch.m              # Batch processing (GUI + CLI)
├── SNAP_prepare.m            # Multi-channel format converter
├── SNAP_classify.m           # SVM classifier training
├── SNAP_train.m              # Programmatic + interactive SVM training from labels
├── compareMaximaWithLabeled.m # Utility for maxima/label comparison
├── README.md                 # This file
│
├── +snap_modules/            # Shared modular engines
│   ├── +signal/              # Shared signal pipeline runner + registry
│   ├── +prepare/             # SNAP_prepare provider/exporter engines
│   └── +plugins/             # Built-in module plugins
│
├── +snap_contrib/            # Example contributor extensions (non-core)
│   └── +svm/                 # SVM feature/expression contribution package
│
├── examples/
│   └── svm_feature_pack/     # Example usage/docs for SVM contribution
│
└── +snap_helpers/            # Core processing functions
    ├── createUI.m            # Build main GUI
    ├── updateControls.m      # Dynamic UI state management
    ├── updateLivePreview.m   # Processing + preview refresh (manual trigger)
    │
    ├── # Image Processing
    ├── loadImage.m           # TIFF stack loading
    ├── processImage.m        # Apply preprocessing pipeline
    ├── preprocessNuclei.m    # Nuclei preprocessing helpers
    ├── preprocessNucleiWithBgCorr.m # Nuclei preprocessing + background correction
    │
    ├── # Nuclei Analysis
    ├── segmentNuclei.m       # Multi-algorithm segmentation
    ├── generateNucleiLabels.m # Consistent nucleus labeling
    ├── computeNucleusMeasurements.m # Morphological measurements
    │
    ├── # Signal Detection
    ├── detectMaxima.m        # Local maxima detection
    ├── fitGaussians.m        # Sub-pixel Gaussian fitting
    ├── applyFitFiltering.m   # Quality-based filtering
    ├── filterMaximaByNuclei.m # Nucleus inclusion/exclusion
    ├── computeSignalMeasurements.m # Signal measurements
    │
    ├── # Export Functions
    ├── exportData.m          # Unified export controller
    ├── exportNucleiDataStandardized.m
    ├── exportChannelDataStandardized.m
    ├── exportNucleiSignalDataStandardized.m
    │
    └── +classification/      # Machine learning module
        ├── trainClassifier.m     # SVM training
        ├── applyClassifier.m     # Apply trained model
        ├── buildFeatureMatrix.m  # Feature extraction
        ├── featureSelectionUI.m  # Feature selection dialog
        ├── evaluateExpression.m  # Custom feature expressions
        ├── saveClassifier.m      # Save model + normalization
        └── loadClassifier.m      # Load trained classifier

external_plugins/
├── signal/                   # External signal-stage plugins
└── prepare/
    ├── readers/              # External library reader providers
    └── exporters/            # External channel exporters
```

### Key Parameters

#### Nuclei Segmentation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `nucSegMainMethod` | Thresholding strategy | `'Absolute'` |
| `nucSegLocalAlgorithm` | Auto Local Threshold algorithm | `'Otsu'` |
| `nucFilterMinSize` | Minimum nucleus size | `100 pixels` |
| `nucFilterMinCircularity` | Minimum circularity (0-1) | `0.3` |
| `nucExcludeEdges` | Remove border-touching nuclei | `true` |

#### Maxima Detection

| Parameter | Description | Default |
|-----------|-------------|---------|
| `maximaMethod` | Detection algorithm | `'Simple Regional'` |
| `maximaNeighborhoodSize` | Search radius (pixels or µm if scaled) | `1` |
| `hMaxValue` | H-maxima threshold | `0` |
| `sigmaValue` | LoG sigma | `2` |

#### Gaussian Fitting

| Parameter | Description | Default |
|-----------|-------------|---------|
| `gaussFitMethod` | Fitting model | `'1D (X,Y,Z) Gaussian'` |
| `gaussFitVoxelWindowSize` | Fitting window size | `7` |
| `gaussFitBgCorrMethod` | Background correction | `'Mean Surrounding Subtraction'` |

#### Fit Filtering

| Parameter | Description | Default |
|-----------|-------------|---------|
| `fitFilterRSquaredMin` | Minimum R² for quality fits | `0.8` |
| `fitFilterSigmaSumMin/Max` | Sigma sum range | `0 - 10` |
| `fitFilterAmplitudeMin/Max` | Amplitude range | `100 - 10000` |

---

### Supported File Formats

#### Direct Loading (SNAP)
- TIFF stacks (`.tif`, `.tiff`, `.ome.tif`) - single or multi-page

#### Via SNAP_prepare (Bio-Formats)
- MetaMorph Stack (`.mvd2`)
- Zeiss CZI (`.czi`)
- Nikon NIS-Elements (`.nd2`)
- Leica Image File (`.lif`)
- Olympus (`.oib`, `.oif`)
- Zeiss LSM (`.lsm`)
- Imaris (`.ims`)
- And many more via Bio-Formats

---

### Batch Processing

SNAP_batch expects input folders organized as:
```
input_directory/
├── Sample_001/
│   ├── nuclei.tif      # or DAPI.tif, hoechst.tif
│   ├── channel1.tif    # or ch1.tif, c01.tif
│   ├── channel2.tif
│   └── dic.tif         # optional
├── Sample_002/
│   └── ...
```

This structure is automatically created by `SNAP_prepare` when converting multi-channel libraries.

#### Command-Line Usage

```matlab
% Basic usage
SNAP_batch('data/', 'params.mat')

% With options
SNAP_batch('data/', 'params.mat', ...
    'OutputDir', 'results/', ...
    'ExportFormat', 'TIFF')
```

---

### Classification Workflow

1. **Export fit results** from SNAP with Gaussian fitting enabled
2. **Launch SNAP_classify**:
   ```matlab
   SNAP_classify
   ```
3. **Load exported data** (MAT file with fit results)
4. **Load channel image** for visual inspection
5. **Select features** to use for classification
6. **Label spots** using keyboard shortcuts:
   - `R` = Real spot
   - `N` = Noise
   - `S` = Skip
   - `←/→` = Navigate
   - `Ctrl+Z` = Undo
7. **Train SVM** (requires ≥5 labels per class)
8. **Export classifier** for use in SNAP or SNAP_batch

---

### Train from Existing Labels (`SNAP_train`)

`SNAP_train` now supports both:
- **GUI mode** (`SNAP_train`) for multi-channel, per-channel SVM training
- **Programmatic mode** (`SNAP_train(exportFiles, labelFiles, outputClassifierPath, ...)`) for scripted training

#### GUI Mode (Recommended for multi-channel workflows)

```matlab
SNAP_train
```

The training UI provides:
- Parameter file loading (`.mat`) to infer the active channel count and per-channel fitting context
- Exactly one potential SVM per detected channel (enable/disable channels with the per-channel train checkbox)
- Per-channel training and validation directory assignment (set independently for each channel; validation is required when sweep mode is enabled)
- Shared output directory with per-channel classifier output file names
- Adjustable **Match Distance** for label-to-candidate assignment
- Two training strategies:
  - Manual hyperparameters (kernel, box constraint, kernel scale, polynomial order, standardization, CV folds)
  - Validation sweep optimization (kernel/grid search with held-out validation files)
- Optional sweep performance reporting (toggle-able log + plot output for hyperparameter comparison)

Training discovers export/label pairs recursively and trains one classifier per selected channel.

#### Programmatic Mode

If you already have labeled spots for many images, you can train directly from file lists:

```matlab
exportFiles = {'image1_ch1_signals.mat', 'image2_ch1_signals.mat'};
labelFiles  = {'image1_labels.csv', 'image2_labels.csv'};

SNAP_train(exportFiles, labelFiles, 'channel1_classifier.mat', ...
    'MatchDistance', 2, ... % voxels
    'ValidationExportFiles', valExportFiles, ...
    'ValidationLabelFiles', valLabelFiles, ...
    'HyperparameterSweep', true);
```

Supported label files:
- CSV/table with coordinate columns (`maxima_y`/`fitted_y`/`y`, `maxima_x`/`fitted_x`/`x`, optional z) plus `label`
- MAT progress files from `SNAP_classify` (`labeledReal`, `labeledNoise`)

For coordinate labels, each manual label is matched to the nearest unassigned candidate within `MatchDistance` voxels (one candidate per manual label). Any unmatched candidates are automatically labeled as noise so every candidate contributes to training. If validation files are provided, `SNAP_train` sweeps SVM parameters (kernel, box constraint, kernel scale / polynomial order) and selects the best model by validation F1 (real class). The training-set real/noise labels are fixed from training volumes and are not relabeled during validation.

The output classifier is saved in the same format as `SNAP_classify` export and can be loaded directly in `SNAP` / `SNAP_batch`.

---


## Output Data Format

### Nuclei Data (CSV)

| Column | Description |
|--------|-------------|
| `image_name` | Source image identifier |
| `nucleus_id` | Unique nucleus ID |
| `centroid_row`, `centroid_col`, `centroid_slice` | 3D centroid position |
| `area_pixels` / `volume_voxels` | Size in pixels/voxels |
| `area_um2` / `volume_um3` | Size in physical units |
| `major_axis_length`, `minor_axis_length` | Ellipse fit |
| `circularity` / `sphericity` | Shape metric |
| `solidity` | Area/ConvexHullArea |
| `mean_intensity`, `integrated_intensity` | Intensity measurements |

### Signal Data (CSV)

| Column | Description |
|--------|-------------|
| `signal_id` | Unique signal ID |
| `maxima_x`, `maxima_y`, `maxima_z` | Original detection coordinates |
| `fitted_x`, `fitted_y`, `fitted_z` | Sub-pixel fitted coordinates |
| `amplitude` or `amplitude_*` | Method-dependent amplitude fields |
| `sigma_x`, `sigma_y`, `sigma_z` | Gaussian widths |
| `r_squared` or `radial_symmetry_score` | Method-dependent quality metric |
| `integrated_intensity` | Total signal |
| `background` | Local background estimate |

### Signal Data (MAT)

Channel and nuclei-signal MAT exports include all CSV-equivalent signal fields plus full fit-window context when fitting is available:
- `rawDataWindow` (full local intensity window used for fitting)
- `fitWindowDimensions`
- `fitWindowOrigin`
- `localMaximaInWindow`
- `originalMaximaCoords`

This enables advanced downstream feature engineering (for example AIC/BIC and residual-model statistics) directly from exported results.
For large datasets, SNAP automatically falls back to MAT `-v7.3` when needed.

---

## Tips & Best Practices

1. **Start with Preview**: Always use "Update Previews" to verify parameters before batch processing

2. **Parameter Optimization**: Use the interactive GUI to find optimal parameters on representative images, then save for batch use

3. **Scaling**: Enable "Scale" checkboxes when working with physical units (µm) rather than pixels

4. **3D vs 2D**: For 3D data, consider:
   - "2D (Slice-by-slice)" for thick samples with independent z-slices
   - "3D" mode for connected structures spanning multiple z-planes
   - "On Z-Projection" for quick 2D analysis of 3D stacks

5. **Classification**: Train classifiers on diverse examples from your dataset for best generalization

6. **Memory**: For very large images, consider processing in tiles or using the batch system with multiple smaller images

---

## Troubleshooting

### Common Issues

**Q: SNAP_prepare can't read my file format**
> Install Bio-Formats Toolbox from MATLAB Add-Ons (Home → Add-Ons → Get Add-Ons → Search "Bio-Formats")

**Q: Processing is slow**
> - Reduce preview image size during parameter tuning
> - Use "On Z-Projection" mode for initial parameter exploration
> - Enable only needed processing steps
> - Use "Abort Processing" if a trial run is clearly too large

**Q: No spots detected**
> - Check that maxima detection is enabled for the channel
> - Lower the neighborhood size or H-maxima threshold
> - Verify preprocessing isn't removing signal (check raw vs processed preview)

**Q: All spots filtered out**
> - Check fit filtering thresholds (R², sigma, amplitude ranges)
> - Verify Gaussian fitting window size is appropriate for spot size
> - Check nuclei inclusion/exclusion settings

---

## Citation

If you use SNAP in your research, please cite:

```bibtex
@software{snap2025,
  title = {SNAP: Spot \& Nuclei Analysis Pipeline},
  author = {Rojas-Cessa, M. A.},
  year = {2025},
  url = {https://github.com/marcorojas-cessa/SNAP}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Auto Local Threshold algorithms based on [ImageJ implementation](https://imagej.net/plugins/auto-local-threshold)
- Bio-Formats library by [Open Microscopy Environment](https://www.openmicroscopy.org/bio-formats/)

---

<p align="center">
  <strong>Developed by the Rothstein Lab</strong><br>
  Columbia University Irving Medical Center
</p>
