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
- Vectorized local thresholding (10-50x faster than pixel-by-pixel)
- Efficient batch processing with progress tracking
- Preview caching for responsive UI

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
% 4. Click "Update Previews" to see results
% 5. Export data when satisfied
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

## Documentation

### Project Structure

```
SNAP/
├── SNAP.m                    # Main interactive GUI
├── SNAP_batch.m              # Batch processing (GUI + CLI)
├── SNAP_prepare.m            # Multi-channel format converter
├── SNAP_classify.m           # SVM classifier training
├── README.md                 # This file
│
└── +snap_helpers/            # Core processing functions
    ├── createUI.m            # Build main GUI
    ├── updateControls.m      # Dynamic UI state management
    ├── updateLivePreview.m   # Real-time preview updates
    │
    ├── # Image Processing
    ├── loadImage.m           # TIFF stack loading
    ├── processImage.m        # Apply preprocessing pipeline
    ├── preprocessNuclei.m    # Nuclei-specific preprocessing
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
- TIFF stacks (`.tif`, `.tiff`) - single or multi-page

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

### Train from Existing Labels (Command Line)

If you already have labeled spots for many images, you can train a classifier without the manual UI:

```matlab
exportFiles = {'image1_ch1_signals.mat', 'image2_ch1_signals.mat'};
labelFiles  = {'image1_labels.csv', 'image2_labels.csv'};

SNAP_trainSVM(exportFiles, labelFiles, 'channel1_classifier.mat', ...
    'MatchDistance', 2, ... % voxels
    'ValidationExportFiles', valExportFiles, ...
    'ValidationLabelFiles', valLabelFiles, ...
    'HyperparameterSweep', true);
```

Or run interactively and provide training/validation directories + output when prompted:

```matlab
SNAP_trainSVM
```

Supported label files:
- CSV/table with coordinate columns (`maxima_y`/`fitted_y`/`y`, `maxima_x`/`fitted_x`/`x`, optional z) plus `label`
- MAT progress files from `SNAP_classify` (`labeledReal`, `labeledNoise`)

For coordinate labels, each manual label is matched to the nearest unassigned candidate within `MatchDistance` voxels (one candidate per manual label). Any unmatched candidates are automatically labeled as noise so every candidate contributes to training. If validation files are provided, `SNAP_trainSVM` sweeps SVM parameters (kernel, box constraint, kernel scale / polynomial order) and selects the best model by validation F1 (real class). The training-set real/noise labels are fixed from training volumes and are not relabeled during validation.

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
| `maxima_coords` | Original detection coordinates |
| `fitted_coords` | Sub-pixel fitted coordinates |
| `amplitude` | Gaussian amplitude |
| `sigma_x`, `sigma_y`, `sigma_z` | Gaussian widths |
| `r_squared` | Fit quality (0-1) |
| `integrated_intensity` | Total signal |
| `background` | Local background estimate |

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
