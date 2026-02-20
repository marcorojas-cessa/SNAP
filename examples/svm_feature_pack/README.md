# SNAP Example Contribution: SVM Feature and Expression Pack

This folder contains an example contributor package for improving classifier robustness in spot/noise discrimination.

## What This Adds

- A **literature-informed expression pack** for `SNAP_train` / `SNAP_classify` workflows.
- Optional **model-selection/statistical features** computed from each fit window:
  - Gaussian-vs-flat AIC/BIC deltas
  - Poisson NLL deltas
  - Residual distribution descriptors (std, MAD, skewness, kurtosis, optional KS p-value)

These are delivered as an example contribution and are not wired into core SNAP UI logic by default.

## Files

- `+snap_contrib/+svm/buildExpressionPack.m`
- `+snap_contrib/+svm/saveExpressionPack.m`
- `+snap_contrib/+svm/augmentFitResultsWithModelStats.m`
- `examples/svm_feature_pack/create_example_expression_pack.m`
- `external_plugins/signal/template_signal_modelstats_plugin.m`

## Scientific Rationale

The expression set emphasizes:

- Dynamic-range stabilization (`log10(...)`) for skewed fluorescence intensity distributions
- Signal-vs-background normalization (SNR-like transforms)
- Shape/compactness terms (`intensity / sigma-product`) as PSF-consistency proxies
- Quality-coupled terms (`r_squared`-weighted intensity metrics)
- Optional model-selection statistics (AIC/BIC, likelihood deltas) for distinguishing structured spot-like signal from local noise models

These are broadly consistent with common practice in fluorescence spot detection and model-selection workflows.

## Basic Usage

```matlab
% Build per-channel expression pack from parameter file
pack = snap_contrib.svm.buildExpressionPack('ParameterFile', 'svm_parameters.mat');

% Save MAT + CSV
snap_contrib.svm.saveExpressionPack(pack, 'snap_svm_expression_pack.mat');
```

Or run:

```matlab
pack = create_example_expression_pack('svm_parameters.mat');
```

For `SNAP_train` UI usage, pick each channel and open **Select Features...**, then copy in that channel’s `selectedFeatures` + `customExpressions` from the saved pack.

## Advanced Optional Features (Model Stats)

If you have fit results containing `rawDataWindow` (from fitting output or exported SNAP MAT `signals`), you can augment them:

```matlab
[fitDataAug, summary] = snap_contrib.svm.augmentFitResultsWithModelStats(fitData, ...
    'ComputeNormalityPValue', false, 'Verbose', true);
```

Then include these added fields as **base selected features** in programmatic training workflows.

## Compartmentalized Extension Path

An optional template plugin is provided at:

- `external_plugins/signal/template_signal_modelstats_plugin.m`

It demonstrates how a contributor can wrap the built-in Gaussian fitting stage and append custom fit statistics using the modular signal pipeline architecture.
