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

Important:
- `'/ABSOLUTE/PATH/to/...'` is placeholder text and will fail if copied literally.
- Use real absolute paths if you are not running from the SNAP repo root, for example:

```matlab
pack = create_example_expression_pack( ...
    '/path/to/SNAP/svm_parameters.mat', ...
    '/path/to/SNAP/examples/svm_feature_pack/snap_svm_expression_pack.mat');
```

For `SNAP_train` UI usage, load your parameter file, then click **Load Expression Pack...** to apply channel-matched `selectedFeatures` + `customExpressions` directly from the saved pack.

Compatibility behavior in `SNAP_train`:
- The loaded pack is checked per-channel against inferred fitting context.
- Incompatible base features/custom expressions are dropped automatically and reported in the training log.
- During training, all-NaN features are automatically pruned and feature extraction is retried.
- If an entire channel feature set becomes incompatible, `SNAP_train` falls back to AUTO base features.

## Advanced Optional Features (Model Stats)

If you have fit results containing `rawDataWindow` (from fitting output or exported SNAP MAT `signals`), you can augment them:

```matlab
[fitDataAug, summary] = snap_contrib.svm.augmentFitResultsWithModelStats(fitData, ...
    'ComputeNormalityPValue', false, 'Verbose', true);
```

Then include these added fields as **base selected features** in programmatic training workflows.

## Creating Your Own Pack (From Scratch)

Minimum required structure:

```matlab
pack = struct();
pack.name = 'My pack';
pack.version = '1.0.0';
pack.channelPacks = struct( ...
    'channelIdx', {1}, ...
    'selectedFeatures', {{'amplitude', 'background', 'r_squared'}}, ...
    'customExpressions', struct('name', {'snr_like'}, 'expression', {'integrated_intensity / background'}));
```

Save with:

```matlab
snap_contrib.svm.saveExpressionPack(pack, 'my_expression_pack.mat');
```

Then load in `SNAP_train` with **Load Expression Pack...**.

## Compartmentalized Extension Path

An optional template plugin is provided at:

- `external_plugins/signal/template_signal_modelstats_plugin.m`

It demonstrates how a contributor can wrap the built-in Gaussian fitting stage and append custom fit statistics using the modular signal pipeline architecture.
