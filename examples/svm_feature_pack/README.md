# SNAP Example: SVM Expression Pack

This folder contains an example feature/expression pack for spot-vs-noise classification.

## One Important Point

Building a pack does **not** train an SVM.

- Pack generation only defines which base features and custom expressions to use.
- SVM training still happens in `SNAP_train`.

## Channel Semantics

- In `SNAP_train`, channel numbers are slot labels for bookkeeping.
- Each SVM is trained independently per channel slot.
- There is no cross-channel pooling during training.
- If you want the same trained SVM on multiple channels, load/inject the same classifier file into each channel slot.

## Quick Start

1. Build a pack:
```matlab
pack = create_example_expression_pack('svm_parameters.mat');
```

2. Train:
```matlab
SNAP_train
```

3. In `SNAP_train`:
- Load parameter file.
- Click **Load Expression Pack...** and choose the generated pack `.mat`.
- Configure training data for each channel slot.
- Train selected channels.

## `SelectionMode`

Controls base-feature selection during pack generation.

- `focused` (default): prioritized compact subset.
- `all_non_position`: all non-position base features.

Used in:
- `snap_contrib.svm.buildExpressionPack(..., 'SelectionMode', ...)`
- `create_example_expression_pack(..., 'SelectionMode', ...)`

## `LintMode`

Validation mode for pack compatibility checks.

- It does **not** train SVMs.
- It only validates/sanitizes pack contents.

Modes:
- `strict`: fail on incompatibilities.
- `permissive`: drop incompatible entries with warnings (and optional `real(...)` safety guard).

Used in:
- `snap_contrib.svm.lintExpressionPack(..., 'Mode', ...)`
- `create_example_expression_pack(..., 'LintMode', ...)`

## Compatibility Behavior

When a pack is loaded in `SNAP_train`:

- schema is normalized
- channel capability checks are run from the loaded parameter file
- incompatible entries are pruned in permissive mode with log messages
- expressions are stress-checked for numeric safety

## Files

- `+snap_contrib/+svm/buildExpressionPack.m`
- `+snap_contrib/+svm/saveExpressionPack.m`
- `+snap_contrib/+svm/lintExpressionPack.m`
- `+snap_contrib/+svm/augmentFitResultsWithModelStats.m`
- `examples/svm_feature_pack/create_example_expression_pack.m`
- `examples/svm_feature_pack/snap_svm_expression_pack_expressions.csv`

## Optional Advanced Stats

If fit results include raw fit windows, you can compute optional model-stat features:

```matlab
[fitDataAug, summary] = snap_contrib.svm.augmentFitResultsWithModelStats(fitData, ...
    'ComputeNormalityPValue', false, ...
    'Verbose', true);
```

## Contributor Check

Before sharing a pack:

```matlab
report = snap_contrib.svm.lintExpressionPack('snap_svm_expression_pack.mat', ...
    'ParameterFile', 'svm_parameters.mat', ...
    'Mode', 'strict');
```
