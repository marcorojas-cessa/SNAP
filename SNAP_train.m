function varargout = SNAP_train(exportFiles, labelFiles, outputClassifierPath, varargin)
% SNAP_train - Train and validate a SNAP-compatible SVM from labeled spot files
%
% UI MODE:
%   SNAP_train
%   Opens a SNAP-style training GUI for channel-aware multi-SVM training.
%   The GUI:
%     - Loads a SNAP parameter file to infer channel count and per-channel defaults
%     - Uses exactly one independent SVM slot per detected channel (toggle on/off)
%     - Supports a parameter file per channel tab (with template file prefill)
%     - Sets match distance per channel (per SVM)
%     - Lets you assign training/validation directories per channel
%     - Lets you choose output directory/classifier filename per channel
%     - Lets you select base features + custom expressions (shared engine with SNAP_classify)
%     - Supports manual SVM hyperparameters OR validation-based optimization
%     - Optionally emits a sweep performance report (log + plot)
%     - Trains one classifier per selected channel without cross-channel pooling
%   NOTE:
%     Channel numbers are slot labels for bookkeeping. Training is always
%     channel-local. Reuse a trained classifier on multiple channels by
%     loading/injecting the same classifier file into each channel slot.
%
% PROGRAMMATIC MODE:
%   SNAP_train(exportFiles, labelFiles, outputClassifierPath, ...)
%
% REQUIRED INPUTS:
%   exportFiles           - char/string/cellstr of SNAP export MAT file(s)
%   labelFiles            - char/string/cellstr of label file(s)
%   outputClassifierPath  - output .mat classifier path
%
% LABEL FILE FORMATS:
%   1) CSV/table with columns:
%      - maxima_y, maxima_x, [maxima_z], and optional label/class/is_real
%      - OR fitted_y, fitted_x, [fitted_z], and optional label/class/is_real
%      - OR y, x, [z or slice], and optional label/class/is_real
%      label values accepted: {1,true,'real','spot','positive'} and
%                             {0,false,'noise','negative','bg'}
%   2) MAT progress file from SNAP_classify containing
%      `labeledReal` and `labeledNoise` indices.
%
% MATCHING + LABELING:
%   - Manual labels are matched to nearest unassigned candidate within
%     MatchDistance voxels (one-to-one).
%   - Any unmatched/unlabeled candidates are automatically labeled noise (0)
%     so all candidates are used for SVM training.
%
% NAME-VALUE OPTIONS:
%   'SelectedFeatures'       - cellstr of base features (default: automatic)
%   'CustomExpressions'      - struct array with fields: name, expression
%   'MatchDistance'          - max coordinate distance in voxels (default: 2)
%   'ConvertFijiCoords'      - convert CSV coords from FIJI [x,y,z] (0-indexed)
%                              to MATLAB array indexing [row,col,slice] via
%                              [y+1, x+1, z+1] (default: false)
%   'ParameterStruct'        - SNAP parameter struct for image->candidate generation
%   'ChannelIndex'           - channel index used with ParameterStruct
%   'FittingMethod'          - override fitting method string
%   'Has3D'                  - override 3D flag
%   'TrainingOptions'        - options struct for trainClassifier
%   'ValidationExportFiles'  - export files for held-out validation
%   'ValidationLabelFiles'   - label files for held-out validation
%   'HyperparameterSweep'    - true/false (default: true if validation set)
%   'SweepKernels'           - kernels to test (default: {'linear','rbf','polynomial'})
%   'SweepBoxConstraints'    - box constraints to test (default: [0.1 1 10])
%   'SweepKernelScales'      - kernel scales for non-linear kernels (default: {'auto',0.5,1,2})
%   'SweepPolynomialOrders'  - polynomial orders (default: [2 3 4])
%   'SweepTieHandling'       - tie policy for equal-best score:
%                              'prefer_simple' (default), 'prompt', 'first'
%   'SweepTiePromptCallback' - optional callback for tie selection in prompt mode
%   'Verbose'                - true/false (default: true)
%   'ProgressCallback'       - optional function handle for step-by-step progress logs
%
% OUTPUT:
%   out - struct with trained model, selected params, train/validation stats

    if nargin == 0
        uiState = runTrainingUI();
        fprintf('SNAP_train initiated.\n');
        if nargout > 0
            varargout{1} = uiState;
        end
        return;
    end

    exportFiles = normalizeFileList(exportFiles);
    labelFiles = normalizeFileList(labelFiles);

    if numel(exportFiles) ~= numel(labelFiles)
        error('exportFiles and labelFiles must have the same number of entries.');
    end

    p = inputParser;
    p.addParameter('SelectedFeatures', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.addParameter('CustomExpressions', struct('name', {}, 'expression', {}), @isstruct);
    p.addParameter('MatchDistance', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    p.addParameter('ConvertFijiCoords', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ParameterStruct', struct(), @isstruct);
    p.addParameter('ChannelIndex', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
    p.addParameter('FittingMethod', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Has3D', [], @(x) islogical(x) || isnumeric(x));
    p.addParameter('TrainingOptions', struct(), @isstruct);
    p.addParameter('ValidationExportFiles', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.addParameter('ValidationLabelFiles', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.addParameter('HyperparameterSweep', [], @(x) islogical(x) || isnumeric(x));
    p.addParameter('SweepKernels', {'linear', 'rbf', 'polynomial'}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.addParameter('SweepBoxConstraints', [0.1 1 10], @(x) isnumeric(x) && isvector(x) && all(x > 0));
    p.addParameter('SweepKernelScales', {'auto', 0.5, 1, 2}, @(x) iscell(x) || isstring(x) || isnumeric(x) || ischar(x));
    p.addParameter('SweepPolynomialOrders', [2 3 4], @(x) isnumeric(x) && isvector(x) && all(x >= 2));
    p.addParameter('SweepTieHandling', 'prefer_simple', @(x) ischar(x) || isstring(x));
    p.addParameter('SweepTiePromptCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ProgressCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.parse(varargin{:});
    opts = p.Results;

    verbose = logical(opts.Verbose);
    progressCb = opts.ProgressCallback;

    datasetBuildOpts = struct();
    datasetBuildOpts.convertFijiCoords = logical(opts.ConvertFijiCoords);
    datasetBuildOpts.parameterStruct = opts.ParameterStruct;
    datasetBuildOpts.channelIndex = opts.ChannelIndex;
    datasetBuildOpts.progressCallback = progressCb;
    datasetBuildOpts.datasetLabel = 'training';

    emitProgress(progressCb, 'Building labeled training dataset (%d pair(s))...', numel(exportFiles));

    [trainFitData, trainLabels, trainSources, inferredFittingMethod, inferredHas3D] = ...
        buildLabeledDataset(exportFiles, labelFiles, opts.MatchDistance, datasetBuildOpts);

    if isempty(trainFitData)
        error('No training candidates found.');
    end

    fittingMethod = string(opts.FittingMethod);
    if strlength(fittingMethod) == 0
        fittingMethod = string(inferredFittingMethod);
    end

    has3D = opts.Has3D;
    if isempty(has3D)
        has3D = inferredHas3D;
    else
        has3D = logical(has3D);
    end

    [~, featureInfo] = snap_helpers.classification.getAvailableFeatures(char(fittingMethod), has3D, false);
    selectedFeatures = normalizeSelectedFeatures(opts.SelectedFeatures, featureInfo);
    customExpressions = normalizeCustomExpressionList(opts.CustomExpressions);

    capability = snap_helpers.classification.resolveCapabilitiesFromContext( ...
        char(fittingMethod), has3D, 'ChannelIndex', 1, 'HasPhysicalSpacing', false);
    preflightPack = struct( ...
        'specVersion', '2.0.0', ...
        'packId', 'snap_train_programmatic_preflight', ...
        'strictModeDefault', false, ...
        'channelPacks', struct( ...
            'channelIdx', 1, ...
            'selectedFeatures', {selectedFeatures}, ...
            'customExpressions', customExpressions, ...
            'requiredFeatures', {{}}, ...
            'requiredCapabilities', struct( ...
                'fittingMethod', capability.fittingMethod, ...
                'has3D', logical(capability.has3D), ...
                'hasPhysicalSpacing', false)));

    [preflightPack, preflightReport] = ...
        snap_helpers.classification.validateExpressionPackAgainstCapabilities( ...
            preflightPack, capability, 'Mode', 'permissive', 'AutoGuardUnsafeExpressions', true);

    if ~isempty(preflightPack.channelPacks)
        selectedFeatures = preflightPack.channelPacks(1).selectedFeatures;
        customExpressions = preflightPack.channelPacks(1).customExpressions;
    end
    if preflightReport.nDroppedBase > 0 || preflightReport.nDroppedCustom > 0
        emitProgress(progressCb, ...
            'Feature preflight pruning: dropped base=%d, dropped custom=%d.', ...
            preflightReport.nDroppedBase, preflightReport.nDroppedCustom);
    end
    if preflightReport.nAutoGuarded > 0
        emitProgress(progressCb, ...
            'Feature preflight safety: auto-guarded %d custom expression(s).', ...
            preflightReport.nAutoGuarded);
    end
    if ~isempty(preflightReport.errors)
        error('Feature preflight validation failed: %s', strjoin(preflightReport.errors, ' | '));
    end
    for i = 1:numel(preflightReport.warnings)
        emitProgress(progressCb, 'Feature preflight warning: %s', preflightReport.warnings{i});
    end

    [XTrain, featureNames, validTrainMask, extractionInfoTrain, selectedFeatures, customExpressions] = ...
        buildTrainingFeatureMatrixWithFallback( ...
            trainFitData, selectedFeatures, featureInfo, customExpressions, progressCb);

    XTrain = XTrain(validTrainMask, :);
    yTrain = trainLabels(validTrainMask);
    trainSources = trainSources(validTrainMask);
    fixedTrainLabels = yTrain; % immutable reference labels used for all sweep candidates
    emitProgress(progressCb, 'Training dataset after filtering: %d samples (real=%d, noise=%d).', ...
        numel(yTrain), sum(yTrain == 1), sum(yTrain == 0));

    if isempty(yTrain)
        detail = '';
        if isstruct(extractionInfoTrain) && isfield(extractionInfoTrain, 'featuresAllNaN') && ...
                ~isempty(extractionInfoTrain.featuresAllNaN)
            bad = extractionInfoTrain.featuresAllNaN;
            bad = bad(1:min(8, numel(bad)));
            detail = sprintf(' Features all-NaN: %s.', strjoin(bad, ', '));
        end
        if isstruct(extractionInfoTrain) && isfield(extractionInfoTrain, 'warnings') && ~isempty(extractionInfoTrain.warnings)
            warnPreview = extractionInfoTrain.warnings;
            warnPreview = warnPreview(1:min(3, numel(warnPreview)));
            detail = sprintf('%s Feature extraction warnings: %s', detail, strjoin(warnPreview, ' | '));
        end
        error(['No valid training samples remained after feature extraction. ', ...
               'This usually means one or more selected base features or custom expressions are incompatible ', ...
               'with the generated fit-result fields for this fitting method.%s'], detail);
    end

    if sum(yTrain == 1) < 5 || sum(yTrain == 0) < 5
        error('Need at least 5 valid samples per class for training. Got real=%d, noise=%d.', ...
            sum(yTrain == 1), sum(yTrain == 0));
    end

    valExport = normalizeOptionalFileList(opts.ValidationExportFiles);
    valLabels = normalizeOptionalFileList(opts.ValidationLabelFiles);

    if ~isempty(valExport) || ~isempty(valLabels)
        if numel(valExport) ~= numel(valLabels)
            error('ValidationExportFiles and ValidationLabelFiles must have the same number of entries.');
        end
    end

    hasValidationSet = ~isempty(valExport);
    doSweep = hasValidationSet;
    if ~isempty(opts.HyperparameterSweep)
        doSweep = logical(opts.HyperparameterSweep);
    end

    validation = struct();
    model = [];
    trainStats = struct();
    normParams = struct();
    bestParams = struct();

    if hasValidationSet
        datasetBuildOptsVal = datasetBuildOpts;
        datasetBuildOptsVal.datasetLabel = 'validation';
        emitProgress(progressCb, 'Building labeled validation dataset (%d pair(s))...', numel(valExport));
        [valFitData, valLabelsRaw] = buildLabeledDataset(valExport, valLabels, opts.MatchDistance, datasetBuildOptsVal);
        [XVal, ~, validValMask, extractionInfoVal] = snap_helpers.classification.buildFeatureMatrix( ...
            valFitData, selectedFeatures, featureInfo, customExpressions);
        emitProgress(progressCb, 'Validation feature matrix built: %d candidate rows, %d feature columns.', size(XVal,1), size(XVal,2));
        if isstruct(extractionInfoVal) && isfield(extractionInfoVal, 'modelStats') && ...
                isstruct(extractionInfoVal.modelStats) && ...
                isfield(extractionInfoVal.modelStats, 'requested') && extractionInfoVal.modelStats.requested
            msVal = extractionInfoVal.modelStats;
            if isfield(msVal, 'augmented') && msVal.augmented && isfield(msVal, 'summary') && isstruct(msVal.summary)
                emitProgress(progressCb, ...
                    ['Validation model-stat augmentation: computed %d/%d windows ', ...
                     '(missingWindow=%d, modelFailures=%d).'], ...
                    getfieldwithdefault(msVal.summary, 'nComputed', 0), ...
                    getfieldwithdefault(msVal.summary, 'nTotal', 0), ...
                    getfieldwithdefault(msVal.summary, 'nMissingWindow', 0), ...
                    getfieldwithdefault(msVal.summary, 'nModelFailures', 0));
            else
                emitProgress(progressCb, 'Validation model-stat features requested: using existing fields (no augmentation run).');
            end
        end
        if isstruct(extractionInfoVal) && isfield(extractionInfoVal, 'featuresAllNaN') && ...
                ~isempty(extractionInfoVal.featuresAllNaN)
            emitProgress(progressCb, 'Validation feature extraction all-NaN feature(s): %s', ...
                strjoin(extractionInfoVal.featuresAllNaN, ', '));
        end
        XVal = XVal(validValMask, :);
        yVal = valLabelsRaw(validValMask);
        emitProgress(progressCb, 'Validation dataset after filtering: %d samples (real=%d, noise=%d).', ...
            numel(yVal), sum(yVal == 1), sum(yVal == 0));
        hasUsableValidation = ~isempty(yVal) && any(yVal == 1) && any(yVal == 0);
        if ~hasUsableValidation
            emitProgress(progressCb, ['Validation set is not usable after filtering (requires both classes). ', ...
                'Will train without sweep and skip validation metrics.']);
            doSweep = false;
        end

        validation.nSamples = numel(yVal);
        validation.nReal = sum(yVal == 1);
        validation.nNoise = sum(yVal == 0);
        validation.extractionInfo = extractionInfoVal;

        if doSweep
            sweepCfg = struct();
            sweepCfg.kernels = normalizeCellstr(opts.SweepKernels);
            sweepCfg.boxConstraints = opts.SweepBoxConstraints(:)';
            sweepCfg.kernelScales = normalizeMixedList(opts.SweepKernelScales);
            sweepCfg.polynomialOrders = opts.SweepPolynomialOrders(:)';
            sweepCfg.tieHandling = normalizeSweepTieHandling(opts.SweepTieHandling);
            sweepCfg.tiePromptCallback = opts.SweepTiePromptCallback;

            trainOptions = applyTrainingDefaults(opts.TrainingOptions);
            emitProgress(progressCb, 'Starting validation sweep across %d kernel(s), %d box value(s), %d scale value(s), %d poly order(s).', ...
                numel(sweepCfg.kernels), numel(sweepCfg.boxConstraints), numel(sweepCfg.kernelScales), numel(sweepCfg.polynomialOrders));
            [model, trainStats, normParams, bestParams, sweepResults] = runHyperparameterSweep( ...
                XTrain, yTrain, fixedTrainLabels, XVal, yVal, trainOptions, sweepCfg, verbose, progressCb);
            validation.sweepResults = sweepResults;
        else
            trainOptions = applyTrainingDefaults(opts.TrainingOptions);
            emitProgress(progressCb, 'Training SVM (no validation sweep)...');
            [model, trainStats, normParams] = snap_helpers.classification.trainClassifier(XTrain, yTrain, trainOptions);
            bestParams = trainOptions;
        end

        if hasUsableValidation
            [validation.metrics, validation.confusionMatrix] = evaluateOnValidation(model, normParams, XVal, yVal);
        else
            validation.metrics = struct('accuracy', NaN, 'precisionReal', NaN, 'recallReal', NaN, 'f1Real', NaN);
            validation.confusionMatrix = zeros(2, 2);
        end
    else
        trainOptions = applyTrainingDefaults(opts.TrainingOptions);
        emitProgress(progressCb, 'Training SVM (no validation set)...');
        [model, trainStats, normParams] = snap_helpers.classification.trainClassifier(XTrain, yTrain, trainOptions);
        bestParams = trainOptions;
    end

    if ~isfield(trainStats, 'success') || ~trainStats.success
        if isfield(trainStats, 'error')
            error('Training failed: %s', trainStats.error);
        else
            error('Training failed with unknown error.');
        end
    end

    metadata = struct();
    metadata.nImages = numel(unique(trainSources));
    metadata.nSamples = numel(yTrain);
    metadata.nRealLabeled = sum(yTrain == 1);
    metadata.nNoiseLabeled = sum(yTrain == 0);
    metadata.trainingSource = 'SNAP_train';
    metadata.exportFiles = exportFiles;
    metadata.labelFiles = labelFiles;
    metadata.matchDistance = opts.MatchDistance;
    metadata.bestParams = bestParams;
    metadata.sweepTieHandling = normalizeSweepTieHandling(opts.SweepTieHandling);
    metadata.trainingLabelsFixedDuringValidation = true;
    if hasValidationSet
        metadata.validation = validation;
    end

    success = snap_helpers.classification.saveClassifier( ...
        outputClassifierPath, model, selectedFeatures, featureInfo, trainStats, ...
        char(fittingMethod), metadata, customExpressions, normParams);

    if ~success
        error('Failed to save classifier to %s', outputClassifierPath);
    end

    out = struct();
    out.model = model;
    out.trainStats = trainStats;
    out.normParams = normParams;
    out.selectedFeatures = selectedFeatures;
    out.featureNames = featureNames;
    out.extractionInfoTrain = extractionInfoTrain;
    out.nSamples = numel(yTrain);
    out.nReal = sum(yTrain == 1);
    out.nNoise = sum(yTrain == 0);
    out.customExpressions = customExpressions;
    out.bestParams = bestParams;
    out.outputClassifierPath = outputClassifierPath;
    if hasValidationSet
        out.validation = validation;
    end

    if verbose
        fprintf('\nSNAP_train complete:\n');
        fprintf('  Train samples: %d (real=%d, noise=%d)\n', out.nSamples, out.nReal, out.nNoise);
        if hasValidationSet
            fprintf('  Validation samples: %d (real=%d, noise=%d)\n', validation.nSamples, validation.nReal, validation.nNoise);
            fprintf('  Validation F1 (real): %.3f\n', validation.metrics.f1Real);
        end
        if isfield(trainStats, 'cvAccuracy')
            fprintf('  CV accuracy: %.2f%%\n', 100 * trainStats.cvAccuracy);
        elseif isfield(trainStats, 'trainAccuracy')
            fprintf('  Train accuracy: %.2f%%\n', 100 * trainStats.trainAccuracy);
        end
        fprintf('  Features: %d\n', numel(selectedFeatures) + numel(customExpressions));
        fprintf('  Classifier: %s\n\n', outputClassifierPath);
    end
    emitProgress(progressCb, 'Classifier saved: %s', outputClassifierPath);

    if nargout > 0
        varargout{1} = out;
    end
end

function out = runTrainingUI()
% runTrainingUI - SNAP-style GUI for channel-aware classifier training

    fig = uifigure('Name', 'SNAP Train - Multi-Channel SVM Training', ...
        'Position', [50 50 1280 820]);

    mainGrid = uigridlayout(fig, [2, 2]);
    mainGrid.RowHeight = {'fit', '1x'};
    mainGrid.ColumnWidth = {'1.45x', '1x'};
    mainGrid.Padding = [10 10 10 10];
    mainGrid.RowSpacing = 10;
    mainGrid.ColumnSpacing = 10;

    setupPanel = uipanel(mainGrid, 'Title', 'Global Setup', 'FontWeight', 'bold');
    setupPanel.Layout.Row = 1;
    setupPanel.Layout.Column = [1 2];
    setupGrid = uigridlayout(setupPanel, [2, 6]);
    setupGrid.RowHeight = {'fit', 'fit'};
    setupGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit', '1x'};
    setupGrid.Padding = [6 6 6 6];
    setupGrid.ColumnSpacing = 8;

    uilabel(setupGrid, 'Text', 'Template Parameter File:', 'FontWeight', 'bold');
    paramPathEdit = uieditfield(setupGrid, 'text', 'Editable', 'off', ...
        'Placeholder', 'Optional: select parameter file to detect channel count');
    paramPathEdit.Layout.Column = 2;
    paramBrowseBtn = uibutton(setupGrid, 'Text', 'Browse');
    paramBrowseBtn.Layout.Column = 3;
    uilabel(setupGrid, 'Text', 'Detected Channel Tabs:', 'FontWeight', 'bold');
    detectedChannelsLabel = uilabel(setupGrid, 'Text', '1');
    detectedChannelsLabel.Layout.Column = 5;
    uilabel(setupGrid, 'Text', 'Configure channel-specific settings in tabs below.', ...
        'FontColor', [0.35 0.35 0.35]);
    channelCountPanel = uipanel(setupGrid, 'Title', 'Channel Tabs');
    channelCountPanel.Layout.Row = 2;
    channelCountPanel.Layout.Column = [1 3];
    channelCountGrid = uigridlayout(channelCountPanel, [1, 4]);
    channelCountGrid.ColumnWidth = {'fit', 'fit', 'fit', 'fit'};
    channelCountGrid.Padding = [4 4 4 4];
    channelCountGrid.ColumnSpacing = 6;
    uilabel(channelCountGrid, 'Text', 'Count:', 'FontWeight', 'bold');
    channelCountValueLabel = uilabel(channelCountGrid, 'Text', '1');
    decreaseChannelBtn = uibutton(channelCountGrid, 'Text', '-');
    increaseChannelBtn = uibutton(channelCountGrid, 'Text', '+');
    manualCountNote = uilabel(setupGrid, ...
        'Text', 'Use this to create tabs manually when no template parameter file is loaded.', ...
        'FontColor', [0.35 0.35 0.35]);
    manualCountNote.Layout.Row = 2;
    manualCountNote.Layout.Column = [4 6];

    channelPanel = uipanel(mainGrid, 'Title', 'Per-Channel Configuration', 'FontWeight', 'bold');
    channelPanel.Layout.Row = 2;
    channelPanel.Layout.Column = 1;
    channelGrid = uigridlayout(channelPanel, [1, 1]);
    channelGrid.Padding = [6 6 6 6];
    channelTabGroup = uitabgroup(channelGrid);

    rightPanel = uipanel(mainGrid, 'Title', 'Training Progress & Logs', 'FontWeight', 'bold');
    rightPanel.Layout.Row = 2;
    rightPanel.Layout.Column = 2;
    rightGrid = uigridlayout(rightPanel, [3, 1]);
    rightGrid.RowHeight = {'fit', '1x', 'fit'};
    rightGrid.RowSpacing = 8;
    rightGrid.Padding = [8 8 8 8];

    progressPanel = uipanel(rightGrid, 'Title', 'Training Progress', 'FontWeight', 'bold');
    progressGrid = uigridlayout(progressPanel, [2, 1]);
    progressGrid.RowHeight = {'fit', 'fit'};
    progressGrid.Padding = [6 6 6 6];
    progressGrid.RowSpacing = 4;
    progressStatusLabel = uilabel(progressGrid, 'Text', 'Ready to begin');
    progressStatusLabel.FontWeight = 'bold';
    progressStatusLabel.FontColor = [0 0 0];
    progressBarContainer = uipanel(progressGrid, 'BorderType', 'line', 'BorderWidth', 1);
    progressBarGrid = uigridlayout(progressBarContainer, [1, 2]);
    progressBarGrid.RowHeight = {20};
    progressBarGrid.ColumnWidth = {0, '1x'};
    progressBarGrid.Padding = [0 0 0 0];
    progressBarGrid.ColumnSpacing = 0;
    progressBarFill = uipanel(progressBarGrid, 'BorderType', 'none', 'BackgroundColor', [0.2 0.8 0.2]);
    progressBarFill.Layout.Row = 1;
    progressBarFill.Layout.Column = 1;
    progressBarEmpty = uipanel(progressBarGrid, 'BorderType', 'none', 'BackgroundColor', [0.9 0.9 0.9]);
    progressBarEmpty.Layout.Row = 1;
    progressBarEmpty.Layout.Column = 2;

    logText = uitextarea(rightGrid, 'Editable', 'off');
    logText.Value = {['SNAP_train ready. Load a template parameter file or set per-channel parameter files, ' ...
        'configure each channel tab, then click Train Selected Channels.']};

    actionPanel = uipanel(rightGrid, 'BorderType', 'none');
    actionGrid = uigridlayout(actionPanel, [1, 2]);
    actionGrid.ColumnWidth = {'1x', 'fit'};
    actionGrid.Padding = [0 0 0 0];
    trainBtn = uibutton(actionGrid, 'Text', 'Train Selected Channels', 'FontWeight', 'bold');
    closeBtn = uibutton(actionGrid, 'Text', 'Close', 'ButtonPushedFcn', @(~,~) close(fig));

    % Hidden table retained as a robust backend store for core per-channel fields.
    channelTable = uitable(fig, 'Visible', 'off');
    channelTable.ColumnName = {'Train', 'Channel Slot', 'Match Distance (voxels)', ...
        'Convert FIJI Coords', 'Parameter File', 'Training Directory', 'Validation Directory', 'Output Classifier'};
    channelTable.ColumnEditable = [true, false, true, true, false, false, false, true];
    channelTable.Data = {true, 'Channel 1', 2, false, '', '', '', 'classifier_ch1.mat'};

    handles = struct();
    handles.fig = fig;
    handles.paramPathEdit = paramPathEdit;
    handles.detectedChannelsLabel = detectedChannelsLabel;
    handles.channelCountValueLabel = channelCountValueLabel;
    handles.detectedNumChannels = 1;
    handles.channelTabGroup = channelTabGroup;
    handles.channelTabs = struct([]);
    handles.channelAdvancedConfigs = repmat(defaultChannelAdvancedConfig(), 1, 1);
    handles.channelFeatureConfigs = repmat(defaultChannelFeatureConfig(), 1, 1);
    handles.channelTable = channelTable;
    handles.selectedChannelRow = 1;
    handles.progressStatusLabel = progressStatusLabel;
    handles.progressBarGrid = progressBarGrid;
    handles.logText = logText;
    handles.loadedParams = struct();
    guidata(fig, handles);

    paramBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseParameterFile(fig);
    decreaseChannelBtn.ButtonPushedFcn = @(~,~) onAdjustChannelCount(fig, -1);
    increaseChannelBtn.ButtonPushedFcn = @(~,~) onAdjustChannelCount(fig, +1);
    trainBtn.ButtonPushedFcn = @(~,~) onTrainSelectedChannels(fig);

    updateChannelTableFromDetectedChannels(fig, 1, true);
    setTrainProgressVisual(fig, 0, 'Ready to begin', [0 0 0]);

    out = struct('status', 'launched', 'figure', fig);
end

function onBrowseParameterFile(fig)
    handles = guidata(fig);
    [file, path] = uigetfile({'*.mat', 'MAT files (*.mat)'}, 'Select SNAP parameter file');
    if isequal(file, 0)
        return;
    end
    paramPath = fullfile(path, file);
    handles.paramPathEdit.Value = paramPath;

    [params, numChannels, errMsg] = loadTrainingParameters(paramPath);
    if ~isempty(errMsg)
        guidata(fig, handles);
        uialert(fig, errMsg, 'Parameter File Error');
        return;
    end

    handles.loadedParams = params;
    if isfield(handles, 'channelCountValueLabel') && isgraphics(handles.channelCountValueLabel)
        handles.channelCountValueLabel.Text = num2str(numChannels);
    end
    guidata(fig, handles);
    updateChannelTableFromDetectedChannels(fig, numChannels, true);
    refreshTrainingFeatureSummary(fig, 1);
    appendTrainLog(fig, sprintf('Loaded template parameter file: %s (detected channels=%d)', file, numChannels));
end

function onAdjustChannelCount(fig, delta)
    handles = guidata(fig);
    currentN = max(1, round(getfieldwithdefault(handles, 'detectedNumChannels', 1)));
    if nargin < 2 || ~isfinite(delta)
        delta = 0;
    end
    n = max(1, min(12, currentN + round(delta)));
    updateChannelTableFromDetectedChannels(fig, n, false);
    appendTrainLog(fig, sprintf('Adjusted channel tabs to %d channel(s).', n));
end

function updateChannelTableFromDetectedChannels(fig, detectedChannels, resetOutputs)
    handles = guidata(fig);
    handles = syncTabsToState(handles);

    n = max(1, round(detectedChannels));
    handles.detectedNumChannels = n;
    handles.detectedChannelsLabel.Text = num2str(n);
    if isfield(handles, 'channelCountValueLabel') && isgraphics(handles.channelCountValueLabel)
        handles.channelCountValueLabel.Text = num2str(n);
    end

    oldFeatureCfgs = repmat(defaultChannelFeatureConfig(), 1, 0);
    if isfield(handles, 'channelFeatureConfigs') && ~isempty(handles.channelFeatureConfigs)
        oldFeatureCfgs = handles.channelFeatureConfigs;
    end

    oldAdvancedCfgs = repmat(defaultChannelAdvancedConfig(), 1, 0);
    if isfield(handles, 'channelAdvancedConfigs') && ~isempty(handles.channelAdvancedConfigs)
        oldAdvancedCfgs = handles.channelAdvancedConfigs;
    end

    oldData = handles.channelTable.Data;
    templateParamPath = strtrim(char(string(handles.paramPathEdit.Value)));
    newData = cell(n, 8);

    for ch = 1:n
        channelName = sprintf('Channel %d', ch);
        matchDistance = 2;
        convertFiji = false;
        paramFile = templateParamPath;
        trainDir = '';
        valDir = '';
        defaultOutput = sprintf('classifier_ch%d.mat', ch);
        trainFlag = true;
        outputName = defaultOutput;

        if ~isempty(oldData)
            rowIdx = find(strcmp(oldData(:, 2), channelName), 1);
            if ~isempty(rowIdx)
                trainFlag = isSelectedChannelFlag(oldData{rowIdx, 1});
                if size(oldData, 2) >= 8
                    matchDistance = normalizeScalarNumeric(oldData{rowIdx, 3}, 2);
                    if ~isfinite(matchDistance) || matchDistance < 0
                        matchDistance = 2;
                    end
                    convertFiji = isSelectedChannelFlag(oldData{rowIdx, 4});
                    paramFile = strtrim(char(string(oldData{rowIdx, 5})));
                    trainDir = char(string(oldData{rowIdx, 6}));
                    valDir = char(string(oldData{rowIdx, 7}));
                    if ~resetOutputs && ~isempty(oldData{rowIdx, 8})
                        outputName = char(string(oldData{rowIdx, 8}));
                    end
                elseif size(oldData, 2) >= 7
                    % Previous layout before per-channel parameter-file field.
                    matchDistance = normalizeScalarNumeric(oldData{rowIdx, 3}, 2);
                    if ~isfinite(matchDistance) || matchDistance < 0
                        matchDistance = 2;
                    end
                    convertFiji = isSelectedChannelFlag(oldData{rowIdx, 4});
                    trainDir = char(string(oldData{rowIdx, 5}));
                    valDir = char(string(oldData{rowIdx, 6}));
                    if ~resetOutputs && ~isempty(oldData{rowIdx, 7})
                        outputName = char(string(oldData{rowIdx, 7}));
                    end
                else
                    % Backward compatibility with older 6-column table layout.
                    if size(oldData, 2) >= 3
                        convertFiji = isSelectedChannelFlag(oldData{rowIdx, 3});
                    end
                    if size(oldData, 2) >= 4
                        trainDir = char(string(oldData{rowIdx, 4}));
                    end
                    if size(oldData, 2) >= 5
                        valDir = char(string(oldData{rowIdx, 5}));
                    end
                    if ~resetOutputs && size(oldData, 2) >= 6 && ~isempty(oldData{rowIdx, 6})
                        outputName = char(string(oldData{rowIdx, 6}));
                    end
                end
            end
        end
        if isempty(strtrim(paramFile))
            paramFile = templateParamPath;
        end
        if resetOutputs && ~isempty(templateParamPath)
            paramFile = templateParamPath;
        end

        newData{ch, 1} = trainFlag;
        newData{ch, 2} = channelName;
        newData{ch, 3} = matchDistance;
        newData{ch, 4} = convertFiji;
        newData{ch, 5} = paramFile;
        newData{ch, 6} = trainDir;
        newData{ch, 7} = valDir;
        newData{ch, 8} = outputName;
    end

    newFeatureCfgs = repmat(defaultChannelFeatureConfig(), 1, n);
    for ch = 1:n
        if ch <= numel(oldFeatureCfgs)
            newFeatureCfgs(ch) = normalizeChannelFeatureConfig(oldFeatureCfgs(ch));
        end
    end

    newAdvancedCfgs = repmat(defaultChannelAdvancedConfig(), 1, n);
    for ch = 1:n
        if ch <= numel(oldAdvancedCfgs)
            newAdvancedCfgs(ch) = normalizeChannelAdvancedConfig(oldAdvancedCfgs(ch));
        end
        if resetOutputs
            newAdvancedCfgs(ch).outputDirectory = '';
        end
    end

    handles.channelFeatureConfigs = newFeatureCfgs;
    handles.channelAdvancedConfigs = newAdvancedCfgs;
    handles.channelTable.Data = newData;
    handles.selectedChannelRow = max(1, min(n, handles.selectedChannelRow));
    guidata(fig, handles);

    rebuildChannelTabs(fig);
    updateOptimizationControlState(fig);
end

function updateOptimizationControlState(fig, channelIdx)
    handles = guidata(fig);
    if ~isfield(handles, 'channelTabs') || isempty(handles.channelTabs)
        return;
    end

    if nargin < 2 || isempty(channelIdx)
        idxList = 1:numel(handles.channelTabs);
    else
        idxList = unique(round(channelIdx(:)'));
    end

    for ch = idxList
        if ch < 1 || ch > numel(handles.channelTabs)
            continue;
        end
        tab = handles.channelTabs(ch);
        if ~isfield(tab, 'optimizeCheck') || ~isgraphics(tab.optimizeCheck)
            continue;
        end

        useOptimization = logical(tab.optimizeCheck.Value);
        if useOptimization
            manualState = 'off';
            sweepState = 'on';
        else
            manualState = 'on';
            sweepState = 'off';
        end

        manualControls = { ...
            tab.kernelDrop, tab.boxConstraintInput, tab.kernelScaleInput, ...
            tab.polyOrderInput, tab.standardizeCheck, tab.crossValidateCheck, ...
            tab.kFoldInput, tab.balanceClassCheck ...
        };
        sweepControls = { ...
            tab.sweepKernelsInput, tab.sweepBoxInput, tab.sweepScaleInput, ...
            tab.sweepPolyInput, tab.sweepTieHandlingDrop, tab.sweepReportCheck ...
        };

        for i = 1:numel(manualControls)
            if isgraphics(manualControls{i})
                manualControls{i}.Enable = manualState;
            end
        end
        for i = 1:numel(sweepControls)
            if isgraphics(sweepControls{i})
                sweepControls{i}.Enable = sweepState;
            end
        end
    end

    guidata(fig, handles);
end

function cfg = defaultChannelAdvancedConfig()
    cfg = struct( ...
        'outputDirectory', '', ...
        'optimizeWithSweep', true, ...
        'includeSweepReport', true, ...
        'kernelFunction', 'rbf', ...
        'boxConstraint', 1, ...
        'kernelScale', 'auto', ...
        'polynomialOrder', 3, ...
        'standardize', true, ...
        'crossValidate', true, ...
        'kFold', 5, ...
        'balanceClassBins', true, ...
        'sweepKernels', 'linear, rbf, polynomial', ...
        'sweepBoxConstraints', '0.1, 1, 10', ...
        'sweepKernelScales', 'auto, 0.5, 1, 2', ...
        'sweepPolynomialOrders', '2, 3, 4', ...
        'sweepTieHandlingLabel', 'Prefer simpler model for ties');
end

function cfg = normalizeChannelAdvancedConfig(rawCfg)
    cfg = defaultChannelAdvancedConfig();
    if ~isstruct(rawCfg)
        return;
    end

    fields = fieldnames(cfg);
    for i = 1:numel(fields)
        fn = fields{i};
        if isfield(rawCfg, fn)
            cfg.(fn) = rawCfg.(fn);
        end
    end

    cfg.outputDirectory = strtrim(char(string(cfg.outputDirectory)));
    cfg.optimizeWithSweep = logical(cfg.optimizeWithSweep);
    cfg.includeSweepReport = logical(cfg.includeSweepReport);
    cfg.kernelFunction = char(string(cfg.kernelFunction));
    cfg.kernelScale = char(string(cfg.kernelScale));
    cfg.polynomialOrder = max(2, round(normalizeScalarNumeric(cfg.polynomialOrder, 3)));
    cfg.boxConstraint = normalizeScalarNumeric(cfg.boxConstraint, 1);
    if ~isfinite(cfg.boxConstraint) || cfg.boxConstraint <= 0
        cfg.boxConstraint = 1;
    end
    cfg.standardize = logical(cfg.standardize);
    cfg.crossValidate = logical(cfg.crossValidate);
    cfg.kFold = max(2, round(normalizeScalarNumeric(cfg.kFold, 5)));
    cfg.balanceClassBins = logical(cfg.balanceClassBins);
    cfg.sweepKernels = char(string(cfg.sweepKernels));
    cfg.sweepBoxConstraints = char(string(cfg.sweepBoxConstraints));
    cfg.sweepKernelScales = char(string(cfg.sweepKernelScales));
    cfg.sweepPolynomialOrders = char(string(cfg.sweepPolynomialOrders));
    cfg.sweepTieHandlingLabel = char(string(cfg.sweepTieHandlingLabel));
    if isempty(cfg.sweepTieHandlingLabel)
        cfg.sweepTieHandlingLabel = 'Prefer simpler model for ties';
    end
end

function cfg = getChannelAdvancedConfig(handles, channelIdx)
    cfg = defaultChannelAdvancedConfig();
    if ~isfield(handles, 'channelAdvancedConfigs') || isempty(handles.channelAdvancedConfigs)
        return;
    end
    idx = max(1, round(channelIdx));
    if idx <= numel(handles.channelAdvancedConfigs)
        cfg = normalizeChannelAdvancedConfig(handles.channelAdvancedConfigs(idx));
    end
end

function handles = setChannelAdvancedConfig(handles, channelIdx, cfg)
    idx = max(1, round(channelIdx));
    cfg = normalizeChannelAdvancedConfig(cfg);
    if ~isfield(handles, 'channelAdvancedConfigs') || isempty(handles.channelAdvancedConfigs)
        handles.channelAdvancedConfigs = repmat(defaultChannelAdvancedConfig(), 1, idx);
    elseif numel(handles.channelAdvancedConfigs) < idx
        nExisting = numel(handles.channelAdvancedConfigs);
        handles.channelAdvancedConfigs(nExisting+1:idx) = ...
            repmat(defaultChannelAdvancedConfig(), 1, idx - nExisting);
    end
    handles.channelAdvancedConfigs(idx) = cfg;
end

function handles = syncTabsToState(handles)
    if ~isfield(handles, 'channelTabs') || isempty(handles.channelTabs) || ...
            ~isfield(handles, 'channelTable') || ~isgraphics(handles.channelTable)
        return;
    end

    data = handles.channelTable.Data;
    if isempty(data)
        return;
    end

    for ch = 1:min(numel(handles.channelTabs), size(data, 1))
        tab = handles.channelTabs(ch);
        if ~isstruct(tab) || ~isfield(tab, 'tab') || ~isgraphics(tab.tab)
            continue;
        end

        data{ch, 1} = logical(tab.trainEnableCheck.Value);
        data{ch, 2} = sprintf('Channel %d', ch);
        data{ch, 3} = tab.matchDistanceInput.Value;
        data{ch, 4} = logical(tab.convertFijiCheck.Value);
        data{ch, 5} = strtrim(char(string(tab.paramFileEdit.Value)));
        data{ch, 6} = strtrim(char(string(tab.trainDirEdit.Value)));
        data{ch, 7} = strtrim(char(string(tab.valDirEdit.Value)));
        outName = strtrim(char(string(tab.outputNameEdit.Value)));
        if isempty(outName)
            outName = sprintf('classifier_ch%d.mat', ch);
        end
        data{ch, 8} = outName;

        adv = getChannelAdvancedConfig(handles, ch);
        adv.outputDirectory = strtrim(char(string(tab.outputDirEdit.Value)));
        adv.optimizeWithSweep = logical(tab.optimizeCheck.Value);
        adv.includeSweepReport = logical(tab.sweepReportCheck.Value);
        adv.kernelFunction = char(string(tab.kernelDrop.Value));
        adv.boxConstraint = tab.boxConstraintInput.Value;
        adv.kernelScale = char(string(tab.kernelScaleInput.Value));
        adv.polynomialOrder = tab.polyOrderInput.Value;
        adv.standardize = logical(tab.standardizeCheck.Value);
        adv.crossValidate = logical(tab.crossValidateCheck.Value);
        adv.kFold = tab.kFoldInput.Value;
        adv.balanceClassBins = logical(tab.balanceClassCheck.Value);
        adv.sweepKernels = char(string(tab.sweepKernelsInput.Value));
        adv.sweepBoxConstraints = char(string(tab.sweepBoxInput.Value));
        adv.sweepKernelScales = char(string(tab.sweepScaleInput.Value));
        adv.sweepPolynomialOrders = char(string(tab.sweepPolyInput.Value));
        adv.sweepTieHandlingLabel = char(string(tab.sweepTieHandlingDrop.Value));
        handles = setChannelAdvancedConfig(handles, ch, adv);
    end

    handles.channelTable.Data = data;
end

function rebuildChannelTabs(fig)
    handles = guidata(fig);
    if ~isfield(handles, 'channelTabGroup') || ~isgraphics(handles.channelTabGroup)
        return;
    end
    if ~isfield(handles, 'channelTable') || ~isgraphics(handles.channelTable)
        return;
    end

    handles = syncTabsToState(handles);
    data = handles.channelTable.Data;
    if isempty(data)
        guidata(fig, handles);
        return;
    end

    selectedChannel = max(1, min(size(data, 1), handles.selectedChannelRow));
    existingTabs = handles.channelTabGroup.Children;
    for i = 1:numel(existingTabs)
        if isgraphics(existingTabs(i))
            delete(existingTabs(i));
        end
    end

    n = size(data, 1);
    channelTabs = repmat(emptyChannelTabRecord(), 1, n);
    hasLoadedParams = isfield(handles, 'loadedParams') && ...
        isstruct(handles.loadedParams) && ~isempty(fieldnames(handles.loadedParams));

    for ch = 1:n
        coreTrain = isSelectedChannelFlag(data{ch, 1});
        coreMatchDist = normalizeScalarNumeric(data{ch, 3}, 2);
        if ~isfinite(coreMatchDist) || coreMatchDist < 0
            coreMatchDist = 2;
        end
        coreConvertFiji = isSelectedChannelFlag(data{ch, 4});
        coreParamPath = '';
        coreTrainDir = '';
        coreValDir = '';
        coreOutName = '';
        if size(data, 2) >= 8
            coreParamPath = char(string(data{ch, 5}));
            coreTrainDir = char(string(data{ch, 6}));
            coreValDir = char(string(data{ch, 7}));
            coreOutName = char(string(data{ch, 8}));
        else
            coreTrainDir = char(string(data{ch, 5}));
            coreValDir = char(string(data{ch, 6}));
            coreOutName = char(string(data{ch, 7}));
            coreParamPath = strtrim(char(string(handles.paramPathEdit.Value)));
        end
        if isempty(strtrim(coreParamPath))
            coreParamPath = strtrim(char(string(handles.paramPathEdit.Value)));
        end
        if isempty(strtrim(coreOutName))
            coreOutName = sprintf('classifier_ch%d.mat', ch);
        end

        adv = getChannelAdvancedConfig(handles, ch);

        tab = uitab(handles.channelTabGroup, 'Title', sprintf('Channel %d', ch));
        if isprop(tab, 'Scrollable')
            tab.Scrollable = 'on';
        end
        tabRootGrid = uigridlayout(tab, [1, 1]);
        tabRootGrid.RowHeight = {'1x'};
        tabRootGrid.ColumnWidth = {'1x'};
        tabRootGrid.Padding = [0 0 0 0];
        tabRootGrid.RowSpacing = 0;
        tabRootGrid.ColumnSpacing = 0;
        tabContainer = uipanel(tabRootGrid, 'BorderType', 'none');
        if isprop(tabContainer, 'Scrollable')
            tabContainer.Scrollable = 'on';
        end
        tabGrid = uigridlayout(tabContainer, [7, 1]);
        tabGrid.RowHeight = {'fit', 'fit', 'fit', 180, 'fit', 'fit', 'fit'};
        tabGrid.Padding = [8 8 8 8];
        tabGrid.RowSpacing = 7;

        generalPanel = uipanel(tabGrid, 'Title', 'Core');
        generalGrid = uigridlayout(generalPanel, [2, 3]);
        generalGrid.ColumnWidth = {'fit', 'fit', 120};
        generalGrid.RowHeight = {'fit', 'fit'};
        generalGrid.Padding = [6 6 6 6];
        trainEnableCheck = uicheckbox(generalGrid, 'Text', 'Train this channel', 'Value', coreTrain);
        trainEnableCheck.Layout.Row = 1;
        trainEnableCheck.Layout.Column = 1;
        matchDistanceLabel = uilabel(generalGrid, 'Text', 'Match Distance (voxels):');
        matchDistanceLabel.Layout.Row = 1;
        matchDistanceLabel.Layout.Column = 2;
        matchDistanceInput = uieditfield(generalGrid, 'numeric', 'Value', coreMatchDist, 'LowerLimit', 0);
        matchDistanceInput.Layout.Row = 1;
        matchDistanceInput.Layout.Column = 3;
        convertFijiCheck = uicheckbox(generalGrid, 'Text', 'Convert FIJI Coords', 'Value', coreConvertFiji);
        convertFijiCheck.Layout.Row = 2;
        convertFijiCheck.Layout.Column = [1 3];

        dataPanel = uipanel(tabGrid, 'Title', 'Channel Inputs');
        dataGrid = uigridlayout(dataPanel, [3, 3]);
        dataGrid.ColumnWidth = {'fit', '1x', 'fit'};
        dataGrid.RowHeight = {'fit', 'fit', 'fit'};
        dataGrid.Padding = [6 6 6 6];
        paramLabel = uilabel(dataGrid, 'Text', 'Parameter File:');
        paramLabel.Layout.Row = 1;
        paramLabel.Layout.Column = 1;
        paramFileEdit = uieditfield(dataGrid, 'text', 'Editable', 'off', 'Value', coreParamPath);
        paramFileEdit.Layout.Row = 1;
        paramFileEdit.Layout.Column = 2;
        paramBrowseBtn = uibutton(dataGrid, 'Text', 'Browse');
        paramBrowseBtn.Layout.Row = 1;
        paramBrowseBtn.Layout.Column = 3;
        trainDirLabel = uilabel(dataGrid, 'Text', 'Training Directory:');
        trainDirLabel.Layout.Row = 2;
        trainDirLabel.Layout.Column = 1;
        trainDirEdit = uieditfield(dataGrid, 'text', 'Editable', 'off', 'Value', coreTrainDir);
        trainDirEdit.Layout.Row = 2;
        trainDirEdit.Layout.Column = 2;
        trainBrowseBtn = uibutton(dataGrid, 'Text', 'Browse');
        trainBrowseBtn.Layout.Row = 2;
        trainBrowseBtn.Layout.Column = 3;
        valDirLabel = uilabel(dataGrid, 'Text', 'Validation Directory:');
        valDirLabel.Layout.Row = 3;
        valDirLabel.Layout.Column = 1;
        valDirEdit = uieditfield(dataGrid, 'text', 'Editable', 'off', 'Value', coreValDir);
        valDirEdit.Layout.Row = 3;
        valDirEdit.Layout.Column = 2;
        valBrowseBtn = uibutton(dataGrid, 'Text', 'Browse');
        valBrowseBtn.Layout.Row = 3;
        valBrowseBtn.Layout.Column = 3;

        outputPanel = uipanel(tabGrid, 'Title', 'Output');
        outputGrid = uigridlayout(outputPanel, [2, 3]);
        outputGrid.ColumnWidth = {'fit', '1x', 'fit'};
        outputGrid.RowHeight = {'fit', 'fit'};
        outputGrid.Padding = [6 6 6 6];
        uilabel(outputGrid, 'Text', 'Output Directory:');
        outputDirEdit = uieditfield(outputGrid, 'text', 'Editable', 'off', 'Value', adv.outputDirectory);
        outputDirEdit.Layout.Column = 2;
        outputDirBrowseBtn = uibutton(outputGrid, 'Text', 'Browse');
        outputDirBrowseBtn.Layout.Column = 3;
        uilabel(outputGrid, 'Text', 'Classifier File:');
        outputNameEdit = uieditfield(outputGrid, 'text', 'Value', coreOutName);
        outputNameEdit.Layout.Column = [2 3];

        featurePanel = uipanel(tabGrid, 'Title', 'Features & Custom Expressions');
        featureGrid = uigridlayout(featurePanel, [4, 1]);
        featureGrid.RowHeight = {'fit', 'fit', 'fit', '1x'};
        featureGrid.Padding = [6 6 6 6];
        featureGrid.RowSpacing = 4;
        hasChannelParamFile = ~isempty(strtrim(coreParamPath)) && exist(coreParamPath, 'file') == 2;
        controlsEnabled = ternaryState(hasLoadedParams || hasChannelParamFile, 'on', 'off');
        selectFeaturesBtn = uibutton(featureGrid, 'Text', sprintf('Select Features (Channel %d)...', ch), ...
            'Enable', controlsEnabled);
        loadFeaturePackBtn = uibutton(featureGrid, 'Text', sprintf('Load Expression Pack (Channel %d)...', ch), ...
            'Enable', controlsEnabled);
        featureCountLabel = uilabel(featureGrid, 'Text', sprintf('Channel %d: AUTO (all non-position features)', ch));
        featureCountLabel.FontColor = [0.45 0.45 0.45];
        featureListArea = uitextarea(featureGrid, 'Editable', 'off', ...
            'Value', {'AUTO: all non-position features', 'Custom expressions: none'});

        optimizePanel = uipanel(tabGrid, 'Title', 'Validation Strategy');
        optimizeGrid = uigridlayout(optimizePanel, [2, 1]);
        optimizeGrid.ColumnWidth = {'1x'};
        optimizeGrid.RowHeight = {'fit', 'fit'};
        optimizeGrid.Padding = [6 6 6 6];
        optimizeCheck = uicheckbox(optimizeGrid, 'Text', 'Optimize with validation sweep', ...
            'Value', logical(adv.optimizeWithSweep));
        sweepReportCheck = uicheckbox(optimizeGrid, 'Text', 'Include sweep report', ...
            'Value', logical(adv.includeSweepReport));
        sweepReportCheck.Layout.Row = 2;

        svmPanel = uipanel(tabGrid, 'Title', 'Manual SVM Hyperparameters');
        svmGrid = uigridlayout(svmPanel, [3, 4]);
        svmGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
        svmGrid.RowHeight = {'fit', 'fit', 'fit'};
        svmGrid.Padding = [6 6 6 6];
        uilabel(svmGrid, 'Text', 'Kernel:');
        kernelDrop = uidropdown(svmGrid, 'Items', {'rbf', 'linear', 'polynomial'}, ...
            'Value', char(string(adv.kernelFunction)));
        kernelDrop.Layout.Column = 2;
        uilabel(svmGrid, 'Text', 'Box Constraint (C):');
        boxConstraintInput = uieditfield(svmGrid, 'numeric', ...
            'Value', normalizeScalarNumeric(adv.boxConstraint, 1), 'LowerLimit', 1e-6);
        boxConstraintInput.Layout.Column = 4;
        uilabel(svmGrid, 'Text', 'Kernel Scale:');
        kernelScaleInput = uieditfield(svmGrid, 'text', 'Value', char(string(adv.kernelScale)));
        kernelScaleInput.Layout.Column = 2;
        uilabel(svmGrid, 'Text', 'Polynomial Order:');
        polyOrderInput = uieditfield(svmGrid, 'numeric', ...
            'Value', max(2, round(normalizeScalarNumeric(adv.polynomialOrder, 3))), ...
            'LowerLimit', 2, 'RoundFractionalValues', true);
        polyOrderInput.Layout.Column = 4;
        standardizeCheck = uicheckbox(svmGrid, 'Text', 'Z-score normalize', 'Value', logical(adv.standardize));
        standardizeCheck.Layout.Row = 3;
        standardizeCheck.Layout.Column = 1;
        crossValidateCheck = uicheckbox(svmGrid, 'Text', 'Cross-validate', 'Value', logical(adv.crossValidate));
        crossValidateCheck.Layout.Row = 3;
        crossValidateCheck.Layout.Column = 2;
        kFoldInput = uieditfield(svmGrid, 'numeric', ...
            'Value', max(2, round(normalizeScalarNumeric(adv.kFold, 5))), ...
            'LowerLimit', 2, 'RoundFractionalValues', true);
        kFoldInput.Layout.Row = 3;
        kFoldInput.Layout.Column = 3;
        balanceClassCheck = uicheckbox(svmGrid, 'Text', 'Balance class bins', 'Value', logical(adv.balanceClassBins));
        balanceClassCheck.Layout.Row = 3;
        balanceClassCheck.Layout.Column = 4;

        sweepPanel = uipanel(tabGrid, 'Title', 'Validation Sweep Grid');
        sweepGrid = uigridlayout(sweepPanel, [5, 2]);
        sweepGrid.ColumnWidth = {'fit', '1x'};
        sweepGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
        sweepGrid.Padding = [6 6 6 6];
        uilabel(sweepGrid, 'Text', 'Kernels:');
        sweepKernelsInput = uieditfield(sweepGrid, 'text', 'Value', char(string(adv.sweepKernels)));
        uilabel(sweepGrid, 'Text', 'Box Constraints:');
        sweepBoxInput = uieditfield(sweepGrid, 'text', 'Value', char(string(adv.sweepBoxConstraints)));
        uilabel(sweepGrid, 'Text', 'Kernel Scales:');
        sweepScaleInput = uieditfield(sweepGrid, 'text', 'Value', char(string(adv.sweepKernelScales)));
        uilabel(sweepGrid, 'Text', 'Polynomial Orders:');
        sweepPolyInput = uieditfield(sweepGrid, 'text', 'Value', char(string(adv.sweepPolynomialOrders)));
        uilabel(sweepGrid, 'Text', 'Tie Handling:');
        sweepTieHandlingDrop = uidropdown(sweepGrid, ...
            'Items', {'Prefer simpler model for ties', 'Prompt me to choose tied best settings', 'Keep first (legacy)'}, ...
            'Value', char(string(adv.sweepTieHandlingLabel)));

        trainEnableCheck.ValueChangedFcn = @(~,~) onChannelTabCoreChanged(fig, ch);
        matchDistanceInput.ValueChangedFcn = @(~,~) onChannelTabCoreChanged(fig, ch);
        convertFijiCheck.ValueChangedFcn = @(~,~) onChannelTabCoreChanged(fig, ch);
        outputNameEdit.ValueChangedFcn = @(~,~) onChannelTabCoreChanged(fig, ch);
        paramFileEdit.ValueChangedFcn = @(~,~) onChannelTabCoreChanged(fig, ch);

        paramBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseChannelDirFromTab(fig, ch, 'param');
        trainBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseChannelDirFromTab(fig, ch, 'train');
        valBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseChannelDirFromTab(fig, ch, 'validation');
        outputDirBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseChannelDirFromTab(fig, ch, 'output');

        optimizeCheck.ValueChangedFcn = @(~,~) onChannelOptimizationModeChanged(fig, ch);
        sweepReportCheck.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        kernelDrop.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        boxConstraintInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        kernelScaleInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        polyOrderInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        standardizeCheck.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        crossValidateCheck.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        kFoldInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        balanceClassCheck.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        sweepKernelsInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        sweepBoxInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        sweepScaleInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        sweepPolyInput.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);
        sweepTieHandlingDrop.ValueChangedFcn = @(~,~) onChannelTabAdvancedChanged(fig, ch);

        selectFeaturesBtn.ButtonPushedFcn = @(~,~) onSelectTrainingFeaturesForChannel(fig, ch);
        loadFeaturePackBtn.ButtonPushedFcn = @(~,~) onLoadExpressionPackForChannel(fig, ch);

        channelTabs(ch).tab = tab;
        channelTabs(ch).trainEnableCheck = trainEnableCheck;
        channelTabs(ch).matchDistanceInput = matchDistanceInput;
        channelTabs(ch).convertFijiCheck = convertFijiCheck;
        channelTabs(ch).paramFileEdit = paramFileEdit;
        channelTabs(ch).trainDirEdit = trainDirEdit;
        channelTabs(ch).valDirEdit = valDirEdit;
        channelTabs(ch).outputDirEdit = outputDirEdit;
        channelTabs(ch).outputNameEdit = outputNameEdit;
        channelTabs(ch).optimizeCheck = optimizeCheck;
        channelTabs(ch).sweepReportCheck = sweepReportCheck;
        channelTabs(ch).kernelDrop = kernelDrop;
        channelTabs(ch).boxConstraintInput = boxConstraintInput;
        channelTabs(ch).kernelScaleInput = kernelScaleInput;
        channelTabs(ch).polyOrderInput = polyOrderInput;
        channelTabs(ch).standardizeCheck = standardizeCheck;
        channelTabs(ch).crossValidateCheck = crossValidateCheck;
        channelTabs(ch).kFoldInput = kFoldInput;
        channelTabs(ch).balanceClassCheck = balanceClassCheck;
        channelTabs(ch).sweepKernelsInput = sweepKernelsInput;
        channelTabs(ch).sweepBoxInput = sweepBoxInput;
        channelTabs(ch).sweepScaleInput = sweepScaleInput;
        channelTabs(ch).sweepPolyInput = sweepPolyInput;
        channelTabs(ch).sweepTieHandlingDrop = sweepTieHandlingDrop;
        channelTabs(ch).selectFeaturesBtn = selectFeaturesBtn;
        channelTabs(ch).loadFeaturePackBtn = loadFeaturePackBtn;
        channelTabs(ch).featureCountLabel = featureCountLabel;
        channelTabs(ch).featureListArea = featureListArea;
    end

    handles.channelTabs = channelTabs;
    handles.channelTabGroup.SelectionChangedFcn = @(~,evt) onChannelTabSelectionChanged(fig, evt);
    guidata(fig, handles);

    if selectedChannel <= numel(channelTabs)
        try
            handles.channelTabGroup.SelectedTab = channelTabs(selectedChannel).tab;
        catch
            % No-op if setting selected tab fails on some MATLAB builds.
        end
    end

    refreshTrainingFeatureSummary(fig);
    updateOptimizationControlState(fig);
end

function tabRec = emptyChannelTabRecord()
    tabRec = struct( ...
        'tab', [], ...
        'trainEnableCheck', [], ...
        'matchDistanceInput', [], ...
        'convertFijiCheck', [], ...
        'paramFileEdit', [], ...
        'trainDirEdit', [], ...
        'valDirEdit', [], ...
        'outputDirEdit', [], ...
        'outputNameEdit', [], ...
        'optimizeCheck', [], ...
        'sweepReportCheck', [], ...
        'kernelDrop', [], ...
        'boxConstraintInput', [], ...
        'kernelScaleInput', [], ...
        'polyOrderInput', [], ...
        'standardizeCheck', [], ...
        'crossValidateCheck', [], ...
        'kFoldInput', [], ...
        'balanceClassCheck', [], ...
        'sweepKernelsInput', [], ...
        'sweepBoxInput', [], ...
        'sweepScaleInput', [], ...
        'sweepPolyInput', [], ...
        'sweepTieHandlingDrop', [], ...
        'selectFeaturesBtn', [], ...
        'loadFeaturePackBtn', [], ...
        'featureCountLabel', [], ...
        'featureListArea', []);
end

function out = ternaryState(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end

function onChannelTabSelectionChanged(fig, eventData)
    handles = guidata(fig);
    if nargin < 2 || isempty(eventData) || ~isfield(handles, 'channelTabs') || isempty(handles.channelTabs)
        return;
    end
    selectedTab = [];
    if isprop(eventData, 'NewValue')
        selectedTab = eventData.NewValue;
    elseif isstruct(eventData) && isfield(eventData, 'NewValue')
        selectedTab = eventData.NewValue;
    end
    if isempty(selectedTab)
        return;
    end
    for ch = 1:numel(handles.channelTabs)
        tab = handles.channelTabs(ch);
        if isfield(tab, 'tab') && isequal(tab.tab, selectedTab)
            handles.selectedChannelRow = ch;
            guidata(fig, handles);
            refreshTrainingFeatureSummary(fig, ch);
            return;
        end
    end
end

function onChannelTabCoreChanged(fig, channelIdx)
    handles = guidata(fig);
    handles = syncTabsToState(handles);
    handles.selectedChannelRow = max(1, round(channelIdx));
    guidata(fig, handles);
end

function onChannelTabAdvancedChanged(fig, channelIdx)
    handles = guidata(fig);
    handles = syncTabsToState(handles);
    handles.selectedChannelRow = max(1, round(channelIdx));
    guidata(fig, handles);
end

function onChannelOptimizationModeChanged(fig, channelIdx)
    onChannelTabAdvancedChanged(fig, channelIdx);
    updateOptimizationControlState(fig, channelIdx);
end

function onBrowseChannelDirFromTab(fig, channelIdx, dirType)
    handles = guidata(fig);
    if ~isfield(handles, 'channelTabs') || channelIdx < 1 || channelIdx > numel(handles.channelTabs)
        return;
    end
    tab = handles.channelTabs(channelIdx);

    dirType = lower(char(string(dirType)));
    switch dirType
        case 'param'
            existingParam = strtrim(char(string(tab.paramFileEdit.Value)));
            startDir = pwd;
            if ~isempty(existingParam)
                if exist(existingParam, 'file') == 2
                    startDir = fileparts(existingParam);
                elseif isfolder(existingParam)
                    startDir = existingParam;
                end
            end
            [file, path] = uigetfile({'*.mat', 'MAT files (*.mat)'}, ...
                sprintf('Select parameter file for Channel %d', channelIdx), startDir);
            if isequal(file, 0)
                return;
            end
            selectedPath = fullfile(path, file);
            [~, ~, errMsg] = loadTrainingParameters(selectedPath);
            if ~isempty(errMsg)
                uialert(fig, errMsg, 'Parameter File Error');
                return;
            end
            tab.paramFileEdit.Value = selectedPath;
            if isfield(tab, 'selectFeaturesBtn') && isgraphics(tab.selectFeaturesBtn)
                tab.selectFeaturesBtn.Enable = 'on';
            end
            if isfield(tab, 'loadFeaturePackBtn') && isgraphics(tab.loadFeaturePackBtn)
                tab.loadFeaturePackBtn.Enable = 'on';
            end
            handles.channelTabs(channelIdx) = tab;
            handles = syncTabsToState(handles);
            guidata(fig, handles);
            refreshTrainingFeatureSummary(fig, channelIdx);
            appendTrainLog(fig, sprintf('[Channel %d] Parameter file set: %s', channelIdx, selectedPath));
            return;
        case 'train'
            startPath = strtrim(char(string(tab.trainDirEdit.Value)));
            dialogTitle = sprintf('Select training directory for Channel %d', channelIdx);
        case 'validation'
            startPath = strtrim(char(string(tab.valDirEdit.Value)));
            dialogTitle = sprintf('Select validation directory for Channel %d', channelIdx);
        case 'output'
            startPath = strtrim(char(string(tab.outputDirEdit.Value)));
            dialogTitle = sprintf('Select output directory for Channel %d', channelIdx);
        otherwise
            return;
    end
    if isempty(startPath) || ~isfolder(startPath)
        startPath = pwd;
    end

    selectedPath = uigetdir(startPath, dialogTitle);
    if isequal(selectedPath, 0)
        return;
    end

    switch dirType
        case 'train'
            tab.trainDirEdit.Value = selectedPath;
        case 'validation'
            tab.valDirEdit.Value = selectedPath;
        case 'output'
            tab.outputDirEdit.Value = selectedPath;
    end
    handles.channelTabs(channelIdx) = tab;
    handles = syncTabsToState(handles);
    guidata(fig, handles);
end

function cfg = defaultChannelFeatureConfig()
    cfg = struct( ...
        'selectedFeatures', {{}}, ...
        'customExpressions', struct('name', {}, 'expression', {}));
end

function cfg = normalizeChannelFeatureConfig(rawCfg)
    cfg = defaultChannelFeatureConfig();
    if ~isstruct(rawCfg)
        return;
    end

    if isfield(rawCfg, 'selectedFeatures')
        sel = rawCfg.selectedFeatures;
        if ischar(sel)
            cfg.selectedFeatures = {sel};
        elseif isstring(sel)
            cfg.selectedFeatures = cellstr(sel(:));
        elseif iscell(sel)
            cfg.selectedFeatures = sel(:)';
        end
    end

    if isfield(rawCfg, 'customExpressions') && isstruct(rawCfg.customExpressions)
        if isempty(rawCfg.customExpressions)
            cfg.customExpressions = struct('name', {}, 'expression', {});
        else
            ce = rawCfg.customExpressions(:);
            if ~(isfield(ce, 'name') && isfield(ce, 'expression'))
                cfg.customExpressions = struct('name', {}, 'expression', {});
            else
                for k = 1:numel(ce)
                    ce(k).name = char(string(ce(k).name));
                    ce(k).expression = char(string(ce(k).expression));
                end
                cfg.customExpressions = normalizeCustomExpressionList(ce);
            end
        end
    end
end

function cfg = getChannelFeatureConfig(handles, channelIdx)
    cfg = defaultChannelFeatureConfig();
    if ~isfield(handles, 'channelFeatureConfigs') || isempty(handles.channelFeatureConfigs)
        return;
    end
    idx = max(1, round(channelIdx));
    if idx <= numel(handles.channelFeatureConfigs)
        cfg = normalizeChannelFeatureConfig(handles.channelFeatureConfigs(idx));
    end
end

function handles = setChannelFeatureConfig(handles, channelIdx, cfg)
    idx = max(1, round(channelIdx));
    cfg = normalizeChannelFeatureConfig(cfg);
    if ~isfield(handles, 'channelFeatureConfigs') || isempty(handles.channelFeatureConfigs)
        handles.channelFeatureConfigs = repmat(defaultChannelFeatureConfig(), 1, idx);
    elseif numel(handles.channelFeatureConfigs) < idx
        nExisting = numel(handles.channelFeatureConfigs);
        handles.channelFeatureConfigs(nExisting+1:idx) = repmat(defaultChannelFeatureConfig(), 1, idx - nExisting);
    end
    handles.channelFeatureConfigs(idx) = cfg;
end

function onSelectTrainingFeatures(fig)
    handles = guidata(fig);
    contextChannel = max(1, handles.selectedChannelRow);
    onSelectTrainingFeaturesForChannel(fig, contextChannel);
end

function onSelectTrainingFeaturesForChannel(fig, channelIdx)
    [handles, channelParams, paramChannelIdx, paramPathUsed, errMsg] = ...
        resolveContextParamsForChannel(fig, channelIdx);
    if ~isempty(errMsg)
        uialert(fig, errMsg, 'Missing Parameter File');
        return;
    end

    [~, fittingMethod, has3D] = getTrainingFeatureSelectionContext(handles, paramChannelIdx, channelParams);
    chanCfg = getChannelFeatureConfig(handles, channelIdx);
    [selected, customExpr, cancelled] = snap_helpers.classification.featureSelectionUI( ...
        fittingMethod, has3D, false, chanCfg.selectedFeatures, chanCfg.customExpressions);

    if cancelled
        return;
    end

    chanCfg.selectedFeatures = selected;
    chanCfg.customExpressions = customExpr;
    handles = setChannelFeatureConfig(handles, channelIdx, chanCfg);
    handles.selectedChannelRow = max(1, round(channelIdx));
    guidata(fig, handles);
    refreshTrainingFeatureSummary(fig, channelIdx);
    appendTrainLog(fig, sprintf('[Channel %d] Updated feature selection using parameter file: %s', ...
        channelIdx, paramPathUsed));
end

function onLoadExpressionPackForChannel(fig, channelIdx)
    if nargin < 2 || ~isfinite(channelIdx)
        channelIdx = [];
    end
    onLoadExpressionPack(fig, channelIdx);
end

function onLoadExpressionPack(fig, targetChannels)
    if nargin < 2
        targetChannels = [];
    end
    handles = guidata(fig);

    [file, path] = uigetfile({'*.mat', 'Expression pack MAT files (*.mat)'}, ...
        'Select SVM expression pack file');
    if isequal(file, 0)
        return;
    end
    packPath = fullfile(path, file);

    [pack, errMsg] = loadExpressionPackForTraining(packPath);
    if ~isempty(errMsg)
        uialert(fig, errMsg, 'Expression Pack Error');
        return;
    end

    tableData = handles.channelTable.Data;
    if isempty(tableData)
        uialert(fig, 'No channels are currently configured in SNAP_train.', 'No Channels');
        return;
    end

    nRows = size(tableData, 1);
    targetChannels = unique(round(targetChannels(:)'));
    targetChannels = targetChannels(isfinite(targetChannels) & targetChannels >= 1 & targetChannels <= nRows);
    rowChannels = nan(1, nRows);
    for r = 1:nRows
        try
            rowChannels(r) = parseChannelFromRowValue(tableData{r, 2});
        catch
            rowChannels(r) = r;
        end
    end

    caps = struct([]);
    paramCtx = repmat(struct('ok', false, 'channel', 1, 'paramPath', '', ...
        'params', struct(), 'paramChannelIdx', 1), 1, nRows);
    capabilityErrors = {};
    for r = 1:nRows
        ch = rowChannels(r);
        if ~isempty(targetChannels) && ~ismember(ch, targetChannels)
            continue;
        end
        [handles, params, paramChannelIdx, paramPathUsed, errMsg] = resolveContextParamsForChannel(fig, ch);
        if ~isempty(errMsg)
            capabilityErrors{end+1} = sprintf('[Channel %d] %s', ch, errMsg); %#ok<AGROW>
            continue;
        end
        try
            capLocal = snap_helpers.classification.resolveChannelCapabilities(params, ...
                'Channels', paramChannelIdx, ...
                'IncludeFeatureInfo', false);
            if isempty(capLocal)
                error('No capability context returned.');
            end
            cap = capLocal(1);
            cap.channelIdx = ch;
            if isempty(caps)
                caps = cap;
            else
                caps(end+1) = cap; %#ok<AGROW>
            end
            paramCtx(r) = struct( ...
                'ok', true, ...
                'channel', ch, ...
                'paramPath', paramPathUsed, ...
                'params', params, ...
                'paramChannelIdx', paramChannelIdx);
        catch ME
            capabilityErrors{end+1} = sprintf('[Channel %d] Failed to resolve capabilities: %s', ...
                ch, ME.message); %#ok<AGROW>
        end
    end
    guidata(fig, handles);

    if isempty(caps)
        uialert(fig, ['Unable to resolve per-channel capabilities for expression-pack validation. ', ...
            'Set a valid parameter file in each channel tab.'], 'Expression Pack Error');
        return;
    end

    validCapChannels = [caps.channelIdx];
    packForValidation = pack;
    if numel(pack.channelPacks) > 1
        keepMask = false(1, numel(pack.channelPacks));
        for i = 1:numel(pack.channelPacks)
            chIdx = inferPackChannelIndex(pack.channelPacks(i), i);
            if ~isfinite(chIdx)
                continue;
            end
            chIdx = round(chIdx);
            if ~isempty(targetChannels) && ~ismember(chIdx, targetChannels)
                continue;
            end
            keepMask(i) = ismember(chIdx, validCapChannels);
        end
        if any(keepMask)
            packForValidation.channelPacks = pack.channelPacks(keepMask);
        end
    end

    [packForValidation, validationReport] = ...
        snap_helpers.classification.validateExpressionPackAgainstCapabilities( ...
            packForValidation, caps, 'Mode', 'permissive', 'AutoGuardUnsafeExpressions', true);

    if ~validationReport.success && ~isempty(validationReport.errors)
        uialert(fig, sprintf('Expression pack is incompatible with configured channel parameter files:\n%s', ...
            strjoin(validationReport.errors, newline)), 'Expression Pack Error');
        return;
    end

    if numel(pack.channelPacks) > 1
        pack = packForValidation;
    end

    appliedChannels = [];
    skippedChannels = [];
    totalBase = 0;
    totalCustom = 0;
    totalDroppedBase = 0;
    totalDroppedCustom = 0;
    compatLogLines = {};

    for i = 1:numel(pack.channelPacks)
        chPack = pack.channelPacks(i);
        channelIdx = inferPackChannelIndex(chPack, i);
        if ~isfinite(channelIdx)
            skippedChannels(end+1) = i; %#ok<AGROW>
            continue;
        end
        channelIdx = round(channelIdx);

        if ~isempty(targetChannels) && ~ismember(channelIdx, targetChannels)
            continue;
        end

        rowMatch = find(rowChannels == channelIdx, 1, 'first');
        if isempty(rowMatch)
            skippedChannels(end+1) = channelIdx; %#ok<AGROW>
            continue;
        end
        if ~paramCtx(rowMatch).ok
            compatLogLines{end+1} = sprintf('[Channel %d] Skipped expression-pack apply: unresolved parameter context.', ...
                channelIdx); %#ok<AGROW>
            continue;
        end

        cfg = normalizeChannelFeatureConfig(chPack);
        [cfg, compatReport] = sanitizeChannelFeatureConfig(cfg, ...
            paramCtx(rowMatch).params, paramCtx(rowMatch).paramChannelIdx);
        handles = setChannelFeatureConfig(handles, channelIdx, cfg);
        appliedChannels(end+1) = channelIdx; %#ok<AGROW>
        totalBase = totalBase + numel(cfg.selectedFeatures);
        totalCustom = totalCustom + numel(cfg.customExpressions);

        if ~isempty(compatReport.droppedBase)
            totalDroppedBase = totalDroppedBase + numel(compatReport.droppedBase);
            compatLogLines{end+1} = sprintf('[Channel %d] Expression pack dropped incompatible base feature(s): %s', ... %#ok<AGROW>
                channelIdx, strjoin(compatReport.droppedBase, ', '));
        end
        if ~isempty(compatReport.droppedCustom)
            totalDroppedCustom = totalDroppedCustom + numel(compatReport.droppedCustom);
            compatLogLines{end+1} = sprintf('[Channel %d] Expression pack dropped incompatible custom expression(s): %s', ... %#ok<AGROW>
                channelIdx, strjoin(compatReport.droppedCustom, ', '));
        end
        if isfield(compatReport, 'autoGuarded') && ~isempty(compatReport.autoGuarded)
            compatLogLines{end+1} = sprintf('[Channel %d] Expression pack auto-guarded expression(s): %s', ... %#ok<AGROW>
                channelIdx, strjoin(compatReport.autoGuarded, ', '));
        end
    end

    if ~isempty(targetChannels)
        missingTargets = setdiff(targetChannels, appliedChannels, 'stable');
        if numel(pack.channelPacks) == 1 && ~isempty(missingTargets)
            cfgSingle = normalizeChannelFeatureConfig(pack.channelPacks(1));
            for mt = missingTargets
                rowMatch = find(rowChannels == mt, 1, 'first');
                if isempty(rowMatch) || ~paramCtx(rowMatch).ok
                    compatLogLines{end+1} = sprintf('[Channel %d] Skipped expression-pack apply: unresolved parameter context.', ...
                        mt); %#ok<AGROW>
                    continue;
                end
                cfg = cfgSingle;
                [cfg, compatReport] = sanitizeChannelFeatureConfig(cfg, ...
                    paramCtx(rowMatch).params, paramCtx(rowMatch).paramChannelIdx);
                handles = setChannelFeatureConfig(handles, mt, cfg);
                appliedChannels(end+1) = mt; %#ok<AGROW>
                totalBase = totalBase + numel(cfg.selectedFeatures);
                totalCustom = totalCustom + numel(cfg.customExpressions);
                if ~isempty(compatReport.droppedBase)
                    totalDroppedBase = totalDroppedBase + numel(compatReport.droppedBase);
                    compatLogLines{end+1} = sprintf('[Channel %d] Expression pack dropped incompatible base feature(s): %s', ... %#ok<AGROW>
                        mt, strjoin(compatReport.droppedBase, ', '));
                end
                if ~isempty(compatReport.droppedCustom)
                    totalDroppedCustom = totalDroppedCustom + numel(compatReport.droppedCustom);
                    compatLogLines{end+1} = sprintf('[Channel %d] Expression pack dropped incompatible custom expression(s): %s', ... %#ok<AGROW>
                        mt, strjoin(compatReport.droppedCustom, ', '));
                end
            end
        end
    end

    for i = 1:numel(validationReport.channelReports)
        chReport = validationReport.channelReports(i);
        if ~isfield(chReport, 'channelIdx') || ~isfinite(chReport.channelIdx) || ...
                ~any(rowChannels == round(chReport.channelIdx))
            continue;
        end
        for w = 1:numel(chReport.warnings)
            compatLogLines{end+1} = chReport.warnings{w}; %#ok<AGROW>
        end
    end

    if isempty(appliedChannels)
        msg = sprintf(['Expression pack loaded (%s), but no channels matched this UI.\n' ...
            'Detected channels in UI: 1..%d'], file, nRows);
        uialert(fig, msg, 'No Matching Channels');
        appendTrainLog(fig, sprintf('Loaded expression pack: %s (no matching channels).', file));
        if ~isempty(skippedChannels)
            appendTrainLog(fig, sprintf('Skipped channel indices: %s', formatChannelIndexList(skippedChannels)));
        end
        return;
    end

    guidata(fig, handles);
    if isempty(targetChannels)
        refreshTrainingFeatureSummary(fig);
    else
        for ch = unique(appliedChannels(:)')
            refreshTrainingFeatureSummary(fig, ch);
        end
    end
    appendTrainLog(fig, sprintf('Loaded expression pack: %s', file));
    appendTrainLog(fig, sprintf(['Applied expression pack to channel(s): %s ' ...
        '(total base features=%d, total custom expressions=%d).'], ...
        formatChannelIndexList(appliedChannels), totalBase, totalCustom));
    appendTrainLog(fig, sprintf('Expression-pack validation mode: %s', validationReport.mode));
    for i = 1:numel(capabilityErrors)
        appendTrainLog(fig, capabilityErrors{i});
    end
    for i = 1:numel(compatLogLines)
        appendTrainLog(fig, compatLogLines{i});
    end
    if totalDroppedBase > 0 || totalDroppedCustom > 0
        appendTrainLog(fig, sprintf('Expression-pack compatibility pruning: dropped base=%d, dropped custom=%d.', ...
            totalDroppedBase, totalDroppedCustom));
    end
    if validationReport.nAutoGuarded > 0
        appendTrainLog(fig, sprintf('Expression-pack safety guard applied to %d expression(s).', ...
            validationReport.nAutoGuarded));
    end
    if ~isempty(skippedChannels)
        appendTrainLog(fig, sprintf('Skipped channel indices (not present in this UI): %s', ...
            formatChannelIndexList(skippedChannels)));
    end
end

function channelIdx = inferPackChannelIndex(chPack, fallbackIdx)
    channelIdx = nan;
    if nargin < 2 || isempty(fallbackIdx)
        fallbackIdx = 1;
    end
    if ~isstruct(chPack)
        return;
    end

    if isfield(chPack, 'channelIdx')
        channelIdx = normalizeScalarNumeric(chPack.channelIdx, nan);
    elseif isfield(chPack, 'channel')
        channelIdx = normalizeScalarNumeric(chPack.channel, nan);
    end

    if ~isfinite(channelIdx)
        channelIdx = fallbackIdx;
    end
end

function txt = formatChannelIndexList(indices)
    vals = unique(round(indices(:)'));
    vals = vals(isfinite(vals));
    if isempty(vals)
        txt = 'none';
        return;
    end
    pieces = arrayfun(@(v) sprintf('%d', v), vals, 'UniformOutput', false);
    txt = strjoin(pieces, ', ');
end

function [pack, errMsg] = loadExpressionPackForTraining(packPath)
    pack = struct();
    errMsg = '';

    if ~(ischar(packPath) || isstring(packPath))
        errMsg = 'Expression pack path must be char or string.';
        return;
    end
    packPath = char(string(packPath));

    if exist(packPath, 'file') ~= 2
        errMsg = sprintf('Expression pack file not found:\n%s', packPath);
        return;
    end

    try
        raw = load(packPath);
    catch ME
        errMsg = sprintf('Failed to load expression pack file:\n%s', ME.message);
        return;
    end

    candidates = {};
    if isfield(raw, 'expressionPack')
        candidates{end+1} = raw.expressionPack; %#ok<AGROW>
    end
    if isfield(raw, 'pack')
        candidates{end+1} = raw.pack; %#ok<AGROW>
    end

    rawFields = fieldnames(raw);
    for i = 1:numel(rawFields)
        value = raw.(rawFields{i});
        if isstruct(value) && (isfield(value, 'channelPacks') || ...
                isfield(value, 'selectedFeatures') || isfield(value, 'customExpressions'))
            candidates{end+1} = value; %#ok<AGROW>
        end
    end

    candidateErrors = {};
    for i = 1:numel(candidates)
        candidate = candidates{i};
        try
            [pack, normalizeReport] = snap_helpers.classification.normalizeExpressionPack(candidate);
            if ~isempty(normalizeReport.warnings)
                warning('SNAP_train:ExpressionPackNormalization', ...
                    'Expression pack normalized with warnings: %s', ...
                    strjoin(normalizeReport.warnings, ' | '));
            end
            return;
        catch ME
            candidateErrors{end+1} = ME.message; %#ok<AGROW>
        end
    end

    errMsg = sprintf(['File does not contain a valid expression pack.\n' ...
        'Expected a pack with channel definitions (selectedFeatures/customExpressions).\n\nFile: %s'], packPath);
    if ~isempty(candidateErrors)
        errMsg = sprintf('%s\n\nNormalization errors:\n- %s', ...
            errMsg, strjoin(candidateErrors, '\n- '));
    end
end

function refreshTrainingFeatureSummary(fig, contextChannel)
    handles = guidata(fig);
    if ~isfield(handles, 'channelTabs') || isempty(handles.channelTabs)
        return;
    end

    if nargin < 2 || isempty(contextChannel)
        contextChannel = max(1, min(numel(handles.channelTabs), handles.selectedChannelRow));
    end
    contextChannel = max(1, min(numel(handles.channelTabs), round(contextChannel)));

    tab = handles.channelTabs(contextChannel);
    if ~isfield(tab, 'featureCountLabel') || ~isgraphics(tab.featureCountLabel) || ...
            ~isfield(tab, 'featureListArea') || ~isgraphics(tab.featureListArea)
        return;
    end

    chanCfg = getChannelFeatureConfig(handles, contextChannel);
    nBase = numel(chanCfg.selectedFeatures);
    nCustom = numel(chanCfg.customExpressions);

    if nBase == 0 && nCustom == 0
        tab.featureCountLabel.Text = sprintf('Channel %d: AUTO (all non-position features)', contextChannel);
        tab.featureCountLabel.FontColor = [0.45 0.45 0.45];
        tab.featureListArea.Value = { ...
            sprintf('Context channel: %d', contextChannel), ...
            'AUTO: all non-position features', ...
            'Custom expressions: none' ...
        };
        handles.channelTabs(contextChannel) = tab;
        guidata(fig, handles);
        return;
    end

    nTotal = nBase + nCustom;
    tab.featureCountLabel.Text = sprintf('Channel %d: %d base + %d custom = %d total', ...
        contextChannel, nBase, nCustom, nTotal);
    tab.featureCountLabel.FontColor = [0.2 0.6 0.2];

    lines = {sprintf('Context channel: %d', contextChannel)};
    if nBase > 0
        for i = 1:nBase
            lines{end+1} = chanCfg.selectedFeatures{i}; %#ok<AGROW>
        end
    else
        lines{end+1} = 'AUTO: all non-position features'; %#ok<AGROW>
    end
    for i = 1:nCustom
        lines{end+1} = sprintf('[EXPR] %s = %s', ...
            char(string(chanCfg.customExpressions(i).name)), ...
            char(string(chanCfg.customExpressions(i).expression))); %#ok<AGROW>
    end
    tab.featureListArea.Value = lines;
    handles.channelTabs(contextChannel) = tab;
    guidata(fig, handles);
end

function [contextChannel, fittingMethod, has3D] = getTrainingFeatureSelectionContext(handles, contextChannel, paramsOverride)
    if nargin < 2 || isempty(contextChannel) || ~isfinite(contextChannel) || contextChannel < 1
        contextChannel = 1;
    end
    fittingMethod = '3D Gaussian';
    has3D = true;

    if nargin >= 2 && ~isempty(contextChannel) && isfinite(contextChannel) && contextChannel >= 1
        contextChannel = round(contextChannel);
    elseif isfield(handles, 'selectedChannelRow') && isfinite(handles.selectedChannelRow)
        contextChannel = max(1, round(handles.selectedChannelRow));
    end

    if isfield(handles, 'channelTable') && isgraphics(handles.channelTable) && ~isempty(handles.channelTable.Data)
        contextChannel = max(1, min(size(handles.channelTable.Data, 1), contextChannel));
    end

    paramsForContext = struct();
    if nargin >= 3 && isstruct(paramsOverride) && ~isempty(fieldnames(paramsOverride))
        paramsForContext = paramsOverride;
    elseif isfield(handles, 'loadedParams') && isstruct(handles.loadedParams) && ~isempty(fieldnames(handles.loadedParams))
        paramsForContext = handles.loadedParams;
    end

    if ~isempty(fieldnames(paramsForContext))
        [fm, h3d] = inferChannelTrainingContext(paramsForContext, contextChannel);
        if ~isempty(fm)
            fittingMethod = char(string(fm));
        end
        if ~isempty(h3d)
            has3D = logical(h3d);
        end
    end
end

function [handles, params, effectiveChannelIdx, paramPathUsed, errMsg] = resolveContextParamsForChannel(fig, channelIdx)
    handles = guidata(fig);
    params = struct();
    effectiveChannelIdx = 1;
    paramPathUsed = '';
    errMsg = '';

    if nargin < 2 || ~isfinite(channelIdx) || channelIdx < 1
        channelIdx = 1;
    end
    channelIdx = round(channelIdx);

    tableData = {};
    if isfield(handles, 'channelTable') && isgraphics(handles.channelTable)
        tableData = handles.channelTable.Data;
    end

    rowIdx = [];
    if ~isempty(tableData)
        for r = 1:size(tableData, 1)
            try
                if parseChannelFromRowValue(tableData{r, 2}) == channelIdx
                    rowIdx = r;
                    break;
                end
            catch
                % Ignore malformed row labels and continue.
            end
        end
        if isempty(rowIdx) && channelIdx <= size(tableData, 1)
            rowIdx = channelIdx;
        end
    end

    templateParamPath = '';
    if isfield(handles, 'paramPathEdit') && isgraphics(handles.paramPathEdit)
        templateParamPath = strtrim(char(string(handles.paramPathEdit.Value)));
    end
    paramPathUsed = resolveChannelParameterPathFromRow(tableData, rowIdx, templateParamPath);

    if isempty(paramPathUsed)
        if isfield(handles, 'loadedParams') && isstruct(handles.loadedParams) && ...
                ~isempty(fieldnames(handles.loadedParams))
            params = handles.loadedParams;
            effectiveChannelIdx = resolveEffectiveParamChannelIndex(params, channelIdx);
            paramPathUsed = '(loaded template parameter struct)';
            return;
        end
        errMsg = sprintf('Channel %d has no parameter file configured.', channelIdx);
        return;
    end

    [params, effectiveChannelIdx, loadErr] = resolveTrainingParams(paramPathUsed, channelIdx);
    if ~isempty(loadErr)
        errMsg = loadErr;
    end
end

function paramPath = resolveChannelParameterPathFromRow(tableData, rowIdx, templateParamPath)
    paramPath = '';
    if nargin < 3 || ~(ischar(templateParamPath) || isstring(templateParamPath))
        templateParamPath = '';
    end
    templateParamPath = strtrim(char(string(templateParamPath)));

    if nargin >= 1 && ~isempty(tableData) && ~isempty(rowIdx) && ...
            isfinite(rowIdx) && rowIdx >= 1 && rowIdx <= size(tableData, 1)
        rowIdx = round(rowIdx);
        if size(tableData, 2) >= 8
            paramPath = strtrim(char(string(tableData{rowIdx, 5})));
        end
    end

    if isempty(paramPath)
        paramPath = templateParamPath;
    end
end

function [params, effectiveChannelIdx, errMsg] = resolveTrainingParams(paramPath, requestedChannelIdx)
    params = struct();
    effectiveChannelIdx = 1;
    errMsg = '';

    paramPath = strtrim(char(string(paramPath)));
    if isempty(paramPath)
        errMsg = 'Parameter file path is empty.';
        return;
    end
    if exist(paramPath, 'file') ~= 2
        errMsg = sprintf('Parameter file not found: %s', paramPath);
        return;
    end

    [params, ~, loadErr] = loadTrainingParameters(paramPath);
    if ~isempty(loadErr)
        errMsg = loadErr;
        return;
    end
    effectiveChannelIdx = resolveEffectiveParamChannelIndex(params, requestedChannelIdx);
end

function idx = resolveEffectiveParamChannelIndex(params, requestedChannelIdx)
    idx = max(1, round(normalizeScalarNumeric(requestedChannelIdx, 1)));
    if nargin < 1 || ~isstruct(params) || isempty(fieldnames(params))
        return;
    end

    nChannels = 1;
    try
        caps = snap_helpers.classification.resolveChannelCapabilities(params, 'IncludeFeatureInfo', false);
        if ~isempty(caps)
            nChannels = numel(caps);
        else
            nChannels = inferNumChannelsFromParameters(params, 1);
        end
    catch
        nChannels = inferNumChannelsFromParameters(params, 1);
    end

    if idx > max(1, nChannels)
        idx = 1;
    end
end

function [cfgOut, report] = sanitizeChannelFeatureConfig(cfgIn, params, channelIdx)
    cfgOut = normalizeChannelFeatureConfig(cfgIn);
    originalSelected = cfgOut.selectedFeatures;
    report = struct( ...
        'droppedBase', {{}}, ...
        'droppedCustom', {{}}, ...
        'autoGuarded', {{}}, ...
        'fittingMethod', '3D Gaussian', ...
        'has3D', true, ...
        'warnings', {{}}, ...
        'errors', {{}});

    if nargin < 2 || ~isstruct(params) || isempty(fieldnames(params))
        return;
    end
    if nargin < 3 || ~isfinite(channelIdx) || channelIdx < 1
        channelIdx = 1;
    end

    effectiveChannelIdx = resolveEffectiveParamChannelIndex(params, channelIdx);

    try
        caps = snap_helpers.classification.resolveChannelCapabilities(params, ...
            'Channels', effectiveChannelIdx, ...
            'IncludeFeatureInfo', false);
    catch ME
        report.errors{end+1} = sprintf('Capability resolution failed: %s', ME.message); %#ok<AGROW>
        return;
    end
    if isempty(caps)
        report.errors{end+1} = 'Capability resolution returned no channel context.'; %#ok<AGROW>
        return;
    end
    cap = caps(1);
    report.fittingMethod = cap.fittingMethod;
    report.has3D = logical(cap.has3D);

    customIn = normalizeCustomExpressionList(cfgOut.customExpressions);
    inputPack = struct( ...
        'specVersion', '2.0.0', ...
        'packId', sprintf('snap_train_sanitize_ch%d', channelIdx), ...
        'strictModeDefault', false, ...
        'channelPacks', struct( ...
            'channelIdx', channelIdx, ...
            'selectedFeatures', {cfgOut.selectedFeatures}, ...
            'customExpressions', customIn, ...
            'requiredFeatures', {{}}, ...
            'requiredCapabilities', struct( ...
                'fittingMethod', cap.fittingMethod, ...
                'has3D', logical(cap.has3D), ...
                'hasPhysicalSpacing', false)));

    [sanitizedPack, validationReport] = ...
        snap_helpers.classification.validateExpressionPackAgainstCapabilities( ...
            inputPack, cap, 'Mode', 'permissive', 'AutoGuardUnsafeExpressions', true);

    if ~isempty(validationReport.channelReports)
        chReport = validationReport.channelReports(1);
        report.warnings = chReport.warnings;
        report.errors = chReport.errors;
        report.autoGuarded = chReport.autoGuarded;
    end

    if ~isempty(sanitizedPack.channelPacks)
        clean = sanitizedPack.channelPacks(1);
        cfgOut.selectedFeatures = clean.selectedFeatures;
        cfgOut.customExpressions = normalizeCustomExpressionList(clean.customExpressions);
    end

    report.droppedBase = setdiff(originalSelected, cfgOut.selectedFeatures, 'stable');

    inNames = {};
    if isstruct(customIn) && ~isempty(customIn)
        inNames = arrayfun(@(s) char(string(s.name)), customIn, 'UniformOutput', false);
    end
    outNames = {};
    if isstruct(cfgOut.customExpressions) && ~isempty(cfgOut.customExpressions)
        outNames = arrayfun(@(s) char(string(s.name)), cfgOut.customExpressions, 'UniformOutput', false);
    end
    report.droppedCustom = setdiff(inNames, outNames, 'stable');
end

function onTrainSelectedChannels(fig)
    handles = guidata(fig);
    handles = syncTabsToState(handles);
    guidata(fig, handles);

    templateParamPath = strtrim(char(string(handles.paramPathEdit.Value)));

    tableData = handles.channelTable.Data;
    if isempty(tableData)
        uialert(fig, 'No channels configured for training.', 'No Channels');
        return;
    end

    selectedRows = find(cellfun(@isSelectedChannelFlag, tableData(:, 1)));
    if isempty(selectedRows)
        uialert(fig, 'Select at least one channel to train.', 'No Channels Selected');
        return;
    end

    defaultOutputDir = fileparts(templateParamPath);
    if isempty(defaultOutputDir) || ~isfolder(defaultOutputDir)
        defaultOutputDir = pwd;
    end

    channelConfigs = struct('channel', {}, 'matchDistance', {}, 'convertFijiCoords', {}, ...
        'parameterFile', {}, 'parameterStruct', {}, 'parameterChannelIndex', {}, ...
        'trainDirectory', {}, 'validationDirectory', {}, 'outputDirectory', {}, ...
        'outputPath', {}, 'useOptimization', {}, 'includeSweepReport', {}, ...
        'trainOptions', {}, 'sweepConfig', {});
    missingParameterChannels = [];
    invalidParameterMessages = {};
    missingTrainChannels = [];
    missingValidationChannels = [];
    invalidMatchDistanceChannels = [];
    invalidOptionMessages = {};

    for idx = 1:numel(selectedRows)
        row = selectedRows(idx);
        ch = parseChannelFromRowValue(tableData{row, 2});
        matchDistance = 2;
        convertFiji = false;

        if size(tableData, 2) >= 8
            matchDistance = normalizeScalarNumeric(tableData{row, 3}, nan);
            if ~isfinite(matchDistance) || matchDistance < 0
                invalidMatchDistanceChannels(end+1) = ch; %#ok<AGROW>
            end
            convertFiji = isSelectedChannelFlag(tableData{row, 4});
        elseif size(tableData, 2) >= 7
            matchDistance = normalizeScalarNumeric(tableData{row, 3}, nan);
            if ~isfinite(matchDistance) || matchDistance < 0
                invalidMatchDistanceChannels(end+1) = ch; %#ok<AGROW>
            end
            convertFiji = isSelectedChannelFlag(tableData{row, 4});
        elseif size(tableData, 2) >= 3
            % Backward compatibility with older layout.
            convertFiji = isSelectedChannelFlag(tableData{row, 3});
        end

        paramPath = resolveChannelParameterPathFromRow(tableData, row, templateParamPath);
        if isempty(paramPath) || exist(paramPath, 'file') ~= 2
            missingParameterChannels(end+1) = ch; %#ok<AGROW>
            continue;
        end
        [channelParams, parameterChannelIndex, paramErr] = resolveTrainingParams(paramPath, ch);
        if ~isempty(paramErr)
            invalidParameterMessages{end+1} = sprintf('Channel %d (%s): %s', ch, paramPath, paramErr); %#ok<AGROW>
            continue;
        end

        trainDir = '';
        if size(tableData, 2) >= 8
            trainDir = strtrim(char(string(tableData{row, 6})));
        elseif size(tableData, 2) >= 7
            trainDir = strtrim(char(string(tableData{row, 5})));
        elseif size(tableData, 2) >= 4
            trainDir = strtrim(char(string(tableData{row, 4})));
        end

        valDir = '';
        if size(tableData, 2) >= 8
            valDir = strtrim(char(string(tableData{row, 7})));
        elseif size(tableData, 2) >= 7
            valDir = strtrim(char(string(tableData{row, 6})));
        elseif size(tableData, 2) >= 5
            valDir = strtrim(char(string(tableData{row, 5})));
        end

        outName = '';
        if size(tableData, 2) >= 8
            outName = strtrim(char(string(tableData{row, 8})));
        elseif size(tableData, 2) >= 7
            outName = strtrim(char(string(tableData{row, 7})));
        elseif size(tableData, 2) >= 6
            outName = strtrim(char(string(tableData{row, 6})));
        end
        if isempty(outName)
            outName = sprintf('classifier_ch%d.mat', ch);
        end
        if ~endsWith(lower(outName), '.mat')
            outName = [outName '.mat'];
        end

        advCfg = getChannelAdvancedConfig(handles, ch);
        outDir = strtrim(char(string(advCfg.outputDirectory)));
        if isempty(outDir)
            outDir = defaultOutputDir;
        end
        useOptimization = logical(advCfg.optimizeWithSweep);
        includeSweepReport = useOptimization && logical(advCfg.includeSweepReport);

        try
            trainOptions = buildManualTrainingOptions(handles, ch);
            if useOptimization
                sweepConfig = buildSweepConfig(handles, ch);
            else
                sweepConfig = struct();
            end
        catch ME
            invalidOptionMessages{end+1} = sprintf('Channel %d: %s', ch, ME.message); %#ok<AGROW>
            continue;
        end

        if isempty(trainDir) || ~isfolder(trainDir)
            missingTrainChannels(end+1) = ch; %#ok<AGROW>
        end
        if useOptimization && (isempty(valDir) || ~isfolder(valDir))
            missingValidationChannels(end+1) = ch; %#ok<AGROW>
        end

        channelConfigs(end+1) = struct( ... %#ok<AGROW>
            'channel', ch, ...
            'matchDistance', matchDistance, ...
            'convertFijiCoords', convertFiji, ...
            'parameterFile', paramPath, ...
            'parameterStruct', channelParams, ...
            'parameterChannelIndex', parameterChannelIndex, ...
            'trainDirectory', trainDir, ...
            'validationDirectory', valDir, ...
            'outputDirectory', outDir, ...
            'outputPath', fullfile(outDir, outName), ...
            'useOptimization', useOptimization, ...
            'includeSweepReport', includeSweepReport, ...
            'trainOptions', trainOptions, ...
            'sweepConfig', sweepConfig);
    end

    if ~isempty(missingParameterChannels)
        msg = sprintf('Set valid per-channel parameter files for: %s', ...
            formatChannelList(missingParameterChannels));
        uialert(fig, msg, 'Missing Parameter Files');
        return;
    end

    if ~isempty(invalidParameterMessages)
        uialert(fig, sprintf('Per-channel parameter file errors:\n%s', ...
            strjoin(invalidParameterMessages, newline)), 'Parameter File Error');
        return;
    end

    if ~isempty(invalidOptionMessages)
        uialert(fig, sprintf('Invalid per-channel training options:\n%s', ...
            strjoin(invalidOptionMessages, newline)), 'Invalid Training Options');
        return;
    end

    if ~isempty(invalidMatchDistanceChannels)
        msg = sprintf('Set non-negative match distance values for: %s', ...
            formatChannelList(invalidMatchDistanceChannels));
        uialert(fig, msg, 'Invalid Match Distance');
        return;
    end

    if ~isempty(missingTrainChannels)
        msg = sprintf('Set valid training directories for: %s', formatChannelList(missingTrainChannels));
        uialert(fig, msg, 'Missing Training Directories');
        return;
    end

    if ~isempty(missingValidationChannels)
        msg = sprintf('Validation sweep is enabled. Set valid validation directories for: %s', ...
            formatChannelList(missingValidationChannels));
        uialert(fig, msg, 'Missing Validation Directories');
        return;
    end

    if isempty(channelConfigs)
        uialert(fig, 'No valid channel configurations remain after validation checks.', 'No Valid Channels');
        return;
    end

    duplicateOutputChannels = localFindDuplicateOutputChannels(channelConfigs);
    if ~isempty(duplicateOutputChannels)
        msg = sprintf(['Each selected channel must write to a unique classifier file. ', ...
            'Duplicate output targets detected for: %s'], ...
            formatChannelList(duplicateOutputChannels));
        uialert(fig, msg, 'Duplicate Output Classifier');
        return;
    end

    appendTrainLog(fig, '------------------------------------------------------------');
    appendTrainLog(fig, sprintf('Starting training (%d selected channels)...', numel(channelConfigs)));
    appendTrainLog(fig, 'Channel numbers are slot labels only; training stays channel-local (no cross-channel pooling).');
    appendTrainLog(fig, 'Parameter mode: PER-CHANNEL (each channel tab can use its own parameter file).');
    appendTrainLog(fig, 'Feature mode: PER-CHANNEL (base features + custom expressions can differ by channel).');
    setTrainProgressVisual(fig, 0, sprintf('Starting training for %d channel(s)...', numel(channelConfigs)), [0 0 0.8]);

    nSuccess = 0;
    nFailed = 0;
    results = struct('channel', {}, 'trainDirectory', {}, 'validationDirectory', {}, ...
        'outputPath', {}, 'success', {}, 'result', {}, 'error', {});

    for idx = 1:numel(channelConfigs)
        cfg = channelConfigs(idx);
        ch = cfg.channel;
        outPath = cfg.outputPath;
        [~, outName, outExt] = fileparts(outPath);
        outName = [outName outExt];
        stageCount = 3 + double(cfg.useOptimization);
        stage = 0;

        if exist(cfg.outputDirectory, 'dir') ~= 7
            try
                mkdir(cfg.outputDirectory);
            catch ME
                nFailed = nFailed + 1;
                results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                    'trainDirectory', cfg.trainDirectory, ...
                    'validationDirectory', cfg.validationDirectory, ...
                    'outputPath', outPath, 'success', false, ...
                    'result', struct(), 'error', sprintf('Failed creating output directory: %s', ME.message));
                appendTrainLog(fig, sprintf('[Channel %d] FAILED: Unable to create output directory "%s" (%s)', ...
                    ch, cfg.outputDirectory, ME.message));
                updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, ...
                    sprintf('Channel %d failed: output directory error.', ch));
                continue;
            end
        end

        stage = stage + 1;
        updateTrainProgressStage(fig, idx, numel(channelConfigs), stage, stageCount, ...
            sprintf('Channel %d: Discovering training pairs (%s)...', ch, outName));
        appendTrainLog(fig, sprintf('[Channel %d] Match distance: %.4g voxels', ch, cfg.matchDistance));
        appendTrainLog(fig, sprintf('[Channel %d] Parameter file: %s (effective parameter-channel index=%d)', ...
            ch, cfg.parameterFile, cfg.parameterChannelIndex));
        appendTrainLog(fig, sprintf('[Channel %d] Training directory: %s', ch, cfg.trainDirectory));
        appendTrainLog(fig, sprintf('[Channel %d] Discovering training pairs for classifier %s...', ch, outName));
        trainDiscoveryCb = @(msg) appendTrainLog(fig, sprintf('[Channel %d] %s', ch, msg));
        try
            [trainExports, trainLabels] = discoverDatasetPairs(cfg.trainDirectory, ch, 'ProgressCallback', trainDiscoveryCb);
        catch ME
            nFailed = nFailed + 1;
            results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                'trainDirectory', cfg.trainDirectory, ...
                'validationDirectory', cfg.validationDirectory, ...
                'outputPath', outPath, 'success', false, ...
                'result', struct(), 'error', sprintf('Training pair discovery failed: %s', ME.message));
            appendTrainLog(fig, sprintf('[Channel %d] FAILED during training pair discovery: %s', ch, ME.message));
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, ...
                sprintf('Channel %d failed: training pair discovery error.', ch));
            continue;
        end
        if isempty(trainExports)
            nFailed = nFailed + 1;
            results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                'trainDirectory', cfg.trainDirectory, ...
                'validationDirectory', cfg.validationDirectory, ...
                'outputPath', outPath, 'success', false, ...
                'result', struct(), 'error', sprintf('No export/label pairs found for channel %d in %s', ch, cfg.trainDirectory));
            appendTrainLog(fig, sprintf('[Channel %d] FAILED: No export/label pairs found.', ch));
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, ...
                sprintf('Channel %d failed: no training pairs found.', ch));
            continue;
        end
        appendTrainLog(fig, sprintf('[Channel %d] Found %d training pair(s).', ch, numel(trainExports)));
        logDiscoveredPairs(fig, ch, 'train', trainExports, trainLabels, 4);

        chFeatureCfg = getChannelFeatureConfig(handles, ch);
        [chFeatureCfg, compatReport] = sanitizeChannelFeatureConfig(chFeatureCfg, cfg.parameterStruct, cfg.parameterChannelIndex);
        if ~isempty(compatReport.droppedBase)
            appendTrainLog(fig, sprintf('[Channel %d] Auto-dropped incompatible base feature(s): %s', ...
                ch, strjoin(compatReport.droppedBase, ', ')));
        end
        if ~isempty(compatReport.droppedCustom)
            appendTrainLog(fig, sprintf('[Channel %d] Auto-dropped incompatible custom expression(s): %s', ...
                ch, strjoin(compatReport.droppedCustom, ', ')));
        end
        handles = setChannelFeatureConfig(handles, ch, chFeatureCfg);
        guidata(fig, handles);
        if isempty(chFeatureCfg.selectedFeatures) && isempty(chFeatureCfg.customExpressions)
            appendTrainLog(fig, sprintf('[Channel %d] Features: AUTO (all non-position base features), custom expressions: none.', ch));
        else
            appendTrainLog(fig, sprintf('[Channel %d] Features: %d base + %d custom expression(s).', ...
                ch, numel(chFeatureCfg.selectedFeatures), numel(chFeatureCfg.customExpressions)));
        end

        trainingProgressCb = @(msg) appendTrainLog(fig, sprintf('[Channel %d] %s', ch, msg));
        args = { ...
            'MatchDistance', cfg.matchDistance, ...
            'ConvertFijiCoords', logical(cfg.convertFijiCoords), ...
            'ParameterStruct', cfg.parameterStruct, ...
            'ChannelIndex', cfg.parameterChannelIndex, ...
            'SelectedFeatures', chFeatureCfg.selectedFeatures, ...
            'CustomExpressions', chFeatureCfg.customExpressions, ...
            'TrainingOptions', cfg.trainOptions, ...
            'HyperparameterSweep', false, ...
            'Verbose', true, ...
            'ProgressCallback', trainingProgressCb ...
        };

        [fittingMethod, has3D] = inferChannelTrainingContext(cfg.parameterStruct, cfg.parameterChannelIndex);
        if ~isempty(fittingMethod)
            args = [args, {'FittingMethod', fittingMethod}]; %#ok<AGROW>
        end
        if ~isempty(has3D)
            args = [args, {'Has3D', has3D}]; %#ok<AGROW>
        end

        if cfg.useOptimization
            stage = stage + 1;
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stage, stageCount, ...
                sprintf('Channel %d: Discovering validation pairs (%s)...', ch, outName));
            appendTrainLog(fig, sprintf('[Channel %d] Validation directory: %s', ch, cfg.validationDirectory));
            appendTrainLog(fig, sprintf('[Channel %d] Discovering validation pairs for classifier %s...', ch, outName));
            valDiscoveryCb = @(msg) appendTrainLog(fig, sprintf('[Channel %d] %s', ch, msg));
            try
                [valExports, valLabels] = discoverDatasetPairs(cfg.validationDirectory, ch, 'ProgressCallback', valDiscoveryCb);
            catch ME
                nFailed = nFailed + 1;
                results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                    'trainDirectory', cfg.trainDirectory, ...
                    'validationDirectory', cfg.validationDirectory, ...
                    'outputPath', outPath, 'success', false, ...
                    'result', struct(), 'error', sprintf('Validation pair discovery failed: %s', ME.message));
                appendTrainLog(fig, sprintf('[Channel %d] FAILED during validation pair discovery: %s', ch, ME.message));
                updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, ...
                    sprintf('Channel %d failed: validation pair discovery error.', ch));
                continue;
            end
            if isempty(valExports)
                nFailed = nFailed + 1;
                results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                    'trainDirectory', cfg.trainDirectory, ...
                    'validationDirectory', cfg.validationDirectory, ...
                    'outputPath', outPath, 'success', false, ...
                    'result', struct(), 'error', sprintf('No validation pairs found for channel %d in %s', ch, cfg.validationDirectory));
                appendTrainLog(fig, sprintf('[Channel %d] FAILED: No validation pairs found.', ch));
                updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, ...
                    sprintf('Channel %d failed: no validation pairs found.', ch));
                continue;
            end
            appendTrainLog(fig, sprintf('[Channel %d] Found %d validation pair(s).', ch, numel(valExports)));
            logDiscoveredPairs(fig, ch, 'validation', valExports, valLabels, 4);

            appendTrainLog(fig, sprintf('[Channel %d] Validation tie handling: %s.', ...
                ch, describeSweepTieHandlingPolicy(cfg.sweepConfig.tieHandling)));
            args = [args, { ...
                'ValidationExportFiles', valExports, ...
                'ValidationLabelFiles', valLabels, ...
                'HyperparameterSweep', true, ...
                'SweepKernels', cfg.sweepConfig.kernels, ...
                'SweepBoxConstraints', cfg.sweepConfig.boxConstraints, ...
                'SweepKernelScales', cfg.sweepConfig.kernelScales, ...
                'SweepPolynomialOrders', cfg.sweepConfig.polynomialOrders, ...
                'SweepTieHandling', cfg.sweepConfig.tieHandling ...
            }]; %#ok<AGROW>
            if isfield(cfg.sweepConfig, 'tiePromptCallback') && isa(cfg.sweepConfig.tiePromptCallback, 'function_handle')
                args = [args, {'SweepTiePromptCallback', cfg.sweepConfig.tiePromptCallback}]; %#ok<AGROW>
            elseif strcmpi(char(string(cfg.sweepConfig.tieHandling)), 'prompt')
                args = [args, {'SweepTiePromptCallback', ...
                    @(tieIdx, sweepResults) promptSweepTieSelection(fig, ch, tieIdx, sweepResults)}]; %#ok<AGROW>
            end
        end

        try
            stage = stage + 1;
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stage, stageCount, ...
                sprintf('Channel %d: Training classifier %s...', ch, outName));
            appendTrainLog(fig, sprintf('[Channel %d] Training classifier: %s', ch, outName));
            appendTrainLog(fig, sprintf('[Channel %d] Step: Feature extraction + SVM training', ch));
            res = SNAP_train(trainExports, trainLabels, outPath, args{:});
            nSuccess = nSuccess + 1;
            results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                'trainDirectory', cfg.trainDirectory, ...
                'validationDirectory', cfg.validationDirectory, ...
                'outputPath', outPath, 'success', true, ...
                'result', res, 'error', '');
            appendTrainLog(fig, sprintf('[Channel %d] SUCCESS', ch));
            if cfg.includeSweepReport
                emitSweepPerformanceReport(fig, ch, res, outPath);
            end
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, sprintf('Channel %d complete.', ch));
        catch ME
            nFailed = nFailed + 1;
            results(end+1) = struct('channel', ch, ... %#ok<AGROW>
                'trainDirectory', cfg.trainDirectory, ...
                'validationDirectory', cfg.validationDirectory, ...
                'outputPath', outPath, 'success', false, ...
                'result', struct(), 'error', ME.message);
            appendTrainLog(fig, sprintf('[Channel %d] FAILED: %s', ch, ME.message));
            updateTrainProgressStage(fig, idx, numel(channelConfigs), stageCount, stageCount, sprintf('Channel %d failed.', ch));
        end
        refreshTrainUiThrottled(0.05);
    end

    summary = struct();
    summary.timestamp = datetime('now');
    summary.templateParameterFile = templateParamPath;
    summary.parameterFilesByChannel = arrayfun(@(c) c.parameterFile, channelConfigs, 'UniformOutput', false);
    summary.parameterChannelIndexByChannel = arrayfun(@(c) c.parameterChannelIndex, channelConfigs);
    summary.outputDirectoriesByChannel = arrayfun(@(c) c.outputDirectory, channelConfigs, 'UniformOutput', false);
    summary.matchDistanceByChannel = arrayfun(@(c) c.matchDistance, channelConfigs);
    summary.optimizedByChannel = arrayfun(@(c) logical(c.useOptimization), channelConfigs);
    summary.sweepReportIncludedByChannel = arrayfun(@(c) logical(c.includeSweepReport), channelConfigs);
    summary.sweepTieHandlingByChannel = arrayfun(@(c) localGetSweepTieHandlingSafe(c), channelConfigs, 'UniformOutput', false);
    summary.optimized = any(summary.optimizedByChannel);
    summary.sweepReportIncluded = any(summary.sweepReportIncludedByChannel);
    summary.channelFeatureConfigs = handles.channelFeatureConfigs;
    summaryChannelConfigs = channelConfigs;
    for iCfg = 1:numel(summaryChannelConfigs)
        summaryChannelConfigs(iCfg).parameterStruct = struct();
    end
    summary.channelConfigs = summaryChannelConfigs;
    summary.results = results;
    summary.successCount = nSuccess;
    summary.failCount = nFailed;

    summaryBaseDir = defaultOutputDir;
    if ~isempty(channelConfigs)
        firstOutDir = strtrim(char(string(channelConfigs(1).outputDirectory)));
        if ~isempty(firstOutDir)
            summaryBaseDir = firstOutDir;
        end
    end
    summaryPath = fullfile(summaryBaseDir, ['SNAP_train_summary_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
    save(summaryPath, '-struct', 'summary');
    appendTrainLog(fig, sprintf('Saved training summary: %s', summaryPath));
    appendTrainLog(fig, sprintf('Training complete: %d success, %d failed.', nSuccess, nFailed));

    if nFailed == 0
        setTrainProgressVisual(fig, 1, sprintf('Training complete: %d/%d channels succeeded.', nSuccess, numel(channelConfigs)), [0 0.6 0]);
        uialert(fig, sprintf('Training complete.\nSaved %d classifier(s).', nSuccess), 'SNAP_train Complete', 'Icon', 'success');
    else
        setTrainProgressVisual(fig, 1, sprintf('Training finished with issues: %d success, %d failed.', nSuccess, nFailed), [0.8 0 0]);
        uialert(fig, sprintf('Training complete with issues.\nSuccess: %d\nFailed: %d', nSuccess, nFailed), ...
            'SNAP_train Complete', 'Icon', 'warning');
    end
end

function updateTrainProgressStage(fig, currentChannel, totalChannels, currentStage, totalStages, message)
    totalChannels = max(1, totalChannels);
    totalStages = max(1, totalStages);
    currentChannel = min(max(1, currentChannel), totalChannels);
    currentStage = min(max(0, currentStage), totalStages);

    fraction = ((currentChannel - 1) + (double(currentStage) / double(totalStages))) / double(totalChannels);
    setTrainProgressVisual(fig, fraction, message, [0 0 0.8]);
end

function setTrainProgressVisual(fig, fraction, message, fontColor)
    try
        handles = guidata(fig);
        if ~isfield(handles, 'progressStatusLabel') || ~isgraphics(handles.progressStatusLabel)
            return;
        end
        if ~isfield(handles, 'progressBarGrid') || ~isgraphics(handles.progressBarGrid)
            return;
        end

        frac = min(max(double(fraction), 0), 1);
        pct = round(100 * frac);
        handles.progressStatusLabel.Text = sprintf('[%d%%] %s', pct, char(string(message)));
        if nargin >= 4 && isnumeric(fontColor) && numel(fontColor) == 3
            handles.progressStatusLabel.FontColor = fontColor;
        end

        if frac <= 0
            handles.progressBarGrid.ColumnWidth = {0, '1x'};
        elseif frac >= 1
            handles.progressBarGrid.ColumnWidth = {'1x', 0};
        else
            greenPct = max(1, pct);
            grayPct = max(1, 100 - pct);
            handles.progressBarGrid.ColumnWidth = {sprintf('%dx', greenPct), sprintf('%dx', grayPct)};
        end

        refreshTrainUiThrottled(0.15);
    catch
        % Never block training due to UI progress updates.
    end
end

function appendTrainLog(fig, msg)
    stamp = datestr(now, 'HH:MM:SS');
    line = sprintf('[%s] %s', stamp, msg);
    fprintf('%s\n', line);

    try
        handles = guidata(fig);
        if ~isfield(handles, 'logText') || ~isgraphics(handles.logText)
            return;
        end
        current = handles.logText.Value;
        if ischar(current)
            current = {current};
        end
        handles.logText.Value = [current; {line}];
        handles.logText.Value = handles.logText.Value(max(1, end-500):end);
        refreshTrainUiThrottled(0.15);
    catch
        % Never block training due to UI log updates.
    end
end

function refreshTrainUiThrottled(minIntervalSec)
    if nargin < 1 || ~isfinite(minIntervalSec) || minIntervalSec <= 0
        minIntervalSec = 0.15;
    end

    persistent lastDrawTic;
    if isempty(lastDrawTic)
        lastDrawTic = tic;
        drawnow limitrate nocallbacks;
        return;
    end

    if toc(lastDrawTic) >= minIntervalSec
        lastDrawTic = tic;
        drawnow limitrate nocallbacks;
    end
end

function txt = formatChannelList(channels)
    channels = unique(channels(:)');
    labels = arrayfun(@(c) sprintf('Channel %d', c), channels, 'UniformOutput', false);
    txt = strjoin(labels, ', ');
end

function tieHandling = localGetSweepTieHandlingSafe(cfg)
    tieHandling = '';
    if ~isstruct(cfg) || ~isfield(cfg, 'useOptimization') || ~logical(cfg.useOptimization)
        return;
    end
    if isfield(cfg, 'sweepConfig') && isstruct(cfg.sweepConfig) && isfield(cfg.sweepConfig, 'tieHandling')
        tieHandling = char(string(cfg.sweepConfig.tieHandling));
    end
end

function channels = localFindDuplicateOutputChannels(channelConfigs)
    channels = [];
    if isempty(channelConfigs)
        return;
    end

    pathKeys = cell(1, numel(channelConfigs));
    for i = 1:numel(channelConfigs)
        pathKeys{i} = localNormalizePathKey(channelConfigs(i).outputPath);
    end

    duplicateChannels = [];
    uniqueKeys = unique(pathKeys, 'stable');
    for i = 1:numel(uniqueKeys)
        memberIdx = find(strcmp(pathKeys, uniqueKeys{i}));
        if numel(memberIdx) > 1
            duplicateChannels = [duplicateChannels, [channelConfigs(memberIdx).channel]]; %#ok<AGROW>
        end
    end

    if isempty(duplicateChannels)
        return;
    end

    channels = unique(duplicateChannels, 'stable');
end

function key = localNormalizePathKey(pathIn)
    key = char(string(pathIn));
    key = strtrim(key);
    key = strrep(key, '/', '\');
    key = lower(key);
end

function logDiscoveredPairs(fig, channelIdx, pairType, exportFiles, labelFiles, maxLines)
    if nargin < 6 || isempty(maxLines)
        maxLines = 4;
    end
    n = min(numel(exportFiles), maxLines);
    for i = 1:n
        [~, eName, eExt] = fileparts(exportFiles{i});
        [~, lName, lExt] = fileparts(labelFiles{i});
        appendTrainLog(fig, sprintf('[Channel %d]   %s pair %d: %s%s <-> %s%s', ...
            channelIdx, pairType, i, eName, eExt, lName, lExt));
    end
    if numel(exportFiles) > n
        appendTrainLog(fig, sprintf('[Channel %d]   ... and %d more %s pair(s)', ...
            channelIdx, numel(exportFiles) - n, pairType));
    end
end

function emitSweepPerformanceReport(fig, channelIdx, resultStruct, outputPath)
    if ~isfield(resultStruct, 'validation') || ~isfield(resultStruct.validation, 'sweepResults')
        appendTrainLog(fig, sprintf('[Channel %d] Sweep report unavailable (no sweep results found).', channelIdx));
        return;
    end

    sweepResults = resultStruct.validation.sweepResults;
    if isempty(sweepResults)
        appendTrainLog(fig, sprintf('[Channel %d] Sweep report unavailable (empty sweep results).', channelIdx));
        return;
    end

    nCombos = numel(sweepResults);
    f1Vals = nan(1, nCombos);
    accVals = nan(1, nCombos);
    for i = 1:nCombos
        if isfield(sweepResults(i), 'f1Real')
            f1Vals(i) = double(sweepResults(i).f1Real);
        end
        if isfield(sweepResults(i), 'accuracy')
            accVals(i) = double(sweepResults(i).accuracy);
        end
    end

    validIdx = find(isfinite(f1Vals));
    if isempty(validIdx)
        appendTrainLog(fig, sprintf('[Channel %d] Sweep report unavailable (no valid scored combinations).', channelIdx));
        return;
    end

    [~, orderRel] = sort(f1Vals(validIdx), 'descend');
    rankedIdx = validIdx(orderRel);
    topN = min(5, numel(rankedIdx));

    appendTrainLog(fig, sprintf('[Channel %d] Sweep performance (%d combinations, top %d by F1-real):', ...
        channelIdx, nCombos, topN));
    for r = 1:topN
        idx = rankedIdx(r);
        params = sweepResults(idx).params;
        appendTrainLog(fig, sprintf('[Channel %d]   #%d: F1=%.4f Acc=%.4f | %s', ...
            channelIdx, r, f1Vals(idx), accVals(idx), formatSweepParameterSummary(params)));
    end

    showSweepPerformanceFigure(channelIdx, sweepResults, f1Vals, accVals, rankedIdx, outputPath);
end

function showSweepPerformanceFigure(channelIdx, sweepResults, f1Vals, accVals, rankedIdx, outputPath)
    reportFig = uifigure('Name', sprintf('SNAP_train Sweep Report - Channel %d', channelIdx), ...
        'Position', [120 120 980 650]);
    reportGrid = uigridlayout(reportFig, [2, 1]);
    reportGrid.RowHeight = {280, '1x'};
    reportGrid.Padding = [10 10 10 10];
    reportGrid.RowSpacing = 8;

    ax = uiaxes(reportGrid);
    comboIdx = 1:numel(sweepResults);
    plot(ax, comboIdx, f1Vals, '-o', 'LineWidth', 1.2, 'Color', [0 0.45 0.74], 'MarkerSize', 4);
    hold(ax, 'on');
    plot(ax, comboIdx, accVals, '-s', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.1], 'MarkerSize', 4);
    if ~isempty(rankedIdx)
        bestIdx = rankedIdx(1);
        scatter(ax, bestIdx, f1Vals(bestIdx), 70, [0.47 0.67 0.19], 'filled');
    end
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'Sweep Combination Index');
    ylabel(ax, 'Metric Value');
    [~, outName, outExt] = fileparts(outputPath);
    title(ax, sprintf('Channel %d Validation Sweep (%s%s)', channelIdx, outName, outExt), 'Interpreter', 'none');
    legend(ax, {'F1 (real)', 'Accuracy', 'Best F1'}, 'Location', 'best');

    topN = min(20, numel(rankedIdx));
    tableData = cell(topN, 7);
    for r = 1:topN
        idx = rankedIdx(r);
        p = sweepResults(idx).params;
        tableData{r, 1} = r;
        tableData{r, 2} = char(string(getfieldwithdefault(p, 'kernelFunction', 'rbf')));
        tableData{r, 3} = getfieldwithdefault(p, 'boxConstraint', NaN);
        tableData{r, 4} = formatSweepScale(getfieldwithdefault(p, 'kernelScale', 'auto'));
        tableData{r, 5} = getfieldwithdefault(p, 'polynomialOrder', NaN);
        tableData{r, 6} = f1Vals(idx);
        tableData{r, 7} = accVals(idx);
    end

    table = uitable(reportGrid);
    table.ColumnName = {'Rank', 'Kernel', 'Box C', 'Scale', 'Poly', 'F1 (real)', 'Accuracy'};
    table.ColumnEditable = false(1, 7);
    table.ColumnWidth = {55, 90, 90, 90, 70, 100, 90};
    table.Data = tableData;
end

function summary = formatSweepParameterSummary(params)
    kernel = char(string(getfieldwithdefault(params, 'kernelFunction', 'rbf')));
    boxC = getfieldwithdefault(params, 'boxConstraint', NaN);
    scale = getfieldwithdefault(params, 'kernelScale', 'auto');
    poly = getfieldwithdefault(params, 'polynomialOrder', NaN);
    summary = sprintf('kernel=%s, C=%.4g, scale=%s, poly=%g', ...
        kernel, boxC, formatSweepScale(scale), poly);
end

function out = formatSweepScale(scale)
    if ischar(scale) || (isstring(scale) && isscalar(scale))
        out = char(string(scale));
    elseif isnumeric(scale) && isscalar(scale) && isfinite(scale)
        out = sprintf('%.4g', scale);
    else
        out = 'auto';
    end
end

function value = getfieldwithdefault(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function [params, numChannels, errMsg] = loadTrainingParameters(paramPath)
    params = struct();
    numChannels = 1;
    errMsg = '';

    if isempty(paramPath) || exist(paramPath, 'file') ~= 2
        errMsg = sprintf('Parameter file not found: %s', paramPath);
        return;
    end

    try
        [params, ~] = snap_helpers.classification.loadParameterStruct(paramPath);
        caps = snap_helpers.classification.resolveChannelCapabilities(params, 'IncludeFeatureInfo', false);
        if isempty(caps)
            numChannels = inferNumChannelsFromParameters(params, 1);
        else
            numChannels = numel(caps);
        end
    catch ME
        loc = '';
        if ~isempty(ME.stack)
            loc = sprintf(' [%s:%d]', ME.stack(1).name, ME.stack(1).line);
        end
        errMsg = sprintf('Failed to load parameter file: %s%s', ME.message, loc);
    end
end

function n = inferNumChannelsFromParameters(params, fallback)
    if ~isstruct(params)
        n = fallback;
        return;
    end

    explicitN = nan;
    if isfield(params, 'numChannels')
        explicitN = normalizeScalarNumeric(params.numChannels, nan);
    elseif isfield(params, 'numChan')
        explicitN = normalizeScalarNumeric(params.numChan, nan);
    elseif isfield(params, 'numChanDrop')
        explicitN = normalizeScalarNumeric(params.numChanDrop, nan);
    elseif isfield(params, 'workflowConfig') && isstruct(params.workflowConfig) && isfield(params.workflowConfig, 'numChannels')
        explicitN = normalizeScalarNumeric(params.workflowConfig.numChannels, nan);
    end

    if isfinite(explicitN) && explicitN >= 1
        n = max(1, round(explicitN));
        return;
    end

    n = fallback;
    channelCountHints = {'gaussFitMethod', 'maximaMode', 'preProcMode', ...
        'channelPath', 'channelPaths', 'classifyClassifierPath', 'classifyEnabled'};

    for i = 1:numel(channelCountHints)
        fieldName = channelCountHints{i};
        if isfield(params, fieldName)
            n = max(n, inferChannelCountFromValue(params.(fieldName)));
        end
    end

    n = max(1, round(n));
end

function n = inferChannelCountFromValue(v)
    n = 1;
    if isempty(v)
        return;
    end
    if iscell(v) || isstring(v)
        n = numel(v);
    elseif islogical(v) || isnumeric(v)
        if isscalar(v)
            n = 1;
        elseif isvector(v)
            n = numel(v);
        end
    elseif ischar(v)
        n = 1;
    end
end

function options = buildManualTrainingOptions(handles, channelIdx)
    options = struct();
    if nargin < 2 || isempty(channelIdx) || ~isfinite(channelIdx)
        channelIdx = max(1, getfieldwithdefault(handles, 'selectedChannelRow', 1));
    end
    channelIdx = max(1, round(channelIdx));

    useTabControls = isfield(handles, 'channelTabs') && channelIdx <= numel(handles.channelTabs) && ...
        isfield(handles.channelTabs(channelIdx), 'kernelDrop') && isgraphics(handles.channelTabs(channelIdx).kernelDrop);

    if useTabControls
        tab = handles.channelTabs(channelIdx);
        kernelFunction = char(string(tab.kernelDrop.Value));
        boxConstraint = tab.boxConstraintInput.Value;
        kernelScaleText = tab.kernelScaleInput.Value;
        polynomialOrder = round(tab.polyOrderInput.Value);
        standardize = logical(tab.standardizeCheck.Value);
        crossValidate = logical(tab.crossValidateCheck.Value);
        kFold = max(2, round(tab.kFoldInput.Value));
        balanceClassBins = logical(tab.balanceClassCheck.Value);
    else
        advCfg = getChannelAdvancedConfig(handles, channelIdx);
        kernelFunction = char(string(advCfg.kernelFunction));
        boxConstraint = advCfg.boxConstraint;
        kernelScaleText = advCfg.kernelScale;
        polynomialOrder = round(advCfg.polynomialOrder);
        standardize = logical(advCfg.standardize);
        crossValidate = logical(advCfg.crossValidate);
        kFold = max(2, round(advCfg.kFold));
        balanceClassBins = logical(advCfg.balanceClassBins);
    end

    options.kernelFunction = kernelFunction;
    options.boxConstraint = boxConstraint;
    options.kernelScale = parseKernelScale(kernelScaleText);
    options.polynomialOrder = polynomialOrder;
    options.standardize = standardize;
    options.crossValidate = crossValidate;
    options.kFold = kFold;
    if balanceClassBins
        options.classWeightMode = 'balanced';
    else
        options.classWeightMode = 'none';
    end
    options.verbose = true;
end

function sweep = buildSweepConfig(handles, channelIdx)
    sweep = struct();
    if nargin < 2 || isempty(channelIdx) || ~isfinite(channelIdx)
        channelIdx = max(1, getfieldwithdefault(handles, 'selectedChannelRow', 1));
    end
    channelIdx = max(1, round(channelIdx));

    useTabControls = isfield(handles, 'channelTabs') && channelIdx <= numel(handles.channelTabs) && ...
        isfield(handles.channelTabs(channelIdx), 'sweepKernelsInput') && isgraphics(handles.channelTabs(channelIdx).sweepKernelsInput);

    if useTabControls
        tab = handles.channelTabs(channelIdx);
        kernelsRaw = tab.sweepKernelsInput.Value;
        boxRaw = tab.sweepBoxInput.Value;
        scaleRaw = tab.sweepScaleInput.Value;
        polyRaw = tab.sweepPolyInput.Value;
        tieRaw = tab.sweepTieHandlingDrop.Value;
    else
        advCfg = getChannelAdvancedConfig(handles, channelIdx);
        kernelsRaw = advCfg.sweepKernels;
        boxRaw = advCfg.sweepBoxConstraints;
        scaleRaw = advCfg.sweepKernelScales;
        polyRaw = advCfg.sweepPolynomialOrders;
        tieRaw = advCfg.sweepTieHandlingLabel;
    end

    sweep.kernels = parseKernelList(kernelsRaw);
    sweep.boxConstraints = parseNumericList(boxRaw, true);
    sweep.kernelScales = parseMixedScaleList(scaleRaw);
    sweep.polynomialOrders = round(parseNumericList(polyRaw, true));
    sweep.tieHandling = parseSweepTieHandlingChoice(tieRaw);
    sweep.tiePromptCallback = [];
end

function scale = parseKernelScale(valueText)
    token = lower(strtrim(char(string(valueText))));
    if isempty(token) || strcmp(token, 'auto')
        scale = 'auto';
        return;
    end
    val = str2double(token);
    if ~isfinite(val) || val <= 0
        error('Kernel scale must be "auto" or a positive number.');
    end
    scale = val;
end

function kernels = parseKernelList(textValue)
    raw = splitListTokens(textValue);
    if isempty(raw)
        kernels = {'linear', 'rbf', 'polynomial'};
        return;
    end
    valid = {'linear', 'rbf', 'polynomial', 'gaussian'};
    kernels = {};
    for i = 1:numel(raw)
        token = lower(strtrim(raw{i}));
        if ~ismember(token, valid)
            error('Invalid kernel "%s". Allowed: linear, rbf, polynomial, gaussian.', token);
        end
        kernels{end+1} = token; %#ok<AGROW>
    end
end

function nums = parseNumericList(textValue, requirePositive)
    if nargin < 2
        requirePositive = false;
    end
    raw = splitListTokens(textValue);
    if isempty(raw)
        nums = [];
        return;
    end
    nums = zeros(1, numel(raw));
    for i = 1:numel(raw)
        val = str2double(raw{i});
        if ~isfinite(val)
            error('Invalid numeric token: "%s"', raw{i});
        end
        if requirePositive && val <= 0
            error('Values must be > 0: "%s"', raw{i});
        end
        nums(i) = val;
    end
end

function mixed = parseMixedScaleList(textValue)
    raw = splitListTokens(textValue);
    if isempty(raw)
        mixed = {'auto', 0.5, 1, 2};
        return;
    end
    mixed = cell(1, numel(raw));
    for i = 1:numel(raw)
        token = lower(strtrim(raw{i}));
        if strcmp(token, 'auto')
            mixed{i} = 'auto';
        else
            val = str2double(token);
            if ~isfinite(val) || val <= 0
                error('Kernel scales must be "auto" or positive numeric values.');
            end
            mixed{i} = val;
        end
    end
end

function tokens = splitListTokens(textValue)
    s = char(string(textValue));
    parts = regexp(s, '[,;\s]+', 'split');
    parts = parts(~cellfun(@isempty, parts));
    tokens = parts(:)';
end

function mode = parseSweepTieHandlingChoice(valueText)
    mode = normalizeSweepTieHandling(valueText);
end

function mode = normalizeSweepTieHandling(valueText)
    raw = lower(strtrim(char(string(valueText))));
    if isempty(raw)
        mode = 'prefer_simple';
        return;
    end

    if strcmp(raw, 'prompt') || contains(raw, 'prompt')
        mode = 'prompt';
    elseif strcmp(raw, 'first') || contains(raw, 'legacy') || contains(raw, 'first')
        mode = 'first';
    elseif strcmp(raw, 'prefer_simple') || contains(raw, 'simple')
        mode = 'prefer_simple';
    else
        error('Invalid SweepTieHandling value: "%s". Use "prefer_simple", "prompt", or "first".', raw);
    end
end

function txt = describeSweepTieHandlingPolicy(mode)
    switch normalizeSweepTieHandling(mode)
        case 'prompt'
            txt = 'Prompt user to choose tied best settings';
        case 'first'
            txt = 'Keep first tied setting (legacy)';
        otherwise
            txt = 'Prefer simpler model among tied best settings';
    end
end

function selectedIdx = promptSweepTieSelection(fig, channelIdx, tiedIdx, sweepResults)
    selectedIdx = [];
    if isempty(tiedIdx)
        return;
    end
    if numel(tiedIdx) == 1
        selectedIdx = tiedIdx(1);
        return;
    end

    listStrings = cell(numel(tiedIdx), 1);
    for i = 1:numel(tiedIdx)
        idx = tiedIdx(i);
        f1 = getfieldwithdefault(sweepResults(idx), 'f1Real', NaN);
        acc = getfieldwithdefault(sweepResults(idx), 'accuracy', NaN);
        params = getfieldwithdefault(sweepResults(idx), 'params', struct());
        listStrings{i} = sprintf('#%d | F1=%.4f Acc=%.4f | %s', ...
            idx, f1, acc, formatSweepParameterSummary(params));
    end

    promptLines = { ...
        sprintf('Channel %d has multiple tied best sweep settings.', channelIdx), ...
        'Choose which hyperparameter set to use:' ...
    };

    try
        [sel, ok] = listdlg( ...
            'PromptString', promptLines, ...
            'SelectionMode', 'single', ...
            'ListString', listStrings, ...
            'InitialValue', 1, ...
            'ListSize', [920 280], ...
            'Name', sprintf('SNAP_train Channel %d Tie Selection', channelIdx));
        if ok && ~isempty(sel)
            selectedIdx = tiedIdx(sel(1));
        end
    catch
        selectedIdx = [];
    end

    if isempty(selectedIdx)
        try
            choice = uiconfirm(fig, ...
                sprintf(['Channel %d has tied best settings.\n' ...
                'Prompt selection was cancelled/unavailable.\n' ...
                'How should SNAP_train continue?'], channelIdx), ...
                'Tie Selection', ...
                'Options', {'Prefer simpler model', 'Use first tied setting'}, ...
                'DefaultOption', 1, ...
                'CancelOption', 1);
            if strcmp(choice, 'Use first tied setting')
                selectedIdx = tiedIdx(1);
            else
                selectedIdx = chooseSimplestSweepCandidate(tiedIdx, sweepResults);
            end
        catch
            selectedIdx = [];
        end
    end
end

function ch = parseChannelFromRowValue(value)
    if isnumeric(value) || islogical(value)
        if isscalar(value) && isfinite(double(value)) && value >= 1
            ch = round(double(value));
            return;
        end
    end

    token = regexp(char(string(value)), '\d+', 'match', 'once');
    if isempty(token)
        error('Could not parse channel index from table row: %s', char(string(value)));
    end
    ch = str2double(token);
    if ~isfinite(ch) || ch < 1
        error('Invalid channel index: %s', char(string(value)));
    end
end

function tf = isSelectedChannelFlag(value)
    tf = false;
    if islogical(value) || isnumeric(value)
        if isscalar(value)
            val = double(value);
            tf = isfinite(val) && (val ~= 0);
        end
    end
end

function [fittingMethod, has3D] = inferChannelTrainingContext(params, ch)
    fittingMethod = '';
    has3D = [];

    if isfield(params, 'gaussFitMethod')
        m = getChannelCellValue(params.gaussFitMethod, ch, '');
        if ~isempty(m)
            fittingMethod = char(string(m));
        end
    end

    if isfield(params, 'maximaMode')
        modeVal = getChannelCellValue(params.maximaMode, ch, '');
        if ~isempty(modeVal)
            has3D = strcmp(char(string(modeVal)), '3D');
        end
    elseif isfield(params, 'preProcMode')
        modeVal = getChannelCellValue(params.preProcMode, ch, '');
        if ~isempty(modeVal)
            has3D = strcmp(char(string(modeVal)), '3D');
        end
    end
end

function v = getChannelCellValue(raw, ch, fallback)
    v = fallback;
    if isempty(raw)
        return;
    end

    if iscell(raw)
        if ch <= numel(raw)
            v = raw{ch};
        end
    elseif isstring(raw)
        if ch <= numel(raw)
            v = char(raw(ch));
        end
    elseif ischar(raw)
        v = raw;
    elseif isnumeric(raw) || islogical(raw)
        if ch <= numel(raw)
            v = raw(ch);
        end
    end
end

function val = normalizeScalarNumeric(raw, fallback)
    val = fallback;
    if isstruct(raw)
        if isfield(raw, 'Value')
            raw = raw.Value;
        else
            return;
        end
    end
    if iscell(raw)
        if isempty(raw)
            return;
        end
        raw = raw{1};
    elseif isstring(raw)
        if isempty(raw)
            return;
        end
        raw = raw(1);
    end

    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
        val = raw;
    elseif ischar(raw) || isstring(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            val = parsed;
        end
    end
end

function [exportFiles, labelFiles] = discoverDatasetPairs(rootDir, channelIdx, varargin)
    if nargin < 2
        channelIdx = [];
    end
    progressCb = [];
    if nargin >= 3 && ~isempty(varargin)
        p = inputParser;
        p.addParameter('ProgressCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
        p.parse(varargin{:});
        progressCb = p.Results.ProgressCallback;
    end

    if ~isfolder(rootDir)
        error('Directory not found: %s', rootDir);
    end

    % Preferred path: image+CSV pairs with matching base names.
    if ~isempty(progressCb)
        progressCb(sprintf('Scanning "%s" for image+CSV pairs...', rootDir));
        if ~isempty(channelIdx)
            progressCb(sprintf('Channel %d pair discovery note: filename channel tags are ignored.', channelIdx));
        end
    end
    [exportFiles, labelFiles] = discoverImageCsvPairs(rootDir, channelIdx, progressCb);
    if ~isempty(exportFiles)
        if ~isempty(progressCb)
            progressCb(sprintf('Image+CSV discovery complete: %d pair(s).', numel(exportFiles)));
        end
        return;
    end

    % Backward-compatible fallback: legacy export MAT + label file pairs.
    if ~isempty(progressCb)
        progressCb('No image+CSV pairs found. Falling back to legacy MAT discovery...');
    end
    [exportFiles, labelFiles] = discoverLegacyMatPairs(rootDir, channelIdx, progressCb);
    if ~isempty(progressCb)
        progressCb(sprintf('Legacy discovery complete: %d pair(s).', numel(exportFiles)));
    end
end

function [exportFiles, labelFiles] = discoverImageCsvPairs(rootDir, ~, progressCb)
    if nargin < 3
        progressCb = [];
    end
    exportFiles = {};
    labelFiles = {};

    imageFiles = [ ...
        collectFilesRecursive(rootDir, '.tif', progressCb, 'images (.tif)'); ...
        collectFilesRecursive(rootDir, '.tiff', progressCb, 'images (.tiff)') ...
    ];
    if ~isempty(progressCb)
        progressCb(sprintf('Found %d candidate image file(s).', numel(imageFiles)));
    end
    if isempty(imageFiles)
        return;
    end

    csvFiles = collectFilesRecursive(rootDir, '.csv', progressCb, 'label files (.csv)');
    if ~isempty(progressCb)
        progressCb(sprintf('Found %d candidate CSV file(s).', numel(csvFiles)));
    end
    if isempty(csvFiles)
        return;
    end

    validCsvMask = cellfun(@isLabelCsvFile, csvFiles);
    csvFiles = csvFiles(validCsvMask);
    if ~isempty(progressCb)
        progressCb(sprintf('%d CSV file(s) passed label-header checks.', numel(csvFiles)));
    end
    if isempty(csvFiles)
        return;
    end

    if ~isempty(progressCb)
        progressCb(sprintf('Matching %d image(s) to %d CSV file(s) (same-name rule)...', ...
            numel(imageFiles), numel(csvFiles)));
    end
    exportFiles = {};
    labelFiles = {};

    for i = 1:numel(imageFiles)
        imagePath = imageFiles{i};
        labelPath = findExactCsvForImageInList(imagePath, csvFiles);
        if isempty(labelPath)
            if ~isempty(progressCb) && (numel(imageFiles) <= 20 || mod(i, 10) == 0)
                [~, imgBase, imgExt] = fileparts(imagePath);
                progressCb(sprintf('No CSV match for image %d/%d: %s%s', i, numel(imageFiles), imgBase, imgExt));
            end
            continue;
        end

        exportFiles{end+1,1} = imagePath; %#ok<AGROW>
        labelFiles{end+1,1} = labelPath; %#ok<AGROW>

        if ~isempty(progressCb) && (numel(imageFiles) <= 20 || mod(i, 10) == 0)
            [~, imgBase, imgExt] = fileparts(imagePath);
            [~, labBase, labExt] = fileparts(labelPath);
            progressCb(sprintf('Matched image %d/%d: %s%s <-> %s%s', ...
                i, numel(imageFiles), imgBase, imgExt, labBase, labExt));
        end
    end
    if ~isempty(progressCb)
        progressCb(sprintf('Pairing loop complete: %d pair(s).', numel(exportFiles)));
    end
end

function labelPath = findExactCsvForImageInList(imagePath, csvFiles)
    labelPath = '';
    [imageDir, imageBase, ~] = fileparts(imagePath);
    imageBaseNoOme = regexprep(imageBase, '\.ome$', '', 'ignorecase');

    sameDirFallbackPath = '';
    globalExactPath = '';
    globalExactCount = 0;
    globalNoOmePath = '';
    globalNoOmeCount = 0;

    for i = 1:numel(csvFiles)
        csvPath = csvFiles{i};
        [csvDir, csvBase, ~] = fileparts(csvPath);

        if strcmpi(csvBase, imageBase)
            globalExactCount = globalExactCount + 1;
            if globalExactCount == 1
                globalExactPath = csvPath;
            end
            if strcmpi(csvDir, imageDir)
                labelPath = csvPath;
                return;
            end
        end

        csvBaseNoOme = regexprep(csvBase, '\.ome$', '', 'ignorecase');
        if strcmpi(csvBaseNoOme, imageBaseNoOme)
            globalNoOmeCount = globalNoOmeCount + 1;
            if globalNoOmeCount == 1
                globalNoOmePath = csvPath;
            end
            if isempty(sameDirFallbackPath) && strcmpi(csvDir, imageDir)
                sameDirFallbackPath = csvPath;
            end
        end
    end

    if ~isempty(sameDirFallbackPath)
        labelPath = sameDirFallbackPath;
        return;
    end

    % If labels are stored in a separate folder, allow unique global base-name matches.
    if globalExactCount == 1
        labelPath = globalExactPath;
        return;
    end
    if globalNoOmeCount == 1
        labelPath = globalNoOmePath;
    end
end

function [exportFiles, labelFiles] = discoverLegacyMatPairs(rootDir, ~, progressCb)
    if nargin < 3
        progressCb = [];
    end
    matFiles = collectFilesRecursive(rootDir, '.mat', progressCb, 'legacy MAT files');
    csvFiles = collectFilesRecursive(rootDir, '.csv', progressCb, 'legacy CSV files');
    if ~isempty(progressCb)
        progressCb(sprintf('Legacy candidates: %d MAT, %d CSV file(s).', numel(matFiles), numel(csvFiles)));
    end
    labelCandidates = [csvFiles; matFiles];
    if ~isempty(labelCandidates)
        validMask = cellfun(@isValidLabelFile, labelCandidates);
        labelCandidates = labelCandidates(validMask);
        if ~isempty(progressCb)
            progressCb(sprintf('Legacy label candidates after validation: %d file(s).', numel(labelCandidates)));
        end
    end
    exportFiles = {};
    labelFiles = {};

    for i = 1:numel(matFiles)
        matPath = matFiles{i};
        if ~isExportFile(matPath)
            continue;
        end

        labelPath = findBestLabelForExport(matPath, labelCandidates, []);
        if ~isempty(labelPath)
            exportFiles{end+1,1} = matPath; %#ok<AGROW>
            labelFiles{end+1,1} = labelPath; %#ok<AGROW>
        end
    end

    if ~isempty(progressCb)
        progressCb(sprintf('Legacy pairing summary: %d pair(s).', numel(exportFiles)));
    end
end

function tf = isExportFile(matPath)
    tf = false;
    try
        info = whos('-file', matPath);
    catch
        return;
    end
    names = string({info.name});
    tf = any(ismember(names, ["signals", "fitResults", "fit_results", "exportData"]));
end

function tf = isLabelMatFile(matPath)
    tf = false;
    if ~endsWith(lower(matPath), '.mat')
        return;
    end
    if isExportFile(matPath)
        return;
    end
    try
        info = whos('-file', matPath);
    catch
        return;
    end
    names = string({info.name});
    tf = any(names == "labeledReal") && any(names == "labeledNoise");
end

function tf = isLabelCsvFile(csvPath)
    tf = false;
    if ~endsWith(lower(csvPath), '.csv')
        return;
    end
    try
        opts = detectImportOptions(csvPath, 'NumHeaderLines', 0);
        names = canonicalizeVarNames(opts.VariableNames);
        hasX = hasAnyCanonicalName(names, ["maximax", "fittedx", "x", "col", "column"]);
        hasY = hasAnyCanonicalName(names, ["maximay", "fittedy", "y", "row"]);
        tf = hasX && hasY;
    catch
        tf = false;
    end
end

function names = canonicalizeVarNames(rawNames)
    names = lower(string(rawNames));
    names = regexprep(names, '[^a-z0-9]', '');
end

function tf = hasAnyCanonicalName(allNames, candidateNames)
    tf = any(ismember(string(allNames), string(candidateNames)));
end

function files = collectFilesRecursive(rootDir, extension, progressCb, labelText)
    if nargin < 3
        progressCb = [];
    end
    if nargin < 4 || isempty(labelText)
        labelText = extension;
    end
    files = {};
    queue = {rootDir};
    extension = lower(extension);
    scannedDirs = 0;
    scannedFiles = 0;

    while ~isempty(queue)
        currentDir = queue{1};
        queue(1) = [];
        scannedDirs = scannedDirs + 1;
        listing = dir(currentDir);

        for i = 1:numel(listing)
            name = listing(i).name;
            if strcmp(name, '.') || strcmp(name, '..')
                continue;
            end
            fullPath = fullfile(currentDir, name);
            if listing(i).isdir
                queue{end+1} = fullPath; %#ok<AGROW>
            else
                scannedFiles = scannedFiles + 1;
                [~, ~, ext] = fileparts(name);
                if strcmpi(ext, extension)
                    files{end+1,1} = fullPath; %#ok<AGROW>
                end
            end
        end

        if ~isempty(progressCb) && (mod(scannedDirs, 50) == 0)
            progressCb(sprintf('Scanning %s: %d folder(s), %d file(s) checked, %d match(es) so far.', ...
                labelText, scannedDirs, scannedFiles, numel(files)));
        end
    end

    if ~isempty(files)
        files = sort(files);
    end
    if ~isempty(progressCb)
        progressCb(sprintf('Finished scanning %s: %d match(es).', labelText, numel(files)));
    end
end

function labelPath = findBestLabelForExport(exportPath, labelCandidates, channelIdx)
    [exportDir, exportBase, ~] = fileparts(exportPath);
    baseKeys = deriveExportBaseKeys(exportBase);
    expectedPaths = buildExpectedLabelPaths(exportDir, baseKeys, channelIdx);

    labelPath = '';
    for i = 1:numel(expectedPaths)
        candidate = expectedPaths{i};
        if strcmp(candidate, exportPath) || exist(candidate, 'file') ~= 2
            continue;
        end
        if isValidLabelFile(candidate)
            labelPath = candidate;
            return;
        end
    end

    if isempty(labelCandidates)
        return;
    end

    sameDirMask = cellfun(@(p) strcmpi(fileparts(p), exportDir), labelCandidates);
    labelPath = chooseBestLabelCandidate(labelCandidates(sameDirMask), exportDir, baseKeys, channelIdx);
    if ~isempty(labelPath)
        return;
    end

    labelPath = chooseBestLabelCandidate(labelCandidates, exportDir, baseKeys, channelIdx);
end

function keys = deriveExportBaseKeys(baseName)
    keys = {baseName};
    keys{end+1} = regexprep(baseName, '_signals?$', '', 'ignorecase'); %#ok<AGROW>
    keys{end+1} = regexprep(baseName, '_export$', '', 'ignorecase'); %#ok<AGROW>
    keys{end+1} = regexprep(baseName, '_ch\d+_signals?$', '', 'ignorecase'); %#ok<AGROW>
    keys{end+1} = regexprep(baseName, '_ch\d+$', '', 'ignorecase'); %#ok<AGROW>
    keys{end+1} = regexprep(baseName, '_channel\d+_signals?$', '', 'ignorecase'); %#ok<AGROW>
    keys{end+1} = regexprep(baseName, '_channel\d+$', '', 'ignorecase'); %#ok<AGROW>

    keepMask = ~cellfun(@isempty, keys);
    keys = keys(keepMask);
    keys = unique(keys, 'stable');
end

function paths = buildExpectedLabelPaths(exportDir, baseKeys, channelIdx)
    expectedNames = {};
    for i = 1:numel(baseKeys)
        k = baseKeys{i};
        expectedNames{end+1} = [k '_labels.csv']; %#ok<AGROW>
        expectedNames{end+1} = [k '_labels.mat']; %#ok<AGROW>
        expectedNames{end+1} = [k '.csv']; %#ok<AGROW>
        expectedNames{end+1} = [k '.mat']; %#ok<AGROW>
        expectedNames{end+1} = [k '_progress.mat']; %#ok<AGROW>
        expectedNames{end+1} = [k '_classification_progress.mat']; %#ok<AGROW>
    end

    if ~isempty(channelIdx)
        expectedNames{end+1} = sprintf('ch%d_labels.csv', channelIdx); %#ok<AGROW>
        expectedNames{end+1} = sprintf('ch%d_labels.mat', channelIdx); %#ok<AGROW>
        expectedNames{end+1} = sprintf('channel%d_labels.csv', channelIdx); %#ok<AGROW>
        expectedNames{end+1} = sprintf('channel%d_labels.mat', channelIdx); %#ok<AGROW>
    end

    expectedNames = unique(expectedNames, 'stable');
    paths = cell(numel(expectedNames), 1);
    for i = 1:numel(expectedNames)
        paths{i} = fullfile(exportDir, expectedNames{i});
    end
end

function pathOut = chooseBestLabelCandidate(candidates, exportDir, baseKeys, channelIdx)
    pathOut = '';
    if isempty(candidates)
        return;
    end

    scores = -inf(numel(candidates), 1);
    for i = 1:numel(candidates)
        candidate = candidates{i};
        scores(i) = scoreLabelCandidate(candidate, exportDir, baseKeys, channelIdx);
    end

    [bestScore, idx] = max(scores);
    if isfinite(bestScore) && bestScore > 0
        pathOut = candidates{idx};
    end
end

function score = scoreLabelCandidate(candidatePath, exportDir, baseKeys, channelIdx)
    [candDir, candBase, candExt] = fileparts(candidatePath);
    candLower = lower(candBase);
    score = 0;

    if strcmpi(candDir, exportDir)
        score = score + 20;
    end

    if strcmpi(candExt, '.csv')
        score = score + 3;
    elseif strcmpi(candExt, '.mat')
        score = score + 1;
    end

    if contains(candLower, 'label')
        score = score + 6;
    end
    if contains(candLower, 'classifier') || contains(candLower, 'summary') || ...
            contains(candLower, 'parameter') || contains(candLower, 'receipt')
        score = score - 50;
    end

    for i = 1:numel(baseKeys)
        keyLower = lower(baseKeys{i});
        if strcmp(candLower, [keyLower '_labels'])
            score = score + 50;
        elseif strcmp(candLower, keyLower)
            score = score + 40;
        elseif startsWith(candLower, [keyLower '_'])
            score = score + 15;
        elseif contains(candLower, keyLower)
            score = score + 8;
        end
    end

    if ~isempty(channelIdx)
        candChannel = parseChannelFromFilename(candidatePath);
        if isfinite(candChannel)
            if candChannel == channelIdx
                score = score + 15;
            else
                score = score - 35;
            end
        elseif channelIdx > 1
            score = score - 5;
        end
    end
end

function tf = isValidLabelFile(filePath)
    tf = false;
    if exist(filePath, 'file') ~= 2
        return;
    end
    [~, ~, ext] = fileparts(filePath);
    ext = lower(ext);
    if strcmp(ext, '.csv')
        tf = isLabelCsvFile(filePath);
    elseif strcmp(ext, '.mat')
        tf = isLabelMatFile(filePath);
    end
end

function ch = parseChannelFromFilename(pathValue)
    [~, base, ~] = fileparts(char(string(pathValue)));
    base = lower(base);
    ch = nan;

    patterns = { ...
        '(?:^|[_\-])ch(?:annel)?[_\-]?(\d+)(?:$|[_\-])', ...
        '(?:^|[_\-])channel[_\-]?(\d+)(?:$|[_\-])', ...
        '(?:^|[_\-])c0*(\d+)(?:$|[_\-])' ...
    };

    for i = 1:numel(patterns)
        tok = regexp(base, patterns{i}, 'tokens', 'once');
        if ~isempty(tok)
            val = str2double(tok{1});
            if isfinite(val) && val >= 1
                ch = val;
                return;
            end
        end
    end
end

function files = normalizeFileList(x)
    if ischar(x) || (isstring(x) && isscalar(x))
        files = {char(x)};
    elseif isstring(x)
        files = cellstr(x(:));
    elseif iscell(x)
        files = x(:);
    else
        error('File input must be char/string/cellstr.');
    end
end

function files = normalizeOptionalFileList(x)
    if isempty(x)
        files = {};
        return;
    end
    files = normalizeFileList(x);
end

function emitProgress(progressCb, msgFmt, varargin)
    if isempty(progressCb) || ~isa(progressCb, 'function_handle')
        return;
    end

    try
        msg = sprintf(msgFmt, varargin{:});
    catch
        msg = char(string(msgFmt));
    end

    try
        progressCb(msg);
    catch
        % Progress callbacks must never stop training.
    end
end

function [allFitData, allLabels, sourceImage, inferredFittingMethod, inferredHas3D] = ...
        buildLabeledDataset(exportFiles, labelFiles, matchDistance, datasetOpts)

    allFitData = struct([]);
    allLabels = [];
    sourceImage = strings(0,1);
    inferredFittingMethod = '';
    inferredHas3D = false;

    if nargin < 4 || isempty(datasetOpts)
        datasetOpts = struct();
    end
    if ~isfield(datasetOpts, 'convertFijiCoords')
        datasetOpts.convertFijiCoords = false;
    end
    if ~isfield(datasetOpts, 'parameterStruct')
        datasetOpts.parameterStruct = struct();
    end
    if ~isfield(datasetOpts, 'channelIndex')
        datasetOpts.channelIndex = [];
    end
    if ~isfield(datasetOpts, 'progressCallback')
        datasetOpts.progressCallback = [];
    end
    if ~isfield(datasetOpts, 'datasetLabel') || isempty(datasetOpts.datasetLabel)
        datasetOpts.datasetLabel = 'dataset';
    end

    progressCb = datasetOpts.progressCallback;
    datasetLabel = char(string(datasetOpts.datasetLabel));
    nPairs = numel(exportFiles);
    emitProgress(progressCb, '[%s] Starting dataset build for %d pair(s).', datasetLabel, nPairs);

    for i = 1:nPairs
        [~, exportName, exportExt] = fileparts(exportFiles{i});
        [~, labelName, labelExt] = fileparts(labelFiles{i});
        emitProgress(progressCb, '[%s] Pair %d/%d: %s%s <-> %s%s', ...
            datasetLabel, i, nPairs, exportName, exportExt, labelName, labelExt);

        fitData = loadFitDataFromTrainingInput( ...
            exportFiles{i}, datasetOpts.parameterStruct, datasetOpts.channelIndex, progressCb);

        if isempty(fitData)
            emitProgress(progressCb, '[%s] Pair %d/%d produced 0 candidates and was skipped.', ...
                datasetLabel, i, nPairs);
            continue;
        end

        labels = loadLabelsForFitData( ...
            labelFiles{i}, fitData, matchDistance, logical(datasetOpts.convertFijiCoords));
        labels(isnan(labels)) = 0;
        emitProgress(progressCb, '[%s] Pair %d/%d labeled candidates: %d (real=%d, noise=%d).', ...
            datasetLabel, i, nPairs, numel(labels), sum(labels == 1), sum(labels == 0));

        if isempty(inferredFittingMethod)
            inferredFittingMethod = inferFittingMethod(fitData);
        end
        inferredHas3D = inferredHas3D || inferHas3D(fitData);

        if isempty(allFitData)
            allFitData = fitData(:);
        else
            allFitData = [allFitData(:); fitData(:)]; %#ok<AGROW>
        end
        allLabels = [allLabels; labels(:)]; %#ok<AGROW>
        sourceImage = [sourceImage; repmat(string(exportFiles{i}), numel(labels), 1)]; %#ok<AGROW>
    end

    emitProgress(progressCb, '[%s] Dataset build complete: %d candidate(s) from %d image file(s).', ...
        datasetLabel, numel(allLabels), numel(unique(sourceImage)));
end

function [model, trainStats, normParams, bestParams, results] = runHyperparameterSweep( ...
        XTrain, yTrain, fixedTrainLabels, XVal, yVal, baseOptions, cfg, verbose, progressCb)
    if nargin < 9
        progressCb = [];
    end

    combos = buildSweepCombinations(cfg);
    results = repmat(struct( ...
        'params', struct(), ...
        'fullOptions', struct(), ...
        'f1Real', nan, ...
        'accuracy', nan, ...
        'score', nan, ...
        'trainStats', struct()), numel(combos), 1);
    emitProgress(progressCb, 'Validation sweep: evaluating %d combination(s).', numel(combos));

    policy = normalizeSweepTieHandling(getfieldwithdefault(cfg, 'tieHandling', 'prefer_simple'));
    tieTol = getfieldwithdefault(cfg, 'tieTolerance', 1e-12);
    if ~isnumeric(tieTol) || ~isscalar(tieTol) || ~isfinite(tieTol) || tieTol < 0
        tieTol = 1e-12;
    end
    tiePromptCallback = getfieldwithdefault(cfg, 'tiePromptCallback', []);
    emitProgress(progressCb, 'Validation sweep tie policy: %s.', describeSweepTieHandlingPolicy(policy));

    model = [];
    trainStats = struct();
    normParams = struct();
    bestParams = baseOptions;
    modelCache = cell(numel(combos), 1);
    normCache = cell(numel(combos), 1);

    for i = 1:numel(combos)
        options = baseOptions;
        fields = fieldnames(combos(i));
        for f = 1:numel(fields)
            options.(fields{f}) = combos(i).(fields{f});
        end
        options.crossValidate = false;
        emitProgress(progressCb, 'Validation sweep %d/%d: %s', ...
            i, numel(combos), formatSweepParameterSummary(combos(i)));

        if ~isequal(yTrain(:), fixedTrainLabels(:))
            error('Training labels changed during validation sweep, which is not allowed.');
        end

        [mdl, stats, norm] = snap_helpers.classification.trainClassifier(XTrain, yTrain, options);
        results(i).params = combos(i);
        results(i).fullOptions = options;
        results(i).trainStats = stats;
        if ~isfield(stats, 'success') || ~stats.success
            emitProgress(progressCb, 'Validation sweep %d/%d failed: %s', ...
                i, numel(combos), getfieldwithdefault(stats, 'error', 'unknown training error'));
            continue;
        end

        [metrics, ~] = evaluateOnValidation(mdl, norm, XVal, yVal);

        results(i).f1Real = metrics.f1Real;
        results(i).accuracy = metrics.accuracy;
        results(i).score = metrics.f1Real + 1e-3 * metrics.accuracy;
        modelCache{i} = mdl;
        normCache{i} = norm;
        emitProgress(progressCb, 'Validation sweep %d/%d metrics: F1(real)=%.4f, accuracy=%.4f.', ...
            i, numel(combos), metrics.f1Real, metrics.accuracy);
    end

    scored = arrayfun(@(r) double(r.score), results);
    validIdx = find(isfinite(scored));
    if isempty(validIdx)
        error('Hyperparameter sweep failed: all parameter combinations failed to train.');
    end

    bestScore = max(scored(validIdx));
    tieIdx = validIdx(abs(scored(validIdx) - bestScore) <= tieTol);
    [bestIdx, tieMsg] = resolveSweepWinnerFromTies(tieIdx, results, policy, tiePromptCallback, progressCb);
    if isempty(bestIdx)
        error('Hyperparameter sweep tie handling failed: no winner selected.');
    end
    if ~isempty(tieMsg)
        emitProgress(progressCb, '%s', tieMsg);
    end

    model = modelCache{bestIdx};
    trainStats = results(bestIdx).trainStats;
    normParams = normCache{bestIdx};
    bestParams = results(bestIdx).fullOptions;

    if isempty(model)
        error('Hyperparameter sweep selected combo #%d but its model is unavailable.', bestIdx);
    end

    if verbose
        fprintf('Hyperparameter sweep complete: %d combinations tested.\n', numel(combos));
        fprintf('Best combo #%d: kernel=%s, C=%.3g, f1_real=%.3f\n', ...
            bestIdx, bestParams.kernelFunction, bestParams.boxConstraint, results(bestIdx).f1Real);
    end
    emitProgress(progressCb, 'Validation sweep complete. Best combo #%d with F1(real)=%.4f.', ...
        bestIdx, results(bestIdx).f1Real);
end

function [bestIdx, message] = resolveSweepWinnerFromTies(tieIdx, results, policy, tiePromptCallback, progressCb)
    bestIdx = [];
    message = '';
    if isempty(tieIdx)
        return;
    end
    if numel(tieIdx) == 1
        bestIdx = tieIdx(1);
        return;
    end

    emitProgress(progressCb, 'Validation sweep tie detected: %d combination(s) share the best score.', numel(tieIdx));

    reason = '';
    switch normalizeSweepTieHandling(policy)
        case 'first'
            bestIdx = tieIdx(1);
            reason = 'first encountered (legacy)';
        case 'prompt'
            if isa(tiePromptCallback, 'function_handle')
                try
                    picked = tiePromptCallback(tieIdx, results);
                    if isnumeric(picked) && isscalar(picked) && any(tieIdx == picked)
                        bestIdx = picked;
                        reason = 'user selection';
                    end
                catch ME
                    emitProgress(progressCb, 'Tie prompt callback failed: %s', ME.message);
                end
            end
            if isempty(bestIdx)
                bestIdx = chooseSimplestSweepCandidate(tieIdx, results);
                reason = 'fallback to simpler model (prompt unavailable/cancelled)';
            end
        otherwise
            bestIdx = chooseSimplestSweepCandidate(tieIdx, results);
            reason = 'prefer simpler model';
    end

    message = sprintf('Validation sweep tie resolved: selected combo #%d via %s.', bestIdx, reason);
end

function bestIdx = chooseSimplestSweepCandidate(tieIdx, results)
    n = numel(tieIdx);
    complexity = inf(n, 5);
    for i = 1:n
        idx = tieIdx(i);
        p = getfieldwithdefault(results(idx), 'params', struct());
        kernel = char(string(getfieldwithdefault(p, 'kernelFunction', 'rbf')));
        complexity(i, 1) = sweepKernelComplexityRank(kernel);

        polyOrder = getfieldwithdefault(p, 'polynomialOrder', 3);
        if ~isnumeric(polyOrder) || ~isscalar(polyOrder) || ~isfinite(polyOrder)
            polyOrder = 3;
        end
        if strcmpi(kernel, 'polynomial')
            complexity(i, 2) = double(polyOrder);
        else
            complexity(i, 2) = 0;
        end

        boxC = getfieldwithdefault(p, 'boxConstraint', 1);
        if ~isnumeric(boxC) || ~isscalar(boxC) || ~isfinite(boxC) || boxC <= 0
            boxC = inf;
        end
        complexity(i, 3) = double(boxC);

        complexity(i, 4) = sweepScaleComplexity(getfieldwithdefault(p, 'kernelScale', 'auto'));
        complexity(i, 5) = idx; % deterministic final tie-breaker
    end

    [~, order] = sortrows(complexity, [1 2 3 4 5]);
    bestIdx = tieIdx(order(1));
end

function rankVal = sweepKernelComplexityRank(kernelName)
    kernel = lower(strtrim(char(string(kernelName))));
    switch kernel
        case 'linear'
            rankVal = 1;
        case {'rbf', 'gaussian'}
            rankVal = 2;
        case 'polynomial'
            rankVal = 3;
        otherwise
            rankVal = 4;
    end
end

function c = sweepScaleComplexity(scaleValue)
    if ischar(scaleValue) || isstring(scaleValue)
        token = lower(strtrim(char(string(scaleValue))));
        if strcmp(token, 'auto')
            c = 0;
            return;
        end
        val = str2double(token);
    else
        val = scaleValue;
    end

    if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
        c = inf;
        return;
    end
    c = abs(log(double(val)));
end

function combos = buildSweepCombinations(cfg)
    kernels = normalizeCellstr(cfg.kernels);
    Cs = cfg.boxConstraints;
    scales = normalizeMixedList(cfg.kernelScales);
    polyOrders = cfg.polynomialOrders;

    combos = struct('kernelFunction', {}, 'boxConstraint', {}, 'kernelScale', {}, 'polynomialOrder', {});

    for i = 1:numel(kernels)
        k = kernels{i};
        for c = 1:numel(Cs)
            if strcmp(k, 'linear')
                combos(end+1) = struct('kernelFunction', k, 'boxConstraint', Cs(c), 'kernelScale', 'auto', 'polynomialOrder', 3); %#ok<AGROW>
            elseif strcmp(k, 'polynomial')
                for po = 1:numel(polyOrders)
                    combos(end+1) = struct('kernelFunction', k, 'boxConstraint', Cs(c), 'kernelScale', 'auto', 'polynomialOrder', polyOrders(po)); %#ok<AGROW>
                end
            else
                for s = 1:numel(scales)
                    combos(end+1) = struct('kernelFunction', k, 'boxConstraint', Cs(c), 'kernelScale', scales{s}, 'polynomialOrder', 3); %#ok<AGROW>
                end
            end
        end
    end
end

function [metrics, cm] = evaluateOnValidation(model, normParams, XVal, yVal)
    if isempty(XVal) || isempty(yVal)
        cm = zeros(2, 2);
        metrics = struct('accuracy', NaN, 'precisionReal', NaN, 'recallReal', NaN, 'f1Real', NaN);
        return;
    end

    XNorm = XVal;
    if isfield(normParams, 'standardized') && normParams.standardized
        XNorm = (XVal - normParams.mu) ./ normParams.sigma;
    end

    yPred = predict(model, XNorm);
    yPred = double(yPred(:));
    yVal = double(yVal(:));

    cm = confusionmat(yVal, yPred, 'Order', [0 1]);
    TN = cm(1,1); FP = cm(1,2);
    FN = cm(2,1); TP = cm(2,2);

    precisionReal = TP / max(TP + FP, 1);
    recallReal = TP / max(TP + FN, 1);
    f1Real = 2 * precisionReal * recallReal / max(precisionReal + recallReal, 1e-6);

    metrics = struct();
    metrics.accuracy = mean(yPred == yVal);
    metrics.precisionReal = precisionReal;
    metrics.recallReal = recallReal;
    metrics.f1Real = f1Real;
end

function c = normalizeCellstr(x)
    if ischar(x)
        c = {x};
    elseif isstring(x)
        c = cellstr(x(:));
    elseif iscell(x)
        c = cellfun(@char, x, 'UniformOutput', false);
    else
        error('Expected char/string/cell input.');
    end
end

function out = normalizeMixedList(x)
    if ischar(x)
        out = {x};
    elseif isnumeric(x)
        out = num2cell(x(:)');
    elseif isstring(x)
        out = cellstr(x(:))';
    elseif iscell(x)
        out = x;
    else
        error('Unsupported list input type.');
    end
end

function fitData = loadFitDataFromTrainingInput(filepath, params, channelIdx, progressCb)
    if nargin < 4
        progressCb = [];
    end
    [~, ~, ext] = fileparts(filepath);
    ext = lower(ext);

    if strcmp(ext, '.mat')
        emitProgress(progressCb, 'Loading candidates from export MAT: %s', filepath);
        fitData = loadFitDataFromExportMat(filepath);
        return;
    end

    if ismember(ext, {'.tif', '.tiff'})
        emitProgress(progressCb, 'Generating candidates from image: %s', filepath);
        fitData = buildFitDataFromImage(filepath, params, channelIdx, progressCb);
        return;
    end

    error('Unsupported training input file type: %s', filepath);
end

function fitData = loadFitDataFromExportMat(filepath)
    s = load(filepath);

    if isfield(s, 'signals')
        fitData = s.signals;
    elseif isfield(s, 'fitResults')
        fitData = s.fitResults;
    elseif isfield(s, 'fit_results')
        fitData = s.fit_results;
    elseif isfield(s, 'exportData') && isfield(s.exportData, 'signals')
        fitData = s.exportData.signals;
    else
        error('Could not find spot data in export MAT file: %s', filepath);
    end

    if istable(fitData)
        fitData = table2struct(fitData);
    end

    fitData = normalizeFitDataLocal(fitData);
end

function fitData = buildFitDataFromImage(imagePath, params, channelIdx, progressCb)
    if nargin < 4
        progressCb = [];
    end
    if nargin < 2 || isempty(params) || ~isstruct(params)
        error('Image-based training requires a valid ParameterStruct.');
    end

    if nargin < 3 || isempty(channelIdx)
        channelIdx = 1;
    end
    channelIdx = max(1, round(channelIdx));
    imageTimer = tic;

    emitProgress(progressCb, 'Image candidate generation started (ch=%d): %s', channelIdx, imagePath);
    rawImage = loadImageVolumeForTraining(imagePath);
    emitProgress(progressCb, 'Loaded image volume (ch=%d): size=%s', channelIdx, mat2str(size(rawImage)));
    handles = createTrainingHandlesFromParams(params, channelIdx);
    maximaEnabled = logical(getChannelParamValue(params, 'maximaEnabled', channelIdx, true));
    if ~maximaEnabled
        emitProgress(progressCb, 'Maxima detection disabled in parameters (ch=%d).', channelIdx);
        fitData = struct([]);
        return;
    end

    fitEnabled = logical(getChannelParamValue(params, 'gaussFitEnabled', channelIdx, true));
    if ~fitEnabled
        error('gaussFitEnabled is false for channel %d; SNAP_train requires fitted candidates for feature extraction.', channelIdx);
    end

    pipelineContext = struct();
    pipelineContext.mode = 'train';
    pipelineContext.channelIdx = channelIdx;
    pipelineContext.params = params;
    pipelineContext.handles = handles;
    pipelineContext.progressCallback = progressCb;
    pipelineContext.enableClassification = false;
    pipelineContext.enableNucleiFiltering = false;
    pipelineContext.fitParams = extractTrainingFitParams(params, channelIdx);
    pipelineContext.fitProgressInterval = 250;
    pipelineContext.fitAbortPollInterval = 5;

    pipelineResult = snap_modules.signal.runPipeline(rawImage, pipelineContext);
    if isempty(pipelineResult.processedImage)
        emitProgress(progressCb, 'Preprocessing returned empty image (ch=%d).', channelIdx);
        fitData = struct([]);
        return;
    end

    maximaCoords = pipelineResult.maximaCoords;
    if isempty(maximaCoords)
        emitProgress(progressCb, 'No maxima detected (ch=%d).', channelIdx);
        fitData = struct([]);
        return;
    end

    fitData = pipelineResult.fitResults;
    fitData = normalizeFitDataLocal(fitData);
    emitProgress(progressCb, 'Image candidate generation complete (ch=%d): %d candidate(s), %.2f s total.', ...
        channelIdx, numel(fitData), toc(imageTimer));
end

function img = loadImageVolumeForTraining(imagePath)
    try
        img = double(tiffreadVolume(imagePath));
    catch
        try
            img = double(imread(imagePath));
        catch ME
            error('Failed to load image file %s: %s', imagePath, ME.message);
        end
    end

    if ndims(img) < 3
        img = reshape(img, size(img, 1), size(img, 2), 1);
    end
end

function handles = createTrainingHandlesFromParams(params, requiredChannelIdx)
    inferredChannels = inferNumChannelsFromParameters(params, requiredChannelIdx);
    handles = struct();
    handles.Nmax = max(requiredChannelIdx, inferredChannels);

    for ch = 1:handles.Nmax
        handles.xySpacingInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'xySpacing', ch, 1));
        handles.zSpacingInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'zSpacing', ch, 1), 'Enable', 'on');

        handles.deconvEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'deconvEnabled', ch, false));
        handles.deconvMethodDrops(ch) = struct('Value', getChannelStringParamValue(params, 'deconvMethod', ch, 'Lucy-Richardson'));
        handles.deconvPSFSourceDrops(ch) = struct('Value', getChannelStringParamValue(params, 'deconvPSFSource', ch, 'Generate'));
        handles.deconvPSFPathTexts(ch) = struct('Value', getChannelStringParamValue(params, 'deconvPSFFilePath', ch, ''));
        handles.deconvPSFSigmaXYInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvPSFSigmaXY', ch, 1));
        handles.deconvPSFSigmaZInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvPSFSigmaZ', ch, 1));
        handles.deconvPSFSizeXYInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvPSFSizeXY', ch, 7));
        handles.deconvPSFSizeZInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvPSFSizeZ', ch, 5));
        handles.deconvLRIterationsInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvLRIterations', ch, 10));
        handles.deconvLRDampingInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvLRDamping', ch, 0));
        handles.deconvWienerNSRInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvWienerNSR', ch, 0.01));
        handles.deconvBlindIterationsInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvBlindIterations', ch, 10));
        handles.deconvBlindUnderRelaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'deconvBlindUnderRelax', ch, 0));

        handles.preprocEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'preprocEnabled', ch, false));
        handles.preprocessModeDrops(ch) = struct('Value', getChannelStringParamValue(params, 'preProcMode', ch, '3D'));
        handles.preprocessScaleChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'preProcScale', ch, false));
        handles.preprocessProjectionDrops(ch) = struct('Value', getChannelStringParamValue(params, 'preProcProjection', ch, 'Max'));
        handles.preprocMethodDrops(ch) = struct('Value', getChannelStringParamValue(params, 'preProcMethod', ch, 'Gaussian'));
        handles.gaussInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'smoothGaussianValues', ch, 1));
        handles.medianInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'smoothMedianValues', ch, 1));
        handles.waveletNameDrops(ch) = struct('Value', getChannelStringParamValue(params, 'waveletName', ch, 'db4'));
        handles.waveletLevelInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'waveletLevel', ch, 2));
        handles.waveletThresholdRuleDrops(ch) = struct('Value', getChannelStringParamValue(params, 'waveletThresholdRule', ch, 'sqtwolog'));
        handles.waveletThresholdMethodDrops(ch) = struct('Value', getChannelStringParamValue(params, 'waveletThresholdMethod', ch, 'soft'));
        handles.nlmFilterStrengthInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'nlmFilterStrength', ch, 10));
        handles.nlmSearchWindowInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'nlmSearchWindow', ch, 21));
        handles.nlmComparisonWindowInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'nlmComparisonWindow', ch, 5));
        handles.preprocClipChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'preprocClipAtZero', ch, false));

        handles.bgCorrEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'bgCorrEnabled', ch, false));
        handles.bgCorrModeDrops(ch) = struct('Value', getChannelStringParamValue(params, 'bgCorrMode', ch, '3D'));
        handles.bgCorrScaleChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'bgCorrScale', ch, false));
        handles.bgCorrProjectionDrops(ch) = struct('Value', getChannelStringParamValue(params, 'bgCorrProjection', ch, 'Max'));
        handles.bgMethodDrops(ch) = struct('Value', getChannelStringParamValue(params, 'bgMethod', ch, 'Gaussian'));
        handles.bgParamInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'bgParam', ch, 5));
        handles.bgCorrClipChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'bgCorrClipAtZero', ch, false));

        handles.maximaEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'maximaEnabled', ch, true));
        handles.maximaModeDrops(ch) = struct('Value', getChannelStringParamValue(params, 'maximaMode', ch, '3D'));
        handles.maximaScaleChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'maximaScale', ch, false));
        handles.maximaProjectionDrops(ch) = struct('Value', getChannelStringParamValue(params, 'maximaProjection', ch, 'Max'));
        handles.maximaMethodDrops(ch) = struct('Value', getChannelStringParamValue(params, 'maximaMethod', ch, 'Simple Regional'));
        handles.maximaNeighborhoodInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'maximaNeighborhoodSize', ch, 3));
        handles.hMaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'hMaxValue', ch, 10));
        handles.logSigmaInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'sigmaValue', ch, 1));
        handles.logThresholdInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'peakThresholdValue', ch, 0));

        handles.gaussFitEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'gaussFitEnabled', ch, true));
        handles.gaussFitVoxelWindowSlider(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitVoxelWindowSize', ch, 7));
        handles.gaussFitBgCorrMethodDrop(ch) = struct('Value', getChannelStringParamValue(params, 'gaussFitBgCorrMethod', ch, 'Mean Surrounding Subtraction'));
        handles.gaussFitBgCorrWidthEdit(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitBgCorrWidth', ch, 2));
        handles.gaussFitPolyDegreeEdit(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitPolyDegree', ch, 2));
        handles.gaussFitMethodDrop(ch) = struct('Value', getChannelStringParamValue(params, 'gaussFitMethod', ch, '1D (X,Y,Z) Gaussian'));
        handles.gaussFitMaxIterationsEdit(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitMaxIterations', ch, 200));
        handles.gaussFitToleranceEdit(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitTolerance', ch, 1e-6));
        handles.gaussFitRadialRadiusEdit(ch) = struct('Value', getChannelNumericParamValue(params, 'gaussFitRadialRadius', ch, 3));

        handles.fitFilterEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'fitFilterEnabled', ch, false));
        handles.fitFilterRSquaredEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'fitFilterRSquaredEnabled', ch, false));
        handles.fitFilterRSquaredMinInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterRSquaredMin', ch, 0));
        handles.fitFilterRSquaredMaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterRSquaredMax', ch, 1));
        handles.fitFilterSigmaSumEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'fitFilterSigmaSumEnabled', ch, false));
        handles.fitFilterSigmaSumMinInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterSigmaSumMin', ch, 0));
        handles.fitFilterSigmaSumMaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterSigmaSumMax', ch, inf));
        handles.fitFilterAmplitudeEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'fitFilterAmplitudeEnabled', ch, false));
        handles.fitFilterAmplitudeMinInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterAmplitudeMin', ch, -inf));
        handles.fitFilterAmplitudeMaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterAmplitudeMax', ch, inf));
        handles.fitFilterIntensityEnabledChecks(ch) = struct('Value', getChannelLogicalParamValue(params, 'fitFilterIntensityEnabled', ch, false));
        handles.fitFilterIntensityMinInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterIntensityMin', ch, -inf));
        handles.fitFilterIntensityMaxInputs(ch) = struct('Value', getChannelNumericParamValue(params, 'fitFilterIntensityMax', ch, inf));
    end
end

function fitParams = extractTrainingFitParams(params, channelIdx)
    fitParams = struct();
    fitParams.gaussFitVoxelWindowSize = getChannelNumericParamValue(params, 'gaussFitVoxelWindowSize', channelIdx, 7);
    fitParams.gaussFitBgCorrMethod = getChannelStringParamValue(params, 'gaussFitBgCorrMethod', channelIdx, 'Mean Surrounding Subtraction');
    fitParams.gaussFitBgCorrWidth = getChannelNumericParamValue(params, 'gaussFitBgCorrWidth', channelIdx, 2);
    fitParams.gaussFitPolyDegree = getChannelNumericParamValue(params, 'gaussFitPolyDegree', channelIdx, 2);
    fitParams.gaussFitMethod = getChannelStringParamValue(params, 'gaussFitMethod', channelIdx, '1D (X,Y,Z) Gaussian');
    fitParams.gaussFitMaxIterations = getChannelNumericParamValue(params, 'gaussFitMaxIterations', channelIdx, 200);
    fitParams.gaussFitTolerance = getChannelNumericParamValue(params, 'gaussFitTolerance', channelIdx, 1e-6);
    fitParams.gaussFitRadialRadius = getChannelNumericParamValue(params, 'gaussFitRadialRadius', channelIdx, 3);
    fitParams.gaussFitPlotCheck = false;
    fitParams.xySpacing = getChannelNumericParamValue(params, 'xySpacing', channelIdx, 1);
    fitParams.zSpacing = getChannelNumericParamValue(params, 'zSpacing', channelIdx, 1);
end

function v = getChannelParamValue(params, fieldName, ch, defaultValue)
    v = defaultValue;
    if ~isstruct(params) || ~isfield(params, fieldName)
        return;
    end
    raw = params.(fieldName);
    if isempty(raw)
        return;
    end

    if iscell(raw)
        if ch <= numel(raw) && ~isempty(raw{ch})
            v = raw{ch};
        elseif ~isempty(raw{1})
            v = raw{1};
        end
    elseif isstring(raw)
        if ch <= numel(raw)
            v = raw(ch);
        else
            v = raw(1);
        end
    elseif isnumeric(raw) || islogical(raw)
        if isscalar(raw)
            v = raw;
        elseif isvector(raw)
            if ch <= numel(raw)
                v = raw(ch);
            else
                v = raw(1);
            end
        else
            v = raw(1);
        end
    else
        v = raw;
    end
end

function v = getChannelNumericParamValue(params, fieldName, ch, defaultValue)
    raw = getChannelParamValue(params, fieldName, ch, defaultValue);
    if iscell(raw) && ~isempty(raw), raw = raw{1}; end
    if isstring(raw) || ischar(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            v = parsed;
            return;
        end
    end
    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
        v = double(raw);
    elseif islogical(raw) && isscalar(raw)
        v = double(raw);
    else
        v = defaultValue;
    end
end

function v = getChannelLogicalParamValue(params, fieldName, ch, defaultValue)
    raw = getChannelParamValue(params, fieldName, ch, defaultValue);
    if islogical(raw) && isscalar(raw)
        v = logical(raw);
    elseif isnumeric(raw) && isscalar(raw) && isfinite(raw)
        v = logical(raw ~= 0);
    elseif ischar(raw) || isstring(raw)
        token = lower(strtrim(char(string(raw))));
        v = ismember(token, {'1', 'true', 'yes', 'on'});
    else
        v = logical(defaultValue);
    end
end

function v = getChannelStringParamValue(params, fieldName, ch, defaultValue)
    raw = getChannelParamValue(params, fieldName, ch, defaultValue);
    if ischar(raw)
        v = raw;
    elseif isstring(raw)
        if isempty(raw)
            v = defaultValue;
        else
            v = char(raw(1));
        end
    elseif isnumeric(raw) || islogical(raw)
        v = char(string(raw));
    else
        v = defaultValue;
    end
    if isempty(v)
        v = defaultValue;
    end
end

function fitData = normalizeFitDataLocal(fitData)
    for i = 1:numel(fitData)
        if ~isfield(fitData(i), 'signal_id') || isempty(fitData(i).signal_id)
            fitData(i).signal_id = i;
        end

        if ~isfield(fitData(i), 'globalFitCenter') || isempty(fitData(i).globalFitCenter)
            if isfield(fitData(i), 'fitted_coords') && ~isempty(fitData(i).fitted_coords)
                center = fitData(i).fitted_coords(:)';
            elseif isfield(fitData(i), 'maxima_coords') && ~isempty(fitData(i).maxima_coords)
                center = fitData(i).maxima_coords(:)';
            else
                center = [nan nan nan];
            end
            fitData(i).globalFitCenter = forceXYZ(center);
        else
            fitData(i).globalFitCenter = forceXYZ(fitData(i).globalFitCenter);
        end
    end
end

function labels = loadLabelsForFitData(labelFile, fitData, matchDistance, convertFijiCoords)
    if nargin < 4
        convertFijiCoords = false;
    end

    [~, ~, ext] = fileparts(labelFile);
    ext = lower(ext);

    switch ext
        case '.mat'
            s = load(labelFile);
            if isfield(s, 'labeledReal') && isfield(s, 'labeledNoise')
                labels = nan(numel(fitData), 1);
                realIdx = s.labeledReal(:);
                noiseIdx = s.labeledNoise(:);
                realIdx = realIdx(realIdx >= 1 & realIdx <= numel(fitData));
                noiseIdx = noiseIdx(noiseIdx >= 1 & noiseIdx <= numel(fitData));
                labels(realIdx) = 1;
                labels(noiseIdx) = 0;
            else
                error('MAT label file must contain labeledReal and labeledNoise: %s', labelFile);
            end
        otherwise
            try
                t = readtable(labelFile, 'VariableNamingRule', 'preserve');
            catch
                t = readtable(labelFile);
            end
            labels = labelsFromTable(t, fitData, matchDistance, convertFijiCoords);
    end
end

function labels = labelsFromTable(t, fitData, matchDistance, convertFijiCoords)
    if nargin < 4
        convertFijiCoords = false;
    end

    rawNames = string(t.Properties.VariableNames);
    vn = canonicalizeVarNames(rawNames);
    labelIdx = find(vn == "label" | vn == "class" | vn == "isreal", 1);
    if isempty(labelIdx)
        y = ones(height(t), 1);
    else
        rawLabel = t.(t.Properties.VariableNames{labelIdx});
        y = normalizeLabelValues(rawLabel);
    end

    labels = nan(numel(fitData), 1);

    yName = find(vn == "maximay" | vn == "fittedy" | vn == "y" | vn == "row", 1);
    xName = find(vn == "maximax" | vn == "fittedx" | vn == "x" | vn == "col" | vn == "column", 1);
    zName = find(vn == "maximaz" | vn == "fittedz" | vn == "z" | vn == "slice", 1);

    if isempty(yName) || isempty(xName)
        error('Label table must include coordinate columns x and y (optional z/slice).');
    end

    yCoord = toNumericCoordColumn(t.(t.Properties.VariableNames{yName}), 'y');
    xCoord = toNumericCoordColumn(t.(t.Properties.VariableNames{xName}), 'x');
    if isempty(zName)
        zCoord = ones(height(t), 1);
    else
        zCoord = toNumericCoordColumn(t.(t.Properties.VariableNames{zName}), 'z');
    end

    if convertFijiCoords
        yCoord = yCoord + 1;
        xCoord = xCoord + 1;
        if ~isempty(zName)
            zCoord = zCoord + 1;
        end
    end

    centers = reshape([fitData.globalFitCenter], 3, [])';
    assignedCandidates = false(numel(fitData), 1);
    matchDistance2 = matchDistance.^2;

    for i = 1:height(t)
        if isnan(y(i))
            continue;
        end

        available = ~assignedCandidates;
        if ~any(available)
            break;
        end

        d2 = (centers(:,1) - yCoord(i)).^2 + (centers(:,2) - xCoord(i)).^2 + (centers(:,3) - zCoord(i)).^2;
        d2(~available) = inf;
        [best, idx] = min(d2);

        if ~isempty(idx) && isfinite(best) && best <= matchDistance2
            labels(idx) = y(i);
            assignedCandidates(idx) = true;
        end
    end
end

function values = toNumericCoordColumn(rawValues, coordName)
    if isnumeric(rawValues) || islogical(rawValues)
        values = double(rawValues(:));
        return;
    end

    if iscell(rawValues)
        rawValues = string(rawValues);
    else
        rawValues = string(rawValues(:));
    end

    values = str2double(rawValues);
    badMask = ~isfinite(values);
    if any(badMask)
        error('Coordinate column "%s" contains non-numeric values.', coordName);
    end
end

function xyz = forceXYZ(coords)
    xyz = coords(:)';
    if numel(xyz) >= 3
        xyz = xyz(1:3);
    elseif numel(xyz) == 2
        xyz = [xyz, 1];
    elseif numel(xyz) == 1
        xyz = [xyz, 1, 1];
    else
        xyz = [nan, nan, nan];
    end
end

function y = normalizeLabelValues(raw)
    if isnumeric(raw) || islogical(raw)
        y = double(raw(:));
        y(y ~= 0) = 1;
        return;
    end

    if iscell(raw)
        raw = string(raw);
    else
        raw = string(raw(:));
    end

    y = nan(numel(raw), 1);
    posSet = ["1", "true", "real", "spot", "positive", "pos"];
    negSet = ["0", "false", "noise", "negative", "bg", "background", "neg"];

    for i = 1:numel(raw)
        token = lower(strtrim(raw(i)));
        if any(token == posSet)
            y(i) = 1;
        elseif any(token == negSet)
            y(i) = 0;
        end
    end
end

function method = inferFittingMethod(fitData)
    method = '3D Gaussian';

    if isempty(fitData)
        return;
    end

    if isfield(fitData, 'fitMethod') && ~isempty(fitData(1).fitMethod)
        method = string(fitData(1).fitMethod);
        return;
    end

    names = fieldnames(fitData);
    if any(strcmp(names, 'radial_symmetry_score'))
        method = 'Radial Symmetry';
    elseif any(strcmp(names, 'amplitude_x'))
        method = '1D (X,Y,Z)';
    elseif any(strcmp(names, 'amplitude_xy'))
        method = '2D (XY) + 1D (Z)';
    elseif any(strcmp(names, 'rho_xy'))
        method = 'Distorted 3D Gaussian';
    end
end

function has3D = inferHas3D(fitData)
    has3D = false;
    if isempty(fitData)
        return;
    end

    if isfield(fitData, 'fitted_coords')
        c = fitData(1).fitted_coords;
        has3D = numel(c) >= 3;
    elseif isfield(fitData, 'maxima_coords')
        c = fitData(1).maxima_coords;
        has3D = numel(c) >= 3;
    elseif isfield(fitData, 'sigma_z')
        has3D = true;
    end
end

function selected = normalizeSelectedFeatures(selectedInput, featureInfo)
    if isempty(selectedInput)
        allFeatures = fieldnames(featureInfo);
        selected = {};
        for i = 1:numel(allFeatures)
            f = allFeatures{i};
            if isfield(featureInfo.(f), 'category') && strcmp(featureInfo.(f).category, 'position')
                continue;
            end
            selected{end+1} = f; %#ok<AGROW>
        end
    else
        if ischar(selectedInput)
            selected = {selectedInput};
        elseif isstring(selectedInput)
            selected = cellstr(selectedInput(:));
        else
            selected = selectedInput;
        end
    end
end

function customExpressions = normalizeCustomExpressionList(customExpressionsInput)
    customExpressions = struct('name', {}, 'expression', {});
    if isempty(customExpressionsInput)
        return;
    end
    if ~isstruct(customExpressionsInput)
        return;
    end

    ce = customExpressionsInput(:);
    keepMask = false(size(ce));
    for i = 1:numel(ce)
        if ~isfield(ce(i), 'name') || ~isfield(ce(i), 'expression')
            continue;
        end
        name = strtrim(char(string(ce(i).name)));
        expr = strtrim(char(string(ce(i).expression)));
        if isempty(name) || isempty(expr)
            continue;
        end
        ce(i).name = name;
        ce(i).expression = expr;
        keepMask(i) = true;
    end
    ce = ce(keepMask);
    if isempty(ce)
        return;
    end

    seenNames = strings(0, 1);
    out = struct('name', {}, 'expression', {});
    for i = 1:numel(ce)
        exprName = string(ce(i).name);
        if any(exprName == seenNames)
            continue;
        end
        out(end+1) = struct('name', ce(i).name, 'expression', ce(i).expression); %#ok<AGROW>
        seenNames(end+1, 1) = exprName; %#ok<AGROW>
    end
    customExpressions = out;
end

function [X, featureNames, validMask, extractionInfo, selectedFeatures, customExpressions] = ...
        buildTrainingFeatureMatrixWithFallback(fitData, selectedFeatures, featureInfo, customExpressions, progressCb)
    maxPasses = max(6, numel(selectedFeatures) + numel(customExpressions) + 2);
    selectedFeatures = selectedFeatures(:)';
    customExpressions = normalizeCustomExpressionList(customExpressions);

    X = [];
    featureNames = {};
    validMask = [];
    extractionInfo = struct();

    for pass = 1:maxPasses
        [X, featureNames, validMask, extractionInfo] = snap_helpers.classification.buildFeatureMatrix( ...
            fitData, selectedFeatures, featureInfo, customExpressions);
        emitProgress(progressCb, 'Training feature matrix built: %d candidate rows, %d feature columns.', ...
            size(X, 1), size(X, 2));
        if isstruct(extractionInfo) && isfield(extractionInfo, 'modelStats') && ...
                isstruct(extractionInfo.modelStats) && ...
                isfield(extractionInfo.modelStats, 'requested') && extractionInfo.modelStats.requested
            ms = extractionInfo.modelStats;
            if isfield(ms, 'augmented') && ms.augmented && isfield(ms, 'summary') && isstruct(ms.summary)
                emitProgress(progressCb, ...
                    ['Model-stat augmentation: computed %d/%d windows ', ...
                     '(missingWindow=%d, modelFailures=%d).'], ...
                    getfieldwithdefault(ms.summary, 'nComputed', 0), ...
                    getfieldwithdefault(ms.summary, 'nTotal', 0), ...
                    getfieldwithdefault(ms.summary, 'nMissingWindow', 0), ...
                    getfieldwithdefault(ms.summary, 'nModelFailures', 0));
            else
                emitProgress(progressCb, 'Model-stat features requested: using existing fields (no augmentation run).');
            end
        end

        allNaNFeatures = {};
        if isstruct(extractionInfo) && isfield(extractionInfo, 'featuresAllNaN') && ~isempty(extractionInfo.featuresAllNaN)
            allNaNFeatures = extractionInfo.featuresAllNaN;
            emitProgress(progressCb, 'Training feature extraction all-NaN feature(s): %s', ...
                strjoin(allNaNFeatures, ', '));
        end

        if isempty(allNaNFeatures) && any(validMask)
            return;
        end

        if ~isempty(allNaNFeatures)
            dropNames = allNaNFeatures;
            dropReason = 'all-NaN';
        else
            [dropNames, dropReason] = chooseFallbackDropFeatures(extractionInfo, featureNames, size(X, 1));
        end

        if isempty(dropNames)
            return;
        end

        [nextSelected, nextCustom, droppedBase, droppedCustom] = ...
            dropSelectedFeaturesByName(selectedFeatures, customExpressions, dropNames);

        if isempty(droppedBase) && isempty(droppedCustom)
            return;
        end

        if ~isempty(droppedBase)
            emitProgress(progressCb, 'Dropping %s base feature(s): %s', dropReason, strjoin(droppedBase, ', '));
        end
        if ~isempty(droppedCustom)
            emitProgress(progressCb, 'Dropping %s custom expression(s): %s', dropReason, strjoin(droppedCustom, ', '));
        end

        selectedFeatures = nextSelected;
        customExpressions = nextCustom;

        if isempty(selectedFeatures) && isempty(customExpressions)
            emitProgress(progressCb, 'All selected features were incompatible. Falling back to AUTO base features.');
            selectedFeatures = normalizeSelectedFeatures({}, featureInfo);
        end

        emitProgress(progressCb, 'Rebuilding training feature matrix after compatibility pruning (pass %d/%d)...', ...
            min(pass + 1, maxPasses), maxPasses);
    end
end

function [selectedOut, customOut, droppedBase, droppedCustom] = ...
        dropSelectedFeaturesByName(selectedIn, customIn, featureNamesToDrop)
    selectedOut = selectedIn(:)';
    customOut = normalizeCustomExpressionList(customIn);
    droppedBase = {};
    droppedCustom = {};

    if isempty(featureNamesToDrop)
        return;
    end

    badSet = string(featureNamesToDrop(:)');

    if ~isempty(selectedOut)
        selStr = string(selectedOut);
        dropMask = ismember(selStr, badSet);
        droppedBase = cellstr(selStr(dropMask));
        selectedOut = selectedOut(~dropMask);
    end

    if isempty(customOut)
        return;
    end

    customNames = strings(1, numel(customOut));
    for i = 1:numel(customOut)
        customNames(i) = string(customOut(i).name);
    end
    dropMask = ismember(customNames, badSet);
    if any(dropMask)
        droppedCustom = cellstr(customNames(dropMask));
        customOut = customOut(~dropMask);
    end
end

function [dropNames, reason] = chooseFallbackDropFeatures(extractionInfo, featureNames, nRows)
    dropNames = {};
    reason = 'incompatible';

    if nargin < 3 || ~isfinite(nRows) || nRows <= 0
        nRows = 0;
    end
    if ~isstruct(extractionInfo) || ~isfield(extractionInfo, 'nanCountByFeature')
        return;
    end

    nanCounts = extractionInfo.nanCountByFeature;
    if isempty(nanCounts)
        return;
    end
    nanCounts = nanCounts(:)';

    names = featureNames;
    if isempty(names) && isfield(extractionInfo, 'featureNames')
        names = extractionInfo.featureNames;
    end
    if isempty(names) || numel(names) ~= numel(nanCounts)
        return;
    end

    hasNaNMask = nanCounts > 0;
    if ~any(hasNaNMask)
        return;
    end

    idxCandidates = find(hasNaNMask);
    [maxNaN, relIdx] = max(nanCounts(idxCandidates));
    worstIdx = idxCandidates(relIdx);
    dropNames = {names{worstIdx}};

    if nRows > 0
        reason = sprintf('highest-NaN (%d/%d rows)', maxNaN, nRows);
    else
        reason = sprintf('highest-NaN (%d rows)', maxNaN);
    end
end

function options = applyTrainingDefaults(options)
    if ~isfield(options, 'kernelFunction'), options.kernelFunction = 'rbf'; end
    if ~isfield(options, 'boxConstraint'), options.boxConstraint = 1.0; end
    if ~isfield(options, 'kernelScale'), options.kernelScale = 'auto'; end
    if ~isfield(options, 'polynomialOrder'), options.polynomialOrder = 3; end
    if ~isfield(options, 'classWeightMode'), options.classWeightMode = 'balanced'; end
    if ~isfield(options, 'standardize'), options.standardize = true; end
    if ~isfield(options, 'crossValidate'), options.crossValidate = true; end
    if ~isfield(options, 'kFold'), options.kFold = 5; end
    if ~isfield(options, 'verbose'), options.verbose = true; end
end
