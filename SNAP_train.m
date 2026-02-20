function varargout = SNAP_train(exportFiles, labelFiles, outputClassifierPath, varargin)
% SNAP_train - Train and validate a SNAP-compatible SVM from labeled spot files
%
% UI MODE:
%   SNAP_train
%   Opens a SNAP-style training GUI for channel-aware multi-SVM training.
%   The GUI:
%     - Loads a SNAP parameter file to infer channel count and per-channel defaults
%     - Uses exactly one SVM slot per detected channel (toggle channels on/off)
%     - Sets match distance per channel (per SVM)
%     - Lets you assign training/validation directories per channel
%     - Lets you browse a shared output directory
%     - Lets you select base features + custom expressions (shared engine with SNAP_classify)
%     - Supports manual SVM hyperparameters OR validation-based optimization
%     - Optionally emits a sweep performance report (log + plot)
%     - Trains one classifier per selected channel
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

    [XTrain, featureNames, validTrainMask, extractionInfoTrain] = snap_helpers.classification.buildFeatureMatrix( ...
        trainFitData, selectedFeatures, featureInfo, opts.CustomExpressions);
    emitProgress(progressCb, 'Training feature matrix built: %d candidate rows, %d feature columns.', size(XTrain,1), size(XTrain,2));

    XTrain = XTrain(validTrainMask, :);
    yTrain = trainLabels(validTrainMask);
    trainSources = trainSources(validTrainMask);
    fixedTrainLabels = yTrain; % immutable reference labels used for all sweep candidates
    emitProgress(progressCb, 'Training dataset after filtering: %d samples (real=%d, noise=%d).', ...
        numel(yTrain), sum(yTrain == 1), sum(yTrain == 0));

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
            valFitData, selectedFeatures, featureInfo, opts.CustomExpressions);
        emitProgress(progressCb, 'Validation feature matrix built: %d candidate rows, %d feature columns.', size(XVal,1), size(XVal,2));
        XVal = XVal(validValMask, :);
        yVal = valLabelsRaw(validValMask);
        emitProgress(progressCb, 'Validation dataset after filtering: %d samples (real=%d, noise=%d).', ...
            numel(yVal), sum(yVal == 1), sum(yVal == 0));

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

        [validation.metrics, validation.confusionMatrix] = evaluateOnValidation(model, normParams, XVal, yVal);
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
    metadata.trainingLabelsFixedDuringValidation = true;
    if hasValidationSet
        metadata.validation = validation;
    end

    success = snap_helpers.classification.saveClassifier( ...
        outputClassifierPath, model, selectedFeatures, featureInfo, trainStats, ...
        char(fittingMethod), metadata, opts.CustomExpressions, normParams);

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
    out.customExpressions = opts.CustomExpressions;
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
        fprintf('  Features: %d\n', numel(selectedFeatures) + numel(opts.CustomExpressions));
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
        'Position', [70 70 1460 860]);

    mainGrid = uigridlayout(fig, [1, 2]);
    mainGrid.ColumnWidth = {520, '1x'};
    mainGrid.Padding = [10 10 10 10];
    mainGrid.ColumnSpacing = 10;

    % ---------------------------------------------------------------------
    % Left panel: configuration
    % ---------------------------------------------------------------------
    cfgPanel = uipanel(mainGrid, 'Title', 'Training Configuration', 'FontWeight', 'bold');
    cfgPanel.Layout.Column = 1;
    cfgGrid = uigridlayout(cfgPanel, [6, 1]);
    cfgGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    cfgGrid.RowSpacing = 8;
    cfgGrid.Padding = [8 8 8 8];

    % Parameter file
    paramPanel = uipanel(cfgGrid, 'Title', 'Parameter File');
    paramGrid = uigridlayout(paramPanel, [3, 3]);
    paramGrid.RowHeight = {'fit', 'fit', 'fit'};
    paramGrid.ColumnWidth = {'fit', '1x', 'fit'};
    paramGrid.RowSpacing = 6;
    paramGrid.Padding = [0 0 0 0];
    uilabel(paramGrid, 'Text', 'SNAP Parameter File:', 'FontWeight', 'bold');
    uilabel(paramGrid, 'Text', '');
    uilabel(paramGrid, 'Text', '');
    paramPathEdit = uieditfield(paramGrid, 'text', 'Editable', 'off', 'Placeholder', 'Select parameter file (.mat)');
    paramPathEdit.Layout.Row = 2;
    paramPathEdit.Layout.Column = [1 2];
    paramBrowseBtn = uibutton(paramGrid, 'Text', 'Browse');
    paramBrowseBtn.Layout.Row = 2;
    paramBrowseBtn.Layout.Column = 3;
    detectedChannelsTextLabel = uilabel(paramGrid, 'Text', 'Detected Active Channels:', 'FontWeight', 'bold');
    detectedChannelsTextLabel.Layout.Row = 3;
    detectedChannelsTextLabel.Layout.Column = 1;
    detectedChannelsLabel = uilabel(paramGrid, 'Text', '1');
    detectedChannelsLabel.Layout.Row = 3;
    detectedChannelsLabel.Layout.Column = 2;

    % Output directory
    outPanel = uipanel(cfgGrid, 'Title', 'Output');
    outGrid = uigridlayout(outPanel, [1, 3]);
    outGrid.ColumnWidth = {'fit', '1x', 'fit'};
    outGrid.Padding = [6 6 6 6];
    outGrid.RowSpacing = 6;
    uilabel(outGrid, 'Text', 'Output Directory:');
    outDirEdit = uieditfield(outGrid, 'text', 'Editable', 'off');
    outDirEdit.Layout.Column = 2;
    outBrowseBtn = uibutton(outGrid, 'Text', 'Browse');
    outBrowseBtn.Layout.Column = 3;

    % Core options
    corePanel = uipanel(cfgGrid, 'Title', 'Core Options');
    coreGrid = uigridlayout(corePanel, [2, 2]);
    coreGrid.ColumnWidth = {'fit', '1x'};
    coreGrid.RowHeight = {'fit', 'fit'};
    coreGrid.Padding = [6 6 6 6];
    optimizeCheck = uicheckbox(coreGrid, 'Text', 'Optimize with validation sweep', 'Value', true);
    optimizeCheck.Layout.Column = [1 2];
    sweepReportCheck = uicheckbox(coreGrid, ...
        'Text', 'Include sweep performance report (log + plot)', ...
        'Value', true);
    sweepReportCheck.Layout.Row = 2;
    sweepReportCheck.Layout.Column = [1 2];

    % Feature / expression selection (shared engine with SNAP_classify)
    featurePanel = uipanel(cfgGrid, 'Title', 'Features & Custom Expressions');
    featureGrid = uigridlayout(featurePanel, [3, 1]);
    featureGrid.RowHeight = {'fit', 'fit', 90};
    featureGrid.Padding = [6 6 6 6];
    featureGrid.RowSpacing = 4;
    selectFeaturesBtn = uibutton(featureGrid, 'Text', 'Select Features...', 'Enable', 'off');
    featureCountLabel = uilabel(featureGrid, 'Text', 'AUTO: all non-position features (per channel)');
    featureCountLabel.FontColor = [0.45 0.45 0.45];
    featureListArea = uitextarea(featureGrid, 'Editable', 'off', ...
        'Value', {'AUTO: all non-position features', 'Custom expressions: none'});

    % Manual SVM options
    svmPanel = uipanel(cfgGrid, 'Title', 'Manual SVM Hyperparameters');
    svmGrid = uigridlayout(svmPanel, [4, 4]);
    svmGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
    svmGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    svmGrid.Padding = [6 6 6 6];
    uilabel(svmGrid, 'Text', 'Kernel:');
    kernelDrop = uidropdown(svmGrid, 'Items', {'rbf', 'linear', 'polynomial'}, 'Value', 'rbf');
    kernelDrop.Layout.Column = 2;
    uilabel(svmGrid, 'Text', 'Box Constraint (C):');
    boxConstraintInput = uieditfield(svmGrid, 'numeric', 'Value', 1, 'LowerLimit', 1e-6);
    boxConstraintInput.Layout.Column = 4;
    uilabel(svmGrid, 'Text', 'Kernel Scale:');
    kernelScaleInput = uieditfield(svmGrid, 'text', 'Value', 'auto');
    kernelScaleInput.Layout.Column = 2;
    uilabel(svmGrid, 'Text', 'Polynomial Order:');
    polyOrderInput = uieditfield(svmGrid, 'numeric', 'Value', 3, 'LowerLimit', 2, 'RoundFractionalValues', true);
    polyOrderInput.Layout.Column = 4;
    standardizeCheck = uicheckbox(svmGrid, 'Text', 'Z-score normalize', 'Value', true);
    standardizeCheck.Layout.Row = 3;
    standardizeCheck.Layout.Column = [1 2];
    crossValidateCheck = uicheckbox(svmGrid, 'Text', 'Cross-validate', 'Value', true);
    crossValidateCheck.Layout.Row = 3;
    crossValidateCheck.Layout.Column = [3 4];
    uilabel(svmGrid, 'Text', 'CV folds (k):');
    kFoldInput = uieditfield(svmGrid, 'numeric', 'Value', 5, 'LowerLimit', 2, 'RoundFractionalValues', true);
    kFoldInput.Layout.Row = 4;
    kFoldInput.Layout.Column = 2;
    balanceClassCheck = uicheckbox(svmGrid, 'Text', 'Balance class bins (Real/Noise)', 'Value', true);
    balanceClassCheck.Layout.Row = 4;
    balanceClassCheck.Layout.Column = [3 4];

    % Sweep options
    sweepPanel = uipanel(cfgGrid, 'Title', 'Validation Sweep Grid');
    sweepGrid = uigridlayout(sweepPanel, [4, 2]);
    sweepGrid.ColumnWidth = {'fit', '1x'};
    sweepGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    sweepGrid.Padding = [6 6 6 6];
    uilabel(sweepGrid, 'Text', 'Kernels:');
    sweepKernelsInput = uieditfield(sweepGrid, 'text', 'Value', 'linear, rbf, polynomial');
    uilabel(sweepGrid, 'Text', 'Box Constraints:');
    sweepBoxInput = uieditfield(sweepGrid, 'text', 'Value', '0.1, 1, 10');
    uilabel(sweepGrid, 'Text', 'Kernel Scales:');
    sweepScaleInput = uieditfield(sweepGrid, 'text', 'Value', 'auto, 0.5, 1, 2');
    uilabel(sweepGrid, 'Text', 'Polynomial Orders:');
    sweepPolyInput = uieditfield(sweepGrid, 'text', 'Value', '2, 3, 4');

    % ---------------------------------------------------------------------
    % Right panel: channel selection and logs
    % ---------------------------------------------------------------------
    rightPanel = uipanel(mainGrid, 'Title', 'Per-Channel Training', 'FontWeight', 'bold');
    rightPanel.Layout.Column = 2;
    rightGrid = uigridlayout(rightPanel, [5, 1]);
    rightGrid.RowHeight = {'fit', 'fit', 240, '1x', 'fit'};
    rightGrid.RowSpacing = 8;
    rightGrid.Padding = [8 8 8 8];

    perChannelDirPanel = uipanel(rightGrid, 'BorderType', 'none');
    perChannelDirGrid = uigridlayout(perChannelDirPanel, [1, 3]);
    perChannelDirGrid.ColumnWidth = {'fit', 'fit', '1x'};
    perChannelDirGrid.Padding = [0 0 0 0];
    browseTrainSelectedBtn = uibutton(perChannelDirGrid, 'Text', 'Browse Train Dir (Selected Row)');
    browseValSelectedBtn = uibutton(perChannelDirGrid, 'Text', 'Browse Validation Dir (Selected Row)');
    uilabel(perChannelDirGrid, 'Text', 'Select a row first, then browse directories for that channel.');

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

    channelTable = uitable(rightGrid);
    channelTable.ColumnName = {'Train', 'Channel', 'Match Distance (voxels)', 'Convert FIJI Coords', 'Training Directory', 'Validation Directory', 'Output Classifier'};
    channelTable.ColumnEditable = [true, false, true, true, false, false, true];
    channelTable.ColumnWidth = {55, 80, 145, 135, '1x', '1x', 200};
    channelTable.Data = {true, 'Channel 1', 2, false, '', '', 'classifier_ch1.mat'};

    logText = uitextarea(rightGrid, 'Editable', 'off');
    logText.Value = {'SNAP_train ready. Load parameter file, configure per-channel train/validation directories, then click Train Selected Channels.'};

    actionPanel = uipanel(rightGrid, 'BorderType', 'none');
    actionGrid = uigridlayout(actionPanel, [1, 2]);
    actionGrid.ColumnWidth = {'1x', 'fit'};
    actionGrid.Padding = [0 0 0 0];
    trainBtn = uibutton(actionGrid, 'Text', 'Train Selected Channels', 'FontWeight', 'bold');
    closeBtn = uibutton(actionGrid, 'Text', 'Close', 'ButtonPushedFcn', @(~,~) close(fig));

    % Store handles/state
    handles = struct();
    handles.fig = fig;
    handles.paramPathEdit = paramPathEdit;
    handles.detectedChannelsLabel = detectedChannelsLabel;
    handles.detectedNumChannels = 1;
    handles.outDirEdit = outDirEdit;
    handles.optimizeCheck = optimizeCheck;
    handles.sweepReportCheck = sweepReportCheck;
    handles.kernelDrop = kernelDrop;
    handles.boxConstraintInput = boxConstraintInput;
    handles.kernelScaleInput = kernelScaleInput;
    handles.polyOrderInput = polyOrderInput;
    handles.standardizeCheck = standardizeCheck;
    handles.crossValidateCheck = crossValidateCheck;
    handles.kFoldInput = kFoldInput;
    handles.balanceClassCheck = balanceClassCheck;
    handles.sweepKernelsInput = sweepKernelsInput;
    handles.sweepBoxInput = sweepBoxInput;
    handles.sweepScaleInput = sweepScaleInput;
    handles.sweepPolyInput = sweepPolyInput;
    handles.selectFeaturesBtn = selectFeaturesBtn;
    handles.featureCountLabel = featureCountLabel;
    handles.featureListArea = featureListArea;
    handles.channelFeatureConfigs = repmat(defaultChannelFeatureConfig(), 1, 1);
    handles.channelTable = channelTable;
    handles.selectedChannelRow = 1;
    handles.progressStatusLabel = progressStatusLabel;
    handles.progressBarGrid = progressBarGrid;
    handles.logText = logText;
    handles.loadedParams = struct();
    guidata(fig, handles);

    % Callbacks
    paramBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseParameterFile(fig);
    outBrowseBtn.ButtonPushedFcn = @(~,~) onBrowseDirectory(fig, 'outDirEdit', 'Select output directory');
    browseTrainSelectedBtn.ButtonPushedFcn = @(~,~) onBrowseSelectedChannelDirectory(fig, 'train');
    browseValSelectedBtn.ButtonPushedFcn = @(~,~) onBrowseSelectedChannelDirectory(fig, 'validation');
    if isprop(channelTable, 'SelectionChangedFcn')
        channelTable.SelectionChangedFcn = @(~,evt) onChannelTableSelectionChanged(fig, evt);
    elseif isprop(channelTable, 'CellSelectionCallback')
        channelTable.CellSelectionCallback = @(~,evt) onChannelTableSelectionChanged(fig, evt);
    end
    selectFeaturesBtn.ButtonPushedFcn = @(~,~) onSelectTrainingFeatures(fig);
    optimizeCheck.ValueChangedFcn = @(~,~) updateOptimizationControlState(fig);
    trainBtn.ButtonPushedFcn = @(~,~) onTrainSelectedChannels(fig);

    updateChannelTableFromDetectedChannels(fig, 1, true);
    refreshTrainingFeatureSummary(fig);
    updateOptimizationControlState(fig);
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
    if isfield(handles, 'selectFeaturesBtn') && isgraphics(handles.selectFeaturesBtn)
        handles.selectFeaturesBtn.Enable = 'on';
    end
    guidata(fig, handles);
    updateChannelTableFromDetectedChannels(fig, numChannels, true);
    refreshTrainingFeatureSummary(fig, 1);
    appendTrainLog(fig, sprintf('Loaded parameter file: %s (numChannels=%d)', file, numChannels));
end

function onBrowseDirectory(fig, fieldName, dialogTitle)
    handles = guidata(fig);
    startPath = pwd;
    currentValue = strtrim(char(string(handles.(fieldName).Value)));
    if ~isempty(currentValue) && isfolder(currentValue)
        startPath = currentValue;
    end
    selectedPath = uigetdir(startPath, dialogTitle);
    if isequal(selectedPath, 0)
        return;
    end
    handles.(fieldName).Value = selectedPath;
    guidata(fig, handles);
end

function onChannelTableSelectionChanged(fig, eventData)
    handles = guidata(fig);
    row = [];

    % Different MATLAB versions expose different selection payloads.
    if nargin >= 2 && ~isempty(eventData)
        row = extractSelectedRowFromTableEvent(eventData);
    end

    if isempty(row) && isprop(handles.channelTable, 'Selection')
        row = extractSelectedRowFromSelectionValue(handles.channelTable.Selection);
    end

    if ~isempty(row) && isfinite(row) && row >= 1
        handles.selectedChannelRow = round(row);
        guidata(fig, handles);
        try
            refreshTrainingFeatureSummary(fig, parseChannelFromRowValue(handles.channelTable.Data{handles.selectedChannelRow, 2}));
        catch
            refreshTrainingFeatureSummary(fig);
        end
    end
end

function row = extractSelectedRowFromTableEvent(eventData)
    row = [];

    if isprop(eventData, 'Indices')
        row = extractSelectedRowFromSelectionValue(eventData.Indices);
    elseif isprop(eventData, 'Selection')
        row = extractSelectedRowFromSelectionValue(eventData.Selection);
    elseif isprop(eventData, 'SelectedCells')
        row = extractSelectedRowFromSelectionValue(eventData.SelectedCells);
    elseif isprop(eventData, 'CurrentSelection')
        row = extractSelectedRowFromSelectionValue(eventData.CurrentSelection);
    elseif isstruct(eventData)
        if isfield(eventData, 'Indices')
            row = extractSelectedRowFromSelectionValue(eventData.Indices);
        elseif isfield(eventData, 'Selection')
            row = extractSelectedRowFromSelectionValue(eventData.Selection);
        elseif isfield(eventData, 'SelectedCells')
            row = extractSelectedRowFromSelectionValue(eventData.SelectedCells);
        elseif isfield(eventData, 'CurrentSelection')
            row = extractSelectedRowFromSelectionValue(eventData.CurrentSelection);
        end
    end
end

function row = extractSelectedRowFromSelectionValue(sel)
    row = [];
    if isempty(sel)
        return;
    end

    if iscell(sel)
        if ~isempty(sel)
            row = extractSelectedRowFromSelectionValue(sel{1});
        end
        return;
    end

    if isnumeric(sel)
        if size(sel, 2) >= 1
            row = sel(1, 1);
        elseif isvector(sel)
            row = sel(1);
        end
        return;
    end

    if isstruct(sel)
        if isfield(sel, 'Indices')
            row = extractSelectedRowFromSelectionValue(sel.Indices);
        elseif isfield(sel, 'Selection')
            row = extractSelectedRowFromSelectionValue(sel.Selection);
        end
        return;
    end

    if isobject(sel)
        if isprop(sel, 'Indices')
            row = extractSelectedRowFromSelectionValue(sel.Indices);
        elseif isprop(sel, 'Selection')
            row = extractSelectedRowFromSelectionValue(sel.Selection);
        end
    end
end

function onBrowseSelectedChannelDirectory(fig, dirType)
    handles = guidata(fig);
    data = handles.channelTable.Data;
    if isempty(data)
        uialert(fig, 'No channels available. Load a valid parameter file first.', 'No Channels');
        return;
    end

    row = handles.selectedChannelRow;
    row = max(1, min(size(data, 1), row));
    channelIdx = parseChannelFromRowValue(data{row, 2});

    switch lower(dirType)
        case 'train'
            colIdx = 5;
            dialogTitle = sprintf('Select training directory for Channel %d', channelIdx);
        case 'validation'
            colIdx = 6;
            dialogTitle = sprintf('Select validation directory for Channel %d', channelIdx);
        otherwise
            error('Unknown channel directory type: %s', dirType);
    end

    startPath = pwd;
    existingValue = strtrim(char(string(data{row, colIdx})));
    if ~isempty(existingValue) && isfolder(existingValue)
        startPath = existingValue;
    end

    selectedPath = uigetdir(startPath, dialogTitle);
    if isequal(selectedPath, 0)
        return;
    end

    data{row, colIdx} = selectedPath;
    handles.channelTable.Data = data;
    guidata(fig, handles);
end

function updateChannelTableFromDetectedChannels(fig, detectedChannels, resetOutputs)
    handles = guidata(fig);
    n = max(1, round(detectedChannels));
    handles.detectedNumChannels = n;
    handles.detectedChannelsLabel.Text = num2str(n);

    oldFeatureCfgs = repmat(defaultChannelFeatureConfig(), 1, 0);
    if isfield(handles, 'channelFeatureConfigs') && ~isempty(handles.channelFeatureConfigs)
        oldFeatureCfgs = handles.channelFeatureConfigs;
    end

    oldData = handles.channelTable.Data;
    newData = cell(n, 7);

    for ch = 1:n
        channelName = sprintf('Channel %d', ch);
        matchDistance = 2;
        convertFiji = false;
        trainDir = '';
        valDir = '';
        defaultOutput = sprintf('classifier_ch%d.mat', ch);
        trainFlag = true;
        outputName = defaultOutput;

        if ~isempty(oldData)
            rowIdx = find(strcmp(oldData(:, 2), channelName), 1);
            if ~isempty(rowIdx)
                trainFlag = isSelectedChannelFlag(oldData{rowIdx, 1});
                if size(oldData, 2) >= 7
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

        newData{ch, 1} = trainFlag;
        newData{ch, 2} = channelName;
        newData{ch, 3} = matchDistance;
        newData{ch, 4} = convertFiji;
        newData{ch, 5} = trainDir;
        newData{ch, 6} = valDir;
        newData{ch, 7} = outputName;
    end

    newFeatureCfgs = repmat(defaultChannelFeatureConfig(), 1, n);
    for ch = 1:n
        if ch <= numel(oldFeatureCfgs)
            newFeatureCfgs(ch) = normalizeChannelFeatureConfig(oldFeatureCfgs(ch));
        end
    end

    handles.channelFeatureConfigs = newFeatureCfgs;
    handles.channelTable.Data = newData;
    handles.selectedChannelRow = max(1, min(n, handles.selectedChannelRow));
    guidata(fig, handles);
end

function updateOptimizationControlState(fig)
    handles = guidata(fig);
    useOptimization = handles.optimizeCheck.Value;

    manualControls = { ...
        handles.kernelDrop, handles.boxConstraintInput, handles.kernelScaleInput, ...
        handles.polyOrderInput, handles.standardizeCheck, handles.crossValidateCheck, ...
        handles.kFoldInput, handles.balanceClassCheck ...
    };
    sweepControls = { ...
        handles.sweepKernelsInput, handles.sweepBoxInput, ...
        handles.sweepScaleInput, handles.sweepPolyInput, handles.sweepReportCheck ...
    };

    if useOptimization
        manualState = 'off';
        sweepState = 'on';
    else
        manualState = 'on';
        sweepState = 'off';
    end

    for i = 1:numel(manualControls)
        manualControls{i}.Enable = manualState;
    end
    for i = 1:numel(sweepControls)
        sweepControls{i}.Enable = sweepState;
    end

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
            ce = rawCfg.customExpressions;
            validMask = isfield(ce, 'name') & isfield(ce, 'expression');
            ce = ce(validMask);
            for k = 1:numel(ce)
                ce(k).name = char(string(ce(k).name));
                ce(k).expression = char(string(ce(k).expression));
            end
            cfg.customExpressions = ce;
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
    if isempty(handles.loadedParams) || ~isstruct(handles.loadedParams) || isempty(fieldnames(handles.loadedParams))
        uialert(fig, 'Load a parameter file first so feature availability can be inferred.', 'Missing Parameter File');
        return;
    end

    [contextChannel, fittingMethod, has3D] = getTrainingFeatureSelectionContext(handles);
    chanCfg = getChannelFeatureConfig(handles, contextChannel);
    [selected, customExpr, cancelled] = snap_helpers.classification.featureSelectionUI( ...
        fittingMethod, has3D, false, chanCfg.selectedFeatures, chanCfg.customExpressions);

    if cancelled
        return;
    end

    chanCfg.selectedFeatures = selected;
    chanCfg.customExpressions = customExpr;
    handles = setChannelFeatureConfig(handles, contextChannel, chanCfg);
    guidata(fig, handles);

    refreshTrainingFeatureSummary(fig, contextChannel);
end

function refreshTrainingFeatureSummary(fig, contextChannel)
    handles = guidata(fig);
    if ~isfield(handles, 'featureCountLabel') || ~isgraphics(handles.featureCountLabel) || ...
            ~isfield(handles, 'featureListArea') || ~isgraphics(handles.featureListArea)
        return;
    end

    if nargin < 2 || isempty(contextChannel)
        contextChannel = 1;
        if isfield(handles, 'channelTable') && isgraphics(handles.channelTable) && ~isempty(handles.channelTable.Data)
            row = max(1, min(size(handles.channelTable.Data, 1), handles.selectedChannelRow));
            contextChannel = parseChannelFromRowValue(handles.channelTable.Data{row, 2});
        end
    end

    chanCfg = getChannelFeatureConfig(handles, contextChannel);
    nBase = numel(chanCfg.selectedFeatures);
    nCustom = numel(chanCfg.customExpressions);

    if nBase == 0 && nCustom == 0
        handles.featureCountLabel.Text = sprintf('Channel %d: AUTO (all non-position features)', contextChannel);
        handles.featureCountLabel.FontColor = [0.45 0.45 0.45];
        handles.featureListArea.Value = { ...
            sprintf('Context channel: %d', contextChannel), ...
            'AUTO: all non-position features', ...
            'Custom expressions: none' ...
        };
        guidata(fig, handles);
        return;
    end

    nTotal = nBase + nCustom;
    handles.featureCountLabel.Text = sprintf('Channel %d: %d base + %d custom = %d total', ...
        contextChannel, nBase, nCustom, nTotal);
    handles.featureCountLabel.FontColor = [0.2 0.6 0.2];

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
    handles.featureListArea.Value = lines;
    guidata(fig, handles);
end

function [contextChannel, fittingMethod, has3D] = getTrainingFeatureSelectionContext(handles)
    contextChannel = 1;
    fittingMethod = '3D Gaussian';
    has3D = true;

    if isfield(handles, 'channelTable') && isgraphics(handles.channelTable) && ...
            ~isempty(handles.channelTable.Data)
        row = max(1, min(size(handles.channelTable.Data, 1), handles.selectedChannelRow));
        contextChannel = parseChannelFromRowValue(handles.channelTable.Data{row, 2});
    end

    if isfield(handles, 'loadedParams') && isstruct(handles.loadedParams) && ~isempty(fieldnames(handles.loadedParams))
        [fm, h3d] = inferChannelTrainingContext(handles.loadedParams, contextChannel);
        if ~isempty(fm)
            fittingMethod = char(string(fm));
        end
        if ~isempty(h3d)
            has3D = logical(h3d);
        end
    end
end

function [isValid, errMsg] = validateTrainingFeatureSelection(handles, channels)
    isValid = true;
    errMsg = '';

    if ~isfield(handles, 'loadedParams') || ~isstruct(handles.loadedParams) || isempty(fieldnames(handles.loadedParams))
        isValid = false;
        errMsg = 'Feature validation requires a loaded parameter file.';
        return;
    end

    channels = unique(round(channels(:)'));
    problems = {};

    for i = 1:numel(channels)
        ch = channels(i);
        chanCfg = getChannelFeatureConfig(handles, ch);
        selectedFeatures = chanCfg.selectedFeatures;
        customExpressions = chanCfg.customExpressions;
        if isempty(selectedFeatures) && isempty(customExpressions)
            continue;
        end

        [fittingMethod, has3D] = inferChannelTrainingContext(handles.loadedParams, ch);
        if isempty(fittingMethod)
            fittingMethod = '3D Gaussian';
        end
        if isempty(has3D)
            has3D = true;
        end

        [~, featureInfo] = snap_helpers.classification.getAvailableFeatures(char(fittingMethod), logical(has3D), false);
        available = fieldnames(featureInfo);

        missingBase = setdiff(selectedFeatures, available);
        if ~isempty(missingBase)
            problems{end+1} = sprintf('Channel %d missing base feature(s): %s', ... %#ok<AGROW>
                ch, strjoin(missingBase, ', '));
            continue;
        end

        if isempty(customExpressions)
            continue;
        end

        dummyData = buildDummyFeatureStructForValidation(available);
        for e = 1:numel(customExpressions)
            exprName = char(string(customExpressions(e).name));
            exprText = char(string(customExpressions(e).expression));
            warnState = warning('query', 'all');
            warning('off', 'all');
            warnCleanup = onCleanup(@() warning(warnState));
            result = snap_helpers.classification.evaluateExpression(exprText, dummyData, available);
            clear warnCleanup;
            if isempty(result) || all(~isfinite(result))
                problems{end+1} = sprintf('Channel %d custom expression "%s" is incompatible with available features.', ... %#ok<AGROW>
                    ch, exprName);
            end
        end
    end

    if ~isempty(problems)
        isValid = false;
        errMsg = strjoin(problems, newline);
    end
end

function data = buildDummyFeatureStructForValidation(featureNames)
    n = 4;
    data = repmat(struct(), n, 1);
    for i = 1:n
        for f = 1:numel(featureNames)
            data(i).(featureNames{f}) = 1 + 0.1 * i;
        end
    end
end

function onTrainSelectedChannels(fig)
    handles = guidata(fig);

    paramPath = strtrim(handles.paramPathEdit.Value);
    outDir = strtrim(handles.outDirEdit.Value);
    useOptimization = handles.optimizeCheck.Value;
    includeSweepReport = useOptimization && logical(handles.sweepReportCheck.Value);

    if isempty(paramPath) || exist(paramPath, 'file') ~= 2
        uialert(fig, 'Select a valid SNAP parameter file before training.', 'Missing Parameter File');
        return;
    end
    if isempty(outDir)
        uialert(fig, 'Select an output directory.', 'Missing Output Directory');
        return;
    end
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    if isempty(handles.loadedParams) || ~isstruct(handles.loadedParams) || isempty(fieldnames(handles.loadedParams))
        [params, numChannels, errMsg] = loadTrainingParameters(paramPath);
        if ~isempty(errMsg)
            uialert(fig, errMsg, 'Parameter File Error');
            return;
        end
        handles.loadedParams = params;
        guidata(fig, handles);
        updateChannelTableFromDetectedChannels(fig, numChannels, false);
        handles = guidata(fig);
    end

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

    channelConfigs = struct('channel', {}, 'matchDistance', {}, 'convertFijiCoords', {}, ...
        'trainDirectory', {}, 'validationDirectory', {}, 'outputPath', {});
    missingTrainChannels = [];
    missingValidationChannels = [];
    invalidMatchDistanceChannels = [];

    for idx = 1:numel(selectedRows)
        row = selectedRows(idx);
        ch = parseChannelFromRowValue(tableData{row, 2});
        matchDistance = 2;
        convertFiji = false;

        if size(tableData, 2) >= 7
            matchDistance = normalizeScalarNumeric(tableData{row, 3}, nan);
            if ~isfinite(matchDistance) || matchDistance < 0
                invalidMatchDistanceChannels(end+1) = ch; %#ok<AGROW>
            end
            convertFiji = isSelectedChannelFlag(tableData{row, 4});
        elseif size(tableData, 2) >= 3
            % Backward compatibility with older layout.
            convertFiji = isSelectedChannelFlag(tableData{row, 3});
        end

        trainDir = '';
        if size(tableData, 2) >= 7
            trainDir = strtrim(char(string(tableData{row, 5})));
        elseif size(tableData, 2) >= 4
            trainDir = strtrim(char(string(tableData{row, 4})));
        end

        valDir = '';
        if size(tableData, 2) >= 7
            valDir = strtrim(char(string(tableData{row, 6})));
        elseif size(tableData, 2) >= 5
            valDir = strtrim(char(string(tableData{row, 5})));
        end

        outName = '';
        if size(tableData, 2) >= 7
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
            'trainDirectory', trainDir, ...
            'validationDirectory', valDir, ...
            'outputPath', fullfile(outDir, outName));
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

    [featureSelectionOk, featureSelectionErr] = validateTrainingFeatureSelection(handles, [channelConfigs.channel]);
    if ~featureSelectionOk
        uialert(fig, featureSelectionErr, 'Feature Selection Incompatible');
        return;
    end

    try
        trainOptions = buildManualTrainingOptions(handles);
        if useOptimization
            sweepConfig = buildSweepConfig(handles);
        else
            sweepConfig = struct();
        end
    catch ME
        uialert(fig, ME.message, 'Invalid Training Options');
        return;
    end

    appendTrainLog(fig, '------------------------------------------------------------');
    appendTrainLog(fig, sprintf('Starting training (%d selected channels)...', numel(selectedRows)));
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
        stageCount = 3 + double(useOptimization);
        stage = 0;

        stage = stage + 1;
        updateTrainProgressStage(fig, idx, numel(channelConfigs), stage, stageCount, ...
            sprintf('Channel %d: Discovering training pairs (%s)...', ch, outName));
        appendTrainLog(fig, sprintf('[Channel %d] Match distance: %.4g voxels', ch, cfg.matchDistance));
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
            'ParameterStruct', handles.loadedParams, ...
            'ChannelIndex', ch, ...
            'SelectedFeatures', chFeatureCfg.selectedFeatures, ...
            'CustomExpressions', chFeatureCfg.customExpressions, ...
            'TrainingOptions', trainOptions, ...
            'HyperparameterSweep', false, ...
            'Verbose', true, ...
            'ProgressCallback', trainingProgressCb ...
        };

        [fittingMethod, has3D] = inferChannelTrainingContext(handles.loadedParams, ch);
        if ~isempty(fittingMethod)
            args = [args, {'FittingMethod', fittingMethod}]; %#ok<AGROW>
        end
        if ~isempty(has3D)
            args = [args, {'Has3D', has3D}]; %#ok<AGROW>
        end

        if useOptimization
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

            args = [args, { ...
                'ValidationExportFiles', valExports, ...
                'ValidationLabelFiles', valLabels, ...
                'HyperparameterSweep', true, ...
                'SweepKernels', sweepConfig.kernels, ...
                'SweepBoxConstraints', sweepConfig.boxConstraints, ...
                'SweepKernelScales', sweepConfig.kernelScales, ...
                'SweepPolynomialOrders', sweepConfig.polynomialOrders ...
            }]; %#ok<AGROW>
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
            if includeSweepReport
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
    summary.parameterFile = paramPath;
    summary.outputDirectory = outDir;
    summary.matchDistanceByChannel = arrayfun(@(c) c.matchDistance, channelConfigs);
    summary.optimized = useOptimization;
    summary.sweepReportIncluded = includeSweepReport;
    summary.channelFeatureConfigs = handles.channelFeatureConfigs;
    summary.channelConfigs = channelConfigs;
    summary.results = results;
    summary.successCount = nSuccess;
    summary.failCount = nFailed;

    summaryPath = fullfile(outDir, ['SNAP_train_summary_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
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
        raw = load(paramPath);
        if isfield(raw, 'batchConfig') && isfield(raw.batchConfig, 'parameters')
            params = raw.batchConfig.parameters;
            if isfield(raw.batchConfig, 'workflowConfig') && isstruct(raw.batchConfig.workflowConfig)
                if ~isfield(params, 'numChannels') && isfield(raw.batchConfig.workflowConfig, 'numChannels')
                    params.numChannels = raw.batchConfig.workflowConfig.numChannels;
                end
            end
        elseif isfield(raw, 'paramData') && isfield(raw.paramData, 'parameters')
            params = raw.paramData.parameters;
            if isfield(raw.paramData, 'workflowConfig') && isstruct(raw.paramData.workflowConfig)
                if ~isfield(params, 'numChannels') && isfield(raw.paramData.workflowConfig, 'numChannels')
                    params.numChannels = raw.paramData.workflowConfig.numChannels;
                end
            end
        elseif isfield(raw, 'paramData') && isfield(raw.paramData, 'workflowConfig')
            params = raw.paramData.workflowConfig;
        elseif isfield(raw, 'workflowConfig')
            params = raw.workflowConfig;
        elseif isfield(raw, 'parameters')
            params = raw.parameters;
        elseif isfield(raw, 'lastUsed')
            params = raw.lastUsed;
        else
            errMsg = ['Could not find parameter struct (expected one of: ', ...
                'batchConfig.parameters, paramData.parameters, paramData.workflowConfig, workflowConfig, parameters, or lastUsed).'];
            return;
        end

        numChannels = inferNumChannelsFromParameters(params, 1);
    catch ME
        errMsg = sprintf('Failed to load parameter file: %s', ME.message);
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

function options = buildManualTrainingOptions(handles)
    options = struct();
    options.kernelFunction = char(string(handles.kernelDrop.Value));
    options.boxConstraint = handles.boxConstraintInput.Value;
    options.kernelScale = parseKernelScale(handles.kernelScaleInput.Value);
    options.polynomialOrder = round(handles.polyOrderInput.Value);
    options.standardize = logical(handles.standardizeCheck.Value);
    options.crossValidate = logical(handles.crossValidateCheck.Value);
    options.kFold = max(2, round(handles.kFoldInput.Value));
    if logical(handles.balanceClassCheck.Value)
        options.classWeightMode = 'balanced';
    else
        options.classWeightMode = 'none';
    end
    options.verbose = true;
end

function sweep = buildSweepConfig(handles)
    sweep = struct();
    sweep.kernels = parseKernelList(handles.sweepKernelsInput.Value);
    sweep.boxConstraints = parseNumericList(handles.sweepBoxInput.Value, true);
    sweep.kernelScales = parseMixedScaleList(handles.sweepScaleInput.Value);
    sweep.polynomialOrders = round(parseNumericList(handles.sweepPolyInput.Value, true));
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

function [exportFiles, labelFiles] = discoverImageCsvPairs(rootDir, channelIdx, progressCb)
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
    primaryExports = {};
    primaryLabels = {};
    fallbackExports = {};
    fallbackLabels = {};
    hasTaggedMatch = false;

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

        imageChannel = parseChannelFromFilename(imagePath);
        if ~isempty(channelIdx)
            if isfinite(imageChannel)
                if imageChannel ~= channelIdx
                    continue;
                end
                hasTaggedMatch = true;
            end
        end

        if ~isempty(channelIdx) && ~isfinite(imageChannel)
            fallbackExports{end+1,1} = imagePath; %#ok<AGROW>
            fallbackLabels{end+1,1} = labelPath; %#ok<AGROW>
        else
            primaryExports{end+1,1} = imagePath; %#ok<AGROW>
            primaryLabels{end+1,1} = labelPath; %#ok<AGROW>
        end

        if ~isempty(progressCb) && (numel(imageFiles) <= 20 || mod(i, 10) == 0)
            [~, imgBase, imgExt] = fileparts(imagePath);
            [~, labBase, labExt] = fileparts(labelPath);
            progressCb(sprintf('Matched image %d/%d: %s%s <-> %s%s', ...
                i, numel(imageFiles), imgBase, imgExt, labBase, labExt));
        end
    end
    if ~isempty(progressCb)
        progressCb(sprintf('Pairing loop complete: %d primary pair(s), %d fallback pair(s).', ...
            numel(primaryExports), numel(fallbackExports)));
    end

    if isempty(channelIdx)
        exportFiles = primaryExports;
        labelFiles = primaryLabels;
        return;
    end

    if hasTaggedMatch && ~isempty(primaryExports)
        exportFiles = primaryExports;
        labelFiles = primaryLabels;
    elseif isempty(primaryExports)
        exportFiles = fallbackExports;
        labelFiles = fallbackLabels;
    else
        exportFiles = primaryExports;
        labelFiles = primaryLabels;
    end
end

function labelPath = findExactCsvForImageInList(imagePath, csvFiles)
    labelPath = '';
    [imageDir, imageBase, ~] = fileparts(imagePath);
    imageBaseNoOme = regexprep(imageBase, '\.ome$', '', 'ignorecase');

    fallbackPath = '';
    for i = 1:numel(csvFiles)
        csvPath = csvFiles{i};
        [csvDir, csvBase, ~] = fileparts(csvPath);
        if ~strcmpi(csvDir, imageDir)
            continue;
        end

        if strcmpi(csvBase, imageBase)
            labelPath = csvPath;
            return;
        end

        csvBaseNoOme = regexprep(csvBase, '\.ome$', '', 'ignorecase');
        if isempty(fallbackPath) && strcmpi(csvBaseNoOme, imageBaseNoOme)
            fallbackPath = csvPath;
        end
    end

    if ~isempty(fallbackPath)
        labelPath = fallbackPath;
    end
end

function [exportFiles, labelFiles] = discoverLegacyMatPairs(rootDir, channelIdx, progressCb)
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
    fallbackExports = {};
    fallbackLabels = {};
    hasChannelTaggedMatch = false;

    for i = 1:numel(matFiles)
        matPath = matFiles{i};
        if ~isExportFile(matPath)
            continue;
        end

        exportChannel = parseChannelFromFilename(matPath);
        if ~isempty(channelIdx)
            if isfinite(exportChannel)
                if exportChannel ~= channelIdx
                    continue;
                end
            end
        end

        labelPath = findBestLabelForExport(matPath, labelCandidates, channelIdx);
        if ~isempty(labelPath)
            if ~isempty(channelIdx) && ~isfinite(exportChannel)
                fallbackExports{end+1,1} = matPath; %#ok<AGROW>
                fallbackLabels{end+1,1} = labelPath; %#ok<AGROW>
            else
                exportFiles{end+1,1} = matPath; %#ok<AGROW>
                labelFiles{end+1,1} = labelPath; %#ok<AGROW>
                if ~isempty(channelIdx)
                    hasChannelTaggedMatch = true;
                end
            end
        end
    end

    if ~isempty(channelIdx) && ~hasChannelTaggedMatch && isempty(exportFiles)
        exportFiles = fallbackExports;
        labelFiles = fallbackLabels;
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
    results = repmat(struct('params', struct(), 'f1Real', nan, 'accuracy', nan, 'trainStats', struct()), numel(combos), 1);
    emitProgress(progressCb, 'Validation sweep: evaluating %d combination(s).', numel(combos));

    bestScore = -inf;
    bestIdx = 1;
    model = [];
    trainStats = struct();
    normParams = struct();
    bestParams = baseOptions;

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
        if ~isfield(stats, 'success') || ~stats.success
            results(i).params = combos(i);
            results(i).trainStats = stats;
            emitProgress(progressCb, 'Validation sweep %d/%d failed: %s', ...
                i, numel(combos), getfieldwithdefault(stats, 'error', 'unknown training error'));
            continue;
        end

        [metrics, ~] = evaluateOnValidation(mdl, norm, XVal, yVal);

        results(i).params = combos(i);
        results(i).f1Real = metrics.f1Real;
        results(i).accuracy = metrics.accuracy;
        results(i).trainStats = stats;
        emitProgress(progressCb, 'Validation sweep %d/%d metrics: F1(real)=%.4f, accuracy=%.4f.', ...
            i, numel(combos), metrics.f1Real, metrics.accuracy);

        score = metrics.f1Real + 1e-3 * metrics.accuracy;
        if score > bestScore
            bestScore = score;
            bestIdx = i;
            model = mdl;
            trainStats = stats;
            normParams = norm;
            bestParams = options;
            emitProgress(progressCb, 'Validation sweep new best at %d/%d: %s', ...
                i, numel(combos), formatSweepParameterSummary(bestParams));
        end
    end

    if isempty(model)
        error('Hyperparameter sweep failed: all parameter combinations failed to train.');
    end

    if verbose
        fprintf('Hyperparameter sweep complete: %d combinations tested.\n', numel(combos));
        fprintf('Best combo #%d: kernel=%s, C=%.3g, f1_real=%.3f\n', ...
            bestIdx, bestParams.kernelFunction, bestParams.boxConstraint, results(bestIdx).f1Real);
    end
    emitProgress(progressCb, 'Validation sweep complete. Best combo #%d with F1(real)=%.4f.', ...
        bestIdx, results(bestIdx).f1Real);
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
