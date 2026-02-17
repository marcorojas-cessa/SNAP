function out = SNAP_trainSVM(exportFiles, labelFiles, outputClassifierPath, varargin)
% SNAP_trainSVM - Train and validate a SNAP-compatible SVM from labeled spot files
%
% INTERACTIVE MODE:
%   SNAP_trainSVM
%   Prompts for:
%     - Match distance (voxels)
%     - Training directory containing exported candidate MAT files + labels
%     - Validation directory containing exported candidate MAT files + labels
%     - Output classifier path
%   Then performs hyperparameter sweep on validation F1 and saves best model.
%   Training real/noise labels are fixed from training volumes and never relabeled during validation.
%
% PROGRAMMATIC MODE:
%   SNAP_trainSVM(exportFiles, labelFiles, outputClassifierPath, ...)
%
% REQUIRED INPUTS:
%   exportFiles           - char/string/cellstr of SNAP export MAT file(s)
%   labelFiles            - char/string/cellstr of label file(s)
%   outputClassifierPath  - output .mat classifier path
%
% LABEL FILE FORMATS:
%   1) CSV/table with columns:
%      - maxima_y, maxima_x, [maxima_z], and label
%      - OR fitted_y, fitted_x, [fitted_z], and label
%      - OR y, x, [z], and label
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
%
% OUTPUT:
%   out - struct with trained model, selected params, train/validation stats

    if nargin == 0
        out = runInteractiveTraining();
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
    p.parse(varargin{:});
    opts = p.Results;

    verbose = logical(opts.Verbose);

    [trainFitData, trainLabels, trainSources, inferredFittingMethod, inferredHas3D] = ...
        buildLabeledDataset(exportFiles, labelFiles, opts.MatchDistance);

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

    XTrain = XTrain(validTrainMask, :);
    yTrain = trainLabels(validTrainMask);
    trainSources = trainSources(validTrainMask);
    fixedTrainLabels = yTrain; % immutable reference labels used for all sweep candidates

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
        [valFitData, valLabelsRaw] = buildLabeledDataset(valExport, valLabels, opts.MatchDistance);
        [XVal, ~, validValMask, extractionInfoVal] = snap_helpers.classification.buildFeatureMatrix( ...
            valFitData, selectedFeatures, featureInfo, opts.CustomExpressions);
        XVal = XVal(validValMask, :);
        yVal = valLabelsRaw(validValMask);

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
            [model, trainStats, normParams, bestParams, sweepResults] = runHyperparameterSweep( ...
                XTrain, yTrain, fixedTrainLabels, XVal, yVal, trainOptions, sweepCfg, verbose);
            validation.sweepResults = sweepResults;
        else
            trainOptions = applyTrainingDefaults(opts.TrainingOptions);
            [model, trainStats, normParams] = snap_helpers.classification.trainClassifier(XTrain, yTrain, trainOptions);
            bestParams = trainOptions;
        end

        [validation.metrics, validation.confusionMatrix] = evaluateOnValidation(model, normParams, XVal, yVal);
    else
        trainOptions = applyTrainingDefaults(opts.TrainingOptions);
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
    metadata.trainingSource = 'SNAP_trainSVM';
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
    out.bestParams = bestParams;
    out.outputClassifierPath = outputClassifierPath;
    if hasValidationSet
        out.validation = validation;
    end

    if verbose
        fprintf('\nSNAP_trainSVM complete:\n');
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
end

function out = runInteractiveTraining()
    fprintf('\nSNAP_trainSVM interactive setup\n');
    fprintf('Provide train/validation directories containing export MAT files and label files.\n');

    matchDistance = str2double(strtrim(input('Match distance in voxels [2]: ', 's')));
    if isnan(matchDistance)
        matchDistance = 2;
    end

    trainDir = strtrim(input('Training directory path: ', 's'));
    valDir = strtrim(input('Validation directory path: ', 's'));
    outputClassifierPath = strtrim(input('Output classifier path (e.g., classifier.mat): ', 's'));

    if isempty(trainDir) || isempty(valDir) || isempty(outputClassifierPath)
        error('Training directory, validation directory, and output path are required in interactive mode.');
    end

    [trainExports, trainLabels] = discoverDatasetPairs(trainDir);
    [valExports, valLabels] = discoverDatasetPairs(valDir);

    fprintf('Discovered %d training pairs and %d validation pairs.\n', numel(trainExports), numel(valExports));

    out = SNAP_trainSVM(trainExports, trainLabels, outputClassifierPath, ...
        'MatchDistance', matchDistance, ...
        'ValidationExportFiles', valExports, ...
        'ValidationLabelFiles', valLabels, ...
        'HyperparameterSweep', true, ...
        'Verbose', true);
end

function [exportFiles, labelFiles] = discoverDatasetPairs(rootDir)
    if ~isfolder(rootDir)
        error('Directory not found: %s', rootDir);
    end

    matFiles = dir(fullfile(rootDir, '*.mat'));
    exportFiles = {};
    labelFiles = {};

    for i = 1:numel(matFiles)
        matPath = fullfile(rootDir, matFiles(i).name);
        if ~isExportFile(matPath)
            continue;
        end

        [~, base, ~] = fileparts(matFiles(i).name);
        key = regexprep(base, '_signals$', '', 'ignorecase');

        candidates = {
            fullfile(rootDir, [key '_labels.csv']), ...
            fullfile(rootDir, [key '_labels.mat']), ...
            fullfile(rootDir, [key '.csv']), ...
            fullfile(rootDir, [key '.mat']) ...
        };

        labelPath = '';
        for j = 1:numel(candidates)
            c = candidates{j};
            if exist(c, 'file') == 2 && ~strcmp(c, matPath)
                if endsWith(lower(c), '.mat') && isExportFile(c)
                    continue;
                end
                labelPath = c;
                break;
            end
        end

        if ~isempty(labelPath)
            exportFiles{end+1,1} = matPath; %#ok<AGROW>
            labelFiles{end+1,1} = labelPath; %#ok<AGROW>
        end
    end

    if isempty(exportFiles)
        error('No export/label pairs discovered in %s', rootDir);
    end
end

function tf = isExportFile(matPath)
    tf = false;
    info = whos('-file', matPath);
    names = string({info.name});
    tf = any(ismember(names, ["signals", "fitResults", "fit_results", "exportData"]));
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

function [allFitData, allLabels, sourceImage, inferredFittingMethod, inferredHas3D] = ...
        buildLabeledDataset(exportFiles, labelFiles, matchDistance)

    allFitData = struct([]);
    allLabels = [];
    sourceImage = strings(0,1);
    inferredFittingMethod = '';
    inferredHas3D = false;

    for i = 1:numel(exportFiles)
        fitData = loadFitDataFromExport(exportFiles{i});
        labels = loadLabelsForFitData(labelFiles{i}, fitData, matchDistance);
        labels(isnan(labels)) = 0;

        if isempty(fitData)
            continue;
        end

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
end

function [model, trainStats, normParams, bestParams, results] = runHyperparameterSweep(XTrain, yTrain, fixedTrainLabels, XVal, yVal, baseOptions, cfg, verbose)
    combos = buildSweepCombinations(cfg);
    results = repmat(struct('params', struct(), 'f1Real', nan, 'accuracy', nan, 'trainStats', struct()), numel(combos), 1);

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

        if ~isequal(yTrain(:), fixedTrainLabels(:))
            error('Training labels changed during validation sweep, which is not allowed.');
        end

        [mdl, stats, norm] = snap_helpers.classification.trainClassifier(XTrain, yTrain, options);
        if ~isfield(stats, 'success') || ~stats.success
            results(i).params = combos(i);
            results(i).trainStats = stats;
            continue;
        end

        [metrics, ~] = evaluateOnValidation(mdl, norm, XVal, yVal);

        results(i).params = combos(i);
        results(i).f1Real = metrics.f1Real;
        results(i).accuracy = metrics.accuracy;
        results(i).trainStats = stats;

        score = metrics.f1Real + 1e-3 * metrics.accuracy;
        if score > bestScore
            bestScore = score;
            bestIdx = i;
            model = mdl;
            trainStats = stats;
            normParams = norm;
            bestParams = options;
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

function fitData = loadFitDataFromExport(filepath)
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
        error('Could not find spot data in export file: %s', filepath);
    end

    if istable(fitData)
        fitData = table2struct(fitData);
    end

    fitData = normalizeFitDataLocal(fitData);
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

function labels = loadLabelsForFitData(labelFile, fitData, matchDistance)
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
            t = readtable(labelFile);
            labels = labelsFromTable(t, fitData, matchDistance);
    end
end

function labels = labelsFromTable(t, fitData, matchDistance)
    vn = lower(string(t.Properties.VariableNames));
    labelIdx = find(vn == "label" | vn == "class" | vn == "is_real", 1);
    if isempty(labelIdx)
        error('Label table must contain one of: label, class, is_real');
    end

    rawLabel = t.(t.Properties.VariableNames{labelIdx});
    y = normalizeLabelValues(rawLabel);

    labels = nan(numel(fitData), 1);

    yName = find(vn == "maxima_y" | vn == "fitted_y" | vn == "y", 1);
    xName = find(vn == "maxima_x" | vn == "fitted_x" | vn == "x", 1);
    zName = find(vn == "maxima_z" | vn == "fitted_z" | vn == "z", 1);

    if isempty(yName) || isempty(xName)
        error('Label table must include coordinate columns (x,y[,z]) plus label.');
    end

    yCoord = t.(t.Properties.VariableNames{yName});
    xCoord = t.(t.Properties.VariableNames{xName});
    if isempty(zName)
        zCoord = ones(height(t), 1);
    else
        zCoord = t.(t.Properties.VariableNames{zName});
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
    if ~isfield(options, 'standardize'), options.standardize = true; end
    if ~isfield(options, 'crossValidate'), options.crossValidate = true; end
    if ~isfield(options, 'kFold'), options.kFold = 5; end
    if ~isfield(options, 'verbose'), options.verbose = true; end
end
