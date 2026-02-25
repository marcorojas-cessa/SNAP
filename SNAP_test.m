function varargout = SNAP_test(varargin)
% SNAP_test - Quick GUI to evaluate a trained SNAP SVM on labeled image data.
%
% This is a lightweight utility for local testing (not part of core release flow).
% It evaluates one trained classifier against a dataset directory containing
% image/CSV pairs (same-name matching, recursive search).
%
% Required UI inputs:
%   - Parameter file (.mat) used to generate candidates from images
%   - Trained classifier file (.mat) from SNAP_train
%   - Dataset directory with .tif/.tiff and matching .csv label files
%   - Channel index (single channel per run)
%   - Match distance (voxel units)
%   - Optional FIJI coordinate conversion
%
% Notes:
%   - Candidate generation uses the SNAP signal pipeline (same as SNAP_train).
%   - CSV labels are matched to nearest unassigned candidate within distance.
%   - Unmatched candidates are labeled as noise (0), consistent with SNAP_train.

    if nargin > 0
        error('SNAP_test does not accept input arguments. Launch it with: SNAP_test');
    end

    ui = createUi();
    if nargout > 0
        varargout{1} = ui;
    end
end

function ui = createUi()
    fig = figure( ...
        'Name', 'SNAP_test - SVM Evaluation Utility', ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'none', ...
        'Color', [0.94 0.94 0.94], ...
        'Units', 'pixels', ...
        'Position', [120 80 1220 760], ...
        'Resize', 'on');

    inputPanel = uipanel( ...
        'Parent', fig, ...
        'Title', 'Inputs', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.73 0.98 0.26]);

    summaryPanel = uipanel( ...
        'Parent', fig, ...
        'Title', 'Summary', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.58 0.98 0.14]);

    tablePanel = uipanel( ...
        'Parent', fig, ...
        'Title', 'Per-Image Results', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.29 0.98 0.28]);

    logPanel = uipanel( ...
        'Parent', fig, ...
        'Title', 'Log', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.27]);

    % Row 1: Parameter file
    uicontrol(inputPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01 0.78 0.14 0.14], ...
        'String', 'Parameter File:', 'HorizontalAlignment', 'left');
    paramEdit = uicontrol(inputPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.15 0.78 0.70 0.15], ...
        'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
    uicontrol(inputPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.86 0.78 0.12 0.15], ...
        'String', 'Browse', ...
        'Callback', @(~,~)onBrowseFile(paramEdit, {'*.mat', 'MAT Files (*.mat)'}, 'Select SNAP parameter file'));

    % Row 2: Classifier file
    uicontrol(inputPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01 0.56 0.14 0.14], ...
        'String', 'Classifier File:', 'HorizontalAlignment', 'left');
    clfEdit = uicontrol(inputPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.15 0.56 0.70 0.15], ...
        'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
    uicontrol(inputPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.86 0.56 0.12 0.15], ...
        'String', 'Browse', ...
        'Callback', @(~,~)onBrowseFile(clfEdit, {'*.mat', 'MAT Files (*.mat)'}, 'Select classifier file'));

    % Row 3: Dataset directory
    uicontrol(inputPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01 0.34 0.14 0.14], ...
        'String', 'Dataset Folder:', 'HorizontalAlignment', 'left');
    datasetEdit = uicontrol(inputPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.15 0.34 0.70 0.15], ...
        'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
    uicontrol(inputPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.86 0.34 0.12 0.15], ...
        'String', 'Browse', ...
        'Callback', @(~,~)onBrowseFolder(datasetEdit, 'Select dataset directory'));

    % Row 4 controls
    uicontrol(inputPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.01 0.13 0.10 0.14], ...
        'String', 'Channel:', 'HorizontalAlignment', 'left');
    channelEdit = uicontrol(inputPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.11 0.12 0.08 0.15], ...
        'String', '1', ...
        'BackgroundColor', 'white', 'HorizontalAlignment', 'center');

    uicontrol(inputPanel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.22 0.13 0.15 0.14], ...
        'String', 'Match Dist (vox):', 'HorizontalAlignment', 'left');
    matchEdit = uicontrol(inputPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.37 0.12 0.08 0.15], ...
        'String', '2', ...
        'BackgroundColor', 'white', 'HorizontalAlignment', 'center');

    fijiCheck = uicontrol(inputPanel, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.48 0.12 0.29 0.15], ...
        'String', 'CSV coordinates are FIJI [x,y,z] (0-indexed)', ...
        'Value', 0, ...
        'HorizontalAlignment', 'left');

    runButton = uicontrol(inputPanel, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.79 0.08 0.19 0.18], ...
        'String', 'Run Evaluation', ...
        'FontWeight', 'bold', ...
        'Callback', @onRunEvaluation);

    summaryEdit = uicontrol(summaryPanel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.01 0.05 0.98 0.90], ...
        'Max', 20, 'Min', 0, ...
        'Enable', 'inactive', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'white', ...
        'String', {'No evaluation run yet.'});

    columnNames = { ...
        'Image', 'Candidates', 'Evaluated', 'Real', 'Noise', ...
        'TP', 'FP', 'FN', 'TN', ...
        'Precision', 'Recall', 'F1', 'Accuracy'};
    resultTable = uitable(tablePanel, ...
        'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.98], ...
        'Data', cell(0, numel(columnNames)), ...
        'ColumnName', columnNames, ...
        'ColumnEditable', false(1, numel(columnNames)));

    logList = uicontrol(logPanel, 'Style', 'listbox', 'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.98], ...
        'BackgroundColor', 'white', ...
        'FontName', 'Consolas', ...
        'String', {'SNAP_test ready.'});

    ui = struct();
    ui.figure = fig;
    ui.controls = struct( ...
        'paramEdit', paramEdit, ...
        'clfEdit', clfEdit, ...
        'datasetEdit', datasetEdit, ...
        'channelEdit', channelEdit, ...
        'matchEdit', matchEdit, ...
        'fijiCheck', fijiCheck, ...
        'runButton', runButton, ...
        'summaryEdit', summaryEdit, ...
        'resultTable', resultTable, ...
        'logList', logList);

    function onBrowseFile(targetEdit, filterSpec, dialogTitle)
        startDir = initialBrowseDir(get(targetEdit, 'String'));
        [fileName, folder] = uigetfile(filterSpec, dialogTitle, startDir);
        if isequal(fileName, 0) || isequal(folder, 0)
            return;
        end
        set(targetEdit, 'String', fullfile(folder, fileName));
    end

    function onBrowseFolder(targetEdit, dialogTitle)
        startDir = initialBrowseDir(get(targetEdit, 'String'));
        selected = uigetdir(startDir, dialogTitle);
        if isequal(selected, 0)
            return;
        end
        set(targetEdit, 'String', selected);
    end

    function onRunEvaluation(~, ~)
        opts = collectOptionsFromUi(ui.controls);
        if isempty(opts)
            return;
        end

        partialTableData = cell(0, numel(columnNames));
        set(ui.controls.runButton, 'Enable', 'off');
        set(ui.controls.resultTable, 'Data', cell(0, numel(columnNames)));
        set(ui.controls.summaryEdit, 'String', {'Running evaluation...'});
        set(ui.controls.logList, 'String', {'SNAP_test run started.'}, 'Value', 1);
        drawnow;

        logProgress('Evaluation started.');
        tStart = tic;
        try
            [report, rows, rowData] = runSvmEvaluation(opts, @logProgress, @onRowEvaluated);
            set(ui.controls.resultTable, 'Data', rowData);
            set(ui.controls.summaryEdit, 'String', reportToSummaryLines(report));
            logProgress('Evaluation completed in %.2f s.', toc(tStart));
            logProgress('Final metrics: accuracy=%.4f, precision=%.4f, recall=%.4f, F1=%.4f', ...
                report.accuracy, report.precisionReal, report.recallReal, report.f1Real);
            if any([rows.nInvalidFeatures] > 0)
                totalInvalid = sum([rows.nInvalidFeatures]);
                logProgress('Note: %d candidate(s) were excluded due to invalid feature rows.', totalInvalid);
            end
        catch ME
            set(ui.controls.summaryEdit, 'String', {'Evaluation failed. See log below.'});
            logProgress('ERROR: %s', ME.message);
            for s = 1:numel(ME.stack)
                logProgress('  at %s (line %d)', ME.stack(s).name, ME.stack(s).line);
            end
            errordlg(ME.message, 'SNAP_test Error', 'modal');
        end

        set(ui.controls.runButton, 'Enable', 'on');
        drawnow;

        function onRowEvaluated(rowIdx, rowStruct, totalRows)
            if nargin < 3 || ~isfinite(totalRows) || totalRows < 1
                totalRows = max(1, rowIdx);
            end
            if isempty(partialTableData) || size(partialTableData, 1) ~= totalRows
                partialTableData = cell(totalRows, numel(columnNames));
            end
            if isfinite(rowIdx) && rowIdx >= 1 && rowIdx <= size(partialTableData, 1)
                partialTableData(rowIdx, :) = rowStructToTableRow(rowStruct);
                set(ui.controls.resultTable, 'Data', partialTableData);
                drawnow limitrate;
            end
        end
    end

    function logProgress(msgFmt, varargin)
        msg = formatMessage(msgFmt, varargin{:});
        timestamp = datestr(now, 'HH:MM:SS');
        entry = sprintf('[%s] %s', timestamp, msg);
        current = get(ui.controls.logList, 'String');
        if ischar(current)
            current = {current};
        end
        current{end+1,1} = entry;
        set(ui.controls.logList, 'String', current, 'Value', numel(current));
        drawnow limitrate;
    end
end

function opts = collectOptionsFromUi(ctrl)
    opts = [];

    paramPath = sanitizePathInput(get(ctrl.paramEdit, 'String'));
    clfPath = sanitizePathInput(get(ctrl.clfEdit, 'String'));
    datasetDir = sanitizePathInput(get(ctrl.datasetEdit, 'String'));

    if isempty(paramPath) || exist(paramPath, 'file') ~= 2
        errordlg('Please select a valid parameter file (.mat).', 'Missing Input', 'modal');
        return;
    end
    if isempty(clfPath) || exist(clfPath, 'file') ~= 2
        errordlg('Please select a valid classifier file (.mat).', 'Missing Input', 'modal');
        return;
    end
    if isempty(datasetDir) || exist(datasetDir, 'dir') ~= 7
        errordlg('Please select a valid dataset directory.', 'Missing Input', 'modal');
        return;
    end

    channelIdx = str2double(strtrim(char(string(get(ctrl.channelEdit, 'String')))));
    if ~isfinite(channelIdx) || channelIdx < 1
        errordlg('Channel must be a positive number.', 'Invalid Input', 'modal');
        return;
    end

    matchDist = str2double(strtrim(char(string(get(ctrl.matchEdit, 'String')))));
    if ~isfinite(matchDist) || matchDist < 0
        errordlg('Match distance must be a non-negative number.', 'Invalid Input', 'modal');
        return;
    end

    opts = struct();
    opts.parameterFile = paramPath;
    opts.classifierFile = clfPath;
    opts.datasetDir = datasetDir;
    opts.channelIndex = max(1, round(channelIdx));
    opts.matchDistance = double(matchDist);
    opts.convertFijiCoords = logical(get(ctrl.fijiCheck, 'Value'));
end

function dirPath = initialBrowseDir(currentText)
    dirPath = pwd;
    if nargin < 1 || isempty(currentText)
        return;
    end
    currentText = char(string(currentText));
    if exist(currentText, 'dir') == 7
        dirPath = currentText;
        return;
    end
    [parentDir, ~, ~] = fileparts(currentText);
    if exist(parentDir, 'dir') == 7
        dirPath = parentDir;
    end
end

function p = sanitizePathInput(raw)
    p = strtrim(char(string(raw)));
    if numel(p) >= 2
        if (p(1) == '"' && p(end) == '"') || (p(1) == '''' && p(end) == '''')
            p = p(2:end-1);
        end
    end
    p = strtrim(p);
end

function [report, rows, rowData] = runSvmEvaluation(opts, progressCb, rowUpdateCb)
    tEval = tic;
    if nargin < 2 || isempty(progressCb)
        progressCb = [];
    end
    if nargin < 3 || isempty(rowUpdateCb)
        rowUpdateCb = [];
    end
    emitProgress(progressCb, 'Loading parameters: %s', opts.parameterFile);
    [params, meta] = snap_helpers.classification.loadParameterStruct(opts.parameterFile);
    if isfield(meta, 'container') && ~isempty(meta.container)
        emitProgress(progressCb, 'Parameter container: %s', meta.container);
    end
    if isfield(meta, 'warnings') && ~isempty(meta.warnings)
        for i = 1:numel(meta.warnings)
            emitProgress(progressCb, 'Parameter warning: %s', meta.warnings{i});
        end
    end

    emitProgress(progressCb, 'Loading classifier: %s', opts.classifierFile);
    [model, baseFeatures, featureInfo, fittingMethod, customExpr, normParams] = ...
        loadClassifierForTesting(opts.classifierFile, progressCb);
    if isempty(baseFeatures)
        error('Classifier file contains no feature list.');
    end
    baseFeatures = normalizeFeatureNameList(baseFeatures);
    customExpr = normalizeCustomExpressionListLocal(customExpr);

    emitProgress(progressCb, 'Discovering image/CSV pairs in: %s', opts.datasetDir);
    [imageFiles, csvFiles] = discoverImageCsvPairs(opts.datasetDir, progressCb);
    nPairs = numel(imageFiles);
    if nPairs == 0
        error('No matching image/CSV pairs found in directory: %s', opts.datasetDir);
    end
    emitProgress(progressCb, 'Found %d pair(s).', nPairs);

    rowTemplate = struct( ...
        'imagePath', '', ...
        'nCandidates', 0, ...
        'nEvaluated', 0, ...
        'nReal', 0, ...
        'nNoise', 0, ...
        'tp', 0, ...
        'fp', 0, ...
        'fn', 0, ...
        'tn', 0, ...
        'precisionReal', NaN, ...
        'recallReal', NaN, ...
        'f1Real', NaN, ...
        'accuracy', NaN, ...
        'nInvalidFeatures', 0);
    rows = repmat(rowTemplate, nPairs, 1);

    totalCandidates = 0;
    totalEvaluated = 0;
    totalReal = 0;
    totalNoise = 0;
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;

    classifierFeatureNames = [reshape(baseFeatures, 1, []), reshape({customExpr.name}, 1, [])];

    for i = 1:nPairs
        imagePath = imageFiles{i};
        csvPath = csvFiles{i};
        [~, imageName, imageExt] = fileparts(imagePath);
        emitProgress(progressCb, 'Pair %d/%d: %s%s', i, nPairs, imageName, imageExt);

        fitData = buildFitDataFromImage(imagePath, params, opts.channelIndex, progressCb);
        nCandidates = numel(fitData);
        rows(i).imagePath = imagePath;
        rows(i).nCandidates = nCandidates;
        totalCandidates = totalCandidates + nCandidates;

        if nCandidates == 0
            if isa(rowUpdateCb, 'function_handle')
                try
                    rowUpdateCb(i, rows(i), nPairs);
                catch ME
                    emitProgress(progressCb, 'Row update callback warning: %s', ME.message);
                end
            end
            emitProgress(progressCb, 'Pair %d/%d: no candidates generated; skipping.', i, nPairs);
            continue;
        end

        labels = loadLabelsForFitData(csvPath, fitData, opts.matchDistance, opts.convertFijiCoords);
        labels(isnan(labels)) = 0;
        rows(i).nReal = sum(labels == 1);
        rows(i).nNoise = sum(labels == 0);
        totalReal = totalReal + rows(i).nReal;
        totalNoise = totalNoise + rows(i).nNoise;

        featureInfoLocal = featureInfo;
        if isempty(featureInfoLocal) || ~isstruct(featureInfoLocal) || isempty(fieldnames(featureInfoLocal))
            method = char(string(fittingMethod));
            if isempty(strtrim(method)) || strcmpi(method, 'Unknown')
                method = inferFittingMethod(fitData);
            end
            has3D = inferHas3D(fitData);
            [~, featureInfoLocal] = snap_helpers.classification.getAvailableFeatures(method, has3D, false);
            emitProgress(progressCb, 'Feature metadata inferred from method="%s", has3D=%d.', method, has3D);
        end

        [X, featureNames, validMask, extractionInfo] = snap_helpers.classification.buildFeatureMatrix( ...
            fitData, baseFeatures, featureInfoLocal, customExpr);

        if isstruct(extractionInfo) && isfield(extractionInfo, 'warnings') && ~isempty(extractionInfo.warnings)
            previewWarnings = extractionInfo.warnings(1:min(3, numel(extractionInfo.warnings)));
            for w = 1:numel(previewWarnings)
                emitProgress(progressCb, 'Feature warning (%s): %s', imageName, previewWarnings{w});
            end
            if numel(extractionInfo.warnings) > numel(previewWarnings)
                emitProgress(progressCb, 'Feature warning (%s): ... and %d more', ...
                    imageName, numel(extractionInfo.warnings) - numel(previewWarnings));
            end
        end

        if isempty(classifierFeatureNames)
            classifierFeatureNames = featureNames;
        end
        [predictions, ~, ~, ~] = snap_helpers.classification.applyClassifier( ...
            model, X, featureNames, classifierFeatureNames, normParams);
        predictions = double(predictions(:));

        validEvalMask = logical(validMask(:)) & isfinite(predictions);
        yTrue = labels(validEvalMask);
        yPred = predictions(validEvalMask);
        m = computeMetrics(yTrue, yPred);

        rows(i).nEvaluated = m.nSamples;
        rows(i).tp = m.tp;
        rows(i).fp = m.fp;
        rows(i).fn = m.fn;
        rows(i).tn = m.tn;
        rows(i).precisionReal = m.precisionReal;
        rows(i).recallReal = m.recallReal;
        rows(i).f1Real = m.f1Real;
        rows(i).accuracy = m.accuracy;
        rows(i).nInvalidFeatures = nCandidates - m.nSamples;

        if isa(rowUpdateCb, 'function_handle')
            try
                rowUpdateCb(i, rows(i), nPairs);
            catch ME
                emitProgress(progressCb, 'Row update callback warning: %s', ME.message);
            end
        end

        totalEvaluated = totalEvaluated + m.nSamples;
        TP = TP + m.tp;
        FP = FP + m.fp;
        FN = FN + m.fn;
        TN = TN + m.tn;

        if m.nSamples > 0
            emitProgress(progressCb, ...
                'Pair %d/%d metrics: n=%d, acc=%.4f, prec=%.4f, rec=%.4f, f1=%.4f', ...
                i, nPairs, m.nSamples, m.accuracy, m.precisionReal, m.recallReal, m.f1Real);
        else
            emitProgress(progressCb, 'Pair %d/%d: no valid feature rows for evaluation.', i, nPairs);
        end
    end

    precisionReal = TP / max(TP + FP, 1);
    recallReal = TP / max(TP + FN, 1);
    f1Real = 2 * precisionReal * recallReal / max(precisionReal + recallReal, 1e-6);
    accuracy = (TP + TN) / max(totalEvaluated, 1);

    report = struct();
    report.datasetDir = opts.datasetDir;
    report.parameterFile = opts.parameterFile;
    report.classifierFile = opts.classifierFile;
    report.channelIndex = opts.channelIndex;
    report.matchDistance = opts.matchDistance;
    report.convertFijiCoords = logical(opts.convertFijiCoords);
    report.numPairs = nPairs;
    report.totalCandidates = totalCandidates;
    report.totalEvaluated = totalEvaluated;
    report.totalReal = totalReal;
    report.totalNoise = totalNoise;
    report.tp = TP;
    report.fp = FP;
    report.fn = FN;
    report.tn = TN;
    report.precisionReal = precisionReal;
    report.recallReal = recallReal;
    report.f1Real = f1Real;
    report.accuracy = accuracy;
    report.elapsedSec = toc(tEval);

    rowData = rowsToTableData(rows);
end

function out = reportToSummaryLines(report)
    out = { ...
        sprintf('Dataset: %s', report.datasetDir), ...
        sprintf('Pairs evaluated: %d', report.numPairs), ...
        sprintf('Candidates total: %d | Evaluated: %d', report.totalCandidates, report.totalEvaluated), ...
        sprintf('Labels total (after matching): real=%d, noise=%d', report.totalReal, report.totalNoise), ...
        sprintf('Confusion: TP=%d, FP=%d, FN=%d, TN=%d', report.tp, report.fp, report.fn, report.tn), ...
        sprintf('Accuracy=%.4f | Precision(real)=%.4f | Recall(real)=%.4f | F1(real)=%.4f', ...
            report.accuracy, report.precisionReal, report.recallReal, report.f1Real), ...
        sprintf('Channel=%d | MatchDistance=%.3g vox | ConvertFIJI=%d | Runtime=%.2f s', ...
            report.channelIndex, report.matchDistance, report.convertFijiCoords, report.elapsedSec)};
end

function rowData = rowsToTableData(rows)
    n = numel(rows);
    rowData = cell(n, 13);
    for i = 1:n
        rowData(i, :) = rowStructToTableRow(rows(i));
    end
end

function rowData = rowStructToTableRow(row)
    rowData = cell(1, 13);
    [~, name, ext] = fileparts(row.imagePath);
    rowData{1,1} = [name, ext];
    rowData{1,2} = row.nCandidates;
    rowData{1,3} = row.nEvaluated;
    rowData{1,4} = row.nReal;
    rowData{1,5} = row.nNoise;
    rowData{1,6} = row.tp;
    rowData{1,7} = row.fp;
    rowData{1,8} = row.fn;
    rowData{1,9} = row.tn;
    rowData{1,10} = row.precisionReal;
    rowData{1,11} = row.recallReal;
    rowData{1,12} = row.f1Real;
    rowData{1,13} = row.accuracy;
end

function m = computeMetrics(yTrue, yPred)
    m = struct('nSamples', 0, 'tp', 0, 'fp', 0, 'fn', 0, 'tn', 0, ...
        'precisionReal', NaN, 'recallReal', NaN, 'f1Real', NaN, 'accuracy', NaN);

    if isempty(yTrue) || isempty(yPred)
        return;
    end

    yTrue = double(yTrue(:));
    yPred = double(yPred(:));
    keep = isfinite(yTrue) & isfinite(yPred);
    yTrue = yTrue(keep);
    yPred = yPred(keep);
    if isempty(yTrue)
        return;
    end

    yTrue(yTrue ~= 0) = 1;
    yPred(yPred ~= 0) = 1;

    cm = confusionmat(yTrue, yPred, 'Order', [0 1]);
    m.tn = cm(1, 1);
    m.fp = cm(1, 2);
    m.fn = cm(2, 1);
    m.tp = cm(2, 2);
    m.nSamples = numel(yTrue);

    m.precisionReal = m.tp / max(m.tp + m.fp, 1);
    m.recallReal = m.tp / max(m.tp + m.fn, 1);
    m.f1Real = 2 * m.precisionReal * m.recallReal / max(m.precisionReal + m.recallReal, 1e-6);
    m.accuracy = mean(yPred == yTrue);
end

function [imageFiles, labelFiles] = discoverImageCsvPairs(rootDir, progressCb)
    if nargin < 2
        progressCb = [];
    end
    if ~isfolder(rootDir)
        error('Directory not found: %s', rootDir);
    end

    imageFiles = [ ...
        collectFilesRecursive(rootDir, '.tif', progressCb, 'images (.tif)'); ...
        collectFilesRecursive(rootDir, '.tiff', progressCb, 'images (.tiff)') ...
    ];
    if isempty(imageFiles)
        labelFiles = {};
        return;
    end

    csvFiles = collectFilesRecursive(rootDir, '.csv', progressCb, 'label files (.csv)');
    validCsvMask = cellfun(@isLabelCsvFile, csvFiles);
    csvFiles = csvFiles(validCsvMask);
    if isempty(csvFiles)
        imageFiles = {};
        labelFiles = {};
        return;
    end

    matchedImages = {};
    matchedLabels = {};
    for i = 1:numel(imageFiles)
        imgPath = imageFiles{i};
        labPath = findExactCsvForImageInList(imgPath, csvFiles);
        if isempty(labPath)
            continue;
        end
        matchedImages{end+1,1} = imgPath; %#ok<AGROW>
        matchedLabels{end+1,1} = labPath; %#ok<AGROW>
    end

    imageFiles = matchedImages;
    labelFiles = matchedLabels;
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
    if globalExactCount == 1
        labelPath = globalExactPath;
        return;
    end
    if globalNoOmeCount == 1
        labelPath = globalNoOmePath;
    end
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

        if ~isempty(progressCb) && mod(scannedDirs, 50) == 0
            emitProgress(progressCb, ...
                'Scanning %s: %d folder(s), %d file(s), %d match(es).', ...
                labelText, scannedDirs, scannedFiles, numel(files));
        end
    end

    if ~isempty(files)
        files = sort(files);
    end
end

function fitData = buildFitDataFromImage(imagePath, params, channelIdx, progressCb)
    if nargin < 4
        progressCb = [];
    end
    if nargin < 2 || isempty(params) || ~isstruct(params)
        error('Image-based evaluation requires a valid parameter struct.');
    end
    if nargin < 3 || isempty(channelIdx)
        channelIdx = 1;
    end
    channelIdx = max(1, round(channelIdx));

    emitProgress(progressCb, 'Generating candidates from image (ch=%d): %s', channelIdx, imagePath);
    rawImage = loadImageVolumeForTraining(imagePath);
    handles = createTrainingHandlesFromParams(params, channelIdx);

    maximaEnabled = logical(getChannelParamValue(params, 'maximaEnabled', channelIdx, true));
    if ~maximaEnabled
        emitProgress(progressCb, 'maximaEnabled is false for channel %d; no candidates generated.', channelIdx);
        fitData = struct([]);
        return;
    end

    fitEnabled = logical(getChannelParamValue(params, 'gaussFitEnabled', channelIdx, true));
    if ~fitEnabled
        error('gaussFitEnabled is false for channel %d. Cannot build classifier features.', channelIdx);
    end

    pipelineContext = struct();
    pipelineContext.mode = 'test';
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
    fitData = pipelineResult.fitResults;
    fitData = normalizeFitDataLocal(fitData);

    emitProgress(progressCb, 'Generated %d candidate(s) from image.', numel(fitData));
end

function img = loadImageVolumeForTraining(imagePath)
    try
        img = double(tiffreadVolume(imagePath));
    catch
        try
            img = double(imread(imagePath));
        catch ME
            error('Failed to load image file "%s": %s', imagePath, ME.message);
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

function n = inferNumChannelsFromParameters(params, fallback)
    if ~isstruct(params)
        n = fallback;
        return;
    end

    explicitN = NaN;
    if isfield(params, 'numChannels')
        explicitN = normalizeScalarNumeric(params.numChannels, NaN);
    elseif isfield(params, 'numChan')
        explicitN = normalizeScalarNumeric(params.numChan, NaN);
    elseif isfield(params, 'numChanDrop')
        explicitN = normalizeScalarNumeric(params.numChanDrop, NaN);
    elseif isfield(params, 'workflowConfig') && isstruct(params.workflowConfig) && isfield(params.workflowConfig, 'numChannels')
        explicitN = normalizeScalarNumeric(params.workflowConfig.numChannels, NaN);
    end

    if isfinite(explicitN) && explicitN >= 1
        n = max(1, round(explicitN));
        return;
    end

    n = max(1, round(fallback));
    hints = {'gaussFitMethod', 'maximaMode', 'preProcMode', ...
        'channelPath', 'channelPaths', 'classifyClassifierPath', 'classifyEnabled'};
    for i = 1:numel(hints)
        fieldName = hints{i};
        if isfield(params, fieldName)
            n = max(n, inferChannelCountFromValue(params.(fieldName)));
        end
    end
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

function val = normalizeScalarNumeric(raw, defaultValue)
    val = defaultValue;
    if nargin < 2
        defaultValue = NaN;
        val = defaultValue;
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
    if iscell(raw) && ~isempty(raw)
        raw = raw{1};
    end
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
    if isempty(fitData)
        fitData = struct([]);
        return;
    end

    if istable(fitData)
        fitData = table2struct(fitData);
    end

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
                center = [NaN NaN NaN];
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
                labels = NaN(numel(fitData), 1);
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

    labels = NaN(numel(fitData), 1);

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
    if any(~isfinite(values))
        error('Coordinate column "%s" contains non-numeric values.', coordName);
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

    y = NaN(numel(raw), 1);
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

function xyz = forceXYZ(coords)
    xyz = coords(:)';
    if numel(xyz) >= 3
        xyz = xyz(1:3);
    elseif numel(xyz) == 2
        xyz = [xyz, 1];
    elseif numel(xyz) == 1
        xyz = [xyz, 1, 1];
    else
        xyz = [NaN, NaN, NaN];
    end
end

function method = inferFittingMethod(fitData)
    method = '3D Gaussian';
    if isempty(fitData)
        return;
    end
    if isfield(fitData, 'fitMethod') && ~isempty(fitData(1).fitMethod)
        method = string(fitData(1).fitMethod);
        method = char(method);
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

function emitProgress(progressCb, msgFmt, varargin)
    if isempty(progressCb) || ~isa(progressCb, 'function_handle')
        return;
    end
    msg = formatMessage(msgFmt, varargin{:});
    try
        progressCb(msg);
    catch
        % Progress callbacks should never halt execution.
    end
end

function msg = formatMessage(msgFmt, varargin)
    if isempty(varargin)
        msg = char(string(msgFmt));
        return;
    end
    try
        msg = sprintf(msgFmt, varargin{:});
    catch
        msg = char(string(msgFmt));
    end
end

function [model, features, featureInfo, fittingMethod, customExpressions, normParams] = ...
        loadClassifierForTesting(classifierPath, progressCb)
    model = [];
    features = {};
    featureInfo = struct();
    fittingMethod = '';
    customExpressions = struct('name', {}, 'expression', {});
    normParams = struct('mu', [], 'sigma', [], 'standardized', false);

    classifierPath = sanitizePathInput(classifierPath);
    if exist(classifierPath, 'file') ~= 2
        error('Classifier file not found: %s', classifierPath);
    end

    % Primary loader (shared SNAP utility).
    warnState = warning('off', 'all');
    warnCleanup = onCleanup(@() warning(warnState)); %#ok<NASGU>
    lastwarn('');
    [model, features, featureInfo, trainStats, fittingMethod, ~, success, customExpressions, normParams] = ...
        snap_helpers.classification.loadClassifier(classifierPath);
    [warnMsg, ~] = lastwarn;
    if success
        if isempty(normParams) || ~isstruct(normParams)
            normParams = inferNormParamsFromTrainStats(trainStats);
        end
        return;
    end
    if ~isempty(warnMsg)
        emitProgress(progressCb, 'Classifier load warning: %s', warnMsg);
    end

    emitProgress(progressCb, 'Primary classifier loader failed; attempting fallback load path.');
    attempts = {@() load(classifierPath), @() load(classifierPath, '-mat'), @() loadClassifierViaMatfile(classifierPath)};
    lastErr = '';
    for a = 1:numel(attempts)
        try
            data = attempts{a}();
            [model, features, featureInfo, fittingMethod, customExpressions, normParams] = ...
                unpackClassifierStruct(data);
            emitProgress(progressCb, 'Fallback classifier load succeeded (attempt %d).', a);
            return;
        catch ME
            lastErr = ME.message;
            emitProgress(progressCb, 'Fallback load attempt %d failed: %s', a, ME.message);
            pause(0.2 * a);
        end
    end

    error('Failed to load classifier file "%s". Last error: %s', classifierPath, lastErr);
end

function data = loadClassifierViaMatfile(classifierPath)
    mf = matfile(classifierPath);
    vars = who(mf);
    vars = string(vars);
    wanted = ["model", "features", "featureInfo", "fittingMethod", ...
              "customExpressions", "normParams", "trainStats"];
    present = wanted(ismember(wanted, vars));
    data = struct();
    for i = 1:numel(present)
        key = char(present(i));
        data.(key) = mf.(key);
    end
end

function [model, features, featureInfo, fittingMethod, customExpressions, normParams] = ...
        unpackClassifierStruct(data)
    model = [];
    features = {};
    featureInfo = struct();
    fittingMethod = '';
    customExpressions = struct('name', {}, 'expression', {});
    normParams = struct('mu', [], 'sigma', [], 'standardized', false);

    if ~isstruct(data)
        error('Classifier file did not load as a struct.');
    end
    if ~isfield(data, 'model') || isempty(data.model)
        error('Classifier file is missing "model".');
    end
    if ~isfield(data, 'features') || isempty(data.features)
        error('Classifier file is missing "features".');
    end

    model = data.model;
    features = normalizeFeatureNameList(data.features);

    if isfield(data, 'featureInfo') && isstruct(data.featureInfo)
        featureInfo = data.featureInfo;
    end
    if isfield(data, 'fittingMethod') && ~isempty(data.fittingMethod)
        fittingMethod = char(string(data.fittingMethod));
    end
    if isfield(data, 'customExpressions')
        customExpressions = normalizeCustomExpressionListLocal(data.customExpressions);
    end
    if isfield(data, 'normParams') && isstruct(data.normParams)
        normParams = data.normParams;
    elseif isfield(data, 'trainStats') && isstruct(data.trainStats)
        normParams = inferNormParamsFromTrainStats(data.trainStats);
    end
end

function normParams = inferNormParamsFromTrainStats(trainStats)
    normParams = struct('mu', [], 'sigma', [], 'standardized', false);
    if ~isstruct(trainStats)
        return;
    end
    if isfield(trainStats, 'featureMeans') && isfield(trainStats, 'featureStds') && ...
            ~isempty(trainStats.featureMeans) && ~isempty(trainStats.featureStds)
        normParams.mu = trainStats.featureMeans;
        normParams.sigma = trainStats.featureStds;
        normParams.standardized = true;
    end
end

function features = normalizeFeatureNameList(rawFeatures)
    if isempty(rawFeatures)
        features = {};
        return;
    end
    if ischar(rawFeatures)
        features = {rawFeatures};
    elseif isstring(rawFeatures)
        features = cellstr(rawFeatures(:))';
    elseif iscell(rawFeatures)
        features = cellfun(@(x) char(string(x)), rawFeatures(:)', 'UniformOutput', false);
    else
        error('Classifier feature list is invalid. Expected char/string/cell.');
    end
end

function ce = normalizeCustomExpressionListLocal(rawCustom)
    ce = struct('name', {}, 'expression', {});
    if isempty(rawCustom) || ~isstruct(rawCustom)
        return;
    end
    rawCustom = rawCustom(:);
    keep = false(size(rawCustom));
    for i = 1:numel(rawCustom)
        if ~isfield(rawCustom(i), 'name') || ~isfield(rawCustom(i), 'expression')
            continue;
        end
        rawCustom(i).name = strtrim(char(string(rawCustom(i).name)));
        rawCustom(i).expression = strtrim(char(string(rawCustom(i).expression)));
        keep(i) = ~isempty(rawCustom(i).name) && ~isempty(rawCustom(i).expression);
    end
    ce = rawCustom(keep);
end
