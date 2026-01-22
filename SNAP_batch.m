function SNAP_batch(varargin)
% SNAP_batch - Batch processing for SNAP (Spot & Nuclei Analysis Pipeline)
%
% GUI MODE (No arguments):
%   SNAP_batch
%   Opens a graphical interface for batch processing configuration
%
% COMMAND-LINE MODE:
%   SNAP_batch(inputDir, paramFile)
%   SNAP_batch(inputDir, paramFile, 'OutputDir', outputDir)
%   SNAP_batch(inputDir, paramFile)
%
% Processes multiple image sets using parameters exported from SNAP GUI.
% Expects a folder containing subfolders, where each subfolder has:
%   - Nuclei image (e.g., nuclei.tif, DAPI.tif, *nuclei*.tif)
%   - Channel images (e.g., channel1.tif, channel2.tif, ch1.tif, ch2.tif)
%   - Optional: DIC image (e.g., dic.tif, brightfield.tif)
%
% NOTE: Fully compatible with SNAP_prepare output structure.
%       SNAP_prepare creates folders with: nuclei.tif, dic.tif, channel1.tif, channel2.tif, etc.
%
% Inputs (Command-line mode):
%   inputDir  - Path to folder containing subfolders with image sets
%   paramFile - Path to exported parameter file from SNAP GUI
%
% Optional Name-Value Pairs (Command-line mode):
%   'OutputDir'     - Output directory (default: inputDir/batch_results_TIMESTAMP)
%                     Note: Required when using GUI mode
%   'ExportFormat'  - Image export format: 'TIFF', 'PNG', 'JPEG' (default: 'TIFF')
%
% Examples:
%   SNAP_batch  % Launch GUI
%   SNAP_batch('experiments/dataset1/', 'optimal_params.mat')
%   SNAP_batch('data/', 'params.mat', 'OutputDir', 'results/')
%   SNAP_batch('data/', 'params.mat', 'ExportFormat', 'PNG')

% Check if GUI mode (no arguments) or command-line mode
if nargin == 0
    % Launch GUI
    createBatchGUI();
    return;
end

% Command-line mode - extract required arguments
inputDir = varargin{1};
paramFile = varargin{2};
extraArgs = varargin(3:end);

% Extract progress handle early for debug logging
progressHandle = [];
for i = 1:2:length(extraArgs)
    if strcmp(extraArgs{i}, 'ProgressHandle')
        progressHandle = extraArgs{i+1};
        break;
    end
end


%% Parse inputs
p = inputParser;
addRequired(p, 'inputDir', @ischar);
addRequired(p, 'paramFile', @ischar);
addParameter(p, 'OutputDir', '', @ischar);
addParameter(p, 'ExportFormat', 'TIFF', @ischar);
addParameter(p, 'ExportVisualizations', false, @islogical); % Export annotated PNG images
addParameter(p, 'ProgressHandle', [], @(x) isempty(x) || isstruct(x)); % For GUI progress updates
addParameter(p, 'ExportOptions', struct(), @isstruct); % Export options from GUI
addParameter(p, 'Classifiers', cell(1, 10), @iscell); % Classifiers per channel
addParameter(p, 'ClassifierFeatures', cell(1, 10), @iscell); % Features per classifier
addParameter(p, 'ClassifierFeatureInfo', cell(1, 10), @iscell); % Feature info per classifier
addParameter(p, 'ClassifierCustomExpressions', cell(1, 10), @iscell); % Custom expressions per classifier
addParameter(p, 'ClassifierNormParams', cell(1, 10), @iscell); % Normalization params per classifier
parse(p, inputDir, paramFile, extraArgs{:});

inputDir = p.Results.inputDir;
paramFile = p.Results.paramFile;
exportFormat = p.Results.ExportFormat;
progressHandle = p.Results.ProgressHandle;
exportOptions = p.Results.ExportOptions;
classifiers = p.Results.Classifiers;
classifierFeatures = p.Results.ClassifierFeatures;
classifierFeatureInfo = p.Results.ClassifierFeatureInfo;
classifierCustomExpressions = p.Results.ClassifierCustomExpressions;
classifierNormParams = p.Results.ClassifierNormParams;

% Set default export options if not provided or empty
if isempty(exportOptions) || ~isstruct(exportOptions) || isempty(fieldnames(exportOptions))
    exportOptions.nucleiData = true;
    exportOptions.channelData = true;
    exportOptions.clusteredData = true;
    exportOptions.visualizations = true;
end

% Initialize struct to track what was actually exported per subfolder
batchResults.exportedData = struct();

% Ensure proper path for package access
current_dir = fileparts(mfilename('fullpath'));
if ~contains(path, current_dir)
    addpath(current_dir);
end

%% Validate inputs
if ~exist(inputDir, 'dir')
    error('Input directory does not exist: %s', inputDir);
end

if ~exist(paramFile, 'file')
    error('Parameter file does not exist: %s', paramFile);
end

%% Setup output directory
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isempty(p.Results.OutputDir)
    outputDir = fullfile(inputDir, sprintf('batch_results_%s', timestamp));
else
    outputDir = p.Results.OutputDir;
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Silent mode - all output goes to GUI Progress Log

% Log progress
if ~isempty(progressHandle)
    appendLog(progressHandle, '→ Loading parameters...');
    drawnow;
end

%% Load parameters
try
    paramData = load(paramFile);
    if isfield(paramData, 'batchConfig')
        params = paramData.batchConfig.parameters;
    elseif isfield(paramData, 'paramData')
        params = paramData.paramData.parameters;
    elseif isfield(paramData, 'lastUsed')
        params = paramData.lastUsed;
    else
        error('Could not find parameters in file. Expected batchConfig, paramData, or lastUsed structure.');
    end
    % Parameters loaded successfully
    if ~isempty(progressHandle)
        appendLog(progressHandle, '✓ Parameters loaded');
        drawnow;
    end
    
    % Add classifiers to params if provided
    if any(~cellfun(@isempty, classifiers))
        params.classifiers = classifiers;
        params.classifierFeatures = classifierFeatures;
        params.classifierFeatureInfo = classifierFeatureInfo;
        params.classifierCustomExpressions = classifierCustomExpressions;
        params.classifierNormParams = classifierNormParams;
        
        % Enable classification for channels that have classifiers
        if ~isfield(params, 'classifyEnabled')
            params.classifyEnabled = cell(1, 10);
            params.classifyFilterNoise = cell(1, 10);
            for c = 1:10
                params.classifyEnabled{c} = false;
                params.classifyFilterNoise{c} = true;
            end
        end
        for c = 1:min(numel(classifiers), 10)
            if ~isempty(classifiers{c})
                params.classifyEnabled{c} = true;
                params.classifyFilterNoise{c} = true;
                if ~isempty(progressHandle)
                    appendLog(progressHandle, sprintf('  → Classifier loaded for Channel %d', c));
                end
            end
        end
    end
catch ME
    error('Failed to load parameter file: %s', ME.message);
end

%% Find all subfolders (each represents one image set)
subfolderInfo = dir(inputDir);
subfolders = subfolderInfo([subfolderInfo.isdir] & ~ismember({subfolderInfo.name}, {'.', '..'}));
numSets = length(subfolders);

if numSets == 0
    error('No subfolders found in input directory: %s', inputDir);
end

if ~isempty(progressHandle)
    appendLog(progressHandle, sprintf('→ Found %d image sets to process', numSets));
    drawnow;
end

%% Initialize batch tracking
batchResults = struct();
batchResults.timestamp = timestamp;
batchResults.inputDir = inputDir;
batchResults.outputDir = outputDir;
batchResults.paramFile = paramFile;
batchResults.numSets = numSets;
batchResults.results = cell(numSets, 1);

startTime = tic;

%% Process each subfolder

for setIdx = 1:numSets
    setName = subfolders(setIdx).name;
    setPath = fullfile(inputDir, setName);
    
    % Update progress bar and log start
    if ~isempty(progressHandle)
        updateProgress(progressHandle, setIdx, numSets, sprintf('Processing: %s', setName));
        appendLog(progressHandle, sprintf('[%d/%d] %s', setIdx, numSets, setName));
        drawnow;
    end
    
    setStartTime = tic;
    
    try
        % Find images in this subfolder
        imageSet = findImagesInFolder(setPath, params);
        
        if isempty(imageSet.nuclei)
            % Create consistent structure for skipped images
            skippedResult = createEmptyResult(setName, params.numChannels);
            skippedResult.status = 'skipped';
            skippedResult.reason = 'No nuclei image';
            batchResults.results{setIdx} = skippedResult;
            
            if ~isempty(progressHandle)
                appendLog(progressHandle, '  ⚠ Skipped: No nuclei image found');
            end
            continue;
        end
        
        % Log detected files
        if ~isempty(progressHandle)
            numFound = sum(~cellfun(@isempty, imageSet.channels));
            appendLog(progressHandle, sprintf('  Files found: Nuclei + %d channels', numFound));
        end
        
        % Process this image set
        try
            results = processImageSet(imageSet, params, setName, progressHandle);
            
            % Export results for this set
            setOutputDir = fullfile(outputDir, setName);
            if ~exist(setOutputDir, 'dir')
                mkdir(setOutputDir);
            end
            % Export results for this set (conditionally based on available data and export options)
            exportedTypes = exportSetResults(results, imageSet, setOutputDir, setName, timestamp, exportOptions, params, progressHandle);
            
            % Track results (fields already initialized in processImageSet)
            setTime = toc(setStartTime);
            results.processingTime = setTime;
            results.status = 'success';
            results.exportedTypes = exportedTypes;
            
            batchResults.results{setIdx} = results;
            
            % Log completion with timing
            if ~isempty(progressHandle)
                appendLog(progressHandle, sprintf('  ✓ Complete: %d nuclei, %d spots (%.1fs)', ...
                    results.numNuclei, results.totalSpots, setTime));
            end
            
        catch processError
            % Inner error - processing failed
            setTime = toc(setStartTime);
            errorResult = createEmptyResult(setName, params.numChannels);
            errorResult.status = 'error';
            errorResult.error_message = processError.message;
            errorResult.processingTime = setTime;
            batchResults.results{setIdx} = errorResult;
            
            if ~isempty(progressHandle)
                appendLog(progressHandle, sprintf('  ✗ Error in %s: %s', setName, processError.message));
            end
        end
        
    catch ME
        errorResult = createEmptyResult(setName, params.numChannels);
        errorResult.status = 'error';
        errorResult.error_message = ME.message;
        errorResult.processingTime = 0;
        batchResults.results{setIdx} = errorResult;
        
        if ~isempty(progressHandle)
            appendLog(progressHandle, sprintf('  ✗ Error in %s: %s', setName, ME.message));
        end
    end
end  % for

totalTime = toc(startTime);

%% Export aggregate results if requested
% Create aggregate for any checked export options
if exportOptions.nucleiData || exportOptions.channelData || exportOptions.clusteredData
    if ~isempty(progressHandle)
        appendLog(progressHandle, '');
        appendLog(progressHandle, 'Creating aggregate data...');
    end
    exportAggregateResults(batchResults, outputDir, timestamp, exportOptions, params);
    if ~isempty(progressHandle)
        appendLog(progressHandle, '✓ Aggregate data complete');
    end
end

% Batch complete - all output in GUI Progress Log

% Log summary to GUI if available
if ~isempty(progressHandle)
    numSuccess = sum(cellfun(@(x) isfield(x, 'status') && strcmp(x.status, 'success'), batchResults.results));
    numFailed = sum(cellfun(@(x) isfield(x, 'status') && strcmp(x.status, 'error'), batchResults.results));
    numSkipped = sum(cellfun(@(x) isfield(x, 'status') && strcmp(x.status, 'skipped'), batchResults.results));
    
    appendLog(progressHandle, '');
    appendLog(progressHandle, '═════════════════════════════════════');
    appendLog(progressHandle, 'BATCH PROCESSING COMPLETE');
    appendLog(progressHandle, '═════════════════════════════════════');
    appendLog(progressHandle, sprintf('Total Time: %.1f minutes', totalTime/60));
    appendLog(progressHandle, sprintf('Success: %d/%d sets', numSuccess, numSets));
    appendLog(progressHandle, sprintf('Failed: %d sets', numFailed));
    appendLog(progressHandle, sprintf('Skipped: %d sets', numSkipped));
    appendLog(progressHandle, sprintf('Results: %s', outputDir));
    appendLog(progressHandle, '═════════════════════════════════════');
end

end

%% Helper Functions

function imageSet = findImagesInFolder(folderPath, params)
% Find nuclei, channel, and DIC images in a folder
% Compatible with SNAP_prepare output format
    imageSet = struct();
    imageSet.folderPath = folderPath;
    imageSet.nuclei = '';
    imageSet.dic = '';
    imageSet.channels = {};
    
    % Get all TIFF files
    tifFiles = dir(fullfile(folderPath, '*.tif*'));
    
    if isempty(tifFiles)
        warning('No TIFF files found in: %s', folderPath);
        return;
    end
    
    % Pattern matching for nuclei (SNAP_prepare creates: nuclei.tif)
    nucleiPatterns = {'nuclei', 'dapi', 'hoechst', 'nuc'};
    for i = 1:length(tifFiles)
        fileName = lower(tifFiles(i).name);
        for pattern = nucleiPatterns
            if contains(fileName, pattern{1}) && ~contains(fileName, 'dic')
                imageSet.nuclei = fullfile(folderPath, tifFiles(i).name);
                break;
            end
        end
        if ~isempty(imageSet.nuclei), break; end
    end
    
    % Pattern matching for DIC (SNAP_prepare creates: dic.tif)
    for i = 1:length(tifFiles)
        fileName = lower(tifFiles(i).name);
        if contains(fileName, 'dic') || contains(fileName, 'brightfield')
            imageSet.dic = fullfile(folderPath, tifFiles(i).name);
            break;
        end
    end
    
    % Pattern matching for channels (SNAP_prepare creates: channel1.tif, channel2.tif, etc.)
    numChannels = params.numChannels;
    
    for ch = 1:numChannels
        for i = 1:length(tifFiles)
            fileName = lower(tifFiles(i).name);
            fullFilePath = fullfile(folderPath, tifFiles(i).name);
            
            % Skip if already identified as nuclei or DIC
            if strcmp(fullFilePath, imageSet.nuclei) || strcmp(fullFilePath, imageSet.dic)
                continue;
            end
            
            % Check for channel patterns
            % Priority order: channel1, ch1, c01, c1
            if contains(fileName, sprintf('channel%d', ch)) || ...
               contains(fileName, sprintf('ch%d', ch)) || ...
               contains(fileName, sprintf('c0%d', ch)) || ...
               (contains(fileName, sprintf('c%d', ch)) && ~contains(fileName, 'nuclei'))
                imageSet.channels{ch} = fullfile(folderPath, tifFiles(i).name);
                break;
            end
        end
        
        % Note if channel not found (silent)
    end
    
    % Files detected (silent)
end

function results = processImageSet(imageSet, params, setName, progressHandle)
% Process a single image set using saved parameters
    if nargin < 4
        progressHandle = [];
    end
    
    if ~isempty(progressHandle)
        appendLog(progressHandle, sprintf('  → Loading %s...', setName));
    end
    
    % Get number of channels first (needed for initialization)
    numChannels = params.numChannels;
    
    % Initialize results structure with ALL fields for consistency
    % Include ALL fields that will be added later to ensure same structure across all image sets
    results = struct();
    results.name = setName;
    results.numNuclei = 0;
    results.totalSpots = 0;
    results.channelSpots = zeros(1, numChannels);
    results.nucleiMask = [];
    results.nucleusLabels = struct();
    results.channels = cell(1, numChannels);
    results.signalComposition = [];
    results.processingTime = 0;  % Will be updated later
    results.status = '';  % Will be updated later
    results.exportedTypes = struct('nucleiData', false, 'channelData', false, 'clusteredData', false, 'visualizations', false);  % Will be updated later
    
    % Load images
    if ~isempty(progressHandle)
        appendLog(progressHandle, '  → Loading images...');
        drawnow;
    end
    
    rawNuclei = double(tiffreadVolume(imageSet.nuclei));
    rawChannels = cell(1, numChannels);
    for ch = 1:numChannels
        if ch <= length(imageSet.channels) && ~isempty(imageSet.channels{ch})
            rawChannels{ch} = double(tiffreadVolume(imageSet.channels{ch}));
        end
    end
    
    if ~isempty(progressHandle)
        appendLog(progressHandle, '    ✓ Images loaded');
        drawnow;
    end
    
    % Process nuclei
    if ~isempty(progressHandle)
        appendLog(progressHandle, '  → Processing nuclei (segmentation)...');
        drawnow;
    end
    
    nucStart = tic;
    [processedNuclei, nucleiMask, nucleusLabels] = processNucleiBatch(rawNuclei, params);
    nucTime = toc(nucStart);
    
    if ~isempty(progressHandle)
        appendLog(progressHandle, sprintf('    ✓ Found %d nuclei (%.2fs)', nucleusLabels.num_nuclei, nucTime));
        drawnow;
    end
    
    % Update nuclei results
    results.numNuclei = nucleusLabels.num_nuclei;
    results.nucleiMask = nucleiMask;
    results.nucleusLabels = nucleusLabels;
    
    % PERFORMANCE: Create handles structure once for all channels
    handles = createMinimalHandles(params);
    
    % Process channels (results.channels and channelSpots already initialized)
    totalSpots = 0;
    
    for ch = 1:numChannels
        if ~isempty(rawChannels{ch})
            chProcStart = tic;
            channelResults = processChannelBatch(rawChannels{ch}, params, ch, nucleiMask, nucleusLabels, progressHandle, handles);
            chProcTime = toc(chProcStart);
            
            if ~isempty(progressHandle)
                appendLog(progressHandle, sprintf('  ✓ Channel %d: %d spots (%.1fs)', ch, channelResults.numSpots, chProcTime));
                drawnow;
            end
            
            results.channels{ch} = channelResults;
            totalSpots = totalSpots + channelResults.numSpots;
            results.channelSpots(ch) = channelResults.numSpots;
        else
            % Create empty channel result to maintain structure consistency
            results.channels{ch} = struct('channelIndex', ch, 'processedImage', [], ...
                'maximaCoords', zeros(0,3), 'fitResults', struct([]), 'numSpots', 0);
        end
    end
    
    results.totalSpots = totalSpots;
    
    % Analyze signal composition if nuclei and spots exist (signalComposition already initialized to [])
    if results.numNuclei > 0 && totalSpots > 0
        results.signalComposition = analyzeSignalCompositionBatch(results, params);
    end
    
    % Completion summary
    if ~isempty(progressHandle)
        appendLog(progressHandle, sprintf('  ✓ %s Complete: %d nuclei, %d total spots', ...
            setName, results.numNuclei, results.totalSpots));
        drawnow;
    end
end

function [processedNuclei, nucleiMask, nucleusLabels] = processNucleiBatch(rawNuclei, params)
% Process nuclei image using saved parameters
    
    % Create minimal handles structure for processing functions
    handles = createMinimalHandles(params);
    handles.rawNuclei = rawNuclei;
    
    % Apply preprocessing and background correction
    processedNuclei = snap_helpers.preprocessNucleiWithBgCorr(handles);
    
    % Segment nuclei
    [nucleiMask, nucleusLabels] = snap_helpers.segmentNuclei(processedNuclei, handles);
end

function channelResults = processChannelBatch(rawChannel, params, channelIdx, nucleiMask, nucleusLabels, progressHandle, handles)
% Process a single channel using saved parameters
    if nargin < 6
        progressHandle = [];
    end
    if nargin < 7
        % Fallback: create handles if not provided (for backwards compatibility)
        handles = createMinimalHandles(params);
    end
    
    channelResults = struct();
    channelResults.channelIndex = channelIdx;
    
    % Set raw channel data
    handles.rawChannel = {rawChannel};
    
    if ~isempty(progressHandle)
        appendLog(progressHandle, sprintf('      → Channel %d: Preprocessing & background correction...', channelIdx));
        drawnow;
    end
    
    % Process image (deconv, preproc, bg correction)
    procStart = tic;
    processedChannel = snap_helpers.processImage(handles, 1);
    procTime = toc(procStart);
    channelResults.processedImage = processedChannel;
    
    if ~isempty(progressHandle)
        appendLog(progressHandle, sprintf('        ✓ Preprocessing complete (%.2fs)', procTime));
        drawnow;
    end
    
    % Detect local maxima if enabled - USE SAME FUNCTION AS SNAP GUI
    maximaCoords = [];
    if params.maximaEnabled{channelIdx}
        if ~isempty(progressHandle)
            appendLog(progressHandle, sprintf('      → Channel %d: Detecting local maxima...', channelIdx));
            drawnow;
        end
        
        maxStart = tic;
        % Use shared maxima detection function (same as SNAP GUI)
        % Note: Use channelIdx (not 1) to access correct parameters from handles
        maximaCoords = snap_helpers.detectMaxima(processedChannel, handles, channelIdx);
        maxTime = toc(maxStart);
        
        if ~isempty(progressHandle)
            appendLog(progressHandle, sprintf('        ✓ Found %d spots (%.2fs)', size(maximaCoords,1), maxTime));
            drawnow;
        end
    end
    
    % Apply Gaussian fitting if enabled
    fitResults = [];
    if params.gaussFitEnabled{channelIdx} && ~isempty(maximaCoords)
        if ~isempty(progressHandle)
            appendLog(progressHandle, sprintf('      → Channel %d: Fitting Gaussians to %d spots...', channelIdx, size(maximaCoords,1)));
            drawnow;
        end
        
        fitParams = extractFitParams(params, channelIdx);
        is3D = size(processedChannel, 3) > 1;
        
        % PERFORMANCE: Disable any verbose output
        fitParams.gaussFitPlotCheck = false;
        fitParams.verbose = false;
        
        fitStart = tic;
        fitResults = snap_helpers.fitGaussians(rawChannel, maximaCoords, fitParams, is3D);
        fitTime = toc(fitStart);
        
        if ~isempty(progressHandle)
            appendLog(progressHandle, sprintf('        ✓ Gaussian fitting complete (%.2fs)', fitTime));
            drawnow;
        end
        
        % Apply fit filtering if enabled - USE SAME FUNCTION AS SNAP GUI
        if params.fitFilterEnabled{channelIdx}
            % Use shared fit filtering function (same as SNAP GUI)
            % This handles radial symmetry properly and matches GUI behavior exactly
            [fitResults, filterMask] = snap_helpers.applyFitFiltering(fitResults, channelIdx, handles);
            maximaCoords = maximaCoords(filterMask, :);
            
            % Verify synchronization
            if length(fitResults) ~= size(maximaCoords, 1)
                warning('Channel %d: Fit filtering sync issue (%d fits vs %d coords)', ...
                    channelIdx, length(fitResults), size(maximaCoords, 1));
            end
        end
        
        % Apply SVM classification filtering if enabled
        if isfield(params, 'classifyEnabled') && length(params.classifyEnabled) >= channelIdx && ...
           params.classifyEnabled{channelIdx} && ~isempty(fitResults)
            
            % Check if classifier is loaded
            if isfield(params, 'classifiers') && length(params.classifiers) >= channelIdx && ...
               ~isempty(params.classifiers{channelIdx})
                
                classifier = params.classifiers{channelIdx};
                features = params.classifierFeatures{channelIdx};
                featureInfo = params.classifierFeatureInfo{channelIdx};
                
                % Get custom expressions
                customExpr = struct('name', {}, 'expression', {});
                if isfield(params, 'classifierCustomExpressions') && ...
                   numel(params.classifierCustomExpressions) >= channelIdx && ...
                   ~isempty(params.classifierCustomExpressions{channelIdx})
                    customExpr = params.classifierCustomExpressions{channelIdx};
                end
                
                % Get normalization parameters
                normParams = struct('mu', [], 'sigma', [], 'standardized', false);
                if isfield(params, 'classifierNormParams') && ...
                   numel(params.classifierNormParams) >= channelIdx && ...
                   ~isempty(params.classifierNormParams{channelIdx})
                    normParams = params.classifierNormParams{channelIdx};
                end
                
                % Build feature matrix with custom expressions
                [X, featureNames, validMask] = snap_helpers.classification.buildFeatureMatrix(...
                    fitResults, features, featureInfo, customExpr);
                
                % Apply classifier with normalization
                [predictions, ~, ~] = snap_helpers.classification.applyClassifier(...
                    classifier, X, featureNames, featureNames, normParams);
                
                % Filter out noise predictions if auto-filter is enabled
                if isfield(params, 'classifyFilterNoise') && length(params.classifyFilterNoise) >= channelIdx && ...
                   params.classifyFilterNoise{channelIdx}
                    
                    % Keep spots that are predicted Real OR have invalid features
                    keepMask = (predictions == 1) | ~validMask;
                    preClassifyCount = numel(fitResults);
                    
                    fitResults = fitResults(keepMask);
                    maximaCoords = maximaCoords(keepMask, :);
                    
                    if ~isempty(progressHandle)
                        nFiltered = preClassifyCount - numel(fitResults);
                        appendLog(progressHandle, sprintf('        → SVM: %d/%d classified as Real (%d filtered as Noise)', ...
                            numel(fitResults), preClassifyCount, nFiltered));
                    end
                end
            end
        end
    end
    
    % Apply nuclei inclusion/exclusion filtering if enabled - USE SAME FUNCTION AS SNAP GUI
    if params.nucInclusionExclusionEnabled && ~isempty(nucleiMask) && ~isempty(maximaCoords)
        % Check if this channel should be filtered based on "Apply to Channels" setting
        applyTo = params.nucInclusionExclusionApplyTo;
        shouldFilter = strcmp(applyTo, 'All Channels') || ...
                       strcmp(applyTo, sprintf('Channel %d', channelIdx));
        
        if shouldFilter
            mode = params.nucInclusionExclusionMode;
            preFilterCount = size(maximaCoords, 1);
            
            % Use shared nuclei filtering function (returns mask to maintain synchronization)
            [maximaCoords, keepMask] = snap_helpers.filterMaximaByNuclei(maximaCoords, nucleiMask, mode);
            
            % CRITICAL: Apply the same mask to fitResults to maintain synchronization
            if ~isempty(fitResults) && length(fitResults) == preFilterCount
                fitResults = fitResults(keepMask);
            end
            
            if ~isempty(progressHandle) && size(maximaCoords, 1) < preFilterCount
                appendLog(progressHandle, sprintf('        → Nuclei filter: %d/%d maxima kept (%s mode)', ...
                    size(maximaCoords, 1), preFilterCount, mode));
            end
        end
    end

    % Ensure maximaCoords is always a double array (even if empty) for consistency
    if isempty(maximaCoords)
        channelResults.maximaCoords = zeros(0, 3);  % Empty Nx3 array
    else
        channelResults.maximaCoords = maximaCoords;
    end
    
    % Ensure fitResults is always a struct array (even if empty) for consistency
    if isempty(fitResults)
        channelResults.fitResults = struct([]);
    else
        channelResults.fitResults = fitResults;
    end
    
    channelResults.numSpots = size(maximaCoords, 1);
end

% Old maxima detection functions removed - now using snap_helpers.detectMaxima()
% to ensure identical processing between SNAP GUI and SNAP_batch

% Old filtering functions removed - now using snap_helpers.filterMaximaByNuclei()
% and snap_helpers.applyFitFiltering() to ensure identical processing

function signalComp = analyzeSignalCompositionBatch(results, params)
% Analyze signal composition with FULL signal data embedded
% USES computeSignalMeasurements for ZERO REDUNDANCY
    
    signalComp = struct();
    numNuclei = results.numNuclei;
    numChannels = params.numChannels;
    labeledMask = results.nucleusLabels.labeled_mask;
    
    % Get active channels
    active_channels = [];
    for ch = 1:numChannels
        if ~isempty(results.channels{ch}) && ~isempty(results.channels{ch}.maximaCoords)
            active_channels(end+1) = ch;
        end
    end

    signalComp.active_channels = active_channels;
    signalComp.num_nuclei = numNuclei;
    
    % Build nuclei_data structure with full signal data
    % Pre-allocate struct array to avoid structure mismatch errors
    nuclei_data = repmat(struct('nucleus_id', 0, 'centroid', [], 'channels', struct()), numNuclei, 1);
    
    % Initialize all channel fields for ALL nuclei first
    for nuc_idx = 1:numNuclei
        for ch = 1:numChannels
            ch_field = sprintf('channel_%d', ch);
            nuclei_data(nuc_idx).channels.(ch_field) = struct('signals', [], 'signal_count', 0);
        end
    end
    
    for nuc_idx = 1:numNuclei
        
        % Update the pre-allocated nucleus with data (don't create new struct)
        nuclei_data(nuc_idx).nucleus_id = nuc_idx;
        nuclei_data(nuc_idx).centroid = results.nucleusLabels.centroids_3d(nuc_idx, :);
        
        % Process each active channel (channels already pre-initialized above)
        for ch = active_channels
            % Get maxima and fits for this channel
            maximaCoords = results.channels{ch}.maximaCoords;
            fitResults = [];
            if isfield(results.channels{ch}, 'fitResults')
                fitResults = results.channels{ch}.fitResults;
            end
    
            % Find signals in this nucleus
            signals_found = [];
            for sig = 1:size(maximaCoords, 1)
                row = round(maximaCoords(sig, 1));
                col = round(maximaCoords(sig, 2));
                slice = round(maximaCoords(sig, 3));
        
                nucID = 0;
                if ndims(labeledMask) == 2
                    if row >= 1 && row <= size(labeledMask, 1) && col >= 1 && col <= size(labeledMask, 2)
                        nucID = labeledMask(row, col);
                    end
                else
                    if row >= 1 && row <= size(labeledMask, 1) && ...
                       col >= 1 && col <= size(labeledMask, 2) && ...
                       slice >= 1 && slice <= size(labeledMask, 3)
                        nucID = labeledMask(row, col, slice);
                    end
                end
        
                if nucID == nuc_idx
                    signals_found(end+1) = sig;
                end
            end
    
            % Use SHARED HELPER to get full signal data (ZERO REDUNDANCY!)
            % Update the pre-initialized channel field in the pre-allocated array
            ch_field = sprintf('channel_%d', ch);
            if ~isempty(signals_found)
                fitMethod = params.gaussFitMethod{ch};
                [all_signal_data, ~] = snap_helpers.computeSignalMeasurements(maximaCoords, fitResults, fitMethod);
                nuclei_data(nuc_idx).channels.(ch_field).signals = all_signal_data(signals_found);
                nuclei_data(nuc_idx).channels.(ch_field).signal_count = length(signals_found);
            end
            % If no signals found, the field keeps its pre-initialized empty values
        end
    end
    
    signalComp.nuclei_data = nuclei_data;
    
    % Build spotsPerNucleus matrix for easy pattern analysis
    % Matrix where each row = nucleus, each column = channel signal count
    spotsPerNucleus = zeros(numNuclei, numChannels);
    for nuc_idx = 1:numNuclei
        for ch = 1:numChannels
            ch_field = sprintf('channel_%d', ch);
            if isfield(nuclei_data(nuc_idx).channels, ch_field)
                spotsPerNucleus(nuc_idx, ch) = nuclei_data(nuc_idx).channels.(ch_field).signal_count;
            end
        end
    end
    signalComp.spotsPerNucleus = spotsPerNucleus;
end

function handles = createMinimalHandles(params)
% Create minimal handles structure needed for processing functions
    % PERFORMANCE CRITICAL: Just copy params, create minimal mocks
    handles = params;
    handles.Nmax = 5;
    
    % Create lightweight UI element mocks - helper functions expect .Value properties
    % Instead of creating hundreds of structs, create them on-demand or make params directly accessible
    
    % Nuclei spacing
    handles.nucXYSpacingInput.Value = params.nucXYSpacing;
    handles.nucZSpacingInput.Value = params.nucZSpacing;
    handles.nucPreprocEnabledCheck.Value = params.nucPreprocEnabled;
    handles.nucPreprocessModeDrop = struct('Value', ensureScalarString(params.nucPreProcMode));
    handles.nucPreprocessScaleCheck = struct('Value', ensureScalar(params.nucPreProcScale));
    handles.nucPreprocessProjectionDrop = struct('Value', ensureScalarString(params.nucPreProcProjection));
    handles.nucPreprocMethodDrop = struct('Value', ensureScalarString(params.nucPreProcMethod));
    handles.nucPreprocClipChecks = struct('Value', ensureScalar(params.nucPreprocClipAtZero));
    
    % Nuclei preprocessing method-specific
    handles.nucPreprocParam1Inputs = struct('Value', ensureScalar(params.nucSmoothGaussianValue));
    handles.nucWaveletNameDrop = struct('Value', ensureScalarString(params.nucWaveletName));
    handles.nucWaveletLevelInput = struct('Value', ensureScalar(params.nucWaveletLevel));
    handles.nucWaveletThresholdRuleDrop = struct('Value', ensureScalarString(params.nucWaveletThresholdRule));
    handles.nucWaveletThresholdMethodDrop = struct('Value', ensureScalarString(params.nucWaveletThresholdMethod));
    handles.nucNlmFilterStrengthInput = struct('Value', ensureScalar(params.nucNlmFilterStrength));
    handles.nucNlmSearchWindowInput = struct('Value', ensureScalar(params.nucNlmSearchWindow));
    handles.nucNlmComparisonWindowInput = struct('Value', ensureScalar(params.nucNlmComparisonWindow));
    
    % Nuclei background correction
    if isfield(params, 'nucBgCorrEnabled')
        handles.nucBgCorrEnabledCheck = struct('Value', ensureScalar(params.nucBgCorrEnabled));
        handles.nucBgCorrModeDrop = struct('Value', ensureScalarString(params.nucBgCorrMode));
        handles.nucBgCorrScaleCheck = struct('Value', ensureScalar(params.nucBgCorrScale));
        handles.nucBgCorrProjectionDrop = struct('Value', ensureScalarString(params.nucBgCorrProjection));
        handles.nucBgMethodDrop = struct('Value', ensureScalarString(params.nucBgMethod));
        handles.nucBgParamInput = struct('Value', ensureScalar(params.nucBgParam));
        handles.nucBgCorrClipChecks = struct('Value', ensureScalar(params.nucBgCorrClipAtZero));
    end
    
    % Nuclei deconvolution
    if isfield(params, 'nucDeconvEnabled')
        handles.nucDeconvEnabledCheck = struct('Value', ensureScalar(params.nucDeconvEnabled));
        handles.nucDeconvMethodDrop = struct('Value', ensureScalarString(params.nucDeconvMethod));
        handles.nucDeconvLRIterationsInput = struct('Value', ensureScalar(params.nucDeconvLRIterations));
        handles.nucDeconvLRDampingInput = struct('Value', ensureScalar(params.nucDeconvLRDamping));
        handles.nucDeconvWienerNSRInput = struct('Value', ensureScalar(params.nucDeconvWienerNSR));
        handles.nucDeconvBlindIterationsInput = struct('Value', ensureScalar(params.nucDeconvBlindIterations));
        handles.nucDeconvBlindUnderRelaxInput = struct('Value', ensureScalar(params.nucDeconvBlindUnderRelax));
        handles.nucDeconvPSFSourceDrop = struct('Value', ensureScalarString(params.nucDeconvPSFSource));
        handles.nucDeconvPSFPathText = struct('Value', ensureScalarString(params.nucDeconvPSFFilePath));
        handles.nucDeconvPSFSigmaXYInput = struct('Value', ensureScalar(params.nucDeconvPSFSigmaXY));
        handles.nucDeconvPSFSigmaZInput = struct('Value', ensureScalar(params.nucDeconvPSFSigmaZ));
        handles.nucDeconvPSFSizeXYInput = struct('Value', ensureScalar(params.nucDeconvPSFSizeXY));
        handles.nucDeconvPSFSizeZInput = struct('Value', ensureScalar(params.nucDeconvPSFSizeZ));
    end
    
    % Nuclei segmentation
    handles.nucSegEnabledCheck = struct('Value', ensureScalar(params.nucSegEnabled));
    handles.nucSegModeDrop = struct('Value', ensureScalarString(params.nucSegMode));
    handles.nucSegProjectionDrop = struct('Value', ensureScalarString(params.nucSegProjection));
    handles.nucSegMainMethodDrop = struct('Value', ensureScalarString(params.nucSegMainMethod));
    handles.nucSegSubMethodDrop = struct('Value', ensureScalarString(params.nucSegSubMethod));
    handles.nucSegAbsoluteInput = struct('Value', ensureScalar(params.nucSegAbsoluteThreshold));
    handles.nucSegStdMultiplierInput = struct('Value', ensureScalar(params.nucSegStdMultiplier));
    handles.nucSegOffsetInput = struct('Value', ensureScalar(params.nucSegOffset));
    handles.nucSegLocalAlgorithmDrop = struct('Value', ensureScalarString(params.nucSegLocalAlgorithm));
    handles.nucSegAlgParamInput = struct('Value', ensureScalar(params.nucSegAlgParam1), 'Visible', 'on');
    handles.nucSegAlgParam2Input = struct('Value', ensureScalar(params.nucSegAlgParam2), 'Visible', 'on');
    
    % Nuclei filtering (ensure scalar values)
    handles.nucFilterEnabledCheck = struct('Value', ensureScalar(params.nucFilterEnabled));
    handles.nucFilterSizeEnabledCheck = struct('Value', ensureScalar(params.nucFilterSizeEnabled));
    handles.nucFilterMinSizeInput = struct('Value', ensureScalar(params.nucFilterMinSize));
    handles.nucFilterSizeUnitDrop = struct('Value', ensureScalarString(params.nucFilterSizeUnit));
    handles.nucFilterCircularityEnabledCheck = struct('Value', ensureScalar(params.nucFilterCircularityEnabled));
    handles.nucFilterMinCircularityInput = struct('Value', ensureScalar(params.nucFilterMinCircularity));
    handles.nucFilterSolidityEnabledCheck = struct('Value', ensureScalar(params.nucFilterSolidityEnabled));
    handles.nucFilterMinSolidityInput = struct('Value', ensureScalar(params.nucFilterMinSolidity));
    handles.nucExcludeEdgesCheck = struct('Value', ensureScalar(params.nucExcludeEdges));
    handles.nucInclusionExclusionEnabledCheck = struct('Value', ensureScalar(params.nucInclusionExclusionEnabled));
    handles.nucInclusionExclusionModeDrop = struct('Value', ensureScalarString(params.nucInclusionExclusionMode));
    handles.nucInclusionExclusionApplyDrop = struct('Value', ensureScalarString(params.nucInclusionExclusionApplyTo));
    
    % Channel parameters (arrays) - ensure all values are scalar
    for ch = 1:handles.Nmax
        handles.xySpacingInputs(ch) = struct('Value', ensureScalar(params.xySpacing{ch}));
        handles.zSpacingInputs(ch) = struct('Value', ensureScalar(params.zSpacing{ch}), 'Enable', 'on');
        
        % Deconvolution
        if isfield(params, 'deconvEnabled')
            handles.deconvEnabledChecks(ch) = struct('Value', ensureScalar(params.deconvEnabled{ch}));
            handles.deconvMethodDrops(ch) = struct('Value', ensureScalarString(params.deconvMethod{ch}));
            handles.deconvPSFSourceDrops(ch) = struct('Value', ensureScalarString(params.deconvPSFSource{ch}));
            handles.deconvPSFPathTexts(ch) = struct('Value', ensureScalarString(params.deconvPSFFilePath{ch}));
            handles.deconvPSFSigmaXYInputs(ch) = struct('Value', ensureScalar(params.deconvPSFSigmaXY{ch}));
            handles.deconvPSFSigmaZInputs(ch) = struct('Value', ensureScalar(params.deconvPSFSigmaZ{ch}));
            handles.deconvPSFSizeXYInputs(ch) = struct('Value', ensureScalar(params.deconvPSFSizeXY{ch}));
            handles.deconvPSFSizeZInputs(ch) = struct('Value', ensureScalar(params.deconvPSFSizeZ{ch}));
            handles.deconvLRIterationsInputs(ch) = struct('Value', ensureScalar(params.deconvLRIterations{ch}));
            handles.deconvLRDampingInputs(ch) = struct('Value', ensureScalar(params.deconvLRDamping{ch}));
            handles.deconvWienerNSRInputs(ch) = struct('Value', ensureScalar(params.deconvWienerNSR{ch}));
            handles.deconvBlindIterationsInputs(ch) = struct('Value', ensureScalar(params.deconvBlindIterations{ch}));
            handles.deconvBlindUnderRelaxInputs(ch) = struct('Value', ensureScalar(params.deconvBlindUnderRelax{ch}));
        end
        
        % Preprocessing
        handles.preprocEnabledChecks(ch) = struct('Value', ensureScalar(params.preprocEnabled{ch}));
        handles.preprocessModeDrops(ch) = struct('Value', ensureScalarString(params.preProcMode{ch}));
        handles.preprocessScaleChecks(ch) = struct('Value', ensureScalar(params.preProcScale{ch}));
        handles.preprocessProjectionDrops(ch) = struct('Value', ensureScalarString(params.preProcProjection{ch}));
        handles.preprocMethodDrops(ch) = struct('Value', ensureScalarString(params.preProcMethod{ch}));
        handles.gaussInputs(ch) = struct('Value', ensureScalar(params.smoothGaussianValues{ch}));
        handles.medianInputs(ch) = struct('Value', ensureScalar(params.smoothMedianValues{ch}));
        handles.waveletNameDrops(ch) = struct('Value', ensureScalarString(params.waveletName{ch}));
        handles.waveletLevelInputs(ch) = struct('Value', ensureScalar(params.waveletLevel{ch}));
        handles.waveletThresholdRuleDrops(ch) = struct('Value', ensureScalarString(params.waveletThresholdRule{ch}));
        handles.waveletThresholdMethodDrops(ch) = struct('Value', ensureScalarString(params.waveletThresholdMethod{ch}));
        handles.nlmFilterStrengthInputs(ch) = struct('Value', ensureScalar(params.nlmFilterStrength{ch}));
        handles.nlmSearchWindowInputs(ch) = struct('Value', ensureScalar(params.nlmSearchWindow{ch}));
        handles.nlmComparisonWindowInputs(ch) = struct('Value', ensureScalar(params.nlmComparisonWindow{ch}));
        handles.preprocClipChecks(ch) = struct('Value', ensureScalar(params.preprocClipAtZero{ch}));
        
        % Background correction
        handles.bgCorrEnabledChecks(ch) = struct('Value', ensureScalar(params.bgCorrEnabled{ch}));
        handles.bgCorrModeDrops(ch) = struct('Value', ensureScalarString(params.bgCorrMode{ch}));
        handles.bgCorrScaleChecks(ch) = struct('Value', ensureScalar(params.bgCorrScale{ch}));
        handles.bgCorrProjectionDrops(ch) = struct('Value', ensureScalarString(params.bgCorrProjection{ch}));
        handles.bgMethodDrops(ch) = struct('Value', ensureScalarString(params.bgMethod{ch}));
        handles.bgParamInputs(ch) = struct('Value', ensureScalar(params.bgParam{ch}));
        handles.bgCorrClipChecks(ch) = struct('Value', ensureScalar(params.bgCorrClipAtZero{ch}));
        
        % Maxima detection
        handles.maximaEnabledChecks(ch) = struct('Value', ensureScalar(params.maximaEnabled{ch}));
        handles.maximaModeDrops(ch) = struct('Value', ensureScalarString(params.maximaMode{ch}));
        handles.maximaScaleChecks(ch) = struct('Value', ensureScalar(params.maximaScale{ch}));
        handles.maximaProjectionDrops(ch) = struct('Value', ensureScalarString(params.maximaProjection{ch}));
        handles.maximaMethodDrops(ch) = struct('Value', ensureScalarString(params.maximaMethod{ch}));
        handles.maximaNeighborhoodInputs(ch) = struct('Value', ensureScalar(params.maximaNeighborhoodSize{ch}));
        handles.hMaxInputs(ch) = struct('Value', ensureScalar(params.hMaxValue{ch}));
        handles.logSigmaInputs(ch) = struct('Value', ensureScalar(params.sigmaValue{ch}));
        handles.logThresholdInputs(ch) = struct('Value', ensureScalar(params.peakThresholdValue{ch}));
        
        % Gaussian fitting
        handles.gaussFitEnabledChecks(ch) = struct('Value', ensureScalar(params.gaussFitEnabled{ch}));
        handles.gaussFitVoxelWindowSlider(ch) = struct('Value', ensureScalar(params.gaussFitVoxelWindowSize{ch}));
        handles.gaussFitBgCorrMethodDrop(ch) = struct('Value', ensureScalarString(params.gaussFitBgCorrMethod{ch}));
        handles.gaussFitBgCorrWidthEdit(ch) = struct('Value', ensureScalar(params.gaussFitBgCorrWidth{ch}));
        handles.gaussFitPolyDegreeEdit(ch) = struct('Value', ensureScalar(params.gaussFitPolyDegree{ch}));
        handles.gaussFitMethodDrop(ch) = struct('Value', ensureScalarString(params.gaussFitMethod{ch}));
        handles.gaussFitMaxIterationsEdit(ch) = struct('Value', ensureScalar(params.gaussFitMaxIterations{ch}));
        handles.gaussFitToleranceEdit(ch) = struct('Value', ensureScalar(params.gaussFitTolerance{ch}));
        handles.gaussFitRadialRadiusEdit(ch) = struct('Value', ensureScalar(params.gaussFitRadialRadius{ch}));
        
        % Fit filtering
        handles.fitFilterEnabledChecks(ch) = struct('Value', ensureScalar(params.fitFilterEnabled{ch}));
        handles.fitFilterRSquaredEnabledChecks(ch) = struct('Value', ensureScalar(params.fitFilterRSquaredEnabled{ch}));
        handles.fitFilterRSquaredMinInputs(ch) = struct('Value', ensureScalar(params.fitFilterRSquaredMin{ch}));
        handles.fitFilterRSquaredMaxInputs(ch) = struct('Value', ensureScalar(params.fitFilterRSquaredMax{ch}));
        handles.fitFilterSigmaSumEnabledChecks(ch) = struct('Value', ensureScalar(params.fitFilterSigmaSumEnabled{ch}));
        handles.fitFilterSigmaSumMinInputs(ch) = struct('Value', ensureScalar(params.fitFilterSigmaSumMin{ch}));
        handles.fitFilterSigmaSumMaxInputs(ch) = struct('Value', ensureScalar(params.fitFilterSigmaSumMax{ch}));
        handles.fitFilterAmplitudeEnabledChecks(ch) = struct('Value', ensureScalar(params.fitFilterAmplitudeEnabled{ch}));
        handles.fitFilterAmplitudeMinInputs(ch) = struct('Value', ensureScalar(params.fitFilterAmplitudeMin{ch}));
        handles.fitFilterAmplitudeMaxInputs(ch) = struct('Value', ensureScalar(params.fitFilterAmplitudeMax{ch}));
        handles.fitFilterIntensityEnabledChecks(ch) = struct('Value', ensureScalar(params.fitFilterIntensityEnabled{ch}));
        handles.fitFilterIntensityMinInputs(ch) = struct('Value', ensureScalar(params.fitFilterIntensityMin{ch}));
        handles.fitFilterIntensityMaxInputs(ch) = struct('Value', ensureScalar(params.fitFilterIntensityMax{ch}));
    end
end

function fitParams = extractFitParams(params, channelIdx)
% Extract fitting parameters for a specific channel
    fitParams = struct();
    fitParams.gaussFitVoxelWindowSize = params.gaussFitVoxelWindowSize{channelIdx};
    fitParams.gaussFitBgCorrMethod = params.gaussFitBgCorrMethod{channelIdx};
    fitParams.gaussFitBgCorrWidth = params.gaussFitBgCorrWidth{channelIdx};
    fitParams.gaussFitPolyDegree = params.gaussFitPolyDegree{channelIdx};
    fitParams.gaussFitMethod = params.gaussFitMethod{channelIdx};
    fitParams.gaussFitMaxIterations = params.gaussFitMaxIterations{channelIdx};
    fitParams.gaussFitTolerance = params.gaussFitTolerance{channelIdx};
    fitParams.gaussFitRadialRadius = params.gaussFitRadialRadius{channelIdx};
    fitParams.gaussFitPlotCheck = false; % No plotting in batch mode
    fitParams.xySpacing = params.xySpacing{channelIdx};
    fitParams.zSpacing = params.zSpacing{channelIdx};
end

function exportedTypes = exportSetResults(results, imageSet, outputDir, setName, timestamp, exportOptions, params, progressHandle)
% Export results for a single image set based on AVAILABLE data and export options
% Returns struct indicating what was actually exported
    
    if nargin < 8
        progressHandle = [];
    end
%
% ============================================================================
% BATCH EXPORT CONSISTENCY SYSTEM (Critical Documentation)
% ============================================================================
%
% This function implements the SAME consistency pattern as SNAP GUI's exportData():
%
% ARCHITECTURE:
%
%   PHASE 1: PRE-COMPUTATION (lines ~900-920)
%   ------------------------------------------
%   - Check if nuclei-related exports requested (nucleiData OR clusteredData)
%   - If yes:
%     → Call computeNucleusMeasurements() ONCE for this image set
%     → Store in sharedMeasurements.nucleiData
%     → Generate unique computation_id
%   - If no → skip (no shared computation needed)
%
%   PHASE 2: CONDITIONAL EXPORTS (lines ~922-1054)
%   -----------------------------------------------
%   - Export nuclei (if requested): Pass sharedMeasurements
%   - Export channels (if requested): Independent (no nucleus data needed)
%   - Export nuclei+signals (if requested): Pass sharedMeasurements
%   - Export visualizations (if requested): Independent
%
% CONSISTENCY GUARANTEE:
%   Within a single image set, all nuclei-related exports (nuclei and nuclei+signal)
%   receive IDENTICAL nucleus measurements from sharedMeasurements.nucleiData.
%   This is verified by matching computation_id in the receipt files.
%
% CONDITIONAL EXPORT LOGIC:
%   - Only exports data that exists AND was requested
%   - Never throws errors if data is missing
%   - Returns exportedTypes struct tracking what was successfully exported
%   - Gracefully handles missing channels, missing signals, etc.
%
% ============================================================================

    exportedTypes = struct();
    exportedTypes.nucleiData = false;
    exportedTypes.channelData = false;
    exportedTypes.clusteredData = false;
    exportedTypes.visualizations = false;
    
    % Check what data is available
    hasNuclei = isfield(results, 'nucleiMask') && ~isempty(results.nucleiMask) && ...
                isfield(results, 'numNuclei') && results.numNuclei > 0;
    hasChannels = isfield(results, 'channels') && ~isempty(results.channels);
    hasSignals = false;
    if hasChannels
        for ch = 1:length(results.channels)
            if ~isempty(results.channels{ch}) && isfield(results.channels{ch}, 'maximaCoords') && ...
               ~isempty(results.channels{ch}.maximaCoords) && size(results.channels{ch}.maximaCoords, 1) > 0
                hasSignals = true;
                break;
            end
        end
    end
    
    
    %% ========================================================================
    %  PRE-COMPUTE SHARED NUCLEUS MEASUREMENTS (Consistency System)
    %  ========================================================================
    %
    %  This implements the SAME "compute once, export many" pattern as SNAP GUI.
    %
    %  LOGIC:
    %    1. Check if ANY nuclei-related exports are requested
    %    2. If yes → compute measurements ONCE before export loop
    %    3. Pass sharedMeasurements to all nuclei-related exports
    %
    %  BENEFIT:
    %    Within a single image set, all exports (nuclei, nuclei+signal) have
    %    GUARANTEED IDENTICAL nucleus measurements. This is verified by the
    %    computation_id that appears in all receipt files.
    %
    %  BATCH CONTEXT:
    %    This happens ONCE per image set (per subfolder).
    %    Aggregate exports combine these consistent measurements.
    %
    %  ========================================================================
    
    sharedMeasurements = struct();
    needsNucleiMeasurements = (exportOptions.nucleiData || exportOptions.clusteredData) && hasNuclei;
    
    if needsNucleiMeasurements
        try
            % Prepare options
            measureOptions = struct();
            measureOptions.is_3d = (ndims(results.nucleiMask) == 3 && size(results.nucleiMask, 3) > 1);
            measureOptions.xy_spacing = params.nucXYSpacing;
            measureOptions.z_spacing = params.nucZSpacing;
            measureOptions.image_name = setName;
            
            % CRITICAL: Compute ONCE - this result will be shared by all nuclei exports
            % The function generates a unique computation_id for this image set
            sharedMeasurements.nucleiData = snap_helpers.computeNucleusMeasurements(results.nucleusLabels, results.nucleiMask, measureOptions);
        catch ME
            warning('Failed to compute shared nucleus measurements for %s: %s', setName, ME.message);
            sharedMeasurements.nucleiData = [];
        end
    end
    
    %% 1. Export Nuclei Data (if requested AND available)
    if exportOptions.nucleiData && hasNuclei
        try
            % Use standardized export function (same as SNAP)
            outputPath = fullfile(outputDir, sprintf('export_%s_nuclei', setName));
            
            % Prepare options
            nucleiExportOptions = struct();
            nucleiExportOptions.is_3d = (ndims(results.nucleiMask) == 3 && size(results.nucleiMask, 3) > 1);
            nucleiExportOptions.xy_spacing = params.nucXYSpacing;
            nucleiExportOptions.z_spacing = params.nucZSpacing;
            nucleiExportOptions.image_name = setName;  % Include image name for batch
            
            % CONSISTENCY: Pass pre-computed measurements if available
            % This ensures nuclei and nuclei+signal exports have IDENTICAL data
            if ~isempty(sharedMeasurements.nucleiData)
                snap_helpers.exportNucleiDataStandardized(results.nucleusLabels, results.nucleiMask, outputPath, nucleiExportOptions, sharedMeasurements.nucleiData);
            else
                % Fallback: let function compute itself (shouldn't happen if pre-computation worked)
            snap_helpers.exportNucleiDataStandardized(results.nucleusLabels, results.nucleiMask, outputPath, nucleiExportOptions);
            end
            
            exportedTypes.nucleiData = true;
        catch ME
            warning('Failed to export nuclei data for %s: %s', setName, ME.message);
        end
    end
    
    %% 2. Export Channel Data (if requested AND available)
    if exportOptions.channelData && hasChannels
        if ~isempty(progressHandle)
            appendLog(progressHandle, '  Exporting channel signal data...');
        end
        
        % Determine if 3D from actual data (same method as nuclei export)
        is_3d_data = (ndims(results.nucleiMask) == 3 && size(results.nucleiMask, 3) > 1);
        if is_3d_data
            imageType = '3D Image Stack';
        else
            imageType = '2D Image';
        end
        
        for ch = 1:length(results.channels)
            % Check if channel has signals (size > 0, not just ~isempty)
            hasSignalsInCh = ~isempty(results.channels{ch}) && ...
                         isfield(results.channels{ch}, 'maximaCoords') && ...
                         ~isempty(results.channels{ch}.maximaCoords) && ...
                         size(results.channels{ch}.maximaCoords, 1) > 0;
            
            if hasSignalsInCh
                try
                    % Build channel results structure
                    channelResults = struct();
                    channelResults.maxima_coords = results.channels{ch}.maximaCoords;
                    if isfield(results.channels{ch}, 'fitResults') && ~isempty(results.channels{ch}.fitResults)
                        channelResults.fit_results = results.channels{ch}.fitResults;
                    end
                    
                    % Build fit parameters structure from params
                    fitParams = struct();
                    fitParams.imageType = imageType;
                    
                    % Gaussian fitting parameters (with defaults for missing fields)
                    if isfield(params, 'gaussFitMode') && ch <= length(params.gaussFitMode)
                        fitParams.gaussFitMode = params.gaussFitMode{ch};
                    else
                        fitParams.gaussFitMode = 'Not specified';
                    end
                    
                    if isfield(params, 'gaussFitMethod') && ch <= length(params.gaussFitMethod)
                        fitParams.gaussFitMethod = params.gaussFitMethod{ch};
                    else
                        fitParams.gaussFitMethod = 'Not specified';
                    end
                    
                    if isfield(params, 'gaussFitBgCorrMethod') && ch <= length(params.gaussFitBgCorrMethod)
                        fitParams.gaussFitBgCorrMethod = params.gaussFitBgCorrMethod{ch};
                    else
                        fitParams.gaussFitBgCorrMethod = 'Mean Surrounding Subtraction';
                    end
                    
                    if isfield(params, 'gaussFitBgCorrWidth') && ch <= length(params.gaussFitBgCorrWidth)
                        fitParams.gaussFitBgCorrWidth = params.gaussFitBgCorrWidth{ch};
                    else
                        fitParams.gaussFitBgCorrWidth = 2;
                    end
                    if isfield(params, 'gaussFitPolyDegree')
                        fitParams.gaussFitPolyDegree = params.gaussFitPolyDegree{ch};
                    else
                        fitParams.gaussFitPolyDegree = 2;
                    end
    
                    % Additional fitting parameters
                    if isfield(params, 'gaussFitVoxelWindowSize')
                        fitParams.gaussFitVoxelWindowSize = params.gaussFitVoxelWindowSize{ch};
                    else
                        fitParams.gaussFitVoxelWindowSize = 7;
                    end
                    
                    if isfield(params, 'gaussFitMaxIterations')
                        fitParams.gaussFitMaxIterations = params.gaussFitMaxIterations{ch};
                    else
                        fitParams.gaussFitMaxIterations = 200;
                    end
                    
                    if isfield(params, 'gaussFitTolerance')
                        fitParams.gaussFitTolerance = params.gaussFitTolerance{ch};
                    else
                        fitParams.gaussFitTolerance = 1e-6;
                    end
                    
                    if isfield(params, 'gaussFitRadialRadius')
                        fitParams.gaussFitRadialRadius = params.gaussFitRadialRadius{ch};
                    else
                        fitParams.gaussFitRadialRadius = 3;
                    end
                    
                    % Build options structure
                    channelOptions = struct();
                    channelOptions.is_3d = is_3d_data;
                    if isfield(params, 'xySpacing') && ch <= length(params.xySpacing)
                        xy_spacing = params.xySpacing{ch};
                        z_spacing = 1.0;
                        if isfield(params, 'zSpacing') && ch <= length(params.zSpacing)
                            z_spacing = params.zSpacing{ch};
                        end
                        channelOptions.spacing = [xy_spacing, xy_spacing, z_spacing];
                    end
                    
                    % Build output path
                    outputPath = fullfile(outputDir, sprintf('export_%s_ch%d', setName, ch));
                    
                    % Call standardized export function
                    snap_helpers.exportChannelDataStandardized(outputPath, setName, channelResults, fitParams, channelOptions);
                    
                    if ~isempty(progressHandle)
                        appendLog(progressHandle, sprintf('    ✓ Channel %d data exported', ch));
                    end
                    
                    exportedTypes.channelData = true;
                catch ME
                    warning('Failed to export channel %d data for %s: %s', ch, setName, ME.message);
                    if ~isempty(progressHandle)
                        appendLog(progressHandle, sprintf('    ✗ Channel %d export failed: %s', ch, ME.message));
                    end
                end
            else
                % Channel has no signals - skip export
                if ~isempty(progressHandle)
                    appendLog(progressHandle, sprintf('    ⊘ Channel %d: No signals to export', ch));
                end
            end
        end
    end
    
    %% 3. Export Nuclei Signal Data (if requested AND nuclei AND signals available)
    if exportOptions.clusteredData && hasNuclei && hasSignals
        try
            % Check if signal composition was already computed
            if ~isfield(results, 'signalComposition') || isempty(results.signalComposition)
                warning('Signal composition not computed for %s. Skipping nuclei signal export.', setName);
            else
                % Build options structure
                nucleiSignalOptions = struct();
                nucleiSignalOptions.is_3d = (ndims(results.nucleiMask) == 3 && size(results.nucleiMask, 3) > 1);
                nucleiSignalOptions.xy_spacing = params.nucXYSpacing;
                nucleiSignalOptions.z_spacing = params.nucZSpacing;
                
                % Build output path
                outputPath = fullfile(outputDir, sprintf('export_%s_nucleisignals', setName));
                
                % CONSISTENCY: Pass SAME pre-computed measurements used in nuclei export above
                % This guarantees nucleus data is IDENTICAL in both export files
                % Verified by matching computation_id in receipt files
                if ~isempty(sharedMeasurements.nucleiData)
                    snap_helpers.exportNucleiSignalDataStandardized(outputPath, setName, results.nucleusLabels, ...
                        results.nucleiMask, results.signalComposition, nucleiSignalOptions, sharedMeasurements.nucleiData);
                else
                    % Fallback: let function compute itself (shouldn't happen if pre-computation worked)
                snap_helpers.exportNucleiSignalDataStandardized(outputPath, setName, results.nucleusLabels, ...
                    results.nucleiMask, results.signalComposition, nucleiSignalOptions);
                end
                
                exportedTypes.clusteredData = true;
            end
        catch ME
            warning('Failed to export nuclei signal data for %s: %s', setName, ME.message);
        end
    end
    
    %% 4. Export Visualization (if requested AND nuclei AND signals available)
    if exportOptions.visualizations && (hasNuclei || hasSignals)
        try
            exportVisualizationOverview(results, imageSet, outputDir, setName, params);
            exportedTypes.visualizations = true;
        catch ME
            warning('Failed to export visualization for %s: %s', setName, ME.message);
        end
    end
end

function oldExportProcessedImages() % PLACEHOLDER - keeping old processed image export for reference
    % Export processed images
    if false  % Disabled for now
        imageFormat = lower(exportOptions.imageFormat);
        
        % Export nuclei image if available
        if isfield(results, 'nucleiProcessed') && ~isempty(results.nucleiProcessed)
            filename = fullfile(outputDir, sprintf('export_%s_nuclei_image.%s', setName, imageFormat));
            if strcmp(imageFormat, 'tiff') || strcmp(imageFormat, 'tif')
                imwrite(results.nucleiProcessed, filename, 'Compression', 'none');
            else
                imwrite(results.nucleiProcessed, filename);
            end
        end
        
        % Export channel images
        for ch = 1:length(results.channels)
            if ~isempty(results.channels{ch}) && isfield(results.channels{ch}, 'imageProcessed')
                filename = fullfile(outputDir, sprintf('export_%s_ch%d_image.%s', setName, ch, imageFormat));
                if strcmp(imageFormat, 'tiff') || strcmp(imageFormat, 'tif')
                    imwrite(results.channels{ch}.imageProcessed, filename, 'Compression', 'none');
                else
                    imwrite(results.channels{ch}.imageProcessed, filename);
                end
            end
        end
    end
end

function exportAggregateResults(batchResults, outputDir, timestamp, exportOptions, params)
% Export aggregate results - Simple concatenation of all CSV files from subfolders
%
% Creates aggregate files in main output folder by concatenating all subfolder exports:
% - export_batch_nuclei.csv: All nuclei from all images
% - export_batch_ch1.csv, export_batch_ch2.csv, etc.: All signals per channel from all images  
% - export_batch_nucleisignals_simple.csv: All nuclei with signal counts from all images
% - export_batch_nucleisignals_expanded.csv: All nuclei with full signal data from all images

% Extract successful results only
successfulResults = {};
for i = 1:length(batchResults.results)
    if ~isempty(batchResults.results{i}) && isfield(batchResults.results{i}, 'status') && ...
       strcmp(batchResults.results{i}.status, 'success')
        successfulResults{end+1} = batchResults.results{i};
    end
end

numSuccessful = length(successfulResults);
if numSuccessful == 0
    fprintf('No successful results to aggregate.\n');
    return;
end

% Silent aggregation

%% 1. AGGREGATE NUCLEI DATA - Simple CSV concatenation
if exportOptions.nucleiData
    
    allCSVData = [];
    totalNuclei = 0;
    
    for i = 1:numSuccessful
        result = successfulResults{i};
        setName = result.name;
        setDir = fullfile(outputDir, setName);
        csvFile = fullfile(setDir, sprintf('export_%s_nuclei.csv', setName));
        
        if exist(csvFile, 'file')
            try
                csvData = readtable(csvFile);
                if isempty(allCSVData)
                    allCSVData = csvData;
                else
                    allCSVData = [allCSVData; csvData];
                end
                totalNuclei = totalNuclei + height(csvData);
            catch ME
                warning('Failed to load nuclei CSV for %s: %s', setName, ME.message);
            end
        end
    end
    
    if ~isempty(allCSVData)
        csvFilename = fullfile(outputDir, 'export_batch_nuclei.csv');
        writetable(allCSVData, csvFilename);
    end
end

%% 2. AGGREGATE CHANNEL DATA - Simple CSV concatenation per channel
if exportOptions.channelData
    
    numChannels = params.numChannels;
    
    for ch = 1:numChannels
        allCSVData = [];
        totalSignals = 0;
        
        for i = 1:numSuccessful
            result = successfulResults{i};
            setName = result.name;
            setDir = fullfile(outputDir, setName);
            csvFile = fullfile(setDir, sprintf('export_%s_ch%d.csv', setName, ch));
            
            if exist(csvFile, 'file')
                try
                    csvData = readtable(csvFile);
                    if isempty(allCSVData)
                        allCSVData = csvData;
                    else
                        allCSVData = [allCSVData; csvData];
                    end
                    totalSignals = totalSignals + height(csvData);
                catch ME
                    warning('Failed to load channel %d CSV for %s: %s', ch, setName, ME.message);
                end
            end
        end
        
        if ~isempty(allCSVData)
            csvFilename = fullfile(outputDir, sprintf('export_batch_ch%d.csv', ch));
            writetable(allCSVData, csvFilename);
        end
    end
end

%% 3. AGGREGATE NUCLEI+SIGNAL DATA - Simple CSV concatenation
if exportOptions.clusteredData
    
    % Aggregate simple CSV (nucleus data + signal counts)
    allSimpleCSV = [];
    totalNuclei_simple = 0;
    
    for i = 1:numSuccessful
        result = successfulResults{i};
        setName = result.name;
        setDir = fullfile(outputDir, setName);
        csvFile = fullfile(setDir, sprintf('export_%s_nucleisignals_simple.csv', setName));
        
        if exist(csvFile, 'file')
            try
                csvData = readtable(csvFile);
                if isempty(allSimpleCSV)
                    allSimpleCSV = csvData;
                else
                    allSimpleCSV = [allSimpleCSV; csvData];
                end
                totalNuclei_simple = totalNuclei_simple + height(csvData);
            catch ME
                warning('Failed to load simple nuclei+signal CSV for %s: %s', setName, ME.message);
            end
        end
    end
    
    if ~isempty(allSimpleCSV)
        csvFilename = fullfile(outputDir, 'export_batch_nucleisignals_simple.csv');
        writetable(allSimpleCSV, csvFilename);
    end
    
    % Aggregate expanded CSV (full nucleus + signal data)
    allExpandedCSV = [];
    totalRows_expanded = 0;
    
    for i = 1:numSuccessful
        result = successfulResults{i};
        setName = result.name;
        setDir = fullfile(outputDir, setName);
        csvFile = fullfile(setDir, sprintf('export_%s_nucleisignals_expanded.csv', setName));
        
        if exist(csvFile, 'file')
            try
                csvData = readtable(csvFile);
                if isempty(allExpandedCSV)
                    allExpandedCSV = csvData;
                else
                    allExpandedCSV = [allExpandedCSV; csvData];
                end
                totalRows_expanded = totalRows_expanded + height(csvData);
            catch ME
                warning('Failed to load expanded nuclei+signal CSV for %s: %s', setName, ME.message);
            end
        end
    end
    
    if ~isempty(allExpandedCSV)
        csvFilename = fullfile(outputDir, 'export_batch_nucleisignals_expanded.csv');
        writetable(allExpandedCSV, csvFilename);
    end
    
    % PATTERN ANALYSIS: Compute signal composition statistics
    
    % Extract signal counts per nucleus from simple CSV
    if ~isempty(allSimpleCSV)
        % Build pattern matrix from the simple CSV data
        % Find columns that contain signal counts (ch1_signals, ch2_signals, etc.)
        varNames = allSimpleCSV.Properties.VariableNames;
        signalCountCols = {};
        for ch = 1:params.numChannels
            colName = sprintf('ch%d_signals', ch);
            if ismember(colName, varNames)
                signalCountCols{end+1} = colName;
            end
        end
        
        if ~isempty(signalCountCols)
            % Extract pattern matrix (each row = nucleus, each col = channel signal count)
            patternMatrix = table2array(allSimpleCSV(:, signalCountCols));
            
            % Compute statistics
            stats = struct();
            stats.totalNuclei = size(patternMatrix, 1);
            stats.numChannels = size(patternMatrix, 2);
            
            % Basic statistics per channel
            stats.signalsPerChannel = sum(patternMatrix, 1);
            stats.meanSignalsPerNucleus = mean(patternMatrix, 1);
            stats.stdSignalsPerNucleus = std(patternMatrix, 0, 1);
            stats.maxSignalsInAnyNucleus = max(patternMatrix, [], 1);
            
            % Nuclei with NO signals (empty nuclei)
            totalSignalsPerNucleus = sum(patternMatrix, 2);
            stats.numEmptyNuclei = sum(totalSignalsPerNucleus == 0);
            stats.percentEmptyNuclei = (stats.numEmptyNuclei / stats.totalNuclei) * 100;
            
            % Nuclei with ANY signals
            stats.numNucleiWithSignals = sum(totalSignalsPerNucleus > 0);
            stats.percentNucleiWithSignals = (stats.numNucleiWithSignals / stats.totalNuclei) * 100;
            
            % Per-channel occupancy
            stats.nucleiPerChannel = sum(patternMatrix > 0, 1);
            stats.percentNucleiPerChannel = (stats.nucleiPerChannel / stats.totalNuclei) * 100;
            
            % Find unique composition patterns
            [uniquePatterns, ~, patternIdx] = unique(patternMatrix, 'rows');
            patternCounts = histcounts(patternIdx, 1:(size(uniquePatterns,1)+1));
            
            % Sort patterns by frequency
            [sortedCounts, sortIdx] = sort(patternCounts, 'descend');
            sortedPatterns = uniquePatterns(sortIdx, :);
            
            stats.numUniquePatterns = size(uniquePatterns, 1);
            
            % Create pattern statistics CSV
            statsCSVFilename = fullfile(outputDir, 'export_batch_pattern_statistics.csv');
            fid = fopen(statsCSVFilename, 'w');
            
            fprintf(fid, '=== Signal Composition Pattern Statistics ===\n');
            fprintf(fid, 'Generated: %s\n', datestr(now));
            fprintf(fid, 'Total Nuclei Analyzed: %d\n', stats.totalNuclei);
            fprintf(fid, 'Number of Image Sets: %d\n\n', numSuccessful);
            
            fprintf(fid, '--- Overall Statistics ---\n');
            fprintf(fid, 'Empty Nuclei (no signals): %d (%.1f%%)\n', stats.numEmptyNuclei, stats.percentEmptyNuclei);
            fprintf(fid, 'Nuclei with Signals: %d (%.1f%%)\n\n', stats.numNucleiWithSignals, stats.percentNucleiWithSignals);
            
            fprintf(fid, '--- Per-Channel Statistics ---\n');
            for ch = 1:stats.numChannels
                fprintf(fid, 'Channel %d:\n', ch);
                fprintf(fid, '  Total Signals: %d\n', stats.signalsPerChannel(ch));
                fprintf(fid, '  Mean per Nucleus: %.2f +/- %.2f\n', stats.meanSignalsPerNucleus(ch), stats.stdSignalsPerNucleus(ch));
                fprintf(fid, '  Max in any Nucleus: %d\n', stats.maxSignalsInAnyNucleus(ch));
                fprintf(fid, '  Nuclei Containing this Channel: %d (%.1f%%)\n\n', stats.nucleiPerChannel(ch), stats.percentNucleiPerChannel(ch));
            end
            
            fprintf(fid, '--- Top Signal Composition Patterns ---\n');
            fprintf(fid, 'Pattern,Count,Percentage\n');
            
            % Report top 20 patterns or all if fewer
            numPatternsToReport = min(20, stats.numUniquePatterns);
            for p = 1:numPatternsToReport
                % Create pattern description (e.g., [2;1;3] for ch1=2, ch2=1, ch3=3)
                patternStr = sprintf('[%d', sortedPatterns(p, 1));
                for ch = 2:size(sortedPatterns, 2)
                    patternStr = [patternStr, sprintf(';%d', sortedPatterns(p, ch))];
                end
                patternStr = [patternStr, ']'];
                
                percentage = (sortedCounts(p) / stats.totalNuclei) * 100;
                fprintf(fid, '%s,%d,%.2f%%\n', patternStr, sortedCounts(p), percentage);
            end
            
            fclose(fid);
        end
    end
end

fprintf('\n=== Aggregate Data Complete ===\n\n');
end

function exportVisualizationOverview(results, imageSet, outputDir, setName, params)
% Export annotated PNG visualization showing nuclei boundaries and signal locations
%
% Creates a comprehensive overview image with:
% - Grayscale base image (max projection of nuclei or first channel)
% - Magenta nucleus boundaries
% - Yellow nucleus labels showing ID and signal composition pattern (e.g., "1\n[2,1,3]")
% - Colored signal markers (crosses) with black borders at fitted/detected locations
% - Legend with channel colors and signal counts

try
    % Create visualizations subdirectory
    vizDir = fullfile(outputDir, 'visualizations');
    if ~exist(vizDir, 'dir'), mkdir(vizDir); end
    
    % Create figure (invisible for batch processing)
    fig = figure('Visible', 'off', 'Position', [100, 100, 1400, 1200], 'Color', 'black');
    ax = axes(fig, 'Position', [0.05 0.05 0.9 0.9]);
    
    %% 1. Prepare base image (nuclei max projection)
    if isfield(results, 'nucleiMask') && ~isempty(results.nucleiMask)
        if ndims(results.nucleiMask) == 3
            baseImg = max(results.nucleiMask, [], 3);
        else
            baseImg = results.nucleiMask;
        end
    else
        % Fallback: use first channel if available
        if ~isempty(results.channels) && ~isempty(results.channels{1}) && ...
           isfield(results.channels{1}, 'processedImage')
            procImg = results.channels{1}.processedImage;
            if ndims(procImg) == 3
                baseImg = max(procImg, [], 3);
            else
                baseImg = procImg;
            end
        else
            warning('No base image available for visualization');
            close(fig);
            return;
        end
    end
    
    % Display base image
    imagesc(ax, baseImg);
    axis(ax, 'image');
    colormap(ax, 'gray');
    hold(ax, 'on');
    
    %% 2. Draw nucleus boundaries (magenta)
    if isfield(results, 'nucleiMask') && ~isempty(results.nucleiMask)
        if ndims(results.nucleiMask) == 3
            maskProj = max(results.nucleiMask, [], 3);
        else
            maskProj = results.nucleiMask;
        end
        boundaries = bwboundaries(maskProj);
        for k = 1:length(boundaries)
            if ~isempty(boundaries{k})
                boundary = boundaries{k};
                % boundary is [row, col]; plot needs [x=col, y=row]
                plot(ax, boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2);
            end
        end
    end
    
    %% 3. Draw nucleus labels with signal composition patterns
    if isfield(results, 'nucleusLabels') && ~isempty(results.nucleusLabels) && ...
       isfield(results.nucleusLabels, 'centroids_3d')
        centroids = results.nucleusLabels.centroids_3d;  % [col, row, slice] (Cartesian)
        numNuclei = size(centroids, 1);
        
        % Check if signal composition data is available
        hasComposition = isfield(results, 'signalComposition') && ~isempty(results.signalComposition) && ...
                         isfield(results.signalComposition, 'spotsPerNucleus');
        
        for i = 1:numNuclei
            % Draw nucleus ID (always shown)
            text(ax, centroids(i,1), centroids(i,2), num2str(i), ...
                 'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            
            % Add signal composition pattern if available (30 pixels below ID)
            if hasComposition && i <= size(results.signalComposition.spotsPerNucleus, 1)
                counts = results.signalComposition.spotsPerNucleus(i, :);
                % Create pattern like [1,2,3]
                patternStr = sprintf('[%d', counts(1));
                for ch = 2:length(counts)
                    patternStr = [patternStr, sprintf(',%d', counts(ch))];
                end
                patternStr = [patternStr, ']'];
                
                % Centroids are [col, row, slice] (Cartesian); text needs [x=col, y=row]
                x_pos = centroids(i,1);
                y_pos = centroids(i,2) + 20;  % 20 pixels below nucleus ID
                
                % Draw black border text first (slightly offset in all directions)
                for dx = -1:1
                    for dy = -1:1
                        if dx ~= 0 || dy ~= 0
                            text(ax, x_pos + dx, y_pos + dy, patternStr, ...
                                 'Color', 'black', 'FontSize', 10, 'FontWeight', 'bold', ...
                                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        end
                    end
                end
                
                % Draw yellow pattern text on top
                text(ax, x_pos, y_pos, patternStr, ...
                     'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end
    
    %% 4. Get channel colors (from parameters or defaults)
    numChannels = length(results.channels);
    channelColors = getChannelColors(params, numChannels);
    
    %% 5. Draw signals from all channels
    legendEntries = {};
    legendColors = {};
    
    for ch = 1:numChannels
        if ~isempty(results.channels{ch}) && isfield(results.channels{ch}, 'maximaCoords') && ...
           ~isempty(results.channels{ch}.maximaCoords)
            
            % Get signal positions - ONLY use the filtered maximaCoords
            % fitResults and maximaCoords should be synchronized after filtering
            signalCoords = results.channels{ch}.maximaCoords;  % This is already filtered!
            
            % If fits are available, replace with fitted coordinates where available
            if isfield(results.channels{ch}, 'fitResults') && ~isempty(results.channels{ch}.fitResults)
                fitResults = results.channels{ch}.fitResults;
                % Iterate only up to the length of maximaCoords (they should match)
                numToPlot = min(size(signalCoords, 1), length(fitResults));
                for i = 1:numToPlot
                    if isfield(fitResults(i), 'globalFitCenter') && ~isempty(fitResults(i).globalFitCenter)
                        signalCoords(i, :) = fitResults(i).globalFitCenter;
                    end
                    % If no globalFitCenter, keep the maxima coord that's already there
                end
            end
    
            if ~isempty(signalCoords)
                markerColor = channelColors{ch};
                
                % Plot each signal marker with black border
                for i = 1:size(signalCoords, 1)
                    % signalCoords are [row, col, slice]; plot needs [x=col, y=row]
                    x_pos = signalCoords(i,2);
                    y_pos = signalCoords(i,1);
                    
                    % Draw black border first (slightly larger)
                    plot(ax, x_pos, y_pos, '+', ...
                         'Color', 'black', 'MarkerSize', 12, 'LineWidth', 3);
                    
                    % Draw colored marker on top
                    plot(ax, x_pos, y_pos, '+', ...
                         'Color', markerColor, 'MarkerSize', 10, 'LineWidth', 2);
                end
                
                % Add to legend
                legendEntries{end+1} = sprintf('Ch%d: %d signals', ch, size(signalCoords, 1));
                legendColors{end+1} = markerColor;
            end
        end
    end
    
    %% 6. Add legend (custom legend with colors)
    if ~isempty(legendEntries)
        % Create custom legend in upper right
        legendX = size(baseImg, 2) * 0.98;
        legendY = size(baseImg, 1) * 0.02;
        yOffset = size(baseImg, 1) * 0.04;
        
        for i = 1:length(legendEntries)
            % Legend marker with black border
            x_leg = legendX - 20;
            y_leg = legendY + (i-1)*yOffset;
            
            % Black border
            plot(ax, x_leg, y_leg, '+', ...
                 'Color', 'black', 'MarkerSize', 12, 'LineWidth', 3);
            % Colored marker
            plot(ax, x_leg, y_leg, '+', ...
                 'Color', legendColors{i}, 'MarkerSize', 10, 'LineWidth', 2);
            
            % Legend text
            text(ax, legendX, legendY + (i-1)*yOffset, legendEntries{i}, ...
                 'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end
    end
    
    %% 7. Add title and metadata
    totalSignals = results.totalSpots;
    title(ax, sprintf('%s | Nuclei: %d | Total Signals: %d', ...
          setName, results.numNuclei, totalSignals), ...
          'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % Style axes
    set(ax, 'Color', 'black', 'XColor', 'white', 'YColor', 'white');
    set(ax, 'XTick', [], 'YTick', []);
    
    %% 8. Save PNG
    filename = fullfile(vizDir, sprintf('export_%s_overview.png', setName));
    exportgraphics(fig, filename, 'Resolution', 300);
    
    % Clean up
    close(fig);
    
catch ME
    warning('Failed to export visualization for %s: %s', setName, ME.message);
    if exist('fig', 'var') && isvalid(fig)
        close(fig);
    end
end
end

function channelColors = getChannelColors(params, numChannels)
% Get channel colors from parameters or use defaults
% Returns cell array of RGB triplets

% Default color order for channels 1-5+
defaultColors = {
    [1, 0.2, 0.2],    % Red
    [0.2, 1, 0.2],    % Green  
    [0.3, 0.7, 1],    % Blue
    [1, 1, 0.2],      % Yellow
    [1, 0.2, 1],      % Magenta
    [0.2, 1, 1]       % Cyan
};

channelColors = cell(1, numChannels);

% Color name to RGB mapping (matches SNAP)
colorMap = containers.Map(...
    {'Red','Green','Blue','Yellow','Magenta','Cyan'}, ...
    defaultColors);

% Try to extract colors from parameters
for ch = 1:numChannels
    if isfield(params, 'maximaColor') && length(params.maximaColor) >= ch
        % Parameters from SNAP have maximaColor field
        colorName = params.maximaColor{ch};
        if isKey(colorMap, colorName)
            channelColors{ch} = colorMap(colorName);
        else
            % Fallback to default
            channelColors{ch} = defaultColors{mod(ch-1, length(defaultColors)) + 1};
        end
    else
        % No color info, use default order
        channelColors{ch} = defaultColors{mod(ch-1, length(defaultColors)) + 1};
    end
end
end

function exportSignalCompositionCSV(signalComp, numNuclei, filename)
% Export signal composition to CSV
    fid = fopen(filename, 'w');
    if fid == -1, return; end
    
    numChannels = size(signalComp.spotsPerNucleus, 2);
    
    % Header
    fprintf(fid, 'Nucleus_ID');
    for ch = 1:numChannels
        fprintf(fid, ',Ch%d_Spots', ch);
    end
    fprintf(fid, ',Total_Spots\n');
    
    % Data
    for nuc = 1:numNuclei
        fprintf(fid, '%d', nuc);
        totalForNuc = 0;
        for ch = 1:numChannels
            spots = signalComp.spotsPerNucleus(nuc, ch);
            fprintf(fid, ',%d', spots);
            totalForNuc = totalForNuc + spots;
        end
        fprintf(fid, ',%d\n', totalForNuc);
    end
    
    fclose(fid);
end

function proj = projectZ(img, type)
% Z-projection helper
    if ndims(img) < 3, proj = img; return; end
    switch type
        case 'Max', proj = max(img, [], 3);
        case 'Min', proj = min(img, [], 3);
        case 'Median', proj = median(img, 3);
        case 'Mean', proj = mean(img, 3);
        otherwise, proj = max(img, [], 3);
    end
end

function scalarVal = ensureScalar(val)
% Ensure value is a proper scalar (handle cells, arrays, etc.)
    if iscell(val)
        scalarVal = val{1};
    elseif numel(val) > 1
        scalarVal = val(1);
    else
        scalarVal = val;
    end
    
    % Ensure it's a scalar double or logical
    if islogical(scalarVal) || isnumeric(scalarVal)
        if ~isscalar(scalarVal)
            scalarVal = scalarVal(1);
        end
    end
end

function scalarStr = ensureScalarString(val)
% Ensure value is a proper char string (handle cells, strings, etc.)
    if iscell(val)
        scalarStr = val{1};
    else
        scalarStr = val;
    end
    
    % Convert to char if needed
    if isstring(scalarStr)
        scalarStr = char(scalarStr);
    elseif ~ischar(scalarStr)
        scalarStr = char(string(scalarStr));
    end
end

%% GUI Functions

function createBatchGUI()
% Create the batch processing GUI
    
    % Create the main figure (wider layout with side-by-side progress + instructions)
    fig = uifigure('Name', 'SNAP Batch Processing', 'Position', [50 30 1200 900]);
    
    % Main grid layout
    mainGrid = uigridlayout(fig, [1, 1]);
    mainGrid.Padding = [15 15 15 15];
    
    % Content panel
    contentPanel = uipanel(mainGrid, 'Title', 'Batch Processing Configuration', 'FontWeight', 'bold', 'FontSize', 12);
    contentGrid = uigridlayout(contentPanel, [9, 1]);
    % Row 8 will have progress log + instructions side by side (500px tall)
    contentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 500, 'fit'};
    contentGrid.Padding = [15 15 15 15];
    contentGrid.RowSpacing = 8;
    
    % Initialize handles structure
    handles = struct();
    handles.fig = fig;
    
    % --- Row 1: Input Directory ---
    inputDirPanel = uipanel(contentGrid, 'BorderType', 'none');
    inputDirPanel.Layout.Row = 1;
    inputDirGrid = uigridlayout(inputDirPanel, [2, 2]);
    inputDirGrid.RowHeight = {'fit', 'fit'};
    inputDirGrid.ColumnWidth = {'1x', 100};
    inputDirGrid.Padding = [0 0 0 0];
    inputDirGrid.RowSpacing = 3;
    
    lbl1 = uilabel(inputDirGrid, 'Text', 'Input Directory:', 'FontWeight', 'bold', ...
        'Tooltip', 'Folder containing subfolders with image sets');
    lbl2 = uilabel(inputDirGrid, 'Text', '');
    
    handles.inputDirEdit = uieditfield(inputDirGrid, 'Value', '', 'Editable', 'off', ...
        'Tooltip', 'Select folder containing experiment subfolders');
    handles.inputDirBrowseBtn = uibutton(inputDirGrid, 'Text', 'Browse...', ...
        'ButtonPushedFcn', @(src,evt) browseInputDir(fig));
    
    % --- Row 2: Parameter File ---
    paramFilePanel = uipanel(contentGrid, 'BorderType', 'none');
    paramFilePanel.Layout.Row = 2;
    paramFileGrid = uigridlayout(paramFilePanel, [2, 2]);
    paramFileGrid.RowHeight = {'fit', 'fit'};
    paramFileGrid.ColumnWidth = {'1x', 100};
    paramFileGrid.Padding = [0 0 0 0];
    paramFileGrid.RowSpacing = 3;
    
    lbl3 = uilabel(paramFileGrid, 'Text', 'Parameter File:', 'FontWeight', 'bold', ...
        'Tooltip', 'Parameter file from SNAP (lastUsedParameters.mat or exported .mat file)');
    lbl4 = uilabel(paramFileGrid, 'Text', '');
    
    handles.paramFileEdit = uieditfield(paramFileGrid, 'Value', '', 'Editable', 'off', ...
        'Tooltip', 'Select .mat parameter file from SNAP');
    handles.paramFileBrowseBtn = uibutton(paramFileGrid, 'Text', 'Browse...', ...
        'ButtonPushedFcn', @(src,evt) browseParamFile(fig));
    
    % --- Row 3: Output Directory ---
    outputDirPanel = uipanel(contentGrid, 'BorderType', 'none');
    outputDirPanel.Layout.Row = 3;
    outputDirGrid = uigridlayout(outputDirPanel, [2, 2]);
    outputDirGrid.RowHeight = {'fit', 'fit'};
    outputDirGrid.ColumnWidth = {'1x', 100};
    outputDirGrid.Padding = [0 0 0 0];
    outputDirGrid.RowSpacing = 3;
    
    lbl5 = uilabel(outputDirGrid, 'Text', 'Output Directory:', 'FontWeight', 'bold', ...
        'Tooltip', 'Directory where batch results will be saved');
    lbl6 = uilabel(outputDirGrid, 'Text', '');
    
    handles.outputDirEdit = uieditfield(outputDirGrid, 'Value', '', 'Editable', 'off', ...
        'Tooltip', 'Select output directory for batch results');
    handles.outputDirBrowseBtn = uibutton(outputDirGrid, 'Text', 'Browse...', ...
        'ButtonPushedFcn', @(src,evt) browseOutputDir(fig));
    
    % --- Row 4: Export Options ---
    exportOptionsPanel = uipanel(contentGrid, 'Title', 'Export Options', 'FontWeight', 'bold');
    exportOptionsPanel.Layout.Row = 4;
    exportOptionsGrid = uigridlayout(exportOptionsPanel, [2, 2]);
    exportOptionsGrid.RowHeight = {'fit', 'fit'};
    exportOptionsGrid.ColumnWidth = {'1x', '1x'};
    exportOptionsGrid.Padding = [10 10 10 10];
    exportOptionsGrid.RowSpacing = 8;
    
    handles.exportNucleiDataCheck = uicheckbox(exportOptionsGrid, 'Text', 'Nuclei Data', ...
        'Value', true, 'Tooltip', 'Export nuclei segmentation masks and properties');
    handles.exportNucleiDataCheck.Layout.Row = 1;
    handles.exportNucleiDataCheck.Layout.Column = 1;
    
    handles.exportChannelDataCheck = uicheckbox(exportOptionsGrid, 'Text', 'Channel Data', ...
        'Value', true, 'Tooltip', 'Export detected signals and fit results for each channel');
    handles.exportChannelDataCheck.Layout.Row = 1;
    handles.exportChannelDataCheck.Layout.Column = 2;
    
    handles.exportClusteredDataCheck = uicheckbox(exportOptionsGrid, 'Text', 'Clustered Signal Data', ...
        'Value', true, 'Tooltip', 'Export comprehensive nucleus-signal associations (requires nuclei segmentation)');
    handles.exportClusteredDataCheck.Layout.Row = 2;
    handles.exportClusteredDataCheck.Layout.Column = 1;
    
    handles.exportVisualizationsCheck = uicheckbox(exportOptionsGrid, 'Text', 'Annotated Visualizations (PNG)', ...
        'Value', true, 'Tooltip', 'Export overview images showing nuclei boundaries, signal locations, and counts');
    handles.exportVisualizationsCheck.Layout.Row = 2;
    handles.exportVisualizationsCheck.Layout.Column = 2;
    
    % --- Row 5: Classifier Loading (optional) ---
    classifierPanel = uipanel(contentGrid, 'Title', 'Classifiers (Optional)', 'FontWeight', 'bold');
    classifierPanel.Layout.Row = 5;
    classifierGrid = uigridlayout(classifierPanel, [2, 4]);
    classifierGrid.RowHeight = {'fit', 'fit'};
    classifierGrid.ColumnWidth = {'fit', '1x', '1x', 'fit'};
    classifierGrid.Padding = [10 10 10 10];
    classifierGrid.RowSpacing = 8;
    classifierGrid.ColumnSpacing = 10;
    
    uilabel(classifierGrid, 'Text', 'Channel 1:', 'Layout', struct('Row', 1, 'Column', 1));
    handles.classifier1Edit = uieditfield(classifierGrid, 'Value', '', 'Editable', 'off', ...
        'Tooltip', 'Classifier for channel 1 (optional)');
    handles.classifier1Edit.Layout.Row = 1; handles.classifier1Edit.Layout.Column = 2;
    handles.classifier1BrowseBtn = uibutton(classifierGrid, 'Text', 'Browse...', ...
        'ButtonPushedFcn', @(src,evt) browseClassifier(fig, 1));
    handles.classifier1BrowseBtn.Layout.Row = 1; handles.classifier1BrowseBtn.Layout.Column = 3;
    handles.classifier1ClearBtn = uibutton(classifierGrid, 'Text', 'Clear', ...
        'ButtonPushedFcn', @(src,evt) clearClassifier(fig, 1));
    handles.classifier1ClearBtn.Layout.Row = 1; handles.classifier1ClearBtn.Layout.Column = 4;
    
    uilabel(classifierGrid, 'Text', 'Channel 2:', 'Layout', struct('Row', 2, 'Column', 1));
    handles.classifier2Edit = uieditfield(classifierGrid, 'Value', '', 'Editable', 'off', ...
        'Tooltip', 'Classifier for channel 2 (optional)');
    handles.classifier2Edit.Layout.Row = 2; handles.classifier2Edit.Layout.Column = 2;
    handles.classifier2BrowseBtn = uibutton(classifierGrid, 'Text', 'Browse...', ...
        'ButtonPushedFcn', @(src,evt) browseClassifier(fig, 2));
    handles.classifier2BrowseBtn.Layout.Row = 2; handles.classifier2BrowseBtn.Layout.Column = 3;
    handles.classifier2ClearBtn = uibutton(classifierGrid, 'Text', 'Clear', ...
        'ButtonPushedFcn', @(src,evt) clearClassifier(fig, 2));
    handles.classifier2ClearBtn.Layout.Row = 2; handles.classifier2ClearBtn.Layout.Column = 4;
    
    % Initialize classifier storage
    handles.classifiers = cell(1, 2);
    handles.classifierFeatures = cell(1, 2);
    handles.classifierFeatureInfo = cell(1, 2);
    
    % --- Row 6: Unused placeholder (was row 5) ---
    placeholderPanel6 = uipanel(contentGrid, 'BorderType', 'none');
    placeholderPanel6.Layout.Row = 6;
    
    % --- Row 7: Image Set Preview ---
    previewPanel = uipanel(contentGrid, 'Title', 'Image Sets Found', 'FontWeight', 'bold');
    previewPanel.Layout.Row = 7;
    previewGrid = uigridlayout(previewPanel, [2, 1]);
    previewGrid.RowHeight = {'fit', 'fit'};
    previewGrid.Padding = [10 10 10 10];
    previewGrid.RowSpacing = 8;
    
    handles.imageSetCountLabel = uilabel(previewGrid, 'Text', 'No directory selected', ...
        'FontWeight', 'normal', 'FontColor', [0.5 0.5 0.5]);
    handles.scanDirBtn = uibutton(previewGrid, 'Text', 'Scan Directory', ...
        'ButtonPushedFcn', @(src,evt) scanDirectory(fig), ...
        'Enable', 'off', 'Tooltip', 'Scan for image sets in selected directory');
    
    % --- Row 8: Progress Display + Instructions (side by side, 500px tall) ---
    progressAndInstructionsPanel = uipanel(contentGrid, 'BorderType', 'none');
    progressAndInstructionsPanel.Layout.Row = 8;
    progressAndInstructionsGrid = uigridlayout(progressAndInstructionsPanel, [1, 2]);
    progressAndInstructionsGrid.RowHeight = {'1x'};
    progressAndInstructionsGrid.ColumnWidth = {'1x', '1x'}; % 50% each
    progressAndInstructionsGrid.Padding = [0 0 0 0];
    progressAndInstructionsGrid.ColumnSpacing = 10;
    
    % Left: Progress Log
    progressPanel = uipanel(progressAndInstructionsGrid, 'Title', 'Progress Log', 'FontWeight', 'bold', 'FontSize', 10);
    progressPanel.Layout.Row = 1;
    progressPanel.Layout.Column = 1;
    progressGrid = uigridlayout(progressPanel, [3, 1]);
    progressGrid.RowHeight = {'fit', 'fit', '1x'};
    progressGrid.Padding = [8 8 8 8];
    progressGrid.RowSpacing = 3;
    
    handles.statusLabel = uilabel(progressGrid, 'Text', 'Ready to begin', ...
        'FontWeight', 'bold', 'FontSize', 10, 'FontColor', [0 0.6 0]);
    
    % Progress bar container (visual loading bar with smooth filling)
    progressBarContainer = uipanel(progressGrid, 'BorderType', 'line', 'BorderWidth', 1);
    handles.progressBarGrid = uigridlayout(progressBarContainer, [1, 2]);
    handles.progressBarGrid.RowHeight = {20};
    handles.progressBarGrid.ColumnWidth = {0, '1x'}; % Start: 0% green, 100% gray
    handles.progressBarGrid.Padding = [0 0 0 0];
    handles.progressBarGrid.ColumnSpacing = 0;
    
    % Green fill panel (grows as processing progresses)
    handles.progressBarFill = uipanel(handles.progressBarGrid, 'BorderType', 'none', ...
        'BackgroundColor', [0.2 0.7 0.3]);
    handles.progressBarFill.Layout.Row = 1;
    handles.progressBarFill.Layout.Column = 1;
    
    % Gray unfilled panel (shrinks as processing progresses)
    progressBarEmpty = uipanel(handles.progressBarGrid, 'BorderType', 'none', ...
        'BackgroundColor', [0.9 0.9 0.9]);
    progressBarEmpty.Layout.Row = 1;
    progressBarEmpty.Layout.Column = 2;
    
    % Log text area (tall to fill vertical space)
    handles.logArea = uitextarea(progressGrid, 'Value', {'Ready to process...'}, 'Editable', 'off', ...
        'FontName', 'Courier New', 'FontSize', 9, ...
        'WordWrap', 'off');  % Keep monospace alignment
    
    % Right: Instructions
    instructionsPanel = uipanel(progressAndInstructionsGrid, 'Title', 'Instructions & Info', ...
        'FontWeight', 'bold', 'FontSize', 10, 'BorderType', 'line', 'BackgroundColor', [0.95 0.95 1]);
    instructionsPanel.Layout.Row = 1;
    instructionsPanel.Layout.Column = 2;
    instructionsGrid = uigridlayout(instructionsPanel, [8, 1]);
    instructionsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    instructionsGrid.Padding = [10 10 10 10];
    instructionsGrid.RowSpacing = 8;
    
    lbl9 = uilabel(instructionsGrid, 'Text', 'How to Use:', 'FontWeight', 'bold', 'FontSize', 10);
    lbl10 = uilabel(instructionsGrid, 'Text', '1. Select input directory containing experiment subfolders (e.g., from SNAP_prepare)', ...
        'FontSize', 9, 'WordWrap', 'on');
    lbl11 = uilabel(instructionsGrid, 'Text', '2. Load parameter file from SNAP (lastUsedParameters.mat)', ...
        'FontSize', 9, 'WordWrap', 'on');
    lbl12 = uilabel(instructionsGrid, 'Text', '3. Select output directory and click "Start Processing"', ...
        'FontSize', 9, 'WordWrap', 'on');
    
    lbl13 = uilabel(instructionsGrid, 'Text', 'What Gets Exported:', 'FontWeight', 'bold', 'FontSize', 10);
    lbl14 = uilabel(instructionsGrid, 'Text', '• Per-subfolder: Data based on export options and available data', ...
        'FontSize', 9, 'WordWrap', 'on');
    lbl15 = uilabel(instructionsGrid, 'Text', '• Aggregate: Combined data from ALL subfolders with composition statistics', ...
        'FontSize', 9, 'WordWrap', 'on');
    
    lblSpacer = uilabel(instructionsGrid, 'Text', ''); % Spacer to push content to top
    
    % --- Row 9: Action Buttons ---
    actionPanel = uipanel(contentGrid, 'BorderType', 'none');
    actionPanel.Layout.Row = 9;
    actionGrid = uigridlayout(actionPanel, [1, 3]);
    actionGrid.ColumnWidth = {'1x', 150, 100};
    actionGrid.Padding = [0 0 0 0];
    actionGrid.ColumnSpacing = 10;
    
    lbl8 = uilabel(actionGrid, 'Text', ''); % Spacer
    
    handles.startBtn = uibutton(actionGrid, 'Text', 'Start Processing', ...
        'FontWeight', 'bold', 'FontSize', 12, ...
        'ButtonPushedFcn', @(src,evt) startBatchProcessing(fig), ...
        'Enable', 'off', 'BackgroundColor', [0.2 0.7 0.3], 'FontColor', [1 1 1]);
    
    handles.closeBtn = uibutton(actionGrid, 'Text', 'Close', ...
        'FontWeight', 'bold', 'FontSize', 12, ...
        'ButtonPushedFcn', @(src,evt) close(fig));
    
    % Store handles in figure
    guidata(fig, handles);
end

%% GUI Callback Functions

function browseInputDir(fig)
    handles = guidata(fig);
    selectedDir = uigetdir('', 'Select Input Directory with Image Subfolders');
    if selectedDir ~= 0
        handles.inputDirEdit.Value = selectedDir;
        guidata(fig, handles);
        updateUIState(fig);
        % Auto-scan if parameter file is also loaded
        if ~isempty(handles.paramFileEdit.Value)
            scanDirectory(fig);
        end
    end
end

function browseParamFile(fig)
    handles = guidata(fig);
    [file, path] = uigetfile('*.mat', 'Select Parameter File from SNAP');
    if file ~= 0
        handles.paramFileEdit.Value = fullfile(path, file);
        guidata(fig, handles);
        updateUIState(fig);
        % Auto-scan if input dir is also loaded
        if ~isempty(handles.inputDirEdit.Value)
            scanDirectory(fig);
        end
    end
end

function browseOutputDir(fig)
    handles = guidata(fig);
    selectedDir = uigetdir('', 'Select Output Directory');
    if selectedDir ~= 0
        handles.outputDirEdit.Value = selectedDir;
        guidata(fig, handles);
        updateUIState(fig);
    end
end

function browseClassifier(fig, channelIdx)
    handles = guidata(fig);
    [file, path] = uigetfile('*.mat', sprintf('Select Classifier for Channel %d', channelIdx));
    if file ~= 0
        filepath = fullfile(path, file);
        
        % Load and validate classifier with custom expressions and normParams
        [model, features, featureInfo, ~, fittingMethod, ~, success, customExpr, normParams] = ...
            snap_helpers.classification.loadClassifier(filepath);
        
        if success
            handles.classifiers{channelIdx} = model;
            handles.classifierFeatures{channelIdx} = features;
            handles.classifierFeatureInfo{channelIdx} = featureInfo;
            
            % Store custom expressions and norm params
            if ~isfield(handles, 'classifierCustomExpressions')
                handles.classifierCustomExpressions = cell(1, 2);
            end
            if ~isfield(handles, 'classifierNormParams')
                handles.classifierNormParams = cell(1, 2);
            end
            handles.classifierCustomExpressions{channelIdx} = customExpr;
            handles.classifierNormParams{channelIdx} = normParams;
            
            % Update UI
            nBase = numel(features);
            nCustom = numel(customExpr);
            if nCustom > 0
                displayStr = sprintf('%s (%d+%d feat, %s)', file, nBase, nCustom, fittingMethod);
            else
                displayStr = sprintf('%s (%d features, %s)', file, nBase, fittingMethod);
            end
            
            if channelIdx == 1
                handles.classifier1Edit.Value = displayStr;
            else
                handles.classifier2Edit.Value = displayStr;
            end
            
            guidata(fig, handles);
            appendLog(handles.progressLog, sprintf('Loaded classifier for Channel %d: %s', channelIdx, file));
        else
            uialert(fig, 'Failed to load classifier', 'Load Error');
        end
    end
end

function clearClassifier(fig, channelIdx)
    handles = guidata(fig);
    handles.classifiers{channelIdx} = [];
    handles.classifierFeatures{channelIdx} = {};
    handles.classifierFeatureInfo{channelIdx} = struct();
    
    if isfield(handles, 'classifierCustomExpressions')
        handles.classifierCustomExpressions{channelIdx} = struct('name', {}, 'expression', {});
    end
    if isfield(handles, 'classifierNormParams')
        handles.classifierNormParams{channelIdx} = struct('mu', [], 'sigma', [], 'standardized', false);
    end
    
    if channelIdx == 1
        handles.classifier1Edit.Value = '';
    else
        handles.classifier2Edit.Value = '';
    end
    
    guidata(fig, handles);
end

function updateUIState(fig)
    % Enable/disable buttons based on input state
    handles = guidata(fig);
    
    hasInputDir = ~isempty(handles.inputDirEdit.Value);
    hasParamFile = ~isempty(handles.paramFileEdit.Value);
    hasOutputDir = ~isempty(handles.outputDirEdit.Value);
    
    % Enable scan button if input dir is selected
    if hasInputDir
        handles.scanDirBtn.Enable = 'on';
    else
        handles.scanDirBtn.Enable = 'off';
    end
    
    % Enable start button if all required fields are selected
    if hasInputDir && hasParamFile && hasOutputDir
        handles.startBtn.Enable = 'on';
    else
        handles.startBtn.Enable = 'off';
    end
    
    guidata(fig, handles);
end

function scanDirectory(fig)
    % Scan the input directory and count image sets
    handles = guidata(fig);
    
    inputDir = handles.inputDirEdit.Value;
    if isempty(inputDir) || ~exist(inputDir, 'dir')
        handles.imageSetCountLabel.Text = 'Invalid directory';
        handles.imageSetCountLabel.FontColor = [0.8 0 0];
        guidata(fig, handles);
        return;
    end
    
    % Find all subfolders
    subfolderInfo = dir(inputDir);
    subfolders = subfolderInfo([subfolderInfo.isdir] & ~ismember({subfolderInfo.name}, {'.', '..'}));
    numSets = length(subfolders);
    
    if numSets == 0
        handles.imageSetCountLabel.Text = 'No subfolders found in directory';
        handles.imageSetCountLabel.FontColor = [0.8 0.5 0];
    else
        % Create a more compact display
        if numSets <= 5
            setNames = strjoin({subfolders.name}, ', ');
        else
            setNames = sprintf('%s, ... and %d more', strjoin({subfolders(1:3).name}, ', '), numSets - 3);
        end
        handles.imageSetCountLabel.Text = sprintf('Found %d image set(s): %s', numSets, setNames);
        handles.imageSetCountLabel.FontColor = [0 0.6 0];
    end
    
    guidata(fig, handles);
end

function startBatchProcessing(fig)
    % Main batch processing function
    handles = guidata(fig);
    
    % Get parameters from UI
    inputDir = handles.inputDirEdit.Value;
    paramFile = handles.paramFileEdit.Value;
    outputDir = handles.outputDirEdit.Value;
    % Get export options
    exportOptions = struct();
    exportOptions.nucleiData = handles.exportNucleiDataCheck.Value;
    exportOptions.channelData = handles.exportChannelDataCheck.Value;
    exportOptions.clusteredData = handles.exportClusteredDataCheck.Value;
    exportOptions.visualizations = handles.exportVisualizationsCheck.Value;
    
    % Validate inputs
    if isempty(inputDir) || ~exist(inputDir, 'dir')
        uialert(fig, 'Please select a valid input directory.', 'Invalid Input');
        return;
    end
    
    if isempty(paramFile) || ~exist(paramFile, 'file')
        uialert(fig, 'Please select a valid parameter file.', 'Invalid Parameter File');
        return;
    end
    
    if isempty(outputDir)
        uialert(fig, 'Please select an output directory.', 'Output Directory Required');
        return;
    end
    
    % Disable UI during processing
    handles.startBtn.Enable = 'off';
    handles.inputDirBrowseBtn.Enable = 'off';
    handles.paramFileBrowseBtn.Enable = 'off';
    handles.outputDirBrowseBtn.Enable = 'off';
    handles.scanDirBtn.Enable = 'off';
    handles.exportNucleiDataCheck.Enable = 'off';
    handles.exportChannelDataCheck.Enable = 'off';
    handles.exportClusteredDataCheck.Enable = 'off';
    handles.exportVisualizationsCheck.Enable = 'off';
    handles.statusLabel.Text = 'Processing...';
    handles.statusLabel.FontColor = [0 0 0.8];
    handles.logArea.Value = {...
        '═════════════════════════════════════', ...
        'SNAP BATCH PROCESSING STARTED', ...
        '═════════════════════════════════════', ...
        sprintf('Input: %s', inputDir), ...
        sprintf('Output: %s', outputDir), ...
        sprintf('Parameters: %s', paramFile), ...
        '═════════════════════════════════════', ...
        ''};
    
    % Reset progress bar to 0%
    if isfield(handles, 'progressBarGrid') && isvalid(handles.progressBarGrid)
        handles.progressBarGrid.ColumnWidth = {0, '1x'}; % 0% green, 100% gray
    end
    
    guidata(fig, handles);
    
    % Force UI update
    drawnow;
    
    % Build optional arguments
    optArgs = {};
    if ~isempty(outputDir)
        optArgs = [optArgs, {'OutputDir', outputDir}];
    end
    
    % Create progress handle for updates
    progressStruct = struct('fig', fig);
    optArgs = [optArgs, {'ProgressHandle', progressStruct}];
    
    % Add export options
    optArgs = [optArgs, {'ExportOptions', exportOptions}];
    
    % Add classifiers if loaded
    if isfield(handles, 'classifiers') && any(~cellfun(@isempty, handles.classifiers))
        optArgs = [optArgs, {'Classifiers', handles.classifiers}];
        optArgs = [optArgs, {'ClassifierFeatures', handles.classifierFeatures}];
        optArgs = [optArgs, {'ClassifierFeatureInfo', handles.classifierFeatureInfo}];
        if isfield(handles, 'classifierCustomExpressions')
            optArgs = [optArgs, {'ClassifierCustomExpressions', handles.classifierCustomExpressions}];
        end
        if isfield(handles, 'classifierNormParams')
            optArgs = [optArgs, {'ClassifierNormParams', handles.classifierNormParams}];
        end
    end
    
    % Run batch processing
    try
        SNAP_batch(inputDir, paramFile, optArgs{:});
        
        % Success - check if GUI still exists
        if isvalid(fig) && isgraphics(fig)
            try
                handles = guidata(fig);
                handles.statusLabel.Text = 'Batch processing completed successfully!';
                handles.statusLabel.FontColor = [0 0.6 0];
                
                % Fill progress bar completely to 100%
                if isfield(handles, 'progressBarGrid') && isvalid(handles.progressBarGrid)
                    handles.progressBarGrid.ColumnWidth = {'1x', 0}; % 100% green, 0% gray
                end
                
                guidata(fig, handles);
                
                % Show completion dialog
                uialert(fig, 'Batch processing completed successfully! Check the log for details.', ...
                    'Success', 'Icon', 'success');
            catch
                % GUI might have been closed during processing
            end
        end
        
    catch ME
        % Error handling - check if GUI still exists
        if isvalid(fig) && isgraphics(fig)
            try
                handles = guidata(fig);
                handles.statusLabel.Text = 'Error during processing';
                handles.statusLabel.FontColor = [0.8 0 0];
                
                logText = handles.logArea.Value;
                logText{end+1} = '';
                logText{end+1} = '=== ERROR ===';
                logText{end+1} = ME.message;
                logText{end+1} = 'Stack trace:';
                for i = 1:length(ME.stack)
                    logText{end+1} = sprintf('  %s (line %d)', ME.stack(i).name, ME.stack(i).line);
                end
                handles.logArea.Value = logText;
                guidata(fig, handles);
                
                uialert(fig, sprintf('Error during batch processing:\n\n%s', ME.message), ...
                    'Processing Error', 'Icon', 'error');
            catch
                % GUI closed
            end
        end
    end
    
    % Re-enable UI - check if GUI still exists
    if isvalid(fig) && isgraphics(fig)
        try
            handles = guidata(fig);
            handles.startBtn.Enable = 'on';
            handles.inputDirBrowseBtn.Enable = 'on';
            handles.paramFileBrowseBtn.Enable = 'on';
            handles.outputDirBrowseBtn.Enable = 'on';
            handles.scanDirBtn.Enable = 'on';
            handles.exportNucleiDataCheck.Enable = 'on';
            handles.exportChannelDataCheck.Enable = 'on';
            handles.exportClusteredDataCheck.Enable = 'on';
            handles.exportVisualizationsCheck.Enable = 'on';
            guidata(fig, handles);
        catch
            % GUI closed
        end
    end
end

function updateProgress(progressHandle, current, total, message)
    % Update progress bar and status
    if isempty(progressHandle) || ~isfield(progressHandle, 'fig')
        return;
    end
    
    % Check if figure is still valid
    if ~isvalid(progressHandle.fig) || ~isgraphics(progressHandle.fig)
        return;
    end
    
    try
        handles = guidata(progressHandle.fig);
        
        % Check if handles components still exist
        if ~isfield(handles, 'statusLabel') || ~isvalid(handles.statusLabel)
            return;
        end
        
        percentage = round((current / total) * 100);
        handles.statusLabel.Text = sprintf('[%d/%d - %d%%] %s', current, total, percentage, message);
        
        % Update progress bar (smooth filling)
        if isfield(handles, 'progressBarGrid') && isvalid(handles.progressBarGrid)
            % Calculate proportions for green (filled) vs gray (unfilled)
            greenPortion = sprintf('%dx', max(1, percentage)); % At least 1x to show some progress
            grayPortion = sprintf('%dx', max(1, 100 - percentage));
            
            if percentage == 0
                handles.progressBarGrid.ColumnWidth = {0, '1x'}; % All gray
            elseif percentage >= 100
                handles.progressBarGrid.ColumnWidth = {'1x', 0}; % All green
            else
                handles.progressBarGrid.ColumnWidth = {greenPortion, grayPortion};
            end
        end
        
        guidata(progressHandle.fig, handles);
        drawnow limitrate;  % Rate-limited for better performance
    catch
        % Silently fail if GUI is closed or invalid
    end
end

function appendLog(progressHandle, message)
    % Append message to log area (MINIMAL overhead version)
    if isempty(progressHandle) || ~isfield(progressHandle, 'fig')
        return;
    end
    
    % Check if figure is still valid
    if ~isvalid(progressHandle.fig) || ~isgraphics(progressHandle.fig)
        return;
    end
    
    try
        handles = guidata(progressHandle.fig);
        
        % Check if handles and logArea still exist
        if ~isfield(handles, 'logArea') || ~isvalid(handles.logArea)
            return;
        end
        
        logText = handles.logArea.Value;
        logText{end+1} = message;
        % Keep only last 100 lines
        if length(logText) > 100
            logText = logText(end-99:end);
        end
        handles.logArea.Value = logText;
        
        % ONLY scroll on major events - NO drawnow for performance
        % GUI will update naturally when MATLAB is idle
    catch
        % Silently fail if GUI is closed or invalid
    end
end

function emptyResult = createEmptyResult(setName, numChannels)
    % Create empty result structure with ALL fields for consistency
    % This ensures error/skipped results have the same structure as successful ones
    emptyResult = struct();
    emptyResult.name = setName;
    emptyResult.numNuclei = 0;
    emptyResult.totalSpots = 0;
    emptyResult.channelSpots = zeros(1, numChannels);
    emptyResult.nucleiMask = [];
    emptyResult.nucleusLabels = struct();
    emptyResult.channels = cell(1, numChannels);
    % Initialize all channel slots with consistent empty structures
    for ch = 1:numChannels
        emptyResult.channels{ch} = struct('channelIndex', ch, 'processedImage', [], ...
            'maximaCoords', zeros(0,3), 'fitResults', struct([]), 'numSpots', 0);
    end
    emptyResult.signalComposition = [];
    emptyResult.processingTime = 0;
    emptyResult.status = '';  % Will be set by caller
    emptyResult.exportedTypes = struct('nucleiData', false, 'channelData', false, 'clusteredData', false, 'visualizations', false);
end

