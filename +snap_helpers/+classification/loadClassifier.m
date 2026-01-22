function [model, features, featureInfo, trainStats, fittingMethod, metadata, success, customExpressions, normParams] = loadClassifier(filepath)
% loadClassifier - Load trained classifier from file with normalization parameters
%
% ============================================================================
% LOAD STANDARDIZED CLASSIFIER FORMAT
% ============================================================================
%
% This function loads a classifier saved by saveClassifier and validates
% its format for use in SNAP, SNAP_batch, or SNAP_classify.
%
% USAGE:
%   [model, features, info, stats, method, meta, ok, expr, norm] = snap_helpers.classification.loadClassifier(filepath)
%
% INPUTS:
%   filepath - Full path to classifier .mat file
%
% OUTPUTS:
%   model             - Trained SVM model (compact form)
%   features          - Cell array of BASE feature names
%   featureInfo       - Feature metadata struct
%   trainStats        - Training statistics struct
%   fittingMethod     - String: fitting method classifier was trained for
%   metadata          - Additional metadata struct
%   success           - Logical: true if load succeeded and format is valid
%   customExpressions - Struct array with custom expression definitions
%   normParams        - Struct with z-score normalization parameters:
%                       .mu    - Feature means from training
%                       .sigma - Feature standard deviations from training
%                       .standardized - Whether normalization was applied
%
% ============================================================================

    % Initialize outputs
    model = [];
    features = {};
    featureInfo = struct();
    trainStats = struct();
    fittingMethod = '';
    metadata = struct();
    success = false;
    customExpressions = struct('name', {}, 'expression', {});
    normParams = struct('mu', [], 'sigma', [], 'standardized', false);
    
    % Check file exists
    if ~isfile(filepath)
        warning('Classifier file not found: %s', filepath);
        return;
    end
    
    % Load file
    try
        data = load(filepath);
    catch ME
        warning('Failed to load classifier file: %s', ME.message);
        return;
    end
    
    % Validate required fields
    requiredFields = {'model', 'features'};
    for i = 1:numel(requiredFields)
        if ~isfield(data, requiredFields{i})
            warning('Classifier file missing required field: %s', requiredFields{i});
            return;
        end
    end
    
    % Extract data
    model = data.model;
    features = data.features;
    
    % Validate model
    if isempty(model)
        warning('Classifier file contains empty model');
        return;
    end
    
    % Validate features
    if isempty(features) || ~iscell(features)
        warning('Classifier file contains invalid features list');
        return;
    end
    
    % Extract optional fields with defaults
    if isfield(data, 'featureInfo')
        featureInfo = data.featureInfo;
    end
    
    if isfield(data, 'trainStats')
        trainStats = data.trainStats;
    end
    
    if isfield(data, 'fittingMethod')
        fittingMethod = data.fittingMethod;
    else
        fittingMethod = 'Unknown';
    end
    
    if isfield(data, 'metadata')
        metadata = data.metadata;
    end
    
    % Extract custom expressions (version 2.0+)
    if isfield(data, 'customExpressions')
        customExpressions = data.customExpressions;
        if isempty(customExpressions)
            customExpressions = struct('name', {}, 'expression', {});
        end
    end
    
    % Extract normalization parameters (version 3.0+)
    if isfield(data, 'normParams')
        normParams = data.normParams;
    else
        % Backwards compatibility: check trainStats for normalization info
        if isfield(trainStats, 'featureMeans') && isfield(trainStats, 'featureStds')
            normParams.mu = trainStats.featureMeans;
            normParams.sigma = trainStats.featureStds;
            normParams.standardized = true;
        end
    end
    
    % Check version compatibility
    if isfield(data, 'version')
        metadata.fileVersion = data.version;
    else
        metadata.fileVersion = '1.0';  % Pre-versioning format
    end
    
    success = true;
    
    % Display info
    fprintf('Loaded classifier from: %s\n', filepath);
    fprintf('  Version: %s\n', metadata.fileVersion);
    fprintf('  Base features: %d', numel(features));
    if ~isempty(features)
        fprintf(' (%s', features{1});
        if numel(features) > 1
            fprintf(', %s', features{2});
        end
        if numel(features) > 2
            fprintf(', ...');
        end
        fprintf(')');
    end
    fprintf('\n');
    
    if numel(customExpressions) > 0
        fprintf('  Custom expressions: %d (', numel(customExpressions));
        for i = 1:min(2, numel(customExpressions))
            if i > 1, fprintf(', '); end
            fprintf('%s', customExpressions(i).name);
        end
        if numel(customExpressions) > 2
            fprintf(', ...');
        end
        fprintf(')\n');
    end
    
    fprintf('  Z-score normalized: %s\n', string(normParams.standardized));
    fprintf('  Fitting method: %s\n', fittingMethod);
    
    if isfield(trainStats, 'cvAccuracy')
        fprintf('  CV Accuracy: %.1f%%\n', trainStats.cvAccuracy * 100);
    elseif isfield(trainStats, 'trainAccuracy')
        fprintf('  Train Accuracy: %.1f%%\n', trainStats.trainAccuracy * 100);
    end
    
    if isfield(trainStats, 'options')
        opts = trainStats.options;
        if isfield(opts, 'kernelFunction')
            fprintf('  Kernel: %s\n', opts.kernelFunction);
        end
        if isfield(opts, 'boxConstraint')
            fprintf('  Box Constraint: %.4g\n', opts.boxConstraint);
        end
    end
    
    if isfield(metadata, 'exportDate')
        fprintf('  Trained: %s\n', string(metadata.exportDate));
    end
end
