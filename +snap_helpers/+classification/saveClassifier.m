function success = saveClassifier(filepath, model, features, featureInfo, trainStats, fittingMethod, metadata, customExpressions, normParams)
% saveClassifier - Save trained classifier to file for SNAP integration
%
% ============================================================================
% STANDARDIZED CLASSIFIER EXPORT FORMAT WITH NORMALIZATION
% ============================================================================
%
% This function saves a trained classifier in a standardized format that
% can be loaded by SNAP, SNAP_batch, or SNAP_classify.
%
% USAGE:
%   success = snap_helpers.classification.saveClassifier(filepath, model, features, featureInfo, trainStats, fittingMethod)
%   success = snap_helpers.classification.saveClassifier(..., metadata, customExpressions, normParams)
%
% INPUTS:
%   filepath          - Full path to save file (should end in .mat)
%   model             - Trained SVM model (from trainClassifier)
%   features          - Cell array of BASE feature names used for training
%   featureInfo       - Struct with feature metadata (from getAvailableFeatures)
%   trainStats        - Struct with training statistics (from trainClassifier)
%   fittingMethod     - String: fitting method this classifier was trained for
%   metadata          - (optional) Additional metadata struct
%   customExpressions - (optional) Struct array with custom expression definitions
%   normParams        - (optional) Struct with z-score normalization parameters:
%                       .mu    - Feature means from training
%                       .sigma - Feature standard deviations from training
%                       .standardized - Whether normalization was applied
%
% OUTPUTS:
%   success - Logical: true if save succeeded
%
% FILE FORMAT:
%   The saved .mat file contains:
%   - model: Trained SVM model (compact form)
%   - features: Cell array of BASE feature names
%   - customExpressions: Struct array of custom expression definitions
%   - normParams: Z-score normalization parameters (CRITICAL for prediction)
%   - featureInfo: Feature metadata
%   - trainStats: Training statistics
%   - fittingMethod: Original fitting method
%   - metadata: Additional info
%   - version: Classifier format version
%
% ============================================================================

    success = false;
    
    % Validate inputs
    if isempty(model)
        warning('Cannot save empty model');
        return;
    end
    
    if isempty(features)
        warning('Cannot save classifier without feature list');
        return;
    end
    
    % Handle optional metadata
    if nargin < 7
        metadata = struct();
    end
    
    % Handle optional customExpressions
    if nargin < 8
        customExpressions = struct('name', {}, 'expression', {});
    end
    
    % Handle optional normParams
    if nargin < 9
        normParams = struct('mu', [], 'sigma', [], 'standardized', false);
    end
    
    % Build classifier export structure
    classifierExport = struct();
    classifierExport.model = model;
    classifierExport.features = features;
    classifierExport.customExpressions = customExpressions;
    classifierExport.normParams = normParams;  % Z-score normalization parameters
    classifierExport.featureInfo = featureInfo;
    classifierExport.trainStats = trainStats;
    classifierExport.fittingMethod = fittingMethod;
    
    % Add metadata
    classifierExport.metadata = metadata;
    classifierExport.metadata.exportDate = datetime('now');
    classifierExport.metadata.matlabVersion = version;
    
    % Add training summary for quick reference
    if isfield(trainStats, 'classCounts') && numel(trainStats.classCounts) >= 2
        classifierExport.metadata.nReal = trainStats.classCounts(end);
        classifierExport.metadata.nNoise = trainStats.classCounts(1);
    end
    if isfield(trainStats, 'cvAccuracy')
        classifierExport.metadata.accuracy = trainStats.cvAccuracy;
    elseif isfield(trainStats, 'trainAccuracy')
        classifierExport.metadata.accuracy = trainStats.trainAccuracy;
    end
    
    % Add feature counts to metadata
    classifierExport.metadata.nBaseFeatures = numel(features);
    classifierExport.metadata.nCustomFeatures = numel(customExpressions);
    classifierExport.metadata.zScoreNormalized = normParams.standardized;
    
    % Version for future compatibility
    classifierExport.version = '3.0';  % Updated for normalization params
    
    % Save
    try
        save(filepath, '-struct', 'classifierExport', '-v7.3');
        success = true;
        fprintf('Classifier saved to: %s\n', filepath);
        fprintf('  Base features: %d\n', numel(features));
        fprintf('  Custom expressions: %d\n', numel(customExpressions));
        fprintf('  Z-score normalized: %s\n', string(normParams.standardized));
        fprintf('  Fitting method: %s\n', fittingMethod);
        if isfield(classifierExport.metadata, 'accuracy')
            fprintf('  Accuracy: %.1f%%\n', classifierExport.metadata.accuracy * 100);
        end
    catch ME
        warning('Failed to save classifier: %s', ME.message);
        success = false;
    end
end

