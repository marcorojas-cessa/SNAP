function [X, featureNames, validMask, extractionInfo] = buildFeatureMatrix(fitResults, selectedFeatures, featureInfo, customExpressions)
% buildFeatureMatrix - Extract feature matrix from fit results for SVM
%
% ============================================================================
% GENERIC FEATURE EXTRACTION WITH CUSTOM EXPRESSION SUPPORT
% ============================================================================
%
% This function extracts a feature matrix from SNAP fit results. It supports:
%   1. Base features (selected from getAvailableFeatures)
%   2. Custom expressions (user-defined mathematical combinations)
%
% USAGE:
%   [X, names, valid] = snap_helpers.classification.buildFeatureMatrix(fitResults, selected, info)
%   [X, names, valid] = snap_helpers.classification.buildFeatureMatrix(fitResults, selected, info, customExpr)
%
% INPUTS:
%   fitResults        - Struct array or table with fit results
%   selectedFeatures  - Cell array of base feature names to include
%   featureInfo       - Struct with feature metadata (from getAvailableFeatures)
%   customExpressions - (optional) Struct array with:
%                       .name       - Name for the custom feature
%                       .expression - Mathematical expression string
%
% OUTPUTS:
%   X               - [N x F] numeric feature matrix (N spots, F features)
%   featureNames    - Cell array of feature names (column headers for X)
%   validMask       - [N x 1] logical array, true for rows without NaN
%   extractionInfo  - Struct with extraction statistics
%
% ============================================================================

    % Handle optional customExpressions
    if nargin < 4
        customExpressions = struct('name', {}, 'expression', {});
    end
    
    % Initialize outputs
    extractionInfo = struct();
    extractionInfo.warnings = {};
    extractionInfo.missingFeatures = {};
    extractionInfo.derivedComputed = {};
    extractionInfo.customComputed = {};
    
    % Handle empty input
    if isempty(fitResults)
        X = [];
        featureNames = selectedFeatures;
        validMask = [];
        extractionInfo.nSpots = 0;
        extractionInfo.nFeatures = numel(selectedFeatures);
        return;
    end
    
    % Also handle empty selectedFeatures
    if isempty(selectedFeatures) && isempty(customExpressions)
        nSpots = getNumSamples(fitResults);
        X = zeros(nSpots, 0);
        featureNames = {};
        validMask = true(nSpots, 1);
        extractionInfo.nSpots = nSpots;
        extractionInfo.nFeatures = 0;
        return;
    end
    
    % Determine input type and convert to consistent format
    if istable(fitResults)
        dataTable = fitResults;
        dataStruct = table2struct(fitResults);
    elseif isstruct(fitResults)
        dataStruct = fitResults(:);  % Ensure column vector
        try
            dataTable = struct2table(dataStruct, 'AsArray', true);
        catch
            dataTable = [];
        end
    else
        error('fitResults must be a struct array or table');
    end
    
    nSpots = numel(dataStruct);
    nBaseFeatures = numel(selectedFeatures);
    nCustomFeatures = numel(customExpressions);
    nTotalFeatures = nBaseFeatures + nCustomFeatures;
    
    extractionInfo.nSpots = nSpots;
    extractionInfo.nBaseFeatures = nBaseFeatures;
    extractionInfo.nCustomFeatures = nCustomFeatures;
    extractionInfo.nFeatures = nTotalFeatures;
    
    % Pre-allocate feature matrix
    X = nan(nSpots, nTotalFeatures);
    featureNames = cell(1, nTotalFeatures);
    
    % Get all available feature names for expression evaluation
    if ~isempty(featureInfo)
        allFeatureNames = fieldnames(featureInfo);
    else
        allFeatureNames = selectedFeatures;
    end
    
    % === PART 1: Extract base features ===
    for f = 1:nBaseFeatures
        fname = selectedFeatures{f};
        featureNames{f} = fname;
        
        % Check if this is a derived feature (needs computation)
        if ~isempty(featureInfo) && isfield(featureInfo, fname) && isfield(featureInfo.(fname), 'compute')
            % Derived feature - compute from raw data
            try
                computeFn = featureInfo.(fname).compute;
                
                if ~isempty(dataTable)
                    computed = computeFn(dataTable);
                else
                    computed = computeFn(dataStruct);
                end
                
                if isscalar(computed)
                    X(:, f) = repmat(computed, nSpots, 1);
                else
                    X(:, f) = computed(:);
                end
                extractionInfo.derivedComputed{end+1} = fname;
                
            catch ME
                extractionInfo.warnings{end+1} = sprintf('Could not compute derived feature %s: %s', fname, ME.message);
                X(:, f) = nan(nSpots, 1);
            end
            
        else
            % Direct feature - extract from data
            extracted = extractDirectFeature(dataStruct, dataTable, fname, nSpots);
            
            if isempty(extracted)
                extractionInfo.missingFeatures{end+1} = fname;
                extractionInfo.warnings{end+1} = sprintf('Feature not found: %s', fname);
                X(:, f) = nan(nSpots, 1);
            else
                X(:, f) = extracted;
            end
        end
    end
    
    % === PART 2: Evaluate custom expressions ===
    for e = 1:nCustomFeatures
        colIdx = nBaseFeatures + e;
        exprStruct = customExpressions(e);
        featureNames{colIdx} = exprStruct.name;
        
        try
            result = snap_helpers.classification.evaluateExpression(...
                exprStruct.expression, dataStruct, allFeatureNames);
            
            if numel(result) == nSpots
                X(:, colIdx) = result;
                extractionInfo.customComputed{end+1} = exprStruct.name;
            else
                extractionInfo.warnings{end+1} = sprintf('Expression "%s" returned wrong size', exprStruct.name);
                X(:, colIdx) = nan(nSpots, 1);
            end
            
        catch ME
            extractionInfo.warnings{end+1} = sprintf('Failed to evaluate "%s": %s', exprStruct.name, ME.message);
            X(:, colIdx) = nan(nSpots, 1);
        end
    end
    
    % Identify valid rows (no NaN values across all selected features)
    validMask = all(~isnan(X), 2);
    
    % Report statistics
    nInvalid = sum(~validMask);
    extractionInfo.nValid = sum(validMask);
    extractionInfo.nInvalid = nInvalid;
    
    if nInvalid > 0
        extractionInfo.warnings{end+1} = sprintf('%d/%d spots have NaN features and will be excluded', nInvalid, nSpots);
    end
end

% ============================================================================
% HELPER FUNCTIONS
% ============================================================================

function n = getNumSamples(data)
    if istable(data)
        n = height(data);
    elseif isstruct(data)
        n = numel(data);
    else
        n = 0;
    end
end

function values = extractDirectFeature(dataStruct, dataTable, fname, nSpots)
    values = [];
    
    % Special handling for coordinate fields
    if strcmp(fname, 'x')
        if isfield(dataStruct, 'globalFitCenter')
            values = arrayfun(@(s) safeGetCoord(s.globalFitCenter, 2), dataStruct);
        elseif isfield(dataStruct, 'fitted_coords')
            values = arrayfun(@(s) safeGetCoord(s.fitted_coords, 2), dataStruct);
        elseif isfield(dataStruct, 'center_y')
            values = [dataStruct.center_y]';
        end
        return;
    elseif strcmp(fname, 'y')
        if isfield(dataStruct, 'globalFitCenter')
            values = arrayfun(@(s) safeGetCoord(s.globalFitCenter, 1), dataStruct);
        elseif isfield(dataStruct, 'fitted_coords')
            values = arrayfun(@(s) safeGetCoord(s.fitted_coords, 1), dataStruct);
        elseif isfield(dataStruct, 'center_x')
            values = [dataStruct.center_x]';
        end
        return;
    elseif strcmp(fname, 'z')
        if isfield(dataStruct, 'globalFitCenter')
            values = arrayfun(@(s) safeGetCoord(s.globalFitCenter, 3), dataStruct);
        elseif isfield(dataStruct, 'fitted_coords')
            values = arrayfun(@(s) safeGetCoord(s.fitted_coords, 3), dataStruct);
        elseif isfield(dataStruct, 'center_z')
            values = [dataStruct.center_z]';
        end
        return;
    end
    
    % Get alternative field names
    altNames = getAlternativeFieldNames(fname);
    allNames = [{fname}, altNames];
    
    % Try to extract from struct
    for i = 1:numel(allNames)
        tryName = allNames{i};
        if isfield(dataStruct, tryName)
            try
                rawValues = [dataStruct.(tryName)];
                if numel(rawValues) == nSpots
                    values = rawValues(:);
                    return;
                elseif numel(rawValues) > nSpots
                    values = nan(nSpots, 1);
                    for j = 1:nSpots
                        val = dataStruct(j).(tryName);
                        if ~isempty(val) && isnumeric(val)
                            values(j) = mean(val(~isnan(val)));
                        end
                    end
                    return;
                end
            catch
            end
        end
    end
    
    % Try to extract from table
    if ~isempty(dataTable)
        for i = 1:numel(allNames)
            tryName = allNames{i};
            if ismember(tryName, dataTable.Properties.VariableNames)
                values = dataTable.(tryName);
                if numel(values) == nSpots
                    return;
                end
            end
        end
    end
    
    values = [];
end

function coord = safeGetCoord(coordArray, idx)
    if isempty(coordArray) || numel(coordArray) < idx
        coord = NaN;
    else
        coord = coordArray(idx);
    end
end

function altNames = getAlternativeFieldNames(fname)
    altNames = {};
    
    switch lower(fname)
        case 'integrated_intensity'
            altNames = {'integratedIntensity', 'IntegratedIntensity', 'intensity'};
        case 'r_squared'
            altNames = {'rSquared', 'rsquared', 'R_squared', 'Rsquared'};
        case 'radial_symmetry_score'
            altNames = {'radialSymmetryScore', 'radialSymmetryQuality', 'symmetry_score'};
        case 'amplitude_over_background'
            altNames = {'amplitudeOverBackground', 'amp_bg_ratio'};
        case 'sigma_xy_ratio'
            altNames = {'sigmaXYRatio', 'ellipticity', 'aspect_ratio'};
        case 'sigma_sum'
            altNames = {'sigmaSum', 'total_sigma'};
        case 'sigma_product'
            altNames = {'sigmaProduct', 'sigma_area'};
        case 'sigma_z_ratio'
            altNames = {'sigmaZRatio', 'z_elongation'};
    end
    
    camelCase = strrep(fname, '_', '');
    if ~strcmp(camelCase, fname)
        altNames{end+1} = camelCase;
    end
    
    camelCaseUpper = [upper(camelCase(1)), camelCase(2:end)];
    altNames{end+1} = camelCaseUpper;
end
