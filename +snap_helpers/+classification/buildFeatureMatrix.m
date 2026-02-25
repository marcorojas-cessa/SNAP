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
    selectedFeatures = normalizeFeatureList(selectedFeatures);
    customExpressions = normalizeCustomExpressionListLocal(customExpressions);
    customExpressions = orderCustomExpressionsByDependencies(customExpressions);

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
    
    [dataStruct, modelStatsInfo] = maybeAugmentModelStatsFeaturesIfNeeded( ...
        dataStruct, selectedFeatures, customExpressions);
    extractionInfo.modelStats = modelStatsInfo;
    if ~isempty(modelStatsInfo.warnings)
        extractionInfo.warnings = [extractionInfo.warnings, modelStatsInfo.warnings];
    end

    if ~isempty(dataStruct)
        try
            dataTable = struct2table(dataStruct, 'AsArray', true);
        catch
            dataTable = [];
        end
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
    
    % Get all available feature names for expression evaluation.
    % Include custom-expression names so already-computed custom values can
    % be reused by downstream expressions (dependency chaining).
    if ~isempty(featureInfo)
        allFeatureNames = fieldnames(featureInfo);
    else
        allFeatureNames = selectedFeatures;
    end
    if ~isempty(customExpressions)
        allFeatureNames = unique([allFeatureNames(:); {customExpressions.name}'], 'stable');
    else
        allFeatureNames = allFeatureNames(:);
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

        % Cache all selected base features (including derived) so custom
        % expressions can reference them reliably.
        dataStruct = cacheCustomFeatureForExpressionReuse(dataStruct, fname, X(:, f));
    end

    [dataStruct, supplementalDerived, supplementalWarnings] = ...
        computeSupplementalDerivedFeaturesForExpressions( ...
            dataStruct, dataTable, featureInfo, selectedFeatures, customExpressions, nSpots);
    if ~isempty(supplementalDerived)
        extractionInfo.derivedComputed = unique([extractionInfo.derivedComputed, supplementalDerived], 'stable');
    end
    if ~isempty(supplementalWarnings)
        extractionInfo.warnings = [extractionInfo.warnings, supplementalWarnings];
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
                dataStruct = cacheCustomFeatureForExpressionReuse(dataStruct, exprStruct.name, result);
            else
                extractionInfo.warnings{end+1} = sprintf('Expression "%s" returned wrong size', exprStruct.name);
                X(:, colIdx) = nan(nSpots, 1);
            end
            
        catch ME
            extractionInfo.warnings{end+1} = sprintf('Failed to evaluate "%s": %s', exprStruct.name, ME.message);
            X(:, colIdx) = nan(nSpots, 1);
        end
    end

    % Complex-value diagnostics and sanitization.
    % SVM training requires real-valued feature matrices.
    if ~isempty(X)
        imagPart = imag(X);
        complexMask = isfinite(imagPart) & (abs(imagPart) > 1e-12);
        complexCountByFeature = sum(complexMask, 1);
        extractionInfo.complexCountByFeature = complexCountByFeature;
        extractionInfo.featuresWithComplex = featureNames(complexCountByFeature > 0);

        if any(complexCountByFeature > 0)
            for i = 1:numel(featureNames)
                if complexCountByFeature(i) > 0
                    extractionInfo.warnings{end+1} = sprintf( ...
                        'Feature "%s" produced %d complex value(s); marking as NaN for compatibility', ...
                        featureNames{i}, complexCountByFeature(i)); %#ok<AGROW>
                end
            end
            X(complexMask) = NaN;
        end
        X = real(X);

        nonFiniteMask = ~isfinite(X);
        if any(nonFiniteMask(:))
            X(nonFiniteMask) = NaN;
        end
    else
        extractionInfo.complexCountByFeature = zeros(1, nTotalFeatures);
        extractionInfo.featuresWithComplex = {};
    end
    
    % Feature-level NaN diagnostics (helps identify incompatible selections)
    if ~isempty(X)
        nanCountByFeature = sum(isnan(X), 1);
    else
        nanCountByFeature = zeros(1, nTotalFeatures);
    end
    extractionInfo.nanCountByFeature = nanCountByFeature;
    extractionInfo.featureNames = featureNames;
    extractionInfo.featuresAllNaN = featureNames(nanCountByFeature == nSpots);
    extractionInfo.featuresWithAnyNaN = featureNames(nanCountByFeature > 0);

    if ~isempty(extractionInfo.featuresAllNaN)
        for i = 1:numel(extractionInfo.featuresAllNaN)
            extractionInfo.warnings{end+1} = sprintf( ...
                'Feature "%s" is NaN for all %d samples', ...
                extractionInfo.featuresAllNaN{i}, nSpots); %#ok<AGROW>
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
        case 'background'
            altNames = {'background_value', 'backgroundValue', 'bg'};
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

function [fitResultsOut, info] = maybeAugmentModelStatsFeaturesIfNeeded(fitResultsIn, selectedFeatures, customExpressions)
    fitResultsOut = fitResultsIn;
    info = struct( ...
        'requested', false, ...
        'requiredFeatures', {{}}, ...
        'augmented', false, ...
        'summary', struct(), ...
        'warnings', {{}});

    if isempty(fitResultsIn) || ~isstruct(fitResultsIn)
        return;
    end

    req = snap_helpers.classification.resolveModelStatsRequirements(selectedFeatures, customExpressions);
    info.requested = req.required;
    info.requiredFeatures = req.requiredFeatures;
    if ~req.required
        return;
    end

    required = req.requiredFeatures(:)';
    if isempty(required)
        return;
    end

    missingMask = ~isfield(fitResultsIn, required);
    if ~any(missingMask)
        return;
    end

    missing = required(missingMask);
    if exist('snap_contrib.svm.augmentFitResultsWithModelStats', 'file') ~= 2
        info.warnings{end+1} = sprintf([ ...
            'Model-stat features were requested but augmentation helper is unavailable. ', ...
            'Missing: %s'], strjoin(missing, ', ')); %#ok<AGROW>
        return;
    end

    try
        [fitResultsOut, summary] = snap_contrib.svm.augmentFitResultsWithModelStats( ...
            fitResultsIn, 'ComputeNormalityPValue', false, 'Verbose', false);
        info.augmented = true;
        info.summary = summary;

        stillMissingMask = ~isfield(fitResultsOut, required);
        if any(stillMissingMask)
            info.warnings{end+1} = sprintf([ ...
                'Model-stat augmentation completed but required fields are still missing: %s'], ...
                strjoin(required(stillMissingMask), ', ')); %#ok<AGROW>
        end

        if isstruct(summary) && isfield(summary, 'nComputed') && summary.nComputed == 0
            info.warnings{end+1} = ['Model-stat augmentation ran but produced zero computed rows. ', ...
                'Check that fit results include raw data windows.']; %#ok<AGROW>
        end
    catch ME
        info.warnings{end+1} = sprintf('Model-stat augmentation failed: %s', ME.message); %#ok<AGROW>
    end
end

function out = orderCustomExpressionsByDependencies(in)
    out = in;
    if isempty(in)
        return;
    end

    nExpr = numel(in);
    exprNames = cell(1, nExpr);
    for i = 1:nExpr
        exprNames{i} = strtrim(char(string(localStructField(in(i), 'name', ''))));
    end

    orderedIdx = zeros(1, nExpr);
    nOrdered = 0;
    added = false(1, nExpr);

    for pass = 1:nExpr
        progressed = false;
        knownNames = exprNames(added);
        for i = 1:nExpr
            if added(i)
                continue;
            end
            deps = getCustomExprDependencies(in(i), exprNames);
            if all(ismember(deps, knownNames))
                nOrdered = nOrdered + 1;
                orderedIdx(nOrdered) = i;
                added(i) = true;
                progressed = true;
            end
        end
        if ~progressed
            break;
        end
    end

    if nOrdered < nExpr
        remaining = find(~added);
        orderedIdx(nOrdered + (1:numel(remaining))) = remaining;
        nOrdered = nExpr;
    end

    out = in(orderedIdx(1:nOrdered));
end

function deps = getCustomExprDependencies(exprStruct, allCustomNames)
    deps = {};

    req = localStructField(exprStruct, 'requiredFeatures', {});
    if ischar(req)
        req = {req};
    elseif isstring(req)
        req = cellstr(req(:)');
    elseif ~iscell(req)
        req = {};
    end
    if ~isempty(req)
        req = cellfun(@(c) strtrim(char(string(c))), req, 'UniformOutput', false);
        deps = intersect(req, allCustomNames, 'stable');
        deps = setdiff(deps, {strtrim(char(string(localStructField(exprStruct, 'name', ''))))}, 'stable');
        return;
    end

    expr = strtrim(char(string(localStructField(exprStruct, 'expression', ''))));
    if isempty(expr)
        return;
    end

    tokens = regexp(expr, '[A-Za-z]\w*', 'match');
    if isempty(tokens)
        return;
    end
    tokens = unique(tokens, 'stable');
    tokens = setdiff(tokens, matlabReservedMathTokens(), 'stable');
    deps = intersect(tokens, allCustomNames, 'stable');
    deps = setdiff(deps, {strtrim(char(string(localStructField(exprStruct, 'name', ''))))}, 'stable');
end

function tokens = matlabReservedMathTokens()
    tokens = { ...
        'abs', 'acos', 'acosh', 'asin', 'asinh', 'atan', 'atan2', 'atanh', ...
        'ceil', 'cos', 'cosh', 'eps', 'exp', 'floor', 'isfinite', ...
        'isinf', 'isnan', 'log', 'log10', 'log2', 'max', 'mean', ...
        'median', 'min', 'mod', 'pi', 'pow2', 'real', 'round', ...
        'sign', 'sin', 'sinh', 'sqrt', 'std', 'sum', 'tan', 'tanh', ...
        'true', 'false', 'workspace'};
end

function dataStructOut = cacheCustomFeatureForExpressionReuse(dataStructIn, name, values)
    dataStructOut = dataStructIn;
    if isempty(dataStructIn) || isempty(name)
        return;
    end

    fieldName = strtrim(char(string(name)));
    if ~isvarname(fieldName)
        return;
    end

    n = numel(dataStructOut);
    vals = values(:);
    if numel(vals) ~= n
        return;
    end

    for i = 1:n
        dataStructOut(i).(fieldName) = vals(i); %#ok<AGROW>
    end
end

function [dataStructOut, computedNames, warningsOut] = computeSupplementalDerivedFeaturesForExpressions( ...
    dataStructIn, dataTable, featureInfo, selectedFeatures, customExpressions, nSpots)
    dataStructOut = dataStructIn;
    computedNames = {};
    warningsOut = {};

    if isempty(featureInfo) || ~isstruct(featureInfo) || isempty(customExpressions)
        return;
    end

    allInfoNames = fieldnames(featureInfo);
    if isempty(allInfoNames)
        return;
    end

    referenced = collectExpressionInputFeatureNames(customExpressions);
    if isempty(referenced)
        return;
    end

    referenced = intersect(referenced, allInfoNames, 'stable');
    if isempty(referenced)
        return;
    end

    target = setdiff(referenced, selectedFeatures, 'stable');
    for i = 1:numel(target)
        fname = target{i};
        if ~isfield(featureInfo, fname) || ~isfield(featureInfo.(fname), 'compute')
            continue;
        end

        try
            computeFn = featureInfo.(fname).compute;
            if ~isempty(dataTable)
                computed = computeFn(dataTable);
            else
                computed = computeFn(dataStructOut);
            end

            vals = normalizeComputedFeatureVector(computed, nSpots);
            dataStructOut = cacheCustomFeatureForExpressionReuse(dataStructOut, fname, vals);
            computedNames{end+1} = fname; %#ok<AGROW>
        catch ME
            warningsOut{end+1} = sprintf('Could not compute expression dependency %s: %s', fname, ME.message); %#ok<AGROW>
        end
    end
end

function vals = normalizeComputedFeatureVector(raw, nSpots)
    if isscalar(raw)
        vals = repmat(double(raw), nSpots, 1);
        return;
    end

    vals = raw(:);
    if numel(vals) < nSpots
        vals(end+1:nSpots, 1) = NaN;
    elseif numel(vals) > nSpots
        vals = vals(1:nSpots);
    end
    vals = double(vals);
    vals(~isfinite(vals)) = NaN;
end

function referenced = collectExpressionInputFeatureNames(customExpressions)
    referenced = {};
    if isempty(customExpressions)
        return;
    end

    for i = 1:numel(customExpressions)
        req = localStructField(customExpressions(i), 'requiredFeatures', {});
        if ischar(req)
            req = {req};
        elseif isstring(req)
            req = cellstr(req(:)');
        elseif ~iscell(req)
            req = {};
        end

        req = cellfun(@(c) strtrim(char(string(c))), req, 'UniformOutput', false);
        req = req(~cellfun(@isempty, req));
        referenced = [referenced, req]; %#ok<AGROW>

        expr = strtrim(char(string(localStructField(customExpressions(i), 'expression', ''))));
        if ~isempty(expr)
            tokens = regexp(expr, '[A-Za-z]\w*', 'match');
            if ~isempty(tokens)
                tokens = setdiff(unique(tokens, 'stable'), matlabReservedMathTokens(), 'stable');
                referenced = [referenced, tokens]; %#ok<AGROW>
            end
        end
    end

    referenced = unique(referenced, 'stable');
end

function out = normalizeFeatureList(in)
    if nargin < 1 || isempty(in)
        out = {};
        return;
    end

    if ischar(in)
        out = {strtrim(char(in))};
    elseif isstring(in)
        out = cellstr(in(:)');
    elseif iscell(in)
        out = in(:)';
    else
        out = {};
    end

    clean = {};
    for i = 1:numel(out)
        if ischar(out{i}) || isstring(out{i})
            name = strtrim(char(string(out{i})));
            if ~isempty(name)
                clean{end+1} = name; %#ok<AGROW>
            end
        end
    end
    out = unique(clean, 'stable');
end

function out = normalizeCustomExpressionListLocal(in)
    out = struct('name', {}, 'expression', {}, 'requiredFeatures', {});
    if nargin < 1 || isempty(in)
        return;
    end

    if istable(in)
        in = table2struct(in);
    end
    if ~isstruct(in)
        return;
    end

    for i = 1:numel(in)
        name = strtrim(char(string(localStructField(in(i), 'name', ''))));
        expr = strtrim(char(string(localStructField(in(i), 'expression', ''))));
        if isempty(name) || isempty(expr)
            continue;
        end
        req = localStructField(in(i), 'requiredFeatures', {});
        if ischar(req)
            req = {req};
        elseif isstring(req)
            req = cellstr(req(:)');
        elseif ~iscell(req)
            req = {};
        end
        out(end+1) = struct( ... %#ok<AGROW>
            'name', name, ...
            'expression', expr, ...
            'requiredFeatures', {req});
    end
end

function value = localStructField(s, fieldName, defaultValue)
    value = defaultValue;
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    end
end
