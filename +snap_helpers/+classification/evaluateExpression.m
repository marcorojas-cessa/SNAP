function result = evaluateExpression(expression, data, availableFeatures)
% evaluateExpression - Evaluate a mathematical expression using feature data
%
% ============================================================================
% EXPRESSION EVALUATOR FOR CUSTOM CLASSIFICATION FEATURES
% ============================================================================
%
% Allows users to define custom features as mathematical expressions
% combining base features. Supports standard MATLAB math operations.
%
% SUPPORTED OPERATIONS:
%   Arithmetic:  +, -, *, /, ^, ()
%   Functions:   log, log10, log2, sqrt, abs, exp, sin, cos, tan
%                min, max, mean, std, sign, floor, ceil, round
%   Comparison:  <, >, <=, >=, ==, ~=  (returns 0 or 1)
%   Logical:     &, |, ~
%
% EXAMPLES:
%   'amplitude / background'
%   'log(integrated_intensity / background)'
%   'sigma_x / sigma_y'
%   'sqrt(sigma_x^2 + sigma_y^2)'
%   '(amplitude - background) / sqrt(background)'
%   'log10(snr + 1)'
%
% USAGE:
%   result = snap_helpers.classification.evaluateExpression(expr, data, features)
%
% INPUTS:
%   expression        - String: mathematical expression using feature names
%   data              - Struct array or table of feature data
%   availableFeatures - Cell array of valid feature names
%
% OUTPUTS:
%   result - [N x 1] numeric array of computed values
%
% ============================================================================

    % Validate inputs
    if isempty(expression) || ~ischar(expression) && ~isstring(expression)
        error('Expression must be a non-empty string');
    end
    expression = char(expression);
    
    % Get number of samples
    if istable(data)
        nSamples = height(data);
    elseif isstruct(data)
        nSamples = numel(data);
    else
        error('Data must be a struct array or table');
    end
    
    % Initialize result
    result = nan(nSamples, 1);
    
    % Extract all feature values and create a workspace
    workspace = struct();
    
    % Add all available base features to workspace (with robust alias handling)
    for i = 1:numel(availableFeatures)
        fname = char(string(availableFeatures{i}));
        safeName = makeValidVarName(fname);

        try
            [values, found] = extractFeatureWithAliases(data, fname, nSamples);
            if found
                workspace.(safeName) = values(:);
            end
        catch
            % Skip features that can't be extracted
        end
    end
    
    % Also add common derived values that might be useful
    % These are computed on-the-fly if the base features exist
    try
        if isfield(workspace, 'integrated_intensity') && isfield(workspace, 'background')
            workspace.snr = workspace.integrated_intensity ./ max(workspace.background, 1);
        end
        if isfield(workspace, 'amplitude') && isfield(workspace, 'background')
            workspace.amp_bg = workspace.amplitude ./ max(workspace.background, 1);
        end
        if isfield(workspace, 'sigma_x') && isfield(workspace, 'sigma_y')
            workspace.sigma_xy_ratio = workspace.sigma_x ./ max(workspace.sigma_y, 0.01);
            workspace.sigma_xy_product = workspace.sigma_x .* workspace.sigma_y;
            workspace.sigma_xy_sum = workspace.sigma_x + workspace.sigma_y;
        end
    catch
        % Ignore errors in derived computations
    end
    
    % Replace feature names in expression with workspace references
    evalExpr = expression;
    
    % Sort feature names by length (longest first) to avoid partial replacements
    allNames = fieldnames(workspace);
    [~, sortIdx] = sort(cellfun(@length, allNames), 'descend');
    allNames = allNames(sortIdx);
    
    % Replace each feature name with workspace.featurename
    for i = 1:numel(allNames)
        fname = allNames{i};
        % Use word boundary matching to avoid partial replacements
        pattern = ['(?<![a-zA-Z0-9_])' fname '(?![a-zA-Z0-9_])'];
        evalExpr = regexprep(evalExpr, pattern, ['workspace.' fname]);
    end
    
    % Evaluate the expression
    try
        result = eval(evalExpr);
        
        % Ensure result is a column vector
        result = result(:);
        
        % Ensure result has correct length
        if numel(result) == 1
            result = repmat(result, nSamples, 1);
        elseif numel(result) ~= nSamples
            error('Expression result size mismatch');
        end
        
    catch ME
        warning('Failed to evaluate expression "%s": %s', expression, ME.message);
        result = nan(nSamples, 1);
    end
end

function values = extractStructField(data, fname)
    % Extract field values from struct array
    nSamples = numel(data);
    values = nan(nSamples, 1);
    if ~isstruct(data) || nSamples == 0
        return;
    end

    for i = 1:nSamples
        if ~isfield(data(i), fname)
            continue;
        end
        raw = data(i).(fname);
        values(i) = scalarizeNumericValue(raw);
    end
end

function safeName = makeValidVarName(name)
    % Convert feature name to valid MATLAB variable name
    safeName = name;
    safeName = strrep(safeName, '-', '_');
    safeName = strrep(safeName, ' ', '_');
    if ~isvarname(safeName)
        safeName = matlab.lang.makeValidName(safeName);
    end
end

function [values, found] = extractFeatureWithAliases(data, fname, nSamples)
    values = nan(nSamples, 1);
    found = false;

    aliases = [{fname}, getFeatureAliases(fname)];

    % Coordinate shortcuts
    if strcmp(fname, 'x') || strcmp(fname, 'y') || strcmp(fname, 'z')
        [coordValues, coordFound] = extractCoordinateFeature(data, fname, nSamples);
        if coordFound
            values = coordValues;
            found = true;
            return;
        end
    end

    if istable(data)
        for i = 1:numel(aliases)
            alias = aliases{i};
            if ismember(alias, data.Properties.VariableNames)
                values = tableColumnToNumericVector(data.(alias), nSamples);
                found = true;
                return;
            end
        end
        return;
    end

    if isstruct(data)
        for i = 1:numel(aliases)
            alias = aliases{i};
            if isfield(data, alias)
                values = extractStructField(data, alias);
                found = true;
                return;
            end
        end
    end
end

function aliases = getFeatureAliases(fname)
    aliases = {};
    switch lower(fname)
        case 'integrated_intensity'
            aliases = {'integratedIntensity', 'IntegratedIntensity', 'intensity'};
        case 'background'
            aliases = {'background_value', 'backgroundValue', 'bg'};
        case 'r_squared'
            aliases = {'rSquared', 'rsquared', 'R_squared', 'Rsquared'};
        case 'radial_symmetry_score'
            aliases = {'radialSymmetryScore', 'radialSymmetryQuality', 'symmetry_score'};
        case 'amplitude_over_background'
            aliases = {'amplitudeOverBackground', 'amp_bg_ratio'};
        case 'sigma_xy_ratio'
            aliases = {'sigmaXYRatio', 'ellipticity', 'aspect_ratio'};
        case 'sigma_sum'
            aliases = {'sigmaSum', 'total_sigma'};
        case 'sigma_product'
            aliases = {'sigmaProduct', 'sigma_area'};
        case 'sigma_z_ratio'
            aliases = {'sigmaZRatio', 'z_elongation'};
    end

    % Also try simple underscore removal / camel variants.
    if contains(fname, '_')
        compact = strrep(fname, '_', '');
        aliases{end+1} = compact; %#ok<AGROW>
        aliases{end+1} = [upper(compact(1)), compact(2:end)]; %#ok<AGROW>
    end

    aliases = unique(aliases, 'stable');
end

function [values, found] = extractCoordinateFeature(data, axisName, nSamples)
    values = nan(nSamples, 1);
    found = false;

    if istable(data)
        coordNames = {'globalFitCenter', 'fitted_coords', 'maxima_coords'};
        for i = 1:numel(coordNames)
            name = coordNames{i};
            if ~ismember(name, data.Properties.VariableNames)
                continue;
            end
            col = data.(name);
            for r = 1:nSamples
                if isnumeric(col) && ismatrix(col) && size(col, 1) >= r && size(col, 2) >= 1
                    raw = col(r, :);
                elseif iscell(col)
                    raw = col{r};
                else
                    raw = col(r);
                end
                values(r) = extractCoordFromValue(raw, axisName);
            end
            found = true;
            return;
        end
        return;
    end

    if isstruct(data)
        for r = 1:nSamples
            if isfield(data(r), 'globalFitCenter')
                values(r) = extractCoordFromValue(data(r).globalFitCenter, axisName);
                found = true;
            elseif isfield(data(r), 'fitted_coords')
                values(r) = extractCoordFromValue(data(r).fitted_coords, axisName);
                found = true;
            elseif isfield(data(r), 'maxima_coords')
                values(r) = extractCoordFromValue(data(r).maxima_coords, axisName);
                found = true;
            end
        end
    end
end

function value = extractCoordFromValue(raw, axisName)
    idx = 1;
    if strcmp(axisName, 'x')
        idx = 2;
    elseif strcmp(axisName, 'z')
        idx = 3;
    end

    if iscell(raw) && ~isempty(raw)
        raw = raw{1};
    end
    if ~isnumeric(raw) || isempty(raw)
        value = NaN;
        return;
    end

    raw = raw(:);
    if numel(raw) < idx
        value = NaN;
        return;
    end
    value = double(raw(idx));
end

function values = tableColumnToNumericVector(col, nSamples)
    values = nan(nSamples, 1);
    if isnumeric(col) || islogical(col)
        arr = double(col);
        if isvector(arr)
            values = arr(:);
        elseif ismatrix(arr) && size(arr, 1) == nSamples
            for i = 1:nSamples
                row = arr(i, :);
                finiteVals = row(isfinite(row));
                if ~isempty(finiteVals)
                    values(i) = mean(finiteVals);
                end
            end
        else
            flat = arr(:);
            values(1:min(nSamples, numel(flat))) = flat(1:min(nSamples, numel(flat)));
        end
        if numel(values) < nSamples
            values(end+1:nSamples,1) = NaN;
        end
        return;
    end

    if isstring(col) || ischar(col)
        parsed = str2double(string(col));
        values = parsed(:);
        if numel(values) > nSamples
            values = values(1:nSamples);
        elseif numel(values) < nSamples
            values(end+1:nSamples,1) = NaN;
        end
        return;
    end

    if iscell(col)
        for i = 1:min(nSamples, numel(col))
            values(i) = scalarizeNumericValue(col{i});
        end
        return;
    end

    for i = 1:nSamples
        try
            values(i) = scalarizeNumericValue(col(i));
        catch
            values(i) = NaN;
        end
    end
end

function v = scalarizeNumericValue(raw)
    v = NaN;
    if isempty(raw)
        return;
    end
    if iscell(raw) && ~isempty(raw)
        raw = raw{1};
    end
    if isstring(raw) || ischar(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            v = double(parsed);
        end
        return;
    end
    if ~(isnumeric(raw) || islogical(raw))
        return;
    end
    raw = double(raw(:));
    finiteVals = raw(isfinite(raw));
    if isempty(finiteVals)
        return;
    end
    if numel(finiteVals) == 1
        v = finiteVals;
    else
        v = mean(finiteVals);
    end
end
