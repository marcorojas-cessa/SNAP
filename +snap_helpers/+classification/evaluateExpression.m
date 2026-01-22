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
    
    % Add all available base features to workspace
    for i = 1:numel(availableFeatures)
        fname = availableFeatures{i};
        safeName = makeValidVarName(fname);
        
        try
            if istable(data)
                if ismember(fname, data.Properties.VariableNames)
                    workspace.(safeName) = data.(fname);
                end
            elseif isstruct(data)
                if isfield(data, fname)
                    workspace.(safeName) = extractStructField(data, fname);
                end
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
    try
        values = [data.(fname)]';
        values = values(:);
    catch
        values = nan(numel(data), 1);
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

