function info = resolveModelStatsRequirements(selectedFeatures, customExpressions)
% resolveModelStatsRequirements - Determine whether model-stat fields are needed.
%
% Usage:
%   info = snap_helpers.classification.resolveModelStatsRequirements(selected, customExpr)
%
% Output fields:
%   required                 - true when any model-stat feature is referenced
%   requiredFeatures         - unique referenced model-stat fields
%   requiredFromSelected     - references coming from selected base features
%   requiredFromCustom       - references coming from custom expressions
%   expressionsNeedingModelStats - custom expression names requiring model stats

    if nargin < 1 || isempty(selectedFeatures)
        selectedFeatures = {};
    end
    if nargin < 2
        customExpressions = struct('name', {}, 'expression', {});
    end

    modelStatsFeatures = string(snap_helpers.classification.getModelStatsFeatureNames());
    selected = string(localNormalizeCellstr(selectedFeatures));

    requiredFromSelected = selected(ismember(selected, modelStatsFeatures));

    requiredFromCustom = strings(0, 1);
    expressionNames = strings(0, 1);
    customExpressions = localNormalizeCustomExpressionList(customExpressions);

    for i = 1:numel(customExpressions)
        exprName = string(customExpressions(i).name);

        reqFeatures = localNormalizeCellstr(localGetField(customExpressions(i), 'requiredFeatures', {}));
        if isempty(reqFeatures)
            exprText = char(string(localGetField(customExpressions(i), 'expression', '')));
            reqFeatures = localInferExpressionRequiredFeatures(exprText);
        end
        reqFeatures = string(reqFeatures(:)');
        hits = reqFeatures(ismember(reqFeatures, modelStatsFeatures));
        if isempty(hits)
            continue;
        end

        requiredFromCustom = [requiredFromCustom; hits(:)]; %#ok<AGROW>
        expressionNames = [expressionNames; repmat(exprName, numel(hits), 1)]; %#ok<AGROW>
    end

    requiredFeatures = unique([requiredFromSelected(:); requiredFromCustom(:)], 'stable');

    info = struct();
    info.required = ~isempty(requiredFeatures);
    info.requiredFeatures = cellstr(requiredFeatures');
    info.requiredFromSelected = cellstr(unique(requiredFromSelected, 'stable')');
    info.requiredFromCustom = cellstr(unique(requiredFromCustom, 'stable')');
    info.expressionsNeedingModelStats = cellstr(unique(expressionNames(expressionNames ~= ""), 'stable')');
end

function out = localNormalizeCellstr(value)
    if isempty(value)
        out = {};
        return;
    end
    if ischar(value)
        out = {char(value)};
    elseif isstring(value)
        out = cellstr(value(:));
    elseif iscell(value)
        out = cell(size(value));
        keep = false(size(value));
        for i = 1:numel(value)
            if ischar(value{i}) || isstring(value{i})
                out{i} = char(string(value{i}));
                keep(i) = ~isempty(strtrim(out{i}));
            end
        end
        out = out(keep);
    else
        out = {};
    end
    out = unique(out, 'stable');
end

function req = localInferExpressionRequiredFeatures(exprText)
    req = {};
    if isempty(exprText)
        return;
    end

    exprText = char(exprText);
    % Strip scientific-notation exponent fragments (e.g., 1e-6) so "e"
    % is not interpreted as a required feature token.
    exprText = regexprep(exprText, '(?<=[0-9])[eEdD][+\-]?[0-9]+', ' ');

    tokens = regexp(exprText, '[A-Za-z_][A-Za-z0-9_]*', 'match');
    if isempty(tokens)
        return;
    end

    reserved = lower({ ...
        'log','log10','log2','sqrt','abs','exp','sin','cos','tan', ...
        'min','max','mean','std','median','sum','prod','sign','floor','ceil','round', ...
        'real','imag','conj','pi','inf','nan','true','false','workspace'});

    keep = true(size(tokens));
    for i = 1:numel(tokens)
        tok = lower(tokens{i});
        if strcmp(tok, 'e') || strcmp(tok, 'd')
            keep(i) = false;
        elseif ismember(tok, reserved)
            keep(i) = false;
        end
    end

    req = unique(tokens(keep), 'stable');
end

function customOut = localNormalizeCustomExpressionList(customIn)
    customOut = struct('name', {}, 'expression', {}, 'requiredFeatures', {});
    if isempty(customIn)
        return;
    end

    if istable(customIn)
        customIn = table2struct(customIn);
    end
    if ~isstruct(customIn)
        return;
    end

    for i = 1:numel(customIn)
        name = strtrim(char(string(localGetField(customIn(i), 'name', ''))));
        expr = strtrim(char(string(localGetField(customIn(i), 'expression', ''))));
        if isempty(name) || isempty(expr)
            continue;
        end
        req = localNormalizeCellstr(localGetField(customIn(i), 'requiredFeatures', {}));
        customOut(end+1) = struct( ... %#ok<AGROW>
            'name', name, ...
            'expression', expr, ...
            'requiredFeatures', {req});
    end
end

function value = localGetField(s, fieldName, defaultValue)
    value = defaultValue;
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    end
end
