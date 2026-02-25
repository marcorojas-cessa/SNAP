function [pack, report] = normalizeExpressionPack(rawPack, varargin)
% normalizeExpressionPack - Canonicalize SVM expression-pack schema.
%
% This function normalizes contributed packs into a strict internal schema so
% SNAP_train/SNAP_classify can validate and consume them safely.
%
% Canonical top-level fields:
%   specVersion, packId, minSnapVersion, strictModeDefault, channelPacks
%
% Canonical channel fields:
%   channelIdx, selectedFeatures, customExpressions, requiredFeatures,
%   requiredCapabilities, fittingMethod, has3D, notes

    p = inputParser;
    p.addParameter('DefaultSpecVersion', '2.0.0', @(x) ischar(x) || isstring(x));
    p.addParameter('DefaultMinSnapVersion', '1.1.0', @(x) ischar(x) || isstring(x));
    p.addParameter('DefaultStrictMode', false, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    report = struct();
    report.warnings = {};
    report.wasLegacy = false;

    if nargin < 1 || isempty(rawPack) || ~isstruct(rawPack)
        error('normalizeExpressionPack:InvalidInput', ...
            'Expression pack must be a non-empty struct.');
    end

    packRoot = rawPack;
    if numel(packRoot) > 1
        report.warnings{end+1} = sprintf( ...
            'Pack input is a struct array (%d entries); using first entry.', numel(packRoot)); %#ok<AGROW>
        packRoot = packRoot(1);
        report.wasLegacy = true;
    end

    if ~isfield(packRoot, 'channelPacks')
        [legacyChannelPacks, legacyOk] = localExtractLegacyChannelPacks(rawPack);
        if ~legacyOk
            error('normalizeExpressionPack:MissingChannelPacks', ...
                'Pack must contain channelPacks or legacy selectedFeatures/customExpressions fields.');
        end
        packRoot = struct('channelPacks', legacyChannelPacks);
        report.wasLegacy = true;
        report.warnings{end+1} = 'Legacy pack format detected and normalized.'; %#ok<AGROW>
    end

    channelPacksRaw = packRoot.channelPacks;
    if isempty(channelPacksRaw)
        error('normalizeExpressionPack:EmptyChannelPacks', ...
            'Expression pack contains no channel definitions.');
    end
    if ~isstruct(channelPacksRaw)
        error('normalizeExpressionPack:InvalidChannelPacks', ...
            'channelPacks must be a struct array.');
    end

    pack = struct();
    pack.specVersion = localNormalizeString( ...
        localGetField(packRoot, 'specVersion', localGetField(packRoot, 'version', opts.DefaultSpecVersion)), ...
        char(string(opts.DefaultSpecVersion)));

    defaultPackId = sprintf('snap_expression_pack_%s', datestr(now, 'yyyymmdd_HHMMSS'));
    pack.packId = localNormalizeString( ...
        localGetField(packRoot, 'packId', localGetField(packRoot, 'name', defaultPackId)), ...
        defaultPackId);

    pack.minSnapVersion = localNormalizeString( ...
        localGetField(packRoot, 'minSnapVersion', opts.DefaultMinSnapVersion), ...
        char(string(opts.DefaultMinSnapVersion)));

    strictModeValue = localGetField(packRoot, 'strictModeDefault', []);
    if isempty(strictModeValue)
        validationMode = lower(localNormalizeString(localGetField(packRoot, 'validationMode', ''), ''));
        if strcmp(validationMode, 'strict')
            strictModeValue = true;
        elseif strcmp(validationMode, 'permissive')
            strictModeValue = false;
        else
            strictModeValue = opts.DefaultStrictMode;
        end
    end
    pack.strictModeDefault = logical(strictModeValue);

    pack.description = localNormalizeString(localGetField(packRoot, 'description', ''), '');
    pack.generatedOn = localNormalizeString(localGetField(packRoot, 'generatedOn', ''), '');

    normalizedChannels = repmat(localEmptyChannelPack(), 1, numel(channelPacksRaw));
    for i = 1:numel(channelPacksRaw)
        chRaw = channelPacksRaw(i);
        chNorm = localEmptyChannelPack();

        chNorm.channelIdx = localNormalizeChannelIndex(chRaw, i);
        chNorm.selectedFeatures = localNormalizeCellstr(localGetField(chRaw, 'selectedFeatures', {}));

        [customExpr, exprWarnings] = localNormalizeCustomExpressions( ...
            localGetField(chRaw, 'customExpressions', struct('name', {}, 'expression', {})), chNorm.channelIdx);
        for w = 1:numel(exprWarnings)
            report.warnings{end+1} = exprWarnings{w}; %#ok<AGROW>
        end
        chNorm.customExpressions = customExpr;

        requiredFeatures = localNormalizeCellstr(localGetField(chRaw, 'requiredFeatures', {}));
        exprReq = localCollectExpressionRequiredFeatures(customExpr);
        chNorm.requiredFeatures = unique([requiredFeatures, chNorm.selectedFeatures, exprReq], 'stable');

        reqCapsRaw = localGetField(chRaw, 'requiredCapabilities', struct());
        chNorm.requiredCapabilities = localNormalizeRequiredCapabilities(reqCapsRaw, chRaw);

        chNorm.fittingMethod = localNormalizeString( ...
            localGetField(chRaw, 'fittingMethod', localGetField(chNorm.requiredCapabilities, 'fittingMethod', '')), ...
            '');

        has3DVal = localGetField(chRaw, 'has3D', localGetField(chNorm.requiredCapabilities, 'has3D', []));
        if isempty(has3DVal)
            chNorm.has3D = [];
        else
            chNorm.has3D = logical(has3DVal);
        end

        chNorm.notes = localNormalizeString(localGetField(chRaw, 'notes', ''), '');
        normalizedChannels(i) = chNorm;
    end

    pack.channelPacks = normalizedChannels;
end

function base = localEmptyChannelPack()
    base = struct( ...
        'channelIdx', 1, ...
        'selectedFeatures', {{}}, ...
        'customExpressions', struct('name', {}, 'expression', {}, 'requiredFeatures', {}, 'requiredCapabilities', {}), ...
        'requiredFeatures', {{}}, ...
        'requiredCapabilities', struct('fittingMethod', '', 'has3D', [], 'hasPhysicalSpacing', []), ...
        'fittingMethod', '', ...
        'has3D', [], ...
        'notes', '');
end

function idx = localNormalizeChannelIndex(chRaw, fallbackIdx)
    idx = fallbackIdx;
    if isfield(chRaw, 'channelIdx')
        idx = localToPositiveInteger(chRaw.channelIdx, fallbackIdx);
        return;
    end
    if isfield(chRaw, 'channel')
        idx = localToPositiveInteger(chRaw.channel, fallbackIdx);
    end
end

function [out, warnings] = localNormalizeCustomExpressions(customRaw, channelIdx)
    out = struct('name', {}, 'expression', {}, 'requiredFeatures', {}, 'requiredCapabilities', {});
    warnings = {};

    if isempty(customRaw)
        out = out;
        return;
    end

    if istable(customRaw)
        customRaw = table2struct(customRaw);
    end

    if ~isstruct(customRaw)
        warnings{end+1} = sprintf('[Channel %d] customExpressions is not a struct/table and was ignored.', channelIdx); %#ok<AGROW>
        out = out;
        return;
    end

    nextIdx = 0;
    for i = 1:numel(customRaw)
        entry = customRaw(i);
        name = localNormalizeString(localGetField(entry, 'name', ''), '');
        expr = localNormalizeString(localGetField(entry, 'expression', ''), '');

        if isempty(expr)
            warnings{end+1} = sprintf('[Channel %d] Dropped custom expression entry %d: missing expression text.', ...
                channelIdx, i); %#ok<AGROW>
            continue;
        end

        if isempty(name)
            name = sprintf('expr_%d', i);
            warnings{end+1} = sprintf('[Channel %d] Custom expression %d had no name; assigned "%s".', ...
                channelIdx, i, name); %#ok<AGROW>
        end

        requiredFeatures = localNormalizeCellstr(localGetField(entry, 'requiredFeatures', {}));
        if isempty(requiredFeatures)
            requiredFeatures = localInferExpressionRequiredFeatures(expr);
        end

        reqCaps = localNormalizeRequiredCapabilities(localGetField(entry, 'requiredCapabilities', struct()), struct());

        nextIdx = nextIdx + 1;
        out(nextIdx) = struct( ... %#ok<AGROW>
            'name', name, ...
            'expression', expr, ...
            'requiredFeatures', {requiredFeatures}, ...
            'requiredCapabilities', reqCaps);
    end

    if nextIdx == 0
        out = struct('name', {}, 'expression', {}, 'requiredFeatures', {}, 'requiredCapabilities', {});
    end
end

function reqCaps = localNormalizeRequiredCapabilities(rawReq, fallbackRaw)
    reqCaps = struct('fittingMethod', '', 'has3D', [], 'hasPhysicalSpacing', []);

    if isstruct(rawReq) && ~isempty(fieldnames(rawReq))
        reqCaps.fittingMethod = localNormalizeString(localGetField(rawReq, 'fittingMethod', ''), '');
        if isfield(rawReq, 'has3D') && ~isempty(rawReq.has3D)
            reqCaps.has3D = logical(rawReq.has3D);
        end
        if isfield(rawReq, 'hasPhysicalSpacing') && ~isempty(rawReq.hasPhysicalSpacing)
            reqCaps.hasPhysicalSpacing = logical(rawReq.hasPhysicalSpacing);
        end
    end

    if isempty(reqCaps.fittingMethod) && isstruct(fallbackRaw) && isfield(fallbackRaw, 'fittingMethod')
        reqCaps.fittingMethod = localNormalizeString(fallbackRaw.fittingMethod, '');
    end

    if isempty(reqCaps.has3D) && isstruct(fallbackRaw) && isfield(fallbackRaw, 'has3D') && ~isempty(fallbackRaw.has3D)
        reqCaps.has3D = logical(fallbackRaw.has3D);
    end
end

function req = localCollectExpressionRequiredFeatures(customExpr)
    req = {};
    for i = 1:numel(customExpr)
        if isfield(customExpr(i), 'requiredFeatures')
            req = [req, localNormalizeCellstr(customExpr(i).requiredFeatures)]; %#ok<AGROW>
        end
    end
    req = unique(req, 'stable');
end

function req = localInferExpressionRequiredFeatures(expr)
    req = {};
    if isempty(expr)
        return;
    end

    exprText = char(expr);
    % Strip scientific-notation exponent fragments (e.g., 1e-6) so "e"
    % is not misclassified as a required feature token.
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

function [channelPacks, ok] = localExtractLegacyChannelPacks(raw)
    channelPacks = struct('channelIdx', {}, 'selectedFeatures', {}, 'customExpressions', {});
    ok = false;

    if ~isstruct(raw)
        return;
    end

    if isfield(raw, 'selectedFeatures') || isfield(raw, 'customExpressions')
        channelPacks = raw;
        ok = true;
        return;
    end

    if numel(raw) == 1
        fns = fieldnames(raw);
        for i = 1:numel(fns)
            candidate = raw.(fns{i});
            if isstruct(candidate) && (isfield(candidate, 'selectedFeatures') || isfield(candidate, 'customExpressions'))
                channelPacks = candidate;
                ok = true;
                return;
            end
        end
    end
end

function out = localNormalizeCellstr(v)
    out = {};
    if isempty(v)
        return;
    end

    if ischar(v)
        out = {strtrim(v)};
        return;
    end

    if isstring(v)
        out = cellstr(v(:));
    elseif iscell(v)
        out = cell(size(v));
        for i = 1:numel(v)
            out{i} = char(string(v{i}));
        end
    else
        out = {char(string(v))};
    end

    out = out(:)';
    out = cellfun(@(s) strtrim(char(string(s))), out, 'UniformOutput', false);
    out = out(~cellfun(@isempty, out));
    out = unique(out, 'stable');
end

function out = localNormalizeString(v, defaultValue)
    if nargin < 2
        defaultValue = '';
    end
    if isempty(v)
        out = defaultValue;
        return;
    end
    out = strtrim(char(string(v)));
    if isempty(out)
        out = defaultValue;
    end
end

function value = localGetField(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function n = localToPositiveInteger(v, fallback)
    n = fallback;
    if iscell(v)
        if isempty(v)
            return;
        end
        v = v{1};
    end
    if isstring(v) || ischar(v)
        v = str2double(string(v));
    end
    if isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1
        n = max(1, round(double(v)));
    end
end
