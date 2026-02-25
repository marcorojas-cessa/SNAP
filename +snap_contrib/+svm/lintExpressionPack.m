function report = lintExpressionPack(packSource, varargin)
% lintExpressionPack - Validate an SVM expression pack against SNAP capabilities.
%
% Usage:
%   report = snap_contrib.svm.lintExpressionPack('/path/pack.mat', ...
%       'ParameterFile', '/path/params.mat', 'Mode', 'strict');
%
% Name-value options:
%   'ParameterFile'                - parameter MAT file for capability context
%   'ParameterStruct'              - already-loaded parameter struct
%   'Mode'                         - validation mode:
%                                      'strict' (default): fail on incompatibilities
%                                      'permissive': prune incompatible entries with warnings
%   'AutoGuardUnsafeExpressions'   - allow real(...) guard in permissive mode
%   'Verbose'                      - print lint summary to command window

    p = inputParser;
    p.addParameter('ParameterFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ParameterStruct', struct(), @isstruct);
    p.addParameter('Mode', 'strict', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoGuardUnsafeExpressions', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    [rawPack, sourceDescription] = localLoadPackSource(packSource);
    [normalizedPack, normalizeReport] = snap_helpers.classification.normalizeExpressionPack(rawPack);

    caps = [];
    capMeta = struct();

    params = opts.ParameterStruct;
    parameterFile = strtrim(char(string(opts.ParameterFile)));

    if isempty(fieldnames(params)) && ~isempty(parameterFile)
        params = snap_helpers.classification.loadParameterStruct(parameterFile);
    end

    if ~isempty(fieldnames(params))
        channels = [normalizedPack.channelPacks.channelIdx];
        [caps, capMeta] = snap_helpers.classification.resolveChannelCapabilities(params, ...
            'Channels', channels, ...
            'IncludeFeatureInfo', false);
    else
        caps = localFallbackCapabilities(normalizedPack);
        capMeta.warning = ['No parameter context provided. Using fallback capabilities from ', ...
            'pack-required fitting context (best-effort).'];
    end

    [sanitizedPack, validationReport] = ...
        snap_helpers.classification.validateExpressionPackAgainstCapabilities( ...
            normalizedPack, caps, ...
            'Mode', opts.Mode, ...
            'AutoGuardUnsafeExpressions', logical(opts.AutoGuardUnsafeExpressions));

    report = struct();
    report.success = validationReport.success;
    report.source = sourceDescription;
    report.mode = validationReport.mode;
    report.normalizeReport = normalizeReport;
    report.validationReport = validationReport;
    report.capabilityMeta = capMeta;
    report.pack = sanitizedPack;

    if logical(opts.Verbose)
        fprintf('Expression-pack lint: %s\n', sourceDescription);
        fprintf('  Mode: %s\n', validationReport.mode);
        fprintf('  Channels checked: %d\n', numel(sanitizedPack.channelPacks));
        fprintf('  Dropped base features: %d\n', validationReport.nDroppedBase);
        fprintf('  Dropped custom expressions: %d\n', validationReport.nDroppedCustom);
        fprintf('  Auto-guarded expressions: %d\n', validationReport.nAutoGuarded);
        fprintf('  Warnings: %d\n', numel(validationReport.warnings));
        fprintf('  Errors: %d\n', numel(validationReport.errors));
        if validationReport.success
            fprintf('  Result: PASS\n');
        else
            fprintf('  Result: FAIL\n');
            for i = 1:numel(validationReport.errors)
                fprintf('    - %s\n', validationReport.errors{i});
            end
        end
    end
end

function [pack, sourceDescription] = localLoadPackSource(packSource)
    if nargin < 1 || isempty(packSource)
        error('lintExpressionPack:MissingPackSource', ...
            'Provide a pack struct or MAT file path.');
    end

    if isstruct(packSource)
        pack = packSource;
        sourceDescription = 'in-memory pack struct';
        return;
    end

    if ~(ischar(packSource) || isstring(packSource))
        error('lintExpressionPack:InvalidPackSource', ...
            'packSource must be a struct or MAT file path.');
    end

    packPath = char(string(packSource));
    if exist(packPath, 'file') ~= 2
        error('lintExpressionPack:PackNotFound', 'Pack file not found: %s', packPath);
    end

    raw = load(packPath);
    candidates = {};

    if isfield(raw, 'expressionPack')
        candidates{end+1} = raw.expressionPack; %#ok<AGROW>
    end
    if isfield(raw, 'pack')
        candidates{end+1} = raw.pack; %#ok<AGROW>
    end

    fns = fieldnames(raw);
    for i = 1:numel(fns)
        value = raw.(fns{i});
        if isstruct(value)
            if isfield(value, 'channelPacks') || isfield(value, 'selectedFeatures') || isfield(value, 'customExpressions')
                candidates{end+1} = value; %#ok<AGROW>
            end
        end
    end

    if isempty(candidates)
        error('lintExpressionPack:NoPackFound', ...
            'Could not find a pack struct in %s', packPath);
    end

    pack = candidates{1};
    sourceDescription = packPath;
end

function caps = localFallbackCapabilities(pack)
    cps = pack.channelPacks;
    caps = struct([]);

    for i = 1:numel(cps)
        ch = cps(i).channelIdx;
        fm = '';
        h3d = [];

        if isfield(cps(i), 'requiredCapabilities') && isstruct(cps(i).requiredCapabilities)
            if isfield(cps(i).requiredCapabilities, 'fittingMethod')
                fm = cps(i).requiredCapabilities.fittingMethod;
            end
            if isfield(cps(i).requiredCapabilities, 'has3D')
                h3d = cps(i).requiredCapabilities.has3D;
            end
        end

        if isempty(fm) && isfield(cps(i), 'fittingMethod')
            fm = cps(i).fittingMethod;
        end
        if isempty(h3d) && isfield(cps(i), 'has3D')
            h3d = cps(i).has3D;
        end

        cap = snap_helpers.classification.resolveCapabilitiesFromContext(fm, h3d, ...
            'ChannelIndex', ch, ...
            'HasPhysicalSpacing', false);
        if isfield(cap, 'featureInfo')
            cap = rmfield(cap, 'featureInfo');
        end
        if isempty(caps)
            caps = cap;
        else
            caps(end + 1) = cap; %#ok<AGROW>
        end
    end
end
