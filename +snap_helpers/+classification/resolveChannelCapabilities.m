function [capabilities, meta] = resolveChannelCapabilities(parameterSource, varargin)
% resolveChannelCapabilities - Resolve per-channel fitting/feature capabilities.
%
% Usage:
%   caps = snap_helpers.classification.resolveChannelCapabilities('/path/params.mat')
%   [caps, meta] = snap_helpers.classification.resolveChannelCapabilities(paramsStruct)
%
% Name-value options:
%   'Channels'            - specific channel indices (default: all inferred)
%   'FallbackNumChannels' - fallback channel count if inference is ambiguous (default: 1)
%   'HasPhysicalSpacing'  - expose physical-unit features (default: false)
%   'IncludeFeatureInfo'  - include large featureInfo struct (default: true)

    p = inputParser;
    p.addParameter('Channels', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
    p.addParameter('FallbackNumChannels', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    p.addParameter('HasPhysicalSpacing', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('IncludeFeatureInfo', true, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    [params, paramMeta] = snap_helpers.classification.loadParameterStruct(parameterSource);
    numChannels = localInferNumChannels(params, max(1, round(opts.FallbackNumChannels)));

    if isempty(opts.Channels)
        channels = 1:numChannels;
    else
        channels = unique(round(opts.Channels(:)'));
        channels = channels(isfinite(channels) & channels >= 1);
        if isempty(channels)
            error('resolveChannelCapabilities:InvalidChannels', ...
                'Channels must contain at least one positive integer index.');
        end
    end

    capabilities = struct([]);
    includeFeatureInfo = logical(opts.IncludeFeatureInfo);

    for i = 1:numel(channels)
        ch = channels(i);
        try
            [fittingMethod, has3D] = localInferChannelContext(params, ch);
            cap = snap_helpers.classification.resolveCapabilitiesFromContext( ...
                fittingMethod, has3D, ...
                'ChannelIndex', ch, ...
                'HasPhysicalSpacing', logical(opts.HasPhysicalSpacing));

            if ~includeFeatureInfo && isfield(cap, 'featureInfo')
                cap = rmfield(cap, 'featureInfo');
            end

            if isempty(capabilities)
                capabilities = cap;
            else
                capabilities(end + 1) = cap; %#ok<AGROW>
            end
        catch ME
            error('resolveChannelCapabilities:ChannelResolveFailed', ...
                'Failed to resolve capabilities for channel %d: %s', ch, ME.message);
        end
    end

    meta = struct();
    meta.parameterMeta = paramMeta;
    meta.inferredNumChannels = numChannels;
    meta.resolvedChannels = channels;
    meta.hasPhysicalSpacing = logical(opts.HasPhysicalSpacing);
end

function numChannels = localInferNumChannels(params, fallback)
    numChannels = fallback;
    if ~isstruct(params)
        return;
    end

    explicitN = nan;
    if isfield(params, 'numChannels')
        explicitN = localToScalarNumeric(params.numChannels, nan);
    elseif isfield(params, 'numChan')
        explicitN = localToScalarNumeric(params.numChan, nan);
    elseif isfield(params, 'numChanDrop')
        explicitN = localToScalarNumeric(params.numChanDrop, nan);
    elseif isfield(params, 'workflowConfig') && isstruct(params.workflowConfig) && isfield(params.workflowConfig, 'numChannels')
        explicitN = localToScalarNumeric(params.workflowConfig.numChannels, nan);
    end

    if isfinite(explicitN) && explicitN >= 1
        numChannels = max(1, round(explicitN));
        return;
    end

    hints = {'gaussFitMethod', 'maximaMode', 'preProcMode', ...
        'channelPath', 'channelPaths', 'classifyClassifierPath', 'classifyEnabled'};
    for i = 1:numel(hints)
        fieldName = hints{i};
        if isfield(params, fieldName)
            numChannels = max(numChannels, localCountChannelsFromValue(params.(fieldName)));
        end
    end

    numChannels = max(1, round(numChannels));
end

function n = localCountChannelsFromValue(v)
    n = 1;
    if isempty(v)
        return;
    end
    if iscell(v) || isstring(v)
        n = numel(v);
    elseif isnumeric(v) || islogical(v)
        if isscalar(v)
            n = 1;
        elseif isvector(v)
            n = numel(v);
        end
    end
end

function [fittingMethod, has3D] = localInferChannelContext(params, ch)
    fittingMethod = '3D Gaussian';
    has3D = true;

    if isfield(params, 'gaussFitMethod')
        fm = localGetChannelValue(params.gaussFitMethod, ch, '');
        if ~isempty(fm)
            fittingMethod = char(string(fm));
        end
    end

    modeVal = '';
    if isfield(params, 'maximaMode')
        modeVal = char(string(localGetChannelValue(params.maximaMode, ch, '')));
    elseif isfield(params, 'preProcMode')
        modeVal = char(string(localGetChannelValue(params.preProcMode, ch, '')));
    end

    if ~isempty(modeVal)
        has3D = strcmpi(modeVal, '3D');
    end
end

function value = localGetChannelValue(raw, ch, fallback)
    value = fallback;
    if isempty(raw)
        return;
    end

    if iscell(raw)
        if ch <= numel(raw) && ~isempty(raw{ch})
            value = raw{ch};
        elseif ~isempty(raw{1})
            value = raw{1};
        end
    elseif isstring(raw)
        if ch <= numel(raw)
            value = char(raw(ch));
        else
            value = char(raw(1));
        end
    elseif ischar(raw)
        value = raw;
    elseif isnumeric(raw) || islogical(raw)
        if isscalar(raw)
            value = raw;
        elseif ch <= numel(raw)
            value = raw(ch);
        else
            value = raw(1);
        end
    end
end

function value = localToScalarNumeric(raw, fallback)
    value = fallback;

    if isstruct(raw)
        if isfield(raw, 'Value')
            raw = raw.Value;
        else
            return;
        end
    end

    if iscell(raw)
        if isempty(raw)
            return;
        end
        raw = raw{1};
    elseif isstring(raw)
        if isempty(raw)
            return;
        end
        raw = raw(1);
    end

    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
        value = double(raw);
    elseif ischar(raw) || isstring(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            value = parsed;
        end
    end
end
