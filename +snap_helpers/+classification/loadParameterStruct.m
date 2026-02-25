function [params, meta] = loadParameterStruct(source)
% loadParameterStruct - Load a SNAP parameter struct from file or struct input.
%
% Usage:
%   [params, meta] = snap_helpers.classification.loadParameterStruct('/path/to/params.mat')
%   [params, meta] = snap_helpers.classification.loadParameterStruct(paramStruct)
%
% Supported parameter containers:
%   - batchConfig.parameters
%   - paramData.parameters
%   - paramData.workflowConfig
%   - workflowConfig
%   - parameters
%   - lastUsed

    meta = struct();
    meta.source = '';
    meta.container = '';
    meta.warnings = {};

    if nargin < 1 || isempty(source)
        error('loadParameterStruct:MissingInput', 'Provide a parameter file path or parameter struct.');
    end

    if ischar(source) || isstring(source)
        filePath = char(string(source));
        if exist(filePath, 'file') ~= 2
            error('loadParameterStruct:FileNotFound', 'Parameter file not found: %s', filePath);
        end
        raw = load(filePath);
        [params, metaInner] = localExtractParameterStruct(raw);
        meta = localMergeMeta(meta, metaInner);
        meta.source = filePath;
        return;
    end

    if ~isstruct(source)
        error('loadParameterStruct:InvalidInputType', ...
            'Input must be a parameter file path or struct.');
    end

    [params, metaInner] = localExtractParameterStruct(source);
    meta = localMergeMeta(meta, metaInner);
    if isempty(meta.source)
        meta.source = 'in-memory struct';
    end
end

function [params, meta] = localExtractParameterStruct(raw)
    meta = struct();
    meta.source = '';
    meta.container = '';
    meta.warnings = {};

    params = struct();

    if ~isstruct(raw)
        error('loadParameterStruct:InvalidStruct', 'Loaded parameter container is not a struct.');
    end

    if isfield(raw, 'batchConfig') && isstruct(raw.batchConfig) && isfield(raw.batchConfig, 'parameters')
        params = raw.batchConfig.parameters;
        meta.container = 'batchConfig.parameters';
        if isfield(raw.batchConfig, 'workflowConfig') && isstruct(raw.batchConfig.workflowConfig)
            params = localMergeWorkflowHints(params, raw.batchConfig.workflowConfig);
        end
        localAssertStruct(params);
        return;
    end

    if isfield(raw, 'paramData') && isstruct(raw.paramData)
        if isfield(raw.paramData, 'parameters')
            params = raw.paramData.parameters;
            meta.container = 'paramData.parameters';
            if isfield(raw.paramData, 'workflowConfig') && isstruct(raw.paramData.workflowConfig)
                params = localMergeWorkflowHints(params, raw.paramData.workflowConfig);
            end
            localAssertStruct(params);
            return;
        end
        if isfield(raw.paramData, 'workflowConfig')
            params = raw.paramData.workflowConfig;
            meta.container = 'paramData.workflowConfig';
            localAssertStruct(params);
            return;
        end
    end

    if isfield(raw, 'workflowConfig')
        params = raw.workflowConfig;
        meta.container = 'workflowConfig';
        localAssertStruct(params);
        return;
    end

    if isfield(raw, 'parameters')
        params = raw.parameters;
        meta.container = 'parameters';
        localAssertStruct(params);
        return;
    end

    if isfield(raw, 'lastUsed')
        params = raw.lastUsed;
        meta.container = 'lastUsed';
        localAssertStruct(params);
        return;
    end

    % If the struct already appears to be a parameter struct, accept it.
    if localLooksLikeParameterStruct(raw)
        params = raw;
        meta.container = 'direct';
        localAssertStruct(params);
        return;
    end

    % Fallback: if there is a single struct field, use it with warning.
    fns = fieldnames(raw);
    if numel(fns) == 1 && isstruct(raw.(fns{1}))
        params = raw.(fns{1});
        meta.container = sprintf('fallback:%s', fns{1});
        meta.warnings{end+1} = sprintf([ ...
            'Used fallback parameter container "%s". ', ...
            'Prefer standard containers (parameters/workflowConfig/lastUsed).'], fns{1}); %#ok<AGROW>
        localAssertStruct(params);
        return;
    end

    error('loadParameterStruct:UnsupportedContainer', [ ...
        'Could not find a SNAP parameter struct. Expected one of: ', ...
        'batchConfig.parameters, paramData.parameters, paramData.workflowConfig, ', ...
        'workflowConfig, parameters, or lastUsed.']);
end

function out = localMergeWorkflowHints(params, workflowConfig)
    out = params;
    if ~isstruct(out)
        return;
    end
    if ~isfield(out, 'numChannels') && isfield(workflowConfig, 'numChannels')
        out.numChannels = workflowConfig.numChannels;
    end
end

function tf = localLooksLikeParameterStruct(s)
    tf = false;
    if ~isstruct(s)
        return;
    end

    hints = {'gaussFitMethod', 'maximaMode', 'preProcMode', 'numChannels', 'numChanDrop', 'classifyEnabled'};
    for i = 1:numel(hints)
        if isfield(s, hints{i})
            tf = true;
            return;
        end
    end
end

function localAssertStruct(v)
    if ~isstruct(v)
        error('loadParameterStruct:InvalidParameterStruct', 'Resolved parameter container is not a struct.');
    end
end

function out = localMergeMeta(base, extra)
    out = base;
    if isfield(extra, 'source') && ~isempty(extra.source)
        out.source = extra.source;
    end
    if isfield(extra, 'container') && ~isempty(extra.container)
        out.container = extra.container;
    end
    if isfield(extra, 'warnings') && ~isempty(extra.warnings)
        out.warnings = [out.warnings, extra.warnings];
    end
end
