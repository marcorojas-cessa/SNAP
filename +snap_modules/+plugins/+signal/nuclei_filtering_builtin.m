function spec = nuclei_filtering_builtin()
% Built-in nuclei inclusion/exclusion filtering plugin.

    spec = struct();
    spec.id = 'builtin.signal.nuclei_filtering';
    spec.stage = 'nuclei_filtering';
    spec.displayName = 'SNAP Nuclei Filtering';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(state, context)
    tf = false;
    if isempty(state.maximaCoords)
        return;
    end
    if ~isfield(context, 'enableNucleiFiltering') || ~logical(context.enableNucleiFiltering)
        return;
    end
    if ~isfield(context, 'nucleiMask') || isempty(context.nucleiMask)
        return;
    end

    if isfield(context, 'params') && isstruct(context.params)
        tf = getScalarLogical(context.params, 'nucInclusionExclusionEnabled', false);
    else
        tf = true;
    end
end

function state = runPlugin(state, context)
    if isempty(state.maximaCoords)
        state.nucleiFilterMask = [];
        return;
    end

    ch = context.channelIdx;

    mode = readMode(context);
    applyTo = readApplyTo(context);
    if ~shouldApplyToChannel(applyTo, ch)
        return;
    end

    nBefore = size(state.maximaCoords, 1);
    [filteredCoords, keepMask] = snap_helpers.filterMaximaByNuclei(state.maximaCoords, context.nucleiMask, mode);

    state.maximaCoords = filteredCoords;
    state.nucleiFilterMask = keepMask;

    if ~isempty(state.fitResults) && numel(state.fitResults) == numel(keepMask)
        state.fitResults = state.fitResults(keepMask);
    end

    snap_modules.signal.emitProgress(context, ...
        '[nuclei_filtering] Channel %d kept %d/%d maxima (%s mode).', ...
        ch, size(filteredCoords, 1), nBefore, mode);
end

function tf = shouldApplyToChannel(applyTo, ch)
    tf = strcmpi(applyTo, 'All Channels') || strcmpi(applyTo, sprintf('Channel %d', ch));
end

function mode = readMode(context)
    mode = 'inside';
    if isfield(context, 'nucleiFilterMode') && ~isempty(context.nucleiFilterMode)
        mode = char(string(context.nucleiFilterMode));
        return;
    end
    if isfield(context, 'params') && isstruct(context.params) && isfield(context.params, 'nucInclusionExclusionMode')
        mode = char(string(context.params.nucInclusionExclusionMode));
    end
end

function applyTo = readApplyTo(context)
    applyTo = 'All Channels';
    if isfield(context, 'nucleiFilterApplyTo') && ~isempty(context.nucleiFilterApplyTo)
        applyTo = char(string(context.nucleiFilterApplyTo));
        return;
    end
    if isfield(context, 'params') && isstruct(context.params) && isfield(context.params, 'nucInclusionExclusionApplyTo')
        applyTo = char(string(context.params.nucInclusionExclusionApplyTo));
    end
end

function value = getScalarLogical(params, fieldName, defaultValue)
    value = defaultValue;
    if ~isstruct(params) || ~isfield(params, fieldName)
        return;
    end
    raw = params.(fieldName);
    try
        if iscell(raw)
            if isempty(raw)
                return;
            end
            value = logical(raw{1});
        else
            value = logical(raw);
        end
    catch
        value = defaultValue;
    end
end
