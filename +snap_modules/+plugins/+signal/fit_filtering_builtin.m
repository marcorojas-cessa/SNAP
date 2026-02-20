function spec = fit_filtering_builtin()
% Built-in fit filtering plugin.

    spec = struct();
    spec.id = 'builtin.signal.fit_filtering';
    spec.stage = 'fit_filtering';
    spec.displayName = 'SNAP Fit Filtering';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(~, context)
    tf = readChannelLogical(context, 'fitFilterEnabledChecks', false);
end

function state = runPlugin(state, context)
    if isempty(state.fitResults)
        state.fitFilterMask = [];
        return;
    end

    ch = context.channelIdx;
    [filteredFits, mask] = snap_helpers.applyFitFiltering(state.fitResults, ch, context.handles);

    if isempty(mask)
        mask = true(1, numel(state.fitResults));
    end

    state.fitFilterMask = logical(mask(:))';
    state.fitResults = filteredFits;

    if ~isempty(state.maximaCoords) && numel(state.fitFilterMask) == size(state.maximaCoords, 1)
        state.maximaCoords = state.maximaCoords(state.fitFilterMask, :);
    end

    snap_modules.signal.emitProgress(context, ...
        '[fit_filtering] Channel %d kept %d/%d fit(s).', ...
        ch, numel(state.fitResults), numel(mask));
end

function tf = readChannelLogical(context, fieldName, defaultValue)
    tf = defaultValue;
    if ~isfield(context, 'handles') || ~isstruct(context.handles)
        return;
    end
    handles = context.handles;
    ch = context.channelIdx;
    if ~isfield(handles, fieldName)
        return;
    end
    try
        raw = handles.(fieldName);
        if numel(raw) < ch
            return;
        end
        entry = raw(ch);
        if isstruct(entry) && isfield(entry, 'Value')
            tf = logical(entry.Value);
        elseif isobject(entry) && isprop(entry, 'Value')
            tf = logical(entry.Value);
        end
    catch
        tf = defaultValue;
    end
end
