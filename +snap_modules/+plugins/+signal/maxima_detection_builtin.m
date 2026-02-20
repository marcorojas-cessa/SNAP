function spec = maxima_detection_builtin()
% Built-in maxima detection plugin.

    spec = struct();
    spec.id = 'builtin.signal.maxima_detection';
    spec.stage = 'maxima_detection';
    spec.displayName = 'SNAP Maxima Detection';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(~, context)
    tf = readChannelLogical(context, 'maximaEnabledChecks', true);
end

function state = runPlugin(state, context)
    if isempty(state.processedImage)
        state.maximaCoords = zeros(0, 3);
        snap_modules.signal.emitProgress(context, '[maxima_detection] No processed image available; skipping maxima detection.');
        return;
    end

    handles = context.handles;
    ch = context.channelIdx;

    t = tic;
    maximaCoords = snap_helpers.detectMaxima(state.processedImage, handles, ch);
    elapsed = toc(t);

    if isempty(maximaCoords)
        maximaCoords = zeros(0, 3);
    elseif size(maximaCoords, 2) == 2
        maximaCoords = [maximaCoords, ones(size(maximaCoords, 1), 1)];
    end

    state.maximaCoords = maximaCoords;
    snap_modules.signal.emitProgress(context, ...
        '[maxima_detection] Channel %d found %d candidate(s) in %.2f s.', ...
        ch, size(maximaCoords, 1), elapsed);
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
