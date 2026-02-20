function spec = gaussian_fitting_builtin()
% Built-in Gaussian fitting plugin.

    spec = struct();
    spec.id = 'builtin.signal.gaussian_fitting';
    spec.stage = 'gaussian_fitting';
    spec.displayName = 'SNAP Gaussian Fitting';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(~, context)
    tf = readChannelLogical(context, 'gaussFitEnabledChecks', true);
end

function state = runPlugin(state, context)
    if isempty(state.maximaCoords)
        state.fitResults = struct([]);
        snap_modules.signal.emitProgress(context, '[gaussian_fitting] No maxima candidates; skipping fitting.');
        return;
    end

    ch = context.channelIdx;

    if isfield(context, 'fitParams') && isstruct(context.fitParams) && ~isempty(context.fitParams)
        fitParams = context.fitParams;
    else
        fitParams = snap_modules.signal.getFitParamsFromHandles(context.handles, ch);
    end

    is3D = ndims(state.rawImage) == 3 && size(state.rawImage, 3) > 1;

    abortCheckFcn = [];
    if isfield(context, 'abortCheckFcn') && isa(context.abortCheckFcn, 'function_handle')
        abortCheckFcn = context.abortCheckFcn;
    end

    t = tic;
    [fitResults, wasAborted] = snap_helpers.fitGaussians( ...
        state.rawImage, state.maximaCoords, fitParams, is3D, ...
        abortCheckFcn, context.fitAbortPollInterval, context.progressCallback, context.fitProgressInterval);
    elapsed = toc(t);

    if isempty(fitResults)
        fitResults = struct([]);
    end

    state.fitResults = fitResults;

    if size(state.maximaCoords, 1) ~= numel(fitResults)
        nKeep = min(size(state.maximaCoords, 1), numel(fitResults));
        if nKeep <= 0
            state.maximaCoords = zeros(0, 3);
            state.fitResults = struct([]);
        else
            state.maximaCoords = state.maximaCoords(1:nKeep, :);
            state.fitResults = state.fitResults(1:nKeep);
        end
    end

    state.aborted = logical(wasAborted);

    snap_modules.signal.emitProgress(context, ...
        '[gaussian_fitting] Channel %d produced %d fit(s) in %.2f s.', ...
        ch, numel(state.fitResults), elapsed);
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
