function spec = processing_builtin()
% Built-in processing plugin: deconvolution + preprocessing + background correction.

    spec = struct();
    spec.id = 'builtin.signal.processing';
    spec.stage = 'signal_processing';
    spec.displayName = 'SNAP Signal Processing';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(~, context)
    tf = isfield(context, 'handles') && ~isempty(context.handles);
end

function state = runPlugin(state, context)
    if ~isfield(context, 'handles') || isempty(context.handles)
        error('processing_builtin requires context.handles.');
    end

    handles = context.handles;
    ch = context.channelIdx;

    if ~isfield(handles, 'rawChannel') || numel(handles.rawChannel) < ch
        rawCell = cell(1, max(ch, 1));
        if isfield(handles, 'rawChannel') && ~isempty(handles.rawChannel)
            rawCell(1:min(numel(handles.rawChannel), numel(rawCell))) = handles.rawChannel(1:min(numel(handles.rawChannel), numel(rawCell)));
        end
        handles.rawChannel = rawCell;
    end
    handles.rawChannel{ch} = state.rawImage;

    t = tic;
    processed = snap_helpers.processImage(handles, ch);
    elapsed = toc(t);

    if isempty(processed)
        state.processedImage = [];
        state.aborted = true;
        snap_modules.signal.emitProgress(context, '[signal_processing] Aborted while processing channel %d.', ch);
        return;
    end

    state.processedImage = processed;
    snap_modules.signal.emitProgress(context, ...
        '[signal_processing] Channel %d processed in %.2f s (size=%s).', ...
        ch, elapsed, mat2str(size(processed)));
end
