function result = runPipeline(rawImage, context, varargin)
% runPipeline - Shared stage-driven signal pipeline for SNAP applications.
%
% This runner executes one selected plugin per stage. Built-in plugins wrap
% existing SNAP algorithms so behavior remains backward-compatible while
% enabling extension via external plugins.
%
% Inputs:
%   rawImage  - 2D/3D numeric channel image
%   context   - struct with optional fields:
%       .channelIdx          (default: 1)
%       .handles             handle-like parameter struct (required by built-ins)
%       .params              exported SNAP parameter struct
%       .mode                label string (e.g., 'batch', 'train', 'snap')
%       .progressCallback    function handle accepting one string message
%       .abortCheckFcn       function handle returning true when abort requested
%       .enableClassification logical flag (default: true)
%       .enableNucleiFiltering logical flag (default: true)
%       .nucleiMask          nuclei mask for nuclei filtering stage
%       .fitProgressInterval maxima interval for fit progress updates (default: 250)
%       .fitAbortPollInterval maxima interval for fit abort checks (default: 5)
%
% Name-value options:
%   'Registry'        - pre-built plugin registry struct array
%   'StageOrder'      - cellstr stage order override
%   'IncludeExternal' - include externally discovered plugins (default: true)
%
% Output:
%   result - state struct with fields:
%       rawImage, processedImage, maximaCoords, fitResults,
%       fitFilterMask, classification, nucleiFilterMask, aborted,
%       stageMetrics, stageOrder, registry

    if nargin < 2 || isempty(context)
        context = struct();
    end

    p = inputParser;
    p.addParameter('Registry', struct([]), @isstruct);
    p.addParameter('StageOrder', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.addParameter('IncludeExternal', true, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    if ~isnumeric(rawImage) || isempty(rawImage)
        error('runPipeline:InvalidInput', 'rawImage must be a non-empty numeric array.');
    end

    context = normalizeContext(context);

    if isempty(opts.Registry)
        registry = getRegistry(logical(opts.IncludeExternal), context);
    else
        registry = applyRegistryDefaults(opts.Registry);
    end

    stageOrder = normalizeStageOrder(opts.StageOrder);

    state = struct();
    state.rawImage = double(rawImage);
    state.processedImage = [];
    state.maximaCoords = zeros(0, 3);
    state.fitResults = struct([]);
    state.fitFilterMask = [];
    state.classification = struct('predictions', [], 'scores', [], 'confidence', [], ...
                                  'classLabels', [], 'validMask', [], 'keepMask', []);
    state.nucleiFilterMask = [];
    state.aborted = false;

    if isfield(context, 'initialState') && isstruct(context.initialState)
        initFields = fieldnames(context.initialState);
        for i = 1:numel(initFields)
            state.(initFields{i}) = context.initialState.(initFields{i});
        end
    end

    stageMetrics = struct('stage', {}, 'pluginId', {}, 'elapsedSec', {});

    for s = 1:numel(stageOrder)
        stageName = stageOrder{s};
        plugin = selectPlugin(registry, stageName, state, context);
        if isempty(plugin)
            continue;
        end

        snap_modules.signal.emitProgress(context, '[%s] Running plugin: %s', stageName, plugin.id);
        stageTic = tic;

        try
            state = plugin.run(state, context);
        catch ME
            error('runPipeline:StageFailed', ...
                'Pipeline stage "%s" failed in plugin "%s": %s', ...
                stageName, plugin.id, ME.message);
        end

        if ~isstruct(state)
            error('runPipeline:InvalidState', ...
                'Plugin "%s" returned a non-struct pipeline state.', plugin.id);
        end
        if ~isfield(state, 'aborted')
            state.aborted = false;
        end

        elapsedSec = toc(stageTic);
        stageMetrics(end+1) = struct('stage', stageName, 'pluginId', plugin.id, 'elapsedSec', elapsedSec); %#ok<AGROW>
        snap_modules.signal.emitProgress(context, '[%s] Completed plugin: %s (%.2f s)', stageName, plugin.id, elapsedSec);

        if state.aborted
            snap_modules.signal.emitProgress(context, '[%s] Pipeline aborted by plugin "%s".', stageName, plugin.id);
            break;
        end

        if localAbortRequested(context)
            state.aborted = true;
            snap_modules.signal.emitProgress(context, '[%s] Pipeline aborted by external abort callback.', stageName);
            break;
        end
    end

    result = state;
    result.stageMetrics = stageMetrics;
    result.stageOrder = stageOrder;
    result.registry = registry;
end

function context = normalizeContext(context)
    if nargin < 1 || isempty(context)
        context = struct();
    end
    if ~isstruct(context)
        error('runPipeline:InvalidContext', 'context must be a struct.');
    end

    if ~isfield(context, 'channelIdx') || isempty(context.channelIdx)
        context.channelIdx = 1;
    end
    context.channelIdx = max(1, round(double(context.channelIdx)));

    if ~isfield(context, 'mode') || isempty(context.mode)
        context.mode = 'generic';
    end

    if ~isfield(context, 'enableClassification') || isempty(context.enableClassification)
        context.enableClassification = true;
    end
    context.enableClassification = logical(context.enableClassification);

    if ~isfield(context, 'enableNucleiFiltering') || isempty(context.enableNucleiFiltering)
        context.enableNucleiFiltering = true;
    end
    context.enableNucleiFiltering = logical(context.enableNucleiFiltering);

    if ~isfield(context, 'fitProgressInterval') || isempty(context.fitProgressInterval)
        context.fitProgressInterval = 250;
    end
    context.fitProgressInterval = max(1, round(double(context.fitProgressInterval)));

    if ~isfield(context, 'fitAbortPollInterval') || isempty(context.fitAbortPollInterval)
        context.fitAbortPollInterval = 5;
    end
    context.fitAbortPollInterval = max(1, round(double(context.fitAbortPollInterval)));

    if ~isfield(context, 'progressCallback') || isempty(context.progressCallback)
        context.progressCallback = [];
    end

    if ~isfield(context, 'abortCheckFcn') || isempty(context.abortCheckFcn)
        context.abortCheckFcn = [];
    end
end

function registry = getRegistry(includeExternal, context)
    registry = snap_modules.signal.registerBuiltinPlugins();

    if includeExternal
        external = discoverExternalPlugins(context);
        if ~isempty(external)
            registry = [registry, external]; %#ok<AGROW>
        end
    end

    registry = applyRegistryDefaults(registry);
end

function registry = applyRegistryDefaults(registry)
    if isempty(registry)
        registry = struct([]);
        return;
    end

    normalized = struct('id', {}, 'stage', {}, 'displayName', {}, 'priority', {}, ...
                        'version', {}, 'source', {}, 'supportsFcn', {}, 'isEnabledFcn', {}, ...
                        'run', {});

    for i = 1:numel(registry)
        spec = registry(i);
        if ~isfield(spec, 'id') || isempty(spec.id)
            spec.id = sprintf('plugin_%d', i);
        end
        if ~isfield(spec, 'displayName') || isempty(spec.displayName)
            spec.displayName = spec.id;
        end
        if ~isfield(spec, 'priority') || isempty(spec.priority)
            spec.priority = 100;
        end
        if ~isfield(spec, 'version') || isempty(spec.version)
            spec.version = '1.0.0';
        end
        if ~isfield(spec, 'source') || isempty(spec.source)
            spec.source = 'builtin';
        end
        if ~isfield(spec, 'supportsFcn') || isempty(spec.supportsFcn)
            spec.supportsFcn = [];
        end
        if ~isfield(spec, 'isEnabledFcn') || isempty(spec.isEnabledFcn)
            spec.isEnabledFcn = [];
        end

        validatePluginSpec(spec);

        normalized(end+1) = struct( ...
            'id', char(string(spec.id)), ...
            'stage', char(string(spec.stage)), ...
            'displayName', char(string(spec.displayName)), ...
            'priority', double(spec.priority), ...
            'version', char(string(spec.version)), ...
            'source', char(string(spec.source)), ...
            'supportsFcn', spec.supportsFcn, ...
            'isEnabledFcn', spec.isEnabledFcn, ...
            'run', spec.run); %#ok<AGROW>
    end

    registry = normalized;
end

function validatePluginSpec(spec)
    required = {'id', 'stage', 'run'};
    for r = 1:numel(required)
        f = required{r};
        if ~isfield(spec, f) || isempty(spec.(f))
            error('runPipeline:InvalidPluginSpec', 'Plugin spec missing required field "%s".', f);
        end
    end

    if ~isa(spec.run, 'function_handle')
        error('runPipeline:InvalidPluginSpec', 'Plugin "%s" has invalid run handle.', char(string(spec.id)));
    end

    if isfield(spec, 'supportsFcn') && ~isempty(spec.supportsFcn) && ~isa(spec.supportsFcn, 'function_handle')
        error('runPipeline:InvalidPluginSpec', 'Plugin "%s" supportsFcn must be a function handle.', char(string(spec.id)));
    end

    if isfield(spec, 'isEnabledFcn') && ~isempty(spec.isEnabledFcn) && ~isa(spec.isEnabledFcn, 'function_handle')
        error('runPipeline:InvalidPluginSpec', 'Plugin "%s" isEnabledFcn must be a function handle.', char(string(spec.id)));
    end
end

function plugin = selectPlugin(registry, stageName, state, context)
    plugin = [];
    if isempty(registry)
        return;
    end

    stageMask = strcmp({registry.stage}, stageName);
    candidates = registry(stageMask);
    if isempty(candidates)
        return;
    end

    usable = false(1, numel(candidates));
    for i = 1:numel(candidates)
        supports = true;
        if ~isempty(candidates(i).supportsFcn)
            supports = evaluatePredicate(candidates(i).supportsFcn, state, context, true);
        end
        if ~supports
            continue;
        end

        enabled = true;
        if ~isempty(candidates(i).isEnabledFcn)
            enabled = evaluatePredicate(candidates(i).isEnabledFcn, state, context, true);
        end

        usable(i) = logical(enabled);
    end

    candidates = candidates(usable);
    if isempty(candidates)
        return;
    end

    priorities = [candidates.priority];
    [~, idx] = max(priorities);
    plugin = candidates(idx(1));
end

function tf = evaluatePredicate(fcn, state, context, fallback)
    tf = fallback;
    try
        tf = logical(fcn(state, context));
    catch
        tf = fallback;
    end
end

function stageOrder = normalizeStageOrder(stageOrder)
    if isempty(stageOrder)
        stageOrder = localDefaultStageOrder();
        return;
    end

    if ischar(stageOrder)
        stageOrder = {stageOrder};
    elseif isstring(stageOrder)
        stageOrder = cellstr(stageOrder(:));
    elseif iscell(stageOrder)
        stageOrder = cellfun(@char, cellstr(string(stageOrder)), 'UniformOutput', false);
    else
        error('runPipeline:InvalidStageOrder', 'StageOrder must be char/string/cellstr.');
    end

    stageOrder = stageOrder(:)';
end

function stageOrder = localDefaultStageOrder()
    stageOrder = { ...
        'signal_processing', ...
        'maxima_detection', ...
        'gaussian_fitting', ...
        'fit_filtering', ...
        'classification', ...
        'nuclei_filtering'};
end

function tf = localAbortRequested(context)
    tf = false;
    if ~isfield(context, 'abortCheckFcn') || isempty(context.abortCheckFcn) || ~isa(context.abortCheckFcn, 'function_handle')
        return;
    end
    try
        tf = logical(context.abortCheckFcn());
    catch
        tf = false;
    end
end

function specs = discoverExternalPlugins(context)
    specs = struct([]);

    pluginDirs = {};

    if isfield(context, 'externalSignalPluginDirs') && ~isempty(context.externalSignalPluginDirs)
        dirs = context.externalSignalPluginDirs;
        if ischar(dirs) || isstring(dirs)
            dirs = cellstr(string(dirs));
        end
        pluginDirs = [pluginDirs, dirs(:)']; %#ok<AGROW>
    end

    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(fileparts(thisFile)));
    pluginDirs{end+1} = fullfile(repoRoot, 'external_plugins', 'signal');

    envPath = getenv('SNAP_EXTERNAL_SIGNAL_PLUGIN_PATH');
    if ~isempty(envPath)
        envDirs = strsplit(envPath, pathsep);
        pluginDirs = [pluginDirs, envDirs]; %#ok<AGROW>
    end

    pluginDirs = unique(pluginDirs, 'stable');

    for d = 1:numel(pluginDirs)
        pluginDir = pluginDirs{d};
        if isempty(pluginDir) || ~isfolder(pluginDir)
            continue;
        end

        if isempty(strfind(path, pluginDir)) %#ok<STREMP>
            addpath(pluginDir, '-begin');
        end

        files = dir(fullfile(pluginDir, '*.m'));
        for f = 1:numel(files)
            [~, funcName] = fileparts(files(f).name);

            lowerName = lower(funcName);
            if startsWith(funcName, '_') || contains(lowerName, 'template')
                continue;
            end

            try
                produced = feval(str2func(funcName));
            catch ME
                warning('External signal plugin "%s" failed to load: %s', funcName, ME.message);
                continue;
            end

            if isempty(produced) || ~isstruct(produced)
                warning('External signal plugin "%s" did not return a struct spec.', funcName);
                continue;
            end

            for k = 1:numel(produced)
                produced(k).source = sprintf('external:%s', pluginDir);
                specs(end+1) = produced(k); %#ok<AGROW>
            end
        end
    end
end
