function result = readLibrary(filePath, varargin)
% readLibrary - Read microscopy library metadata/images via pluggable providers.
%
% Inputs:
%   filePath - source library file path
%
% Name-value options:
%   'ProgressCallback' - function handle accepting one string
%   'ProviderRegistry' - struct array of provider specs
%   'IncludeExternal'  - include external providers (default: true)
%   'ExternalDirs'     - additional provider directories
%
% Provider spec fields:
%   id, displayName, priority, canReadFcn, readFcn

    p = inputParser;
    p.addParameter('ProgressCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('ProviderRegistry', struct([]), @isstruct);
    p.addParameter('IncludeExternal', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ExternalDirs', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.parse(varargin{:});
    opts = p.Results;

    if ~(ischar(filePath) || isstring(filePath))
        error('readLibrary:InvalidInput', 'filePath must be char or string.');
    end
    filePath = char(string(filePath));

    if ~isfile(filePath)
        error('readLibrary:FileNotFound', 'Library file not found: %s', filePath);
    end

    providers = opts.ProviderRegistry;
    if isempty(providers)
        providers = [snap_modules.prepare.registerBuiltinProviders(), ...
                     discoverExternalProviders(logical(opts.IncludeExternal), opts.ExternalDirs)];
    end
    providers = normalizeProviders(providers);

    if isempty(providers)
        error('readLibrary:NoProviders', 'No library reader providers are available.');
    end

    canReadMask = false(1, numel(providers));
    for i = 1:numel(providers)
        canReadMask(i) = evaluateCanRead(providers(i), filePath);
    end

    providers = providers(canReadMask);
    if isempty(providers)
        error('readLibrary:NoMatchingProvider', 'No reader provider can handle file: %s', filePath);
    end

    [~, order] = sort([providers.priority], 'descend');
    providers = providers(order);

    errors = strings(0, 1);
    for i = 1:numel(providers)
        provider = providers(i);
        emit(opts.ProgressCallback, 'Library read provider "%s" selected.', provider.id);
        try
            out = provider.readFcn(filePath, opts.ProgressCallback);
            out = validateReadResult(out, provider.id);
            out.providerId = provider.id;
            out.providerDisplayName = provider.displayName;
            out.providerVersion = provider.version;
            result = out;
            emit(opts.ProgressCallback, 'Library read completed with provider "%s".', provider.id);
            return;
        catch ME
            errors(end+1) = sprintf('%s: %s', provider.id, ME.message); %#ok<AGROW>
            emit(opts.ProgressCallback, 'Provider "%s" failed: %s', provider.id, ME.message);
        end
    end

    error('readLibrary:AllProvidersFailed', 'All providers failed to read file %s:\n%s', ...
        filePath, strjoin(cellstr(errors), '\n'));
end

function providers = discoverExternalProviders(includeExternal, externalDirs)
    providers = struct([]);
    if ~includeExternal
        return;
    end

    dirs = {};
    if ~isempty(externalDirs)
        if ischar(externalDirs) || isstring(externalDirs)
            dirs = [dirs, cellstr(string(externalDirs))]; %#ok<AGROW>
        else
            dirs = [dirs, cellstr(string(externalDirs(:)'))]; %#ok<AGROW>
        end
    end

    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(fileparts(thisFile)));
    dirs{end+1} = fullfile(repoRoot, 'external_plugins', 'prepare', 'readers');

    envDirs = getenv('SNAP_EXTERNAL_PREPARE_PROVIDER_PATH');
    if ~isempty(envDirs)
        dirs = [dirs, strsplit(envDirs, pathsep)]; %#ok<AGROW>
    end

    dirs = unique(dirs, 'stable');

    for d = 1:numel(dirs)
        pluginDir = dirs{d};
        if isempty(pluginDir) || ~isfolder(pluginDir)
            continue;
        end

        if isempty(strfind(path, pluginDir)) %#ok<STREMP>
            addpath(pluginDir, '-begin');
        end

        files = dir(fullfile(pluginDir, '*.m'));
        for f = 1:numel(files)
            [~, name] = fileparts(files(f).name);
            lname = lower(name);
            if startsWith(name, '_') || contains(lname, 'template')
                continue;
            end

            try
                produced = feval(str2func(name));
            catch ME
                warning('External prepare provider "%s" failed to load: %s', name, ME.message);
                continue;
            end

            if ~isstruct(produced) || isempty(produced)
                warning('External prepare provider "%s" did not return a struct spec.', name);
                continue;
            end

            for k = 1:numel(produced)
                produced(k).source = sprintf('external:%s', pluginDir);
                providers(end+1) = produced(k); %#ok<AGROW>
            end
        end
    end
end

function providers = normalizeProviders(providers)
    if isempty(providers)
        providers = struct([]);
        return;
    end

    normalized = struct('id', {}, 'displayName', {}, 'priority', {}, 'version', {}, ...
                        'source', {}, 'canReadFcn', {}, 'readFcn', {});

    for i = 1:numel(providers)
        spec = providers(i);

        if ~isfield(spec, 'id') || isempty(spec.id)
            spec.id = sprintf('provider_%d', i);
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
        if ~isfield(spec, 'canReadFcn') || ~isa(spec.canReadFcn, 'function_handle')
            error('readLibrary:InvalidProviderSpec', 'Provider "%s" missing valid canReadFcn.', char(string(spec.id)));
        end
        if ~isfield(spec, 'readFcn') || ~isa(spec.readFcn, 'function_handle')
            error('readLibrary:InvalidProviderSpec', 'Provider "%s" missing valid readFcn.', char(string(spec.id)));
        end

        normalized(end+1) = struct( ...
            'id', char(string(spec.id)), ...
            'displayName', char(string(spec.displayName)), ...
            'priority', double(spec.priority), ...
            'version', char(string(spec.version)), ...
            'source', char(string(spec.source)), ...
            'canReadFcn', spec.canReadFcn, ...
            'readFcn', spec.readFcn); %#ok<AGROW>
    end

    providers = normalized;
end

function tf = evaluateCanRead(provider, filePath)
    tf = false;
    try
        tf = logical(provider.canReadFcn(filePath));
    catch
        tf = false;
    end
end

function out = validateReadResult(out, providerId)
    if ~isstruct(out)
        error('Provider "%s" returned non-struct output.', providerId);
    end

    required = {'imageData', 'tableData', 'numChannels'};
    for i = 1:numel(required)
        f = required{i};
        if ~isfield(out, f)
            error('Provider "%s" output missing field "%s".', providerId, f);
        end
    end

    if ~iscell(out.imageData)
        error('Provider "%s" returned invalid imageData (must be cell).', providerId);
    end
    if ~iscell(out.tableData)
        error('Provider "%s" returned invalid tableData (must be cell).', providerId);
    end
    if ~isnumeric(out.numChannels) || ~isscalar(out.numChannels) || ~isfinite(out.numChannels) || out.numChannels < 0
        error('Provider "%s" returned invalid numChannels.', providerId);
    end
end

function emit(progressCb, msgFmt, varargin)
    if isempty(progressCb) || ~isa(progressCb, 'function_handle')
        return;
    end
    try
        progressCb(sprintf(msgFmt, varargin{:}));
    catch
    end
end
