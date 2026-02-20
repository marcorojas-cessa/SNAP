function exportMappedChannels(imageInfo, mapping, outputPath, varargin)
% exportMappedChannels - Export mapped channels via pluggable exporters.
%
% Inputs:
%   imageInfo  - metadata entry for a single image/series
%   mapping    - struct with fields: dic, nuclei, fluorescence{1..N}
%   outputPath - destination folder path
%
% Name-value options:
%   'ProgressCallback' - function handle accepting one string
%   'ExporterRegistry' - struct array of exporter specs
%   'IncludeExternal'  - include external exporters (default: true)
%   'ExternalDirs'     - additional exporter directories

    p = inputParser;
    p.addParameter('ProgressCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('ExporterRegistry', struct([]), @isstruct);
    p.addParameter('IncludeExternal', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ExternalDirs', {}, @(x) iscell(x) || isstring(x) || ischar(x));
    p.parse(varargin{:});
    opts = p.Results;

    if ~isstruct(imageInfo)
        error('exportMappedChannels:InvalidInput', 'imageInfo must be a struct.');
    end
    if ~isstruct(mapping)
        error('exportMappedChannels:InvalidInput', 'mapping must be a struct.');
    end
    if ~(ischar(outputPath) || isstring(outputPath))
        error('exportMappedChannels:InvalidInput', 'outputPath must be char or string.');
    end

    outputPath = char(string(outputPath));
    if ~isfolder(outputPath)
        mkdir(outputPath);
    end

    exporters = opts.ExporterRegistry;
    if isempty(exporters)
        exporters = [snap_modules.prepare.registerBuiltinExporters(), ...
                     discoverExternalExporters(logical(opts.IncludeExternal), opts.ExternalDirs)];
    end
    exporters = normalizeExporters(exporters);

    if isempty(exporters)
        error('exportMappedChannels:NoExporters', 'No channel exporters are available.');
    end

    canExportMask = false(1, numel(exporters));
    for i = 1:numel(exporters)
        canExportMask(i) = evaluateCanExport(exporters(i), imageInfo, mapping, outputPath);
    end

    exporters = exporters(canExportMask);
    if isempty(exporters)
        error('exportMappedChannels:NoMatchingExporter', 'No exporter can handle this image/mapping.');
    end

    [~, order] = sort([exporters.priority], 'descend');
    exporters = exporters(order);

    errors = strings(0, 1);
    for i = 1:numel(exporters)
        exporter = exporters(i);
        emit(opts.ProgressCallback, 'Channel exporter "%s" selected.', exporter.id);
        try
            exporter.exportFcn(imageInfo, mapping, outputPath, opts.ProgressCallback);
            emit(opts.ProgressCallback, 'Channel export completed with "%s".', exporter.id);
            return;
        catch ME
            errors(end+1) = sprintf('%s: %s', exporter.id, ME.message); %#ok<AGROW>
            emit(opts.ProgressCallback, 'Exporter "%s" failed: %s', exporter.id, ME.message);
        end
    end

    error('exportMappedChannels:AllExportersFailed', 'All exporters failed for output path %s:\n%s', ...
        outputPath, strjoin(cellstr(errors), '\n'));
end

function exporters = discoverExternalExporters(includeExternal, externalDirs)
    exporters = struct([]);
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
    dirs{end+1} = fullfile(repoRoot, 'plugins', 'prepare', 'exporters');

    envDirs = getenv('SNAP_PREPARE_EXPORTER_PATH');
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
                warning('External prepare exporter "%s" failed to load: %s', name, ME.message);
                continue;
            end

            if ~isstruct(produced) || isempty(produced)
                warning('External prepare exporter "%s" did not return a struct spec.', name);
                continue;
            end

            for k = 1:numel(produced)
                produced(k).source = sprintf('external:%s', pluginDir);
                exporters(end+1) = produced(k); %#ok<AGROW>
            end
        end
    end
end

function exporters = normalizeExporters(exporters)
    if isempty(exporters)
        exporters = struct([]);
        return;
    end

    normalized = struct('id', {}, 'displayName', {}, 'priority', {}, 'version', {}, ...
                        'source', {}, 'canExportFcn', {}, 'exportFcn', {});

    for i = 1:numel(exporters)
        spec = exporters(i);

        if ~isfield(spec, 'id') || isempty(spec.id)
            spec.id = sprintf('exporter_%d', i);
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
        if ~isfield(spec, 'canExportFcn') || ~isa(spec.canExportFcn, 'function_handle')
            error('exportMappedChannels:InvalidExporterSpec', 'Exporter "%s" missing valid canExportFcn.', char(string(spec.id)));
        end
        if ~isfield(spec, 'exportFcn') || ~isa(spec.exportFcn, 'function_handle')
            error('exportMappedChannels:InvalidExporterSpec', 'Exporter "%s" missing valid exportFcn.', char(string(spec.id)));
        end

        normalized(end+1) = struct( ...
            'id', char(string(spec.id)), ...
            'displayName', char(string(spec.displayName)), ...
            'priority', double(spec.priority), ...
            'version', char(string(spec.version)), ...
            'source', char(string(spec.source)), ...
            'canExportFcn', spec.canExportFcn, ...
            'exportFcn', spec.exportFcn); %#ok<AGROW>
    end

    exporters = normalized;
end

function tf = evaluateCanExport(exporter, imageInfo, mapping, outputPath)
    tf = false;
    try
        tf = logical(exporter.canExportFcn(imageInfo, mapping, outputPath));
    catch
        tf = false;
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
