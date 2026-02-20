function spec = bioformats_exporter_builtin()
% Built-in Bio-Formats channel exporter for SNAP_prepare.

    spec = struct();
    spec.id = 'builtin.prepare.bioformats_exporter';
    spec.displayName = 'Bio-Formats Channel Exporter';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.canExportFcn = @canExport;
    spec.exportFcn = @exportChannels;
end

function tf = canExport(imageInfo, ~, ~)
    tf = false;
    if ~isstruct(imageInfo)
        return;
    end

    required = {'libraryPath', 'seriesIndex', 'sizeZ', 'sizeT', 'channelNames'};
    for i = 1:numel(required)
        if ~isfield(imageInfo, required{i})
            return;
        end
    end

    if exist('bfGetReader', 'file') ~= 2 || exist('bfGetPlane', 'file') ~= 2
        return;
    end

    tf = isfile(imageInfo.libraryPath);
end

function exportChannels(imageInfo, mapping, outputPath, progressCb)
    channelNames = imageInfo.channelNames;
    sizeZ = imageInfo.sizeZ;
    sizeT = imageInfo.sizeT;

    getChannelIndex = @(name) findChannelIndex(channelNames, name);

    if ~strcmp(mapping.dic, 'None')
        idx = getChannelIndex(mapping.dic);
        if ~isempty(idx)
            emit(progressCb, 'Exporting DIC channel (%s)...', mapping.dic);
            exportSingleChannelBF(imageInfo.libraryPath, imageInfo.seriesIndex, idx, sizeZ, sizeT, fullfile(outputPath, 'dic.tif'));
        end
    end

    if ~strcmp(mapping.nuclei, 'None')
        idx = getChannelIndex(mapping.nuclei);
        if ~isempty(idx)
            emit(progressCb, 'Exporting nuclei channel (%s)...', mapping.nuclei);
            exportSingleChannelBF(imageInfo.libraryPath, imageInfo.seriesIndex, idx, sizeZ, sizeT, fullfile(outputPath, 'nuclei.tif'));
        end
    end

    fluorescentCount = 0;
    nFluo = numel(mapping.fluorescence);
    for i = 1:nFluo
        selected = mapping.fluorescence{i};
        if strcmp(selected, 'None')
            continue;
        end
        idx = getChannelIndex(selected);
        if isempty(idx)
            continue;
        end

        fluorescentCount = fluorescentCount + 1;
        outName = sprintf('channel%d.tif', fluorescentCount);
        emit(progressCb, 'Exporting fluorescence channel %d (%s)...', fluorescentCount, selected);
        exportSingleChannelBF(imageInfo.libraryPath, imageInfo.seriesIndex, idx, sizeZ, sizeT, fullfile(outputPath, outName));
    end
end

function idx = findChannelIndex(channelNames, selectedName)
    idx = [];
    if isempty(selectedName) || strcmp(selectedName, 'None')
        return;
    end

    idx = find(strcmp(channelNames, selectedName), 1);
    if isempty(idx)
        idx = find(strcmpi(channelNames, selectedName), 1);
    end
end

function exportSingleChannelBF(libraryPath, seriesIdx, channelIdx, sizeZ, sizeT, outputPath)
    reader = bfGetReader(libraryPath);
    cleanupObj = onCleanup(@() safeClose(reader)); %#ok<NASGU>

    reader.setSeries(seriesIdx - 1);

    for t = 1:sizeT
        for z = 1:sizeZ
            iPlane = reader.getIndex(z - 1, channelIdx - 1, t - 1) + 1;
            img = bfGetPlane(reader, iPlane);

            if z == 1 && t == 1
                imwrite(img, outputPath, 'tif', 'Compression', 'none');
            else
                imwrite(img, outputPath, 'tif', 'WriteMode', 'append', 'Compression', 'none');
            end
        end
    end
end

function safeClose(reader)
    try
        if ~isempty(reader)
            reader.close();
        end
    catch
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
