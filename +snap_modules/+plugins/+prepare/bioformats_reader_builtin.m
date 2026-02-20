function spec = bioformats_reader_builtin()
% Built-in Bio-Formats library reader provider for SNAP_prepare.

    spec = struct();
    spec.id = 'builtin.prepare.bioformats_reader';
    spec.displayName = 'Bio-Formats Reader';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.canReadFcn = @canRead;
    spec.readFcn = @readLibrary;
end

function tf = canRead(filePath)
    tf = false;
    if ~(ischar(filePath) || isstring(filePath))
        return;
    end
    filePath = char(string(filePath));
    if ~isfile(filePath)
        return;
    end
    tf = exist('bfopen', 'file') == 2;
end

function out = readLibrary(filePath, progressCb)
    emit(progressCb, 'Loading image library: %s', filePath);
    data = bfopen(filePath);
    numSeries = size(data, 1);
    emit(progressCb, 'Found %d series in library.', numSeries);

    imageData = cell(numSeries, 1);
    tableData = cell(numSeries, 4);

    for i = 1:numSeries
        emit(progressCb, 'Reading series %d/%d...', i, numSeries);

        seriesData = data{i, 1};
        seriesMetadata = data{i, 4};
        seriesLabel = data{i, 2};

        metadata = struct();
        metadata.seriesIndex = i;
        metadata.libraryPath = filePath;

        metadata.name = extractSeriesName(seriesLabel, i);

        numPlanes = size(seriesData, 1);
        if numPlanes > 0
            firstImage = seriesData{1, 1};
            metadata.sizeX = size(firstImage, 2);
            metadata.sizeY = size(firstImage, 1);
        else
            metadata.sizeX = 0;
            metadata.sizeY = 0;
        end

        try
            metadata.sizeZ = double(seriesMetadata.getPixelsSizeZ(i - 1).getValue());
            metadata.sizeC = double(seriesMetadata.getPixelsSizeC(i - 1).getValue());
            metadata.sizeT = double(seriesMetadata.getPixelsSizeT(i - 1).getValue());
        catch
            metadata.sizeZ = 1;
            metadata.sizeC = numPlanes;
            metadata.sizeT = 1;
        end

        metadata.channelNames = cell(1, metadata.sizeC);
        for c = 1:metadata.sizeC
            try
                channelName = char(seriesMetadata.getChannelName(i - 1, c - 1));
                if ~isempty(channelName) && ~strcmp(strtrim(channelName), '')
                    metadata.channelNames{c} = char(channelName);
                else
                    metadata.channelNames{c} = sprintf('Channel %d', c);
                end
            catch
                metadata.channelNames{c} = sprintf('Channel %d', c);
            end
        end

        metadata.rawData = seriesData;
        metadata.rawMetadata = seriesMetadata;

        imageData{i} = metadata;

        tableData{i, 1} = false;
        tableData{i, 2} = char(metadata.name);
        tableData{i, 3} = sprintf('%dx%dx%d', metadata.sizeX, metadata.sizeY, metadata.sizeZ);
        tableData{i, 4} = sprintf('%d', metadata.sizeC);
    end

    numChannels = 0;
    if numSeries > 0
        numChannels = imageData{1}.sizeC;
    end

    out = struct();
    out.imageData = imageData;
    out.tableData = tableData;
    out.numChannels = numChannels;
end

function name = extractSeriesName(seriesLabel, idx)
    if ischar(seriesLabel)
        rawName = seriesLabel;
    elseif isstring(seriesLabel)
        rawName = char(seriesLabel);
    else
        rawName = char(string(seriesLabel));
    end

    name = extractNameFromString(rawName);
    if isempty(name) || strcmp(strtrim(name), '')
        name = sprintf('Image_%03d', idx);
    end
end

function extractedName = extractNameFromString(rawString)
    if isempty(rawString)
        extractedName = '';
        return;
    end

    patterns = {
        'name\s*=\s*"([^"]+)"', ...
        'name\s*=\s*''([^'']+)''', ...
        'Name\s*=\s*"([^"]+)"', ...
        'Name\s*=\s*''([^'']+)''', ...
        'NAME\s*=\s*"([^"]+)"', ...
        'NAME\s*=\s*''([^'']+)''', ...
        'name\s*=\s*([^,;\s]+)', ...
        'Name\s*=\s*([^,;\s]+)', ...
        'NAME\s*=\s*([^,;\s]+)'};

    for p = 1:numel(patterns)
        tokens = regexp(rawString, patterns{p}, 'tokens', 'once', 'ignorecase');
        if ~isempty(tokens)
            extractedName = strtrim(tokens{1});
            return;
        end
    end

    extractedName = rawString;
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
