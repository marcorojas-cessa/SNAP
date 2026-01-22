function SNAP_prepare
% SNAP_prepare - Prepare multi-channel image libraries for SNAP batch processing
%
% Opens a GUI to:
%   1. Load Bio-Formats image libraries (see supported formats below)
%   2. View available multi-channel images
%   3. Assign channel identities (DIC, Nuclei, Fluorescence 1-6)
%   4. Select which images to export
%   5. Create folder structure compatible with SNAP_batch
%
% Supported File Formats (via Bio-Formats):
%   - MetaMorph Stack (.mvd2)
%   - Zeiss CZI (.czi)
%   - Nikon NIS-Elements (.nd2)
%   - Leica Image File (.lif)
%   - Olympus (.oib, .oif)
%   - Zeiss LSM (.lsm)
%   - Imaris (.ims)
%   - Olympus VSI (.vsi)
%   - Leica SCN (.scn)
%   - TIFF Stack (.tif, .tiff)
%   - And many more supported by Bio-Formats
%
% Output Structure:
%   outputDir/
%   ├── Image_0001/ (or original_name/)
%   │   ├── dic.tif (if DIC channel assigned)
%   │   ├── nuclei.tif (if nuclei channel assigned)
%   │   ├── channel1.tif (first fluorescence channel)
%   │   ├── channel2.tif (second fluorescence channel)
%   │   └── ...
%   ├── Image_0002/
%   │   └── ...
%
% Requires: Bio-Formats Toolbox (available in MATLAB Add-On Explorer)
%   Install from: Home > Add-Ons > Get Add-Ons > Search "Bio-Formats"

% Check for Bio-Formats (will show dialog if not found, but won't prevent GUI from opening)
checkBioFormats();

% Create the main figure
fig = uifigure('Name', 'SNAP Prepare - Image Library Preparation', ...
    'Position', [50 50 1100 700]);

% Main grid layout
mainGrid = uigridlayout(fig, [1, 2]);
mainGrid.ColumnWidth = {'1x', 400};
mainGrid.Padding = [10 10 10 10];

% Initialize handles
handles = struct();
handles.fig = fig;
handles.imageData = []; % Will store loaded image information
handles.numChannels = 0; % Number of channels in library images

% --- Left Panel: Image Library View ---
leftPanel = uipanel(mainGrid, 'Title', 'Image Library', 'FontWeight', 'bold');
leftPanel.Layout.Column = 1;
leftGrid = uigridlayout(leftPanel, [4, 1]);
leftGrid.RowHeight = {'fit', 'fit', '1x', 'fit'};
leftGrid.Padding = [10 10 10 10];
leftGrid.RowSpacing = 10;

% File selection
filePanel = uipanel(leftGrid, 'BorderType', 'none');
fileGrid = uigridlayout(filePanel, [2, 3]);
fileGrid.RowHeight = {'fit', 'fit'};
fileGrid.ColumnWidth = {'1x', 100, 100};
fileGrid.Padding = [0 0 0 0];
fileGrid.ColumnSpacing = 5;

uilabel(fileGrid, 'Text', 'Image Library File:', 'FontWeight', 'bold');
uilabel(fileGrid, 'Text', '');
uilabel(fileGrid, 'Text', '');
handles.filePathEdit = uieditfield(fileGrid, 'Value', '', 'Editable', 'off');
handles.browseFileBtn = uibutton(fileGrid, 'Text', 'Browse...', ...
    'ButtonPushedFcn', @(src,evt) browseImageFile(fig));
handles.readLibraryBtn = uibutton(fileGrid, 'Text', 'Read', ...
    'ButtonPushedFcn', @(src,evt) readImageLibrary(fig), ...
    'Enable', 'off', 'Tooltip', 'Read selected image library');

% Library info
handles.libraryInfoLabel = uilabel(leftGrid, 'Text', 'No file loaded', ...
    'FontColor', [0.5 0.5 0.5]);

% Image list table
handles.imageTable = uitable(leftGrid);
handles.imageTable.ColumnName = {'Select', 'Image Name', 'Size (XYZ)', 'Channels'};
handles.imageTable.ColumnWidth = {50, 'auto', 100, 80};
handles.imageTable.ColumnEditable = [true false false false];
handles.imageTable.Data = {};
handles.imageTable.CellEditCallback = @(src,evt) updateSelectionCount(fig);

% Selection controls
selectionPanel = uipanel(leftGrid, 'BorderType', 'none');
selectionGrid = uigridlayout(selectionPanel, [1, 3]);
selectionGrid.ColumnWidth = {'1x', 100, 100};
selectionGrid.Padding = [0 0 0 0];

handles.selectionLabel = uilabel(selectionGrid, 'Text', 'Selected: 0 images', ...
    'FontWeight', 'bold');
handles.selectAllBtn = uibutton(selectionGrid, 'Text', 'Select All', ...
    'ButtonPushedFcn', @(src,evt) selectAllImages(fig, true));
handles.deselectAllBtn = uibutton(selectionGrid, 'Text', 'Deselect All', ...
    'ButtonPushedFcn', @(src,evt) selectAllImages(fig, false));

% --- Right Panel: Channel Configuration ---
rightPanel = uipanel(mainGrid, 'Title', 'Channel Configuration (Applies to All Images)', 'FontWeight', 'bold');
rightPanel.Layout.Column = 2;
rightGrid = uigridlayout(rightPanel, [6, 1]);
rightGrid.RowHeight = {'fit', 'fit', '1x', 'fit', 'fit', 'fit'};
rightGrid.Padding = [10 10 10 10];
rightGrid.RowSpacing = 10;

% Instructions
instructionsLabel = uilabel(rightGrid, ...
    'Text', 'Assign channel identities (same order for all images):', ...
    'FontWeight', 'bold', 'WordWrap', 'on');

% Folder naming option
namingPanel = uipanel(rightGrid, 'BorderType', 'none');
namingGrid = uigridlayout(namingPanel, [2, 1]);
namingGrid.RowHeight = {'fit', 'fit'};
namingGrid.Padding = [0 0 0 0];
namingGrid.RowSpacing = 3;

uilabel(namingGrid, 'Text', 'Folder Naming:', 'FontWeight', 'bold');
handles.folderNamingDrop = uidropdown(namingGrid, ...
    'Items', {'Use Original Names (sanitized)', 'Use Simple Indices (Image_0001, Image_0002, ...)'}, ...
    'Value', 'Use Original Names (sanitized)', ...
    'Tooltip', 'Choose how to name exported folders');

% Channel mapping panel
channelMapPanel = uipanel(rightGrid, 'Title', 'Channel Assignments', 'FontWeight', 'bold');
channelMapGrid = uigridlayout(channelMapPanel, [8, 1]); % Max: DIC + Nuclei + 6 fluorescence
channelMapGrid.RowHeight = repmat({'fit'}, 1, 8);
channelMapGrid.Padding = [5 5 5 5];
channelMapGrid.RowSpacing = 5;

% Define channel functions (static labels)
channelFunctions = {'DIC', 'Nuclei', 'Fluorescence 1', 'Fluorescence 2', ...
                    'Fluorescence 3', 'Fluorescence 4', 'Fluorescence 5', 'Fluorescence 6'};

% Pre-allocate channel mapping controls
handles.channelFunctionPanels = gobjects(1, 8);
handles.channelFunctionDropdowns = gobjects(1, 8);
handles.channelFunctionLabels = gobjects(1, 8);

for i = 1:8
    panel = uipanel(channelMapGrid, 'BorderType', 'none', 'Visible', 'on');
    panel.Layout.Row = i;
    grid = uigridlayout(panel, [1, 2]);
    grid.ColumnWidth = {120, '1x'};
    grid.Padding = [0 0 0 0];
    
    handles.channelFunctionLabels(i) = uilabel(grid, 'Text', [channelFunctions{i} ':'], ...
        'FontWeight', 'bold');
    handles.channelFunctionDropdowns(i) = uidropdown(grid, ...
        'Items', {'None'}, ...
        'Value', 'None', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(src,evt) updateExportButtonState(fig), ...
        'Tooltip', sprintf('Select which library channel to use for %s', channelFunctions{i}));
    
    handles.channelFunctionPanels(i) = panel;
end

% Output directory
outputPanel = uipanel(rightGrid, 'BorderType', 'none');
outputGrid = uigridlayout(outputPanel, [2, 2]);
outputGrid.RowHeight = {'fit', 'fit'};
outputGrid.ColumnWidth = {'1x', 100};
outputGrid.Padding = [0 0 0 0];

uilabel(outputGrid, 'Text', 'Output Directory:', 'FontWeight', 'bold');
uilabel(outputGrid, 'Text', '');
handles.outputDirEdit = uieditfield(outputGrid, 'Value', '', 'Editable', 'off');
handles.browseOutputBtn = uibutton(outputGrid, 'Text', 'Browse...', ...
    'ButtonPushedFcn', @(src,evt) browseOutputDir(fig));

% Action buttons
actionPanel = uipanel(rightGrid, 'BorderType', 'none');
actionGrid = uigridlayout(actionPanel, [1, 2]);
actionGrid.ColumnWidth = {'1x', '1x'};
actionGrid.Padding = [0 0 0 0];
actionGrid.ColumnSpacing = 10;

handles.exportBtn = uibutton(actionGrid, 'Text', 'Export Selected', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'ButtonPushedFcn', @(src,evt) exportImages(fig), ...
    'Enable', 'off', 'BackgroundColor', [0.2 0.7 0.3], 'FontColor', [1 1 1]);

handles.closeBtn = uibutton(actionGrid, 'Text', 'Close', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'ButtonPushedFcn', @(src,evt) close(fig));

% Store handles
guidata(fig, handles);
end

%% Callback Functions

function browseImageFile(fig)
    handles = guidata(fig);
    
    % File dialog for Bio-Formats compatible files
    [file, path] = uigetfile({...
        '*.mvd2', 'MetaMorph Stack (*.mvd2)'; ...
        '*.czi', 'Zeiss CZI (*.czi)'; ...
        '*.nd2', 'Nikon NIS-Elements (*.nd2)'; ...
        '*.lif', 'Leica Image File (*.lif)'; ...
        '*.oib;*.oif', 'Olympus Image (*.oib, *.oif)'; ...
        '*.lsm', 'Zeiss LSM (*.lsm)'; ...
        '*.ims', 'Imaris (*.ims)'; ...
        '*.vsi', 'Olympus VSI (*.vsi)'; ...
        '*.scn', 'Leica SCN (*.scn)'; ...
        '*.tif;*.tiff', 'TIFF Stack (*.tif, *.tiff)'; ...
        '*.mvd2;*.czi;*.nd2;*.lif;*.oib;*.oif;*.lsm;*.ims;*.vsi;*.scn;*.tif;*.tiff', 'All Supported Formats'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select Image Library File');
    
    if file == 0
        return;
    end
    
    fullPath = fullfile(path, file);
    handles.filePathEdit.Value = fullPath;
    
    % Enable the Read button
    handles.readLibraryBtn.Enable = 'on';
    
    % Update status
    handles.libraryInfoLabel.Text = 'File selected. Click "Read" to load image library.';
    handles.libraryInfoLabel.FontColor = [0 0 0.8];
    
    % Clear any previous data
    handles.imageTable.Data = {};
    handles.imageData = [];
    handles.numChannels = 0;
    
    % Reset all channel function dropdowns
    for i = 1:8
        handles.channelFunctionDropdowns(i).Items = {'None'};
        handles.channelFunctionDropdowns(i).Value = 'None';
        handles.channelFunctionDropdowns(i).Enable = 'off';
    end
    
    guidata(fig, handles);
end

function readImageLibrary(fig)
    % Read the selected image library file
    handles = guidata(fig);
    
    fullPath = handles.filePathEdit.Value;
    
    if isempty(fullPath)
        uialert(fig, 'Please select an image library file first.', 'No File Selected');
        return;
    end
    
    % Disable Read button during loading
    handles.readLibraryBtn.Enable = 'off';
    handles.browseFileBtn.Enable = 'off';
    guidata(fig, handles);
    
    try
        % Load image library using Bio-Formats
        loadImageLibrary(fig, fullPath);
        
        % Check if fig is still valid (user might have closed during load)
        if ~isvalid(fig)
            return;
        end
        
        handles = guidata(fig);
        numImages = size(handles.imageTable.Data, 1);
        numChannels = handles.numChannels;
        handles.libraryInfoLabel.Text = sprintf('✓ Loaded: %d images, %d channels each', numImages, numChannels);
        handles.libraryInfoLabel.FontColor = [0 0.6 0];
        
    catch ME
        % Check if fig is still valid
        if ~isvalid(fig)
            return;
        end
        
        handles = guidata(fig);
        handles.libraryInfoLabel.Text = 'Error loading file (see console)';
        handles.libraryInfoLabel.FontColor = [0.8 0 0];
        
        uialert(fig, sprintf('Failed to load image library:\n\n%s', ME.message), ...
            'Load Error', 'Icon', 'error');
    end
    
    % Re-enable buttons
    handles = guidata(fig);
    handles.browseFileBtn.Enable = 'on';
    % Keep Read button enabled so user can try again
    handles.readLibraryBtn.Enable = 'on';
    guidata(fig, handles);
end

function loadImageLibrary(fig, filePath)
    handles = guidata(fig);
    
    % Use MATLAB's Bio-Formats reader
    try
        % Update status
        handles.libraryInfoLabel.Text = 'Loading image library...';
        handles.libraryInfoLabel.FontColor = [0 0 0.8];
        guidata(fig, handles);
        drawnow;
        
        % Try using bfopen (from MATLAB Bio-Formats toolbox)
        fprintf('Loading image library: %s\n', filePath);
        data = bfopen(filePath);
        numSeries = size(data, 1);
        fprintf('Found %d series in library\n', numSeries);
        
        % Extract metadata from bfopen format
        imageData = cell(numSeries, 1);
        tableData = cell(numSeries, 4);
        
        for i = 1:numSeries
            % Update progress
            fprintf('Reading series %d/%d...', i, numSeries);
            
            % Update GUI if still valid
            if isvalid(fig)
                handles = guidata(fig);
                handles.libraryInfoLabel.Text = sprintf('Reading series %d/%d...', i, numSeries);
                guidata(fig, handles);
                drawnow;
            end
            
            seriesData = data{i, 1}; % Image planes
            seriesMetadata = data{i, 4}; % OME metadata
            seriesLabel = data{i, 2}; % Series name
            
            % Get metadata
            metadata = struct();
            metadata.seriesIndex = i;
            metadata.libraryPath = filePath; % Store path for later reopening
            
            % Convert series label to proper MATLAB char
            if ischar(seriesLabel)
                rawName = seriesLabel;
            elseif isstring(seriesLabel)
                rawName = char(seriesLabel);
            else
                rawName = char(string(seriesLabel));
            end
            
            % Extract name from patterns like "name=..." or "Name=..."
            metadata.name = extractNameFromString(rawName);
            
            % Ensure name is valid
            if isempty(metadata.name) || strcmp(strtrim(metadata.name), '')
                metadata.name = sprintf('Image_%03d', i);
            end
            
            % Parse dimensions from seriesData
            % Each row in seriesData is: {image, label}
            % We need to figure out dimensions from the structure
            numPlanes = size(seriesData, 1);
            
            if numPlanes > 0
                firstImage = seriesData{1, 1};
                metadata.sizeX = size(firstImage, 2);
                metadata.sizeY = size(firstImage, 1);
            else
                metadata.sizeX = 0;
                metadata.sizeY = 0;
            end
            
            % Try to extract Z, C, T from metadata
            try
                metadata.sizeZ = double(seriesMetadata.getPixelsSizeZ(i-1).getValue());
                metadata.sizeC = double(seriesMetadata.getPixelsSizeC(i-1).getValue());
                metadata.sizeT = double(seriesMetadata.getPixelsSizeT(i-1).getValue());
            catch
                % Fallback: assume single Z, multiple channels
                metadata.sizeZ = 1;
                metadata.sizeC = numPlanes;
                metadata.sizeT = 1;
            end
            
            % Get channel names if available
            metadata.channelNames = cell(1, metadata.sizeC);
            for c = 1:metadata.sizeC
                try
                    channelName = char(seriesMetadata.getChannelName(i-1, c-1));
                    if ~isempty(channelName) && ~strcmp(strtrim(channelName), '')
                        metadata.channelNames{c} = char(channelName);
                    else
                        metadata.channelNames{c} = sprintf('Channel %d', c);
                    end
                catch
                    metadata.channelNames{c} = sprintf('Channel %d', c);
                end
            end
            
            % Store raw data for later export
            metadata.rawData = seriesData;
            metadata.rawMetadata = seriesMetadata;
            
            imageData{i} = metadata;
            
            % Populate table row - ENSURE ALL VALUES ARE PROPER MATLAB TYPES
            tableData{i, 1} = false; % Logical
            tableData{i, 2} = char(metadata.name); % Ensure char, not string or Java object
            tableData{i, 3} = sprintf('%d×%d×%d', metadata.sizeX, metadata.sizeY, metadata.sizeZ); % char
            tableData{i, 4} = sprintf('%d', metadata.sizeC); % char
            
            fprintf(' Done.\n');
        end
        
        fprintf('All series loaded successfully.\n');
        
        % Validate table data before setting
        fprintf('Validating table data...\n');
        for i = 1:size(tableData, 1)
            for j = 1:size(tableData, 2)
                val = tableData{i, j};
                if ~(islogical(val) || isnumeric(val) || ischar(val))
                    warning('Invalid table entry at row %d, col %d: class = %s', i, j, class(val));
                    % Try to convert to char
                    if isstring(val)
                        tableData{i, j} = char(val);
                    else
                        tableData{i, j} = char(string(val));
                    end
                    fprintf('  Converted to: %s\n', tableData{i, j});
                end
            end
        end
        fprintf('Table data validated.\n');
        
    catch ME
        % If bfopen fails, throw error with helpful message
        error('Failed to load image library: %s\n\nMake sure you have the Bio-Formats Toolbox installed from MATLAB Add-Ons.', ME.message);
    end
    
    % Update handles
    handles = guidata(fig);
    handles.imageData = imageData;
    
    % Set table data with error handling
    try
        handles.imageTable.Data = tableData;
    catch ME2
        fprintf('Error setting table data: %s\n', ME2.message);
        fprintf('Attempting to diagnose issue...\n');
        for i = 1:min(3, size(tableData, 1))
            fprintf('Row %d: [%s, %s, %s, %s]\n', i, ...
                class(tableData{i,1}), class(tableData{i,2}), class(tableData{i,3}), class(tableData{i,4}));
        end
        rethrow(ME2);
    end
    
    handles.filePath = filePath;
    
    % Determine number of channels (use first image as reference)
    if numSeries > 0
        numChannels = imageData{1}.sizeC;
        handles.numChannels = numChannels;
        
        fprintf('Configuring channel assignments for %d channels...\n', numChannels);
        
        % Build dropdown items: 'None' + all channel names
        dropdownItems = {'None'};
        for i = 1:numChannels
            channelName = imageData{1}.channelNames{i};
            if ~isempty(channelName) && ~strcmp(strtrim(channelName), '')
                dropdownItems{end+1} = char(channelName);
            else
                dropdownItems{end+1} = sprintf('Channel %d', i);
            end
            fprintf('  Library channel %d: %s\n', i, dropdownItems{end});
        end
        
        % Update all function dropdowns with these items
        for i = 1:8
            handles.channelFunctionDropdowns(i).Items = dropdownItems;
            handles.channelFunctionDropdowns(i).Value = 'None';
            handles.channelFunctionDropdowns(i).Enable = 'on';
        end
        
        fprintf('Channel configuration ready. Assign library channels to functions on the right panel.\n');
    else
        handles.numChannels = 0;
    end
    
    guidata(fig, handles);
    updateSelectionCount(fig);
    updateExportButtonState(fig);
end

function updateSelectionCount(fig)
    handles = guidata(fig);
    
    if isempty(handles.imageTable.Data)
        handles.selectionLabel.Text = 'Selected: 0 images';
        return;
    end
    
    tableData = handles.imageTable.Data;
    numSelected = sum(cell2mat(tableData(:, 1)));
    handles.selectionLabel.Text = sprintf('Selected: %d images', numSelected);
    
    guidata(fig, handles);
    updateExportButtonState(fig);
end

function selectAllImages(fig, selectState)
    handles = guidata(fig);
    
    if isempty(handles.imageTable.Data)
        return;
    end
    
    tableData = handles.imageTable.Data;
    tableData(:, 1) = {selectState};
    handles.imageTable.Data = tableData;
    
    guidata(fig, handles);
    updateSelectionCount(fig);
end

function browseOutputDir(fig)
    handles = guidata(fig);
    selectedDir = uigetdir('', 'Select Output Directory');
    if selectedDir ~= 0
        handles.outputDirEdit.Value = selectedDir;
        guidata(fig, handles);
        updateExportButtonState(fig);
    end
end

function updateExportButtonState(fig)
    handles = guidata(fig);
    
    % Enable export if: images selected, output dir set, at least one non-None mapping
    hasSelection = false;
    hasMapping = false;
    hasOutputDir = ~isempty(handles.outputDirEdit.Value);
    
    if ~isempty(handles.imageTable.Data)
        tableData = handles.imageTable.Data;
        selectedIndices = find(cell2mat(tableData(:, 1)));
        hasSelection = ~isempty(selectedIndices);
    end
    
    % Check if any function has a channel assigned (not 'None')
    for i = 1:8
        if ~strcmp(handles.channelFunctionDropdowns(i).Value, 'None')
            hasMapping = true;
            break;
        end
    end
    
    if hasSelection && hasMapping && hasOutputDir
        handles.exportBtn.Enable = 'on';
    else
        handles.exportBtn.Enable = 'off';
    end
    
    guidata(fig, handles);
end

function exportImages(fig)
    handles = guidata(fig);
    
    % Get selected images
    tableData = handles.imageTable.Data;
    selectedIndices = find(cell2mat(tableData(:, 1)));
    
    if isempty(selectedIndices)
        uialert(fig, 'No images selected for export.', 'No Selection');
        return;
    end
    
    outputDir = handles.outputDirEdit.Value;
    if isempty(outputDir)
        uialert(fig, 'Please select an output directory.', 'No Output Directory');
        return;
    end
    
    % Get folder naming preference
    useSimpleNames = strcmp(handles.folderNamingDrop.Value, 'Use Simple Indices (Image_0001, Image_0002, ...)');
    
    % Build channel mapping from dropdowns
    % For each function (DIC, Nuclei, Fluorescence 1-6), get which channel is assigned
    channelMapping = struct();
    channelMapping.dic = handles.channelFunctionDropdowns(1).Value;
    channelMapping.nuclei = handles.channelFunctionDropdowns(2).Value;
    channelMapping.fluorescence = cell(1, 6);
    for i = 1:6
        channelMapping.fluorescence{i} = handles.channelFunctionDropdowns(i+2).Value;
    end
    
    % Check if any channels are assigned
    hasAssignment = ~strcmp(channelMapping.dic, 'None') || ...
                    ~strcmp(channelMapping.nuclei, 'None') || ...
                    any(~strcmp(channelMapping.fluorescence, 'None'));
    
    if ~hasAssignment
        uialert(fig, 'No channels assigned. Please assign at least one channel.', ...
            'No Channels Assigned');
        return;
    end
    
    % Disable UI during export
    handles.exportBtn.Enable = 'off';
    handles.browseFileBtn.Enable = 'off';
    handles.browseOutputBtn.Enable = 'off';
    handles.selectAllBtn.Enable = 'off';
    handles.deselectAllBtn.Enable = 'off';
    handles.readLibraryBtn.Enable = 'off';
    guidata(fig, handles);
    
    % Create progress dialog
    progressDlg = uiprogressdlg(fig, 'Title', 'Exporting Images', ...
        'Message', 'Initializing...', 'Cancelable', true);
    
    try
        numExported = 0;
        
        for i = 1:length(selectedIndices)
            if progressDlg.CancelRequested
                break;
            end
            
            idx = selectedIndices(i);
            imageInfo = handles.imageData{idx};
            
            % Update progress
            progressDlg.Value = i / length(selectedIndices);
            progressDlg.Message = sprintf('Processing %d/%d: %s', i, length(selectedIndices), imageInfo.name);
            
            % Create folder for this image
            if useSimpleNames
                % Use simple index-based naming
                safeName = sprintf('Image_%04d', i); % Use loop index for sequential numbering
            else
                % Use original name with sanitization
                safeName = createSafeFolderName(imageInfo.name, idx);
            end
            
            imageFolderPath = fullfile(outputDir, safeName);
            
            % Double-check path length (filesystem limitation)
            if length(imageFolderPath) > 200
                fprintf('  Warning: Path too long (%d chars), using simplified name\n', length(imageFolderPath));
                safeName = sprintf('Image_%04d', i);
                imageFolderPath = fullfile(outputDir, safeName);
            end
            
            % Log folder creation
            if ~useSimpleNames && ~strcmp(safeName, imageInfo.name)
                fprintf('  Creating folder: %s → %s\n', imageInfo.name, safeName);
            else
                fprintf('  Creating folder: %s\n', safeName);
            end
            
            if ~exist(imageFolderPath, 'dir')
                mkdir(imageFolderPath);
            end
            
            % Export each assigned channel using channel mapping
            exportImageChannels(imageInfo, channelMapping, imageFolderPath);
            numExported = numExported + 1;
        end
        
        close(progressDlg);
        
        % Success message
        uialert(fig, sprintf('Export complete!\n\n%d images exported to:\n%s', ...
            numExported, outputDir), 'Export Complete', 'Icon', 'success');
        
    catch ME
        if isvalid(progressDlg)
            close(progressDlg);
        end
        uialert(fig, sprintf('Error during export:\n\n%s', ME.message), ...
            'Export Error', 'Icon', 'error');
    end
    
    % Re-enable UI
    handles = guidata(fig);
    handles.exportBtn.Enable = 'on';
    handles.browseFileBtn.Enable = 'on';
    handles.browseOutputBtn.Enable = 'on';
    handles.selectAllBtn.Enable = 'on';
    handles.deselectAllBtn.Enable = 'on';
    handles.readLibraryBtn.Enable = 'on';
    guidata(fig, handles);
end

function exportImageChannels(imageInfo, mapping, outputPath)
    % Export channels according to inverted mapping
    % mapping is a struct with: .dic, .nuclei, .fluorescence{1:6}
    % Each field contains the channel name from library (or 'None')
    
    rawData = imageInfo.rawData;
    channelNames = imageInfo.channelNames;
    sizeZ = imageInfo.sizeZ;
    sizeC = imageInfo.sizeC;
    sizeT = imageInfo.sizeT;
    
    % Helper function to find channel index by name
    findChannelIndex = @(name) find(strcmp(channelNames, name), 1);
    
    % Use Bio-Formats reader for proper channel extraction
    libraryPath = imageInfo.libraryPath;
    seriesIdx = imageInfo.seriesIndex;
    
    % Export DIC if assigned
    if ~strcmp(mapping.dic, 'None')
        channelIdx = findChannelIndex(mapping.dic);
        if ~isempty(channelIdx)
            exportSingleChannelBF(libraryPath, seriesIdx, channelIdx, sizeZ, sizeT, ...
                fullfile(outputPath, 'dic.tif'));
        end
    end
    
    % Export Nuclei if assigned
    if ~strcmp(mapping.nuclei, 'None')
        channelIdx = findChannelIndex(mapping.nuclei);
        if ~isempty(channelIdx)
            exportSingleChannelBF(libraryPath, seriesIdx, channelIdx, sizeZ, sizeT, ...
                fullfile(outputPath, 'nuclei.tif'));
        end
    end
    
    % Export Fluorescence channels if assigned
    fluorescentCount = 0;
    for i = 1:6
        if ~strcmp(mapping.fluorescence{i}, 'None')
            channelIdx = findChannelIndex(mapping.fluorescence{i});
            if ~isempty(channelIdx)
                fluorescentCount = fluorescentCount + 1;
                filename = sprintf('channel%d.tif', fluorescentCount);
                exportSingleChannelBF(libraryPath, seriesIdx, channelIdx, sizeZ, sizeT, ...
                    fullfile(outputPath, filename));
            end
        end
    end
end

function exportSingleChannelBF(libraryPath, seriesIdx, channelIdx, sizeZ, sizeT, outputPath)
    % Export a single channel using Bio-Formats native reader
    % This ensures proper Z-stack extraction without plane interleaving issues
    
    % Create a Bio-Formats reader
    reader = bfGetReader(libraryPath);
    reader.setSeries(seriesIdx - 1); % Bio-Formats uses 0-based indexing
    
    try
        % Export all Z slices for this channel
        for t = 1:sizeT
            for z = 1:sizeZ
                % Bio-Formats getIndex: converts Z, C, T to linear plane index
                % This handles dimension ordering automatically!
                iPlane = reader.getIndex(z - 1, channelIdx - 1, t - 1) + 1; % Convert to 1-based
                
                % Read the plane using bfGetPlane
                img = bfGetPlane(reader, iPlane);
                
                % Write to TIFF (append mode for z-stack)
                if z == 1 && t == 1
                    imwrite(img, outputPath, 'tif', 'Compression', 'none');
                else
                    imwrite(img, outputPath, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                end
            end
        end
    catch ME
        reader.close();
        rethrow(ME);
    end
    
    % Clean up
    reader.close();
end

function exportSingleChannel(rawData, channelIdx, sizeZ, sizeC, sizeT, outputPath)
    % OLD FUNCTION - kept for reference but no longer used
    % Export a single channel to TIFF file
    % rawData: cell array from bfopen
    % channelIdx: which channel to export (1-based)
    
    % Bio-Formats plane ordering: Parse from labels to determine dimension order
    % Labels are typically: "Z:1/10; C:1/3; T:1/1" or similar
    
    % Detect dimension order from first few planes
    dimensionOrder = detectDimensionOrder(rawData, sizeZ, sizeC, sizeT);
    fprintf('  Detected dimension order: %s\n', dimensionOrder);
    
    for t = 1:sizeT
        for z = 1:sizeZ
            % Calculate plane index based on detected dimension order
            iPlane = calculatePlaneIndex(z, channelIdx, t, sizeZ, sizeC, sizeT, dimensionOrder);
            
            if iPlane < 1 || iPlane > size(rawData, 1)
                warning('Plane index %d out of bounds (1-%d)', iPlane, size(rawData, 1));
                continue;
            end
            
            img = rawData{iPlane, 1}; % Get image from {image, label} pair
            
            % Write to TIFF (append mode for z-stack)
            if z == 1 && t == 1
                imwrite(img, outputPath, 'tif', 'Compression', 'none');
            else
                imwrite(img, outputPath, 'tif', 'WriteMode', 'append', 'Compression', 'none');
            end
        end
    end
end

function dimOrder = detectDimensionOrder(rawData, sizeZ, sizeC, sizeT)
    % Detect dimension order from Bio-Formats plane labels
    % Returns 'XYZCT', 'XYZTC', 'XYCZT', 'XYCTZ', 'XYTCZ', or 'XYTZC'
    
    % Default assumption
    dimOrder = 'XYZCT'; % Z fastest, then C, then T
    
    % Check first few planes to detect which dimension varies first
    if size(rawData, 1) >= 2
        % Parse first two plane labels
        label1 = rawData{1, 2};
        label2 = rawData{2, 2};
        
        % Extract Z, C, T values from labels (format: "Z:1/10; C:1/3; T:1/1")
        z1 = extractDimValue(label1, 'Z');
        c1 = extractDimValue(label1, 'C');
        t1 = extractDimValue(label1, 'T');
        
        z2 = extractDimValue(label2, 'Z');
        c2 = extractDimValue(label2, 'C');
        t2 = extractDimValue(label2, 'T');
        
        fprintf('  Plane 1: Z=%d, C=%d, T=%d\n', z1, c1, t1);
        fprintf('  Plane 2: Z=%d, C=%d, T=%d\n', z2, c2, t2);
        
        % Determine which dimension changed
        if z2 ~= z1
            % Z varies fastest
            if sizeC > 1 && size(rawData, 1) > sizeZ
                % Check if C or T varies second
                label_next = rawData{min(sizeZ+1, size(rawData,1)), 2};
                c_next = extractDimValue(label_next, 'C');
                t_next = extractDimValue(label_next, 'T');
                if c_next ~= c1
                    dimOrder = 'XYZCT'; % Z, then C, then T
                else
                    dimOrder = 'XYZTC'; % Z, then T, then C
                end
            end
        elseif c2 ~= c1
            % C varies fastest
            if sizeZ > 1 && size(rawData, 1) > sizeC
                label_next = rawData{min(sizeC+1, size(rawData,1)), 2};
                z_next = extractDimValue(label_next, 'Z');
                t_next = extractDimValue(label_next, 'T');
                if z_next ~= z1
                    dimOrder = 'XYCZT'; % C, then Z, then T
                else
                    dimOrder = 'XYCTZ'; % C, then T, then Z
                end
            end
        elseif t2 ~= t1
            % T varies fastest (rare)
            dimOrder = 'XYTCZ'; % T, then C, then Z
        end
    end
end

function val = extractDimValue(label, dim)
    % Extract dimension value from Bio-Formats label
    % Example: "Z:5/10; C:2/3; T:1/1" -> extractDimValue(label, 'Z') = 5
    
    pattern = sprintf('%s:(\\d+)/', dim);
    tokens = regexp(label, pattern, 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        val = str2double(tokens{1}{1});
    else
        val = 1; % Default
    end
end

function iPlane = calculatePlaneIndex(z, c, t, sizeZ, sizeC, sizeT, dimOrder)
    % Calculate plane index based on dimension order
    % All indices are 1-based
    
    switch dimOrder
        case 'XYZCT' % Z fastest, then C, then T
            iPlane = (t-1)*sizeZ*sizeC + (c-1)*sizeZ + z;
        case 'XYZTC' % Z fastest, then T, then C
            iPlane = (c-1)*sizeZ*sizeT + (t-1)*sizeZ + z;
        case 'XYCZT' % C fastest, then Z, then T
            iPlane = (t-1)*sizeC*sizeZ + (z-1)*sizeC + c;
        case 'XYCTZ' % C fastest, then T, then Z
            iPlane = (z-1)*sizeC*sizeT + (t-1)*sizeC + c;
        case 'XYTCZ' % T fastest, then C, then Z
            iPlane = (z-1)*sizeT*sizeC + (c-1)*sizeT + t;
        case 'XYTZC' % T fastest, then Z, then C
            iPlane = (c-1)*sizeT*sizeZ + (z-1)*sizeT + t;
        otherwise
            % Default to XYZCT
            iPlane = (t-1)*sizeZ*sizeC + (c-1)*sizeZ + z;
    end
end

%% Utility Functions

function extractedName = extractNameFromString(rawString)
    % Extract name from patterns like "name=..." or "Name=..." in metadata string
    % Falls back to full string if pattern not found
    
    if isempty(rawString)
        extractedName = '';
        return;
    end
    
    % Try to find "name=" or "Name=" pattern (various formats)
    patterns = {
        'name\s*=\s*"([^"]+)"',      % name="value"
        'name\s*=\s*''([^'']+)''',   % name='value'
        'Name\s*=\s*"([^"]+)"',      % Name="value"
        'Name\s*=\s*''([^'']+)''',   % Name='value'
        'NAME\s*=\s*"([^"]+)"',      % NAME="value"
        'NAME\s*=\s*''([^'']+)''',   % NAME='value'
        'name\s*=\s*([^,;\s]+)',     % name=value (no quotes, up to delimiter)
        'Name\s*=\s*([^,;\s]+)',     % Name=value (no quotes, up to delimiter)
        'NAME\s*=\s*([^,;\s]+)'      % NAME=value (no quotes, up to delimiter)
    };
    
    for p = 1:length(patterns)
        tokens = regexp(rawString, patterns{p}, 'tokens', 'once', 'ignorecase');
        if ~isempty(tokens)
            extractedName = strtrim(tokens{1});
            fprintf('    Extracted name: "%s" from "%s"\n', extractedName, rawString);
            return;
        end
    end
    
    % No pattern found, use the full string
    extractedName = rawString;
end

function safeName = createSafeFolderName(originalName, index)
    % Create a safe folder name from original image name
    % Handles long names, special characters, and filesystem limitations
    
    % Remove or replace invalid characters for folder names
    % Invalid: / \ : * ? " < > |
    safeName = originalName;
    invalidChars = '/\:*?"<>|';
    for i = 1:length(invalidChars)
        safeName = strrep(safeName, invalidChars(i), '_');
    end
    
    % Remove leading/trailing spaces and dots
    safeName = strtrim(safeName);
    safeName = regexprep(safeName, '^\.+|\.+$', '');
    
    % Replace multiple spaces/underscores with single underscore
    safeName = regexprep(safeName, '[\s_]+', '_');
    
    % Truncate to maximum length (leave room for path and filenames)
    maxLength = 100; % Conservative limit for folder names
    if length(safeName) > maxLength
        % Keep first part and add index
        safeName = [safeName(1:maxLength-10) sprintf('_%04d', index)];
    end
    
    % If name is empty or invalid, use index-based name
    if isempty(safeName) || strcmp(safeName, '_')
        safeName = sprintf('Image_%04d', index);
    end
    
    % Final validation: ensure it's a valid MATLAB variable name structure
    % (helps avoid weird edge cases)
    if ~isvarname(['x' safeName]) % Prepend 'x' since varnames can't start with numbers
        safeName = matlab.lang.makeValidName(safeName);
    end
    
    % If all else fails, use simple index
    if isempty(safeName)
        safeName = sprintf('Image_%04d', index);
    end
end

function hasIt = checkBioFormats()
    % Check if MATLAB's Bio-Formats Toolbox is available
    
    hasIt = false;
    
    fprintf('SNAP_prepare: Checking for Bio-Formats Toolbox...\n');
    
    % Check if bfopen exists (main function from MATLAB's Bio-Formats Toolbox)
    if exist('bfopen', 'file')
        fprintf('  ✓ Bio-Formats Toolbox found (bfopen detected)\n');
        hasIt = true;
        return;
    end
    
    fprintf('  ✗ Bio-Formats Toolbox not found.\n');
    fprintf('     Please install from MATLAB Add-On Explorer.\n');
    
    % Show installation dialog
    showBioFormatsDialog();
end

function showBioFormatsDialog()
    % Show dialog to help user install Bio-Formats from Add-On Explorer
    
    f = uifigure('Name', 'Bio-Formats Toolbox Required', 'Position', [100 100 650 400]);
    grid = uigridlayout(f, [3, 1]);
    grid.RowHeight = {'1x', 'fit', 'fit'};
    grid.Padding = [20 20 20 20];
    grid.RowSpacing = 15;
    
    % Instructions
    msg = ['Bio-Formats Toolbox is required but not installed.\n\n' ...
           'To install Bio-Formats Toolbox:\n\n' ...
           '1. In MATLAB, go to: Home → Add-Ons → Get Add-Ons\n' ...
           '2. Search for "Bio-Formats"\n' ...
           '3. Install "Bio-Formats Toolbox" by OME Consortium\n' ...
           '4. Restart MATLAB\n' ...
           '5. Run SNAP_prepare again\n\n' ...
           'Alternative: You can also install from:\n' ...
           'File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/\n' ...
           'Search for "Bio-Formats Toolbox"\n\n' ...
           'Note: This is different from the standalone bfmatlab package.\n' ...
           'SNAP_prepare requires the MATLAB Add-On version.'];
    
    textArea = uitextarea(grid, 'Value', sprintf(msg), 'Editable', 'off', ...
        'FontSize', 11, 'FontName', 'Helvetica');
    textArea.Layout.Row = 1;
    
    % Info label
    infoLabel = uilabel(grid, 'Text', 'After installing, restart MATLAB and run SNAP_prepare again.', ...
        'FontSize', 10, 'FontColor', [0 0.5 0.8], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    infoLabel.Layout.Row = 2;
    
    % Buttons
    btnPanel = uipanel(grid, 'BorderType', 'none');
    btnPanel.Layout.Row = 3;
    btnGrid = uigridlayout(btnPanel, [1, 3]);
    btnGrid.ColumnWidth = {'1x', 180, 100};
    btnGrid.Padding = [0 0 0 0];
    btnGrid.ColumnSpacing = 10;
    
    uilabel(btnGrid, 'Text', ''); % Spacer
    
    openAddOnsBtn = uibutton(btnGrid, 'Text', 'Open Add-On Explorer', ...
        'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(src,evt) openAddOnExplorer(f));
    
    okBtn = uibutton(btnGrid, 'Text', 'OK', ...
        'ButtonPushedFcn', @(src,evt) close(f));
end

function openAddOnExplorer(parentFig)
    % Open MATLAB's Add-On Explorer
    try
        matlab.addons.supportpackageinstaller;
    catch
        try
            % Alternative method
            web('https://www.mathworks.com/matlabcentral/fileexchange/', '-browser');
        catch
            uialert(parentFig, 'Could not open Add-On Explorer. Please open it manually from the Home tab.', ...
                'Unable to Open');
        end
    end
    close(parentFig);
end


