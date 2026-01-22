function exportData(fig, exportType)
% Export various types of data based on user selection
% exportType: 'selected', 'browse'
%
% ============================================================================
% EXPORT CONSISTENCY ARCHITECTURE (Critical for Future Maintenance)
% ============================================================================
%
% This function implements a "COMPUTE ONCE, EXPORT MANY" pattern to guarantee
% perfect consistency between nuclei, channel, and nuclei+signal exports.
%
% KEY PRINCIPLE:
%   When multiple nuclei-related exports are selected (e.g., both "Nuclei Data"
%   and "Signal Composition"), nucleus measurements are computed EXACTLY ONCE
%   and the SAME data structure is passed to all export functions.
%
% ARCHITECTURE:
%   
%   PHASE 1: PRE-COMPUTATION (lines ~72-123)
%   ----------------------------------------
%   - Check which exports are selected
%   - If nuclei_data OR signal_composition selected:
%     → Call computeNucleusMeasurements() ONCE
%     → Store result in sharedMeasurements.nucleiData
%     → This creates a unique computation_id for verification
%   
%   PHASE 2: EXPORT LOOP (lines ~125-158)
%   --------------------------------------
%   - Pass sharedMeasurements to nuclei-related export functions
%   - Each function receives PRE-COMPUTED data (never recomputes)
%   - Result: Identical measurements across all exports
%
% WHY THIS MATTERS:
%   Without this system, each export could compute measurements independently,
%   leading to:
%   - Potential discrepancies (e.g., rounding differences, algorithm changes)
%   - User confusion (why do nuclei and nuclei+signal exports differ?)
%   - Data integrity issues in downstream analysis
%
% VERIFICATION:
%   All exports include a computation_id field. If two exports have the same
%   computation_id, their nucleus measurements are GUARANTEED identical.
%
% BACKWARD COMPATIBILITY:
%   The sharedMeasurements parameter is optional in all export functions.
%   If not provided, functions compute measurements themselves (old behavior).
%
% ============================================================================

handles = guidata(fig);

% Handle browse case
if strcmp(exportType, 'browse')
    if isfield(handles, 'exportDirectory') && ~isempty(handles.exportDirectory)
        exportDir = handles.exportDirectory;
    else
        exportDir = pwd;
    end
    selectedDir = uigetdir(exportDir, 'Select Export Directory');
    if selectedDir ~= 0
        handles.exportDirectory = selectedDir;
        guidata(fig, handles);
        fprintf('Export directory set to: %s\n', selectedDir);
    end
    return;
end

% Get export directory
if isfield(handles, 'exportDirectory') && ~isempty(handles.exportDirectory)
    exportDir = handles.exportDirectory;
else
    exportDir = pwd;
end

% Ensure export directory exists
if ~exist(exportDir, 'dir')
    mkdir(exportDir);
end

% Generate timestamp for unique filenames
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

try
    % Export selected items based on adaptive checklist
    exportCount = 0;
    
    % Ensure export checklist is initialized
    if ~isfield(handles, 'numExportItems') || handles.numExportItems == 0
        % Try to initialize the export checklist
        snap_helpers.updateExportChecklist(fig);
        handles = guidata(fig);
        
        % Check again after initialization
        if ~isfield(handles, 'numExportItems') || handles.numExportItems == 0
            uialert(handles.fig, 'No export items available. Please load some data and run "Update Previews" first.', 'Export Error');
            return;
        end
    end
    
    % Reload handles to get latest checkbox states (user may have clicked checkboxes)
    handles = guidata(fig);
    
    % Check if any items are selected
    hasSelectedItems = false;
    for i = 1:handles.numExportItems
        if handles.exportItemChecks(i).Value
            hasSelectedItems = true;
            break;
        end
    end
    
    if ~hasSelectedItems
        uialert(handles.fig, 'No items selected for export. Please check at least one export option.', 'Export Error');
        return;
    end
    
    %% ========================================================================
    %  PHASE 1: PRE-COMPUTATION OF SHARED MEASUREMENTS
    %  ========================================================================
    %  
    %  PURPOSE: Compute nucleus measurements ONCE and share across all exports
    %           to guarantee perfect consistency.
    %
    %  LOGIC:
    %    1. Scan selected exports to see if any are nuclei-related
    %    2. If yes → compute measurements once, store in sharedMeasurements
    %    3. If no → skip (no shared computation needed)
    %
    %  CRITICAL: This prevents the following BAD scenario:
    %    - exportNucleiData() computes measurements → gets result A
    %    - exportClusteredData() computes measurements → gets result B
    %    - User sees different nucleus sizes in different exports (INCONSISTENT!)
    %
    %  GOOD scenario (with this code):
    %    - Measurements computed ONCE here → stored
    %    - exportNucleiData() receives stored measurements
    %    - exportClusteredData() receives SAME stored measurements
    %    - User sees IDENTICAL nucleus data (CONSISTENT!)
    %
    %  ========================================================================
    
    sharedMeasurements = struct();
    needsNucleiMeasurements = false;
    
    % Check if any nuclei-related exports are selected
    for i = 1:handles.numExportItems
        if handles.exportItemChecks(i).Value
            itemType = handles.exportItemTypes{i};
            if strcmp(itemType, 'nuclei_data') || strcmp(itemType, 'signal_composition')
                needsNucleiMeasurements = true;
                break;
            end
        end
    end
    
    % Compute nuclei measurements ONCE if needed
    if needsNucleiMeasurements
        % Get nuclei data from cache
        nucleus_labels = [];
        nuclei_mask = [];
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei')
            if isfield(handles.previewCache.nuclei, 'labels')
                nucleus_labels = handles.previewCache.nuclei.labels;
            end
            if isfield(handles.previewCache.nuclei, 'mask')
                nuclei_mask = handles.previewCache.nuclei.mask;
            end
        end
        
        if ~isempty(nucleus_labels) && ~isempty(nuclei_mask) && nucleus_labels.num_nuclei > 0
            % Compute measurements ONCE using the SINGLE SOURCE OF TRUTH function
            fprintf('Computing shared nucleus measurements (ensures consistency)...\n');
            
            % Determine if 3D from actual image data
            is_3d = false;
            if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
                is_3d = (ndims(handles.rawNuclei) == 3 && size(handles.rawNuclei, 3) > 1);
            end
            
            % Build options for measurement computation
            measureOptions = struct();
            measureOptions.is_3d = is_3d;
            measureOptions.xy_spacing = handles.nucXYSpacingInput.Value;
            measureOptions.z_spacing = handles.nucZSpacingInput.Value;
            measureOptions.image_name = '';
            
            % CRITICAL: Call computeNucleusMeasurements() ONCE here
            % This result will be passed to ALL nuclei-related exports below
            % The function generates a unique computation_id for verification
            sharedMeasurements.nucleiData = snap_helpers.computeNucleusMeasurements(nucleus_labels, nuclei_mask, measureOptions);
            fprintf('  ✓ Shared measurements computed (%d nuclei)\n', sharedMeasurements.nucleiData.num_nuclei);
        end
    end
    
    %% ========================================================================
    %  PHASE 2: EXPORT LOOP - All nuclei exports use shared measurements
    %  ========================================================================
    %
    %  CRITICAL: sharedMeasurements is passed to nuclei-related exports
    %            This ensures they all use THE SAME nucleus measurements
    %
    %  ========================================================================
    
    for i = 1:handles.numExportItems
        if handles.exportItemChecks(i).Value
            itemType = handles.exportItemTypes{i};
            
            % Parse the item type and export accordingly
            if strcmp(itemType, 'nuclei_data')
                % CONSISTENCY: Pass sharedMeasurements to ensure identical data
                exportNucleiData(handles, exportDir, timestamp, sharedMeasurements);
                exportCount = exportCount + 1;
                
            elseif strcmp(itemType, 'channel_data')
                % INDEPENDENT: Channel export doesn't need nucleus measurements
                exportChannelData(handles, exportDir, timestamp);
                exportCount = exportCount + 1;
                
            elseif strcmp(itemType, 'signal_composition')
                % CONSISTENCY: Pass sharedMeasurements to ensure identical nucleus data as nuclei export
                exportClusteredData(handles, exportDir, timestamp, sharedMeasurements);
                exportCount = exportCount + 1;
                
            elseif strcmp(itemType, 'parameters')
                exportParameters(handles, exportDir, timestamp);
                exportCount = exportCount + 1;
                
            elseif strcmp(itemType, 'image_nuclei')
                exportNucleiImage(handles, exportDir, timestamp);
                exportCount = exportCount + 1;
                
            elseif startsWith(itemType, 'image_channel_')
                % Extract channel number
                ch = str2double(itemType(15:end));
                exportChannelImage(handles, exportDir, timestamp, ch);
                exportCount = exportCount + 1;
            end
        end
    end
    
    if exportCount > 0
        msgbox(sprintf('Successfully exported %d item(s) to:\n%s', exportCount, exportDir), 'Export Complete', 'help');
    else
        msgbox('No items selected for export. Please check at least one item.', 'Export', 'warn');
    end
    
catch ME
    errordlg(['Export failed: ' ME.message], 'Export Error');
end
end

function exportNucleiData(handles, exportDir, timestamp, sharedMeasurements)
% Export nuclei segmentation data using standardized format
%
% INPUTS:
%   handles            - SNAP handles structure
%   exportDir          - Export directory
%   timestamp          - Export timestamp
%   sharedMeasurements - (OPTIONAL) Pre-computed measurements struct with:
%                        .nucleiData - from computeNucleusMeasurements()
%                        Ensures consistency with other exports

fprintf('Exporting nuclei segmentation data...\n');

% Try to get nuclei data from cache first, then compute if needed
nuclei_mask = [];
nucleus_labels = [];

% Check if preview cache exists and has nuclei data
if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei')
    if isfield(handles.previewCache.nuclei, 'mask')
        nuclei_mask = handles.previewCache.nuclei.mask;
    end
    if isfield(handles.previewCache.nuclei, 'labels')
        nucleus_labels = handles.previewCache.nuclei.labels;
    end
end

% If no cached data, try to compute it
if isempty(nuclei_mask) && isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
    try
        fprintf('Computing nuclei segmentation for export...\n');
        processed_nuclei = snap_helpers.preprocessNucleiWithBgCorr(handles);
        [nuclei_mask, nucleus_labels] = snap_helpers.segmentNuclei(processed_nuclei, handles);
    catch ME
        warning('Failed to compute nuclei segmentation: %s', ME.message);
        return;
    end
end

% Check if we have valid data
if isempty(nuclei_mask) || isempty(nucleus_labels)
    warning('No nuclei segmentation data available for export. Please run segmentation first.');
    return;
end

% Use standardized export function (no timestamp in filename)
outputPath = fullfile(exportDir, 'export_nuclei');

% Prepare options
exportOptions = struct();
exportOptions.is_3d = strcmp(handles.nucSegModeDrop.Value, 'On 3D Volume');
exportOptions.xy_spacing = handles.nucXYSpacingInput.Value;
exportOptions.z_spacing = handles.nucZSpacingInput.Value;
exportOptions.image_name = '';  % No image name for single SNAP export

% Call standardized export function with pre-computed measurements if available
if nargin >= 4 && ~isempty(sharedMeasurements) && isfield(sharedMeasurements, 'nucleiData')
    snap_helpers.exportNucleiDataStandardized(nucleus_labels, nuclei_mask, outputPath, exportOptions, sharedMeasurements.nucleiData);
else
snap_helpers.exportNucleiDataStandardized(nucleus_labels, nuclei_mask, outputPath, exportOptions);
end

% OLD CODE BELOW - Now handled by exportNucleiDataStandardized
return;

%%%% LEGACY CODE (for reference, will be removed) %%%%
filename = fullfile(exportDir, sprintf('nuclei_data_%s_OLD.mat', timestamp));

% Prepare nuclei data structure
nucleiData = struct();
nucleiData.exportTimestamp = timestamp;
nucleiData.exportDate = datestr(now);

% Binary mask (what's shown as overlay)
nucleiData.binary_mask = nuclei_mask;

% Labeled mask with nucleus IDs
if ~isempty(nucleus_labels) && isfield(nucleus_labels, 'labeled_mask')
    nucleiData.labeled_mask = nucleus_labels.labeled_mask;
    nucleiData.nucleus_labels = nucleus_labels;
    
    % Add nucleus information
    nucleiData.num_nuclei = nucleus_labels.num_nuclei;
    nucleiData.centroids_3d = nucleus_labels.centroids_3d;
    
    % Create ROI information for mapping signals to nuclei
    nucleiData.roi_info = struct();
    nucleiData.roi_info.description = 'Each nucleus ROI is labeled with its ID number for signal mapping';
    nucleiData.roi_info.nucleus_ids = 1:nucleus_labels.num_nuclei;
    
    if isfield(nucleus_labels, 'connected_components')
        nucleiData.roi_info.pixel_lists = nucleus_labels.connected_components.PixelIdxList;
    end
    
    fprintf('Exported %d labeled nuclei ROIs\n', nucleus_labels.num_nuclei);
else
    nucleiData.labeled_mask = [];
    nucleiData.nucleus_labels = [];
    fprintf('Warning: No nucleus labeling information available\n');
end

% Add image dimensions and spacing information
if isfield(handles, 'rawNuclei')
    nucleiData.image_size = size(handles.rawNuclei);
end

if isfield(handles, 'nucXYSpacingInput') && isvalid(handles.nucXYSpacingInput)
    nucleiData.xy_spacing = handles.nucXYSpacingInput.Value;
end

if isfield(handles, 'nucZSpacingInput') && isvalid(handles.nucZSpacingInput)
    nucleiData.z_spacing = handles.nucZSpacingInput.Value;
end

% Add segmentation parameters used
nucleiData.segmentation_parameters = struct();
if isfield(handles, 'nucSegModeDrop') && isvalid(handles.nucSegModeDrop)
    nucleiData.segmentation_parameters.mode = handles.nucSegModeDrop.Value;
end
if isfield(handles, 'nucSegMainMethodDrop') && isvalid(handles.nucSegMainMethodDrop)
    nucleiData.segmentation_parameters.method = handles.nucSegMainMethodDrop.Value;
end

% Add metadata for cross-referencing
nucleiData.metadata = struct();
nucleiData.metadata.description = 'Nucleus segmentation data with binary and labeled masks';
nucleiData.metadata.companion_exports = 'For signal analysis, see Channel Data export; for nucleus-signal associations, see Clustered Signal Data export';
nucleiData.metadata.note = 'This export focuses on nucleus morphology. Combine with Channel/Clustered exports for complete analysis.';

save(filename, 'nucleiData');
fprintf('Nuclei data exported to: %s\n', filename);

% Also create a summary text file
txtFilename = fullfile(exportDir, sprintf('nuclei_summary_%s.txt', timestamp));
exportNucleiSummary(nucleiData, txtFilename);
end

function exportNucleiSummary(nucleiData, filename)
% Export a human-readable summary of nuclei data
fid = fopen(filename, 'w');
if fid == -1
    warning('Could not create nuclei summary file: %s', filename);
    return;
end

fprintf(fid, 'Nuclei Segmentation Export Summary\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

if ~isempty(nucleiData.binary_mask)
    fprintf(fid, 'Image Dimensions: %s\n', mat2str(size(nucleiData.binary_mask)));
end

if isfield(nucleiData, 'num_nuclei') && ~isempty(nucleiData.num_nuclei)
    fprintf(fid, 'Number of Nuclei: %d\n', nucleiData.num_nuclei);
else
    fprintf(fid, 'Number of Nuclei: Not available (binary mask only)\n');
end

if isfield(nucleiData, 'xy_spacing')
    fprintf(fid, 'XY Spacing: %.3f µm/pixel\n', nucleiData.xy_spacing);
end

if isfield(nucleiData, 'z_spacing')
    fprintf(fid, 'Z Spacing: %.3f µm/slice\n', nucleiData.z_spacing);
end

if isfield(nucleiData, 'segmentation_parameters')
    fprintf(fid, '\nSegmentation Parameters:\n');
    params = nucleiData.segmentation_parameters;
    fields = fieldnames(params);
    for i = 1:length(fields)
        fprintf(fid, '  %s: %s\n', fields{i}, char(params.(fields{i})));
    end
end

if isfield(nucleiData, 'roi_info') && ~isempty(nucleiData.roi_info)
    fprintf(fid, '\nROI Information:\n');
    fprintf(fid, '  %s\n', nucleiData.roi_info.description);
    fprintf(fid, '  Nucleus IDs: %s\n', mat2str(nucleiData.roi_info.nucleus_ids));
end

fprintf(fid, '\nData Files:\n');
fprintf(fid, '  binary_mask: Binary segmentation mask (logical array)\n');
fprintf(fid, '  labeled_mask: Labeled mask with nucleus IDs (integer array)\n');
fprintf(fid, '  nucleus_labels: Complete labeling structure with centroids and metadata\n');
fprintf(fid, '  roi_info: ROI information for signal-to-nucleus mapping\n');

fclose(fid);
fprintf('Nuclei summary exported to: %s\n', filename);
end

function exportChannelData(handles, exportDir, timestamp)
% Export channel signal data using standardized export function
fprintf('Exporting channel signal data...\n');

% Get number of channels
numChannels = str2double(handles.numChanDrop.Value);

% Determine if data is 3D from actual image data
is_3d = false;
if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
    is_3d = (ndims(handles.rawNuclei) == 3 && size(handles.rawNuclei, 3) > 1);
elseif isfield(handles, 'rawChannel') && ~isempty(handles.rawChannel)
    for ch = 1:length(handles.rawChannel)
        if ~isempty(handles.rawChannel{ch})
            is_3d = (ndims(handles.rawChannel{ch}) == 3 && size(handles.rawChannel{ch}, 3) > 1);
            break;
        end
    end
end

% Determine image type string
if is_3d
    imageType = '3D Image Stack';
else
    imageType = '2D Image';
end

% Export each channel separately using standardized function
for ch = 1:numChannels
    fprintf('  Exporting channel %d...\n', ch);
    
    % Get maxima coordinates and fit results
    maxima_coords = [];
    fit_results = [];
    
    % Get from handles (already processed during Update Previews)
    if isfield(handles, 'maximaCoords') && ch <= length(handles.maximaCoords)
        maxima_coords = handles.maximaCoords{ch};
    end
    if isfield(handles, 'gaussFitResults') && ch <= length(handles.gaussFitResults)
        fit_results = handles.gaussFitResults{ch};
    end
    
    % Skip if no data for this channel
    if isempty(maxima_coords)
        fprintf('  No data for channel %d, skipping...\n', ch);
        continue;
    end
    
    % Build channel results structure
    channelResults = struct();
    channelResults.maxima_coords = maxima_coords;
    if ~isempty(fit_results)
        channelResults.fit_results = fit_results;
    end
    
    % Build fit parameters structure
    fitParams = struct();
    fitParams.imageType = imageType;
    
    if isfield(handles, 'gaussFitModeDrop') && length(handles.gaussFitModeDrop) >= ch
        fitParams.gaussFitMode = handles.gaussFitModeDrop(ch).Value;
    else
        fitParams.gaussFitMode = 'Not specified';
    end
    
    if isfield(handles, 'gaussFitMethodDrop') && length(handles.gaussFitMethodDrop) >= ch
        fitParams.gaussFitMethod = handles.gaussFitMethodDrop(ch).Value;
    else
        fitParams.gaussFitMethod = 'Not specified';
    end
    
    if isfield(handles, 'gaussFitBgCorrMethodDrop') && length(handles.gaussFitBgCorrMethodDrop) >= ch
        fitParams.gaussFitBgCorrMethod = handles.gaussFitBgCorrMethodDrop(ch).Value;
    else
        fitParams.gaussFitBgCorrMethod = 'Mean Surrounding Subtraction';
    end
    
    if isfield(handles, 'gaussFitBgCorrWidthInputs') && length(handles.gaussFitBgCorrWidthInputs) >= ch
        fitParams.gaussFitBgCorrWidth = str2double(handles.gaussFitBgCorrWidthInputs(ch).Value);
    else
        fitParams.gaussFitBgCorrWidth = 2;
    end
    
    if isfield(handles, 'gaussFitPolyDegreeInputs') && length(handles.gaussFitPolyDegreeInputs) >= ch
        fitParams.gaussFitPolyDegree = str2double(handles.gaussFitPolyDegreeInputs(ch).Value);
    else
        fitParams.gaussFitPolyDegree = 2;
    end
    
    % Additional fitting parameters
    if isfield(handles, 'gaussFitVoxelWindowSlider') && length(handles.gaussFitVoxelWindowSlider) >= ch
        fitParams.gaussFitVoxelWindowSize = handles.gaussFitVoxelWindowSlider(ch).Value;
    else
        fitParams.gaussFitVoxelWindowSize = 7;
    end
    
    if isfield(handles, 'gaussFitMaxIterationsEdit') && length(handles.gaussFitMaxIterationsEdit) >= ch
        fitParams.gaussFitMaxIterations = handles.gaussFitMaxIterationsEdit(ch).Value;
    else
        fitParams.gaussFitMaxIterations = 200;
    end
    
    if isfield(handles, 'gaussFitToleranceEdit') && length(handles.gaussFitToleranceEdit) >= ch
        fitParams.gaussFitTolerance = handles.gaussFitToleranceEdit(ch).Value;
    else
        fitParams.gaussFitTolerance = 1e-6;
    end
    
    if isfield(handles, 'gaussFitRadialRadiusEdit') && length(handles.gaussFitRadialRadiusEdit) >= ch
        fitParams.gaussFitRadialRadius = handles.gaussFitRadialRadiusEdit(ch).Value;
    else
        fitParams.gaussFitRadialRadius = 3;
    end
    
    % Build options structure
    options = struct();
    options.is_3d = is_3d;
    if isfield(handles, 'xySpacingInputs') && length(handles.xySpacingInputs) >= ch
        xy_spacing = str2double(handles.xySpacingInputs(ch).Value);
        z_spacing = 1.0;
        if isfield(handles, 'zSpacingInputs') && length(handles.zSpacingInputs) >= ch
            z_spacing = str2double(handles.zSpacingInputs(ch).Value);
        end
        options.spacing = [xy_spacing, xy_spacing, z_spacing];  % [dy, dx, dz]
    end
    
    % Build output path (no timestamp in filename)
    outputPath = fullfile(exportDir, sprintf('export_ch%d', ch));
    
    % Call standardized export function
    % NOTE: Pass empty string for imageName (SNAP GUI exports don't need image_name column)
    % Only SNAP_batch provides image names (subfolder names)
    snap_helpers.exportChannelDataStandardized(outputPath, '', channelResults, fitParams, options);
end

fprintf('Channel signal data export complete.\n');
end

function exportClusteredData(handles, exportDir, timestamp, sharedMeasurements)
% Export nuclei signal composition data using standardized export function
%
% INPUTS:
%   handles            - SNAP handles structure
%   exportDir          - Export directory
%   timestamp          - Export timestamp
%   sharedMeasurements - (OPTIONAL) Pre-computed measurements struct with:
%                        .nucleiData - from computeNucleusMeasurements()
%                        Ensures consistency with other exports

fprintf('Exporting nuclei signal composition data...\n');

% Get nuclei labels and mask from cache
nucleus_labels = [];
nuclei_mask = [];
if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei')
    if isfield(handles.previewCache.nuclei, 'labels')
        nucleus_labels = handles.previewCache.nuclei.labels;
    end
    if isfield(handles.previewCache.nuclei, 'mask')
        nuclei_mask = handles.previewCache.nuclei.mask;
    end
end

if isempty(nucleus_labels) || isempty(nuclei_mask)
    warning('No nuclei segmentation data available. Run "Update Previews" with nuclei segmentation enabled first.');
    return;
end

% Get or compute signal composition (with full signal data embedded)
signalComposition = [];
if isfield(handles, 'signalComposition') && ~isempty(handles.signalComposition)
    signalComposition = handles.signalComposition;
else
    % Compute it fresh if not cached
    fprintf('  Computing signal composition...\n');
    signalComposition = snap_helpers.analyzeNucleiSignalComposition(handles);
end

if isempty(signalComposition)
    warning('No signal composition data available.');
    return;
end

% Determine if 3D from actual image data
is_3d = false;
if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
    is_3d = (ndims(handles.rawNuclei) == 3 && size(handles.rawNuclei, 3) > 1);
end

% Build options structure
options = struct();
options.is_3d = is_3d;
options.xy_spacing = handles.nucXYSpacingInput.Value;
options.z_spacing = handles.nucZSpacingInput.Value;

% Build output path (no timestamp in filename)
outputPath = fullfile(exportDir, 'export_nucleisignals');

% Call standardized export function with pre-computed measurements if available
% NOTE: Pass empty string for imageName (SNAP GUI exports don't need image_name column)
% Only SNAP_batch provides image names (subfolder names)
if nargin >= 4 && ~isempty(sharedMeasurements) && isfield(sharedMeasurements, 'nucleiData')
    snap_helpers.exportNucleiSignalDataStandardized(outputPath, '', nucleus_labels, nuclei_mask, signalComposition, options, sharedMeasurements.nucleiData);
else
    snap_helpers.exportNucleiSignalDataStandardized(outputPath, '', nucleus_labels, nuclei_mask, signalComposition, options);
end

fprintf('Nuclei signal composition export complete.\n');
end

function signal_indices = findSignalsInNucleus(maxima_coords, nucleus_pixel_indices, nuclei_mask, nuc_idx)
% LEGACY FUNCTION - Kept for backward compatibility
% NEW CODE SHOULD USE getSignalsInNucleus from analyzeNucleiSignalComposition.m
%
% Find which signals (by index) fall within this nucleus
% COORDINATE CONVENTION: Uses ARRAY CONVENTION
%   maxima_coords are [row, col, slice]
%   Array access: nucleiMask(row, col, slice) - DIRECT indexing

    signal_indices = [];
    
    if isempty(maxima_coords)
        return;
    end
    
    % Create a binary mask for this specific nucleus
    nucleus_mask = false(size(nuclei_mask));
    nucleus_mask(nucleus_pixel_indices) = true;
    
    % Check each signal's coordinates (ARRAY CONVENTION)
    for i = 1:size(maxima_coords, 1)
        row = round(maxima_coords(i, 1));
        col = round(maxima_coords(i, 2));
        slice = round(maxima_coords(i, 3));
        
        % Handle 2D/3D compatibility
        if ndims(nucleus_mask) == 2
            % 2D nucleus mask - only check row, col
            if row >= 1 && row <= size(nucleus_mask, 1) && col >= 1 && col <= size(nucleus_mask, 2)
                % DIRECT array access with ARRAY CONVENTION
                if nucleus_mask(row, col)
                    signal_indices(end+1) = i;
                end
            end
        else
            % 3D nucleus mask - check row, col, slice
            if row >= 1 && row <= size(nucleus_mask, 1) && ...
               col >= 1 && col <= size(nucleus_mask, 2) && ...
               slice >= 1 && slice <= size(nucleus_mask, 3)
                % DIRECT array access with ARRAY CONVENTION
                if nucleus_mask(row, col, slice)
                    signal_indices(end+1) = i;
                end
            end
        end
    end
end

function exportParameters(handles, exportDir, timestamp)
% Export current processing parameters - COLLECTS FROM UI, NOT CACHED lastUsed
fprintf('Exporting current parameters from UI...\n');

% First, collect CURRENT parameter values from UI (like closeAndSave does)
lastUsed = collectCurrentParameters(handles);

% Save MAT file only
filename = fullfile(exportDir, 'export_parameters.mat');
save(filename, 'lastUsed');
fprintf('Parameters exported to: %s\n', filename);
end

function exportParametersToText(handles, filename)
% Export parameters in human-readable text format
fid = fopen(filename, 'w');
if fid == -1
    warning('Could not create text parameter file: %s', filename);
    return;
end

fprintf(fid, 'Processing Parameters Export\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

if isfield(handles, 'lastUsed')
    params = handles.lastUsed;
    fields = fieldnames(params);
    
    for i = 1:length(fields)
        fieldName = fields{i};
        value = params.(fieldName);
        
        if isnumeric(value)
            if length(value) == 1
                fprintf(fid, '%s: %g\n', fieldName, value);
            else
                fprintf(fid, '%s: [%s]\n', fieldName, num2str(value));
            end
        elseif ischar(value) || isstring(value)
            fprintf(fid, '%s: %s\n', fieldName, char(value));
        elseif islogical(value)
            if length(value) == 1
                fprintf(fid, '%s: %s\n', fieldName, logical2str(value));
            else
                fprintf(fid, '%s: [%s]\n', fieldName, strjoin(arrayfun(@logical2str, value, 'UniformOutput', false), ', '));
            end
        elseif iscell(value)
            try
                % Handle cell arrays more carefully
                cellStrs = cell(size(value));
                for i = 1:numel(value)
                    if ischar(value{i}) || isstring(value{i})
                        cellStrs{i} = char(value{i});
                    elseif isnumeric(value{i})
                        cellStrs{i} = num2str(value{i});
                    elseif islogical(value{i})
                        cellStrs{i} = logical2str(value{i});
                    else
                        cellStrs{i} = '[complex]';
                    end
                end
                fprintf(fid, '%s: {%s}\n', fieldName, strjoin(cellStrs, ', '));
            catch
                fprintf(fid, '%s: [cell array - complex data]\n', fieldName);
            end
        else
            fprintf(fid, '%s: [complex data type]\n', fieldName);
        end
    end
end

fclose(fid);
fprintf('Human-readable parameters exported to: %s\n', filename);
end

function str = logical2str(logicalValue)
% Convert logical value to readable string
if logicalValue
    str = 'true';
else
    str = 'false';
end
end

function checksum = generateParameterChecksum(parameters)
% Generate a simple checksum for parameter validation
try
    paramStr = mat2str(struct2cell(parameters));
    checksum = sum(double(paramStr));
catch
    checksum = 0;
end
end

function createBatchInstructions(paramData, filename)
% Create batch processing instructions file
fid = fopen(filename, 'w');
if fid == -1
    warning('Could not create batch instructions file: %s', filename);
    return;
end

fprintf(fid, 'Batch Processing Instructions\n');
fprintf(fid, '============================\n\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Version: %s\n\n', paramData.version);

fprintf(fid, 'Workflow Configuration:\n');
fprintf(fid, '  Number of Channels: %d\n', paramData.workflowConfig.numChannels);
fprintf(fid, '  3D Data Processing: %s\n', logical2str(paramData.workflowConfig.has3DData));

if isfield(paramData.workflowConfig, 'nucleiProcessing')
    fprintf(fid, '\nNuclei Processing:\n');
    nuc = paramData.workflowConfig.nucleiProcessing;
    if isfield(nuc, 'segmentationEnabled')
        fprintf(fid, '  Segmentation: %s\n', logical2str(nuc.segmentationEnabled));
    end
    if isfield(nuc, 'filteringEnabled')
        fprintf(fid, '  Filtering: %s\n', logical2str(nuc.filteringEnabled));
    end
end

fprintf(fid, '\nChannel Processing:\n');
for k = 1:paramData.workflowConfig.numChannels
    channelField = ['channel_' num2str(k)];
    if isfield(paramData.workflowConfig.channelProcessing, channelField)
        ch = paramData.workflowConfig.channelProcessing.(channelField);
        fprintf(fid, '  Channel %d:\n', k);
        if isfield(ch, 'preprocessingEnabled')
            fprintf(fid, '    Preprocessing: %s\n', logical2str(ch.preprocessingEnabled));
        end
        if isfield(ch, 'backgroundCorrectionEnabled')
            fprintf(fid, '    Background Correction: %s\n', logical2str(ch.backgroundCorrectionEnabled));
        end
        if isfield(ch, 'maximaDetectionEnabled')
            fprintf(fid, '    Maxima Detection: %s\n', logical2str(ch.maximaDetectionEnabled));
        end
        if isfield(ch, 'gaussianFittingEnabled')
            fprintf(fid, '    Gaussian Fitting: %s\n', logical2str(ch.gaussianFittingEnabled));
        end
    end
end

fprintf(fid, '\nExported Files:\n');
fprintf(fid, '  processing_parameters_*.mat - Complete parameter export\n');
fprintf(fid, '  processing_parameters_*.txt - Human-readable parameters\n');

fprintf(fid, '\nUsage:\n');
fprintf(fid, '  1. Load processing_parameters_*.mat in SNAP or SNAP_batch\n');
fprintf(fid, '  2. Apply parameters to new image sets\n');
fprintf(fid, '  3. Ensure image file paths are updated for new datasets\n');
fprintf(fid, '  4. Verify spacing parameters match new image properties\n');

fclose(fid);
fprintf('Batch processing instructions exported to: %s\n', filename);
end

function exportProcessedImages(handles, exportDir, timestamp)
% Export processed images that match what's displayed in previews
% Uses cached processed images from the preview system - much simpler and faster!
fprintf('Exporting processed images (from preview cache)...\n');

% Create subdirectory for processed images
processedImagesDir = fullfile(exportDir, sprintf('processed_images_%s', timestamp));
if ~exist(processedImagesDir, 'dir')
    mkdir(processedImagesDir);
end

% Get export format preference (default to TIFF)
exportFormat = 'tiff';
if isfield(handles, 'exportImageFormatDrop') && isvalid(handles.exportImageFormatDrop)
    exportFormat = lower(handles.exportImageFormatDrop.Value);
end

% Prepare metadata structure
imageMetadata = struct();
imageMetadata.exportTimestamp = timestamp;
imageMetadata.exportDate = datestr(now);
imageMetadata.description = 'Processed images exported from preview cache - exactly what you see in the preview panels';
imageMetadata.export_format = exportFormat;

images_exported = 0;

% Check if preview cache exists
if ~isfield(handles, 'previewCache')
    fprintf('Warning: No preview cache found. Please run "Update Previews" first to process images.\n');
    fprintf('No images exported.\n');
    return;
end

% Export processed nuclei image if available in cache
if isfield(handles.previewCache, 'nuclei') && isfield(handles.previewCache.nuclei, 'processed') && ~isempty(handles.previewCache.nuclei.processed)
    try
        fprintf('Exporting cached nuclei image...\n');
        
        processed_nuclei = handles.previewCache.nuclei.processed;
        nucleiFilename = fullfile(processedImagesDir, sprintf('nuclei_image.%s', exportFormat));
        saveProcessedImage(processed_nuclei, nucleiFilename, exportFormat);
        
        % Store metadata
        imageMetadata.nuclei = struct();
        imageMetadata.nuclei.filename = sprintf('nuclei_image.%s', exportFormat);
        imageMetadata.nuclei.size = size(processed_nuclei);
        imageMetadata.nuclei.processing_applied = getNucleiProcessingDescription(handles);
        
        % Add spacing information
        if isfield(handles, 'nucXYSpacingInput') && isvalid(handles.nucXYSpacingInput)
            imageMetadata.nuclei.xy_spacing = handles.nucXYSpacingInput.Value;
        end
        if isfield(handles, 'nucZSpacingInput') && isvalid(handles.nucZSpacingInput)
            imageMetadata.nuclei.z_spacing = handles.nucZSpacingInput.Value;
        end
        
        fprintf('✓ Nuclei image exported: %s\n', nucleiFilename);
        images_exported = images_exported + 1;
        
    catch ME
        fprintf('Warning: Could not export nuclei image: %s\n', ME.message);
        imageMetadata.nuclei = struct('error', ME.message);
    end
else
    fprintf('No processed nuclei image found in cache. Run "Update Previews" first.\n');
end

% Export processed channel images from cache
numChannels = str2double(handles.numChanDrop.Value);
imageMetadata.channels = struct();

if isfield(handles.previewCache, 'channels')
    for ch = 1:numChannels
        if length(handles.previewCache.channels) >= ch && ...
           isfield(handles.previewCache.channels{ch}, 'processed') && ...
           ~isempty(handles.previewCache.channels{ch}.processed)
            
            try
                fprintf('Exporting cached channel %d image...\n', ch);
                
                processed_channel = handles.previewCache.channels{ch}.processed;
                channelFilename = fullfile(processedImagesDir, sprintf('ch%d_image.%s', ch, exportFormat));
                saveProcessedImage(processed_channel, channelFilename, exportFormat);
                
                % Store metadata for this channel
                channelMetadata = struct();
                channelMetadata.filename = sprintf('ch%d_image.%s', ch, exportFormat);
                channelMetadata.size = size(processed_channel);
                channelMetadata.processing_applied = getChannelProcessingDescription(handles, ch);
                
                % Add spacing information
                if isfield(handles, 'xySpacingInputs') && length(handles.xySpacingInputs) >= ch
                    channelMetadata.xy_spacing = handles.xySpacingInputs(ch).Value;
                end
                if isfield(handles, 'zSpacingInputs') && length(handles.zSpacingInputs) >= ch
                    channelMetadata.z_spacing = handles.zSpacingInputs(ch).Value;
                end
                
                imageMetadata.channels.(sprintf('channel_%d', ch)) = channelMetadata;
                
                fprintf('✓ Channel %d image exported: %s\n', ch, channelFilename);
                images_exported = images_exported + 1;
                
            catch ME
                fprintf('Warning: Could not export channel %d image: %s\n', ch, ME.message);
                imageMetadata.channels.(sprintf('channel_%d', ch)) = struct('error', ME.message);
            end
        else
            fprintf('No processed image found for channel %d in cache. Run "Update Previews" first.\n', ch);
            imageMetadata.channels.(sprintf('channel_%d', ch)) = struct('status', 'not_processed_in_cache');
        end
    end
else
    fprintf('No processed channel images found in cache. Run "Update Previews" first.\n');
end

% Export DIC image if available (raw image)
if isfield(handles, 'rawDIC') && ~isempty(handles.rawDIC)
    try
        fprintf('Exporting DIC image...\n');
        
        dicFilename = fullfile(processedImagesDir, sprintf('dic_image_%s.%s', timestamp, exportFormat));
        saveProcessedImage(handles.rawDIC, dicFilename, exportFormat);
        
        % Store metadata
        imageMetadata.dic = struct();
        imageMetadata.dic.filename = sprintf('dic_image_%s.%s', timestamp, exportFormat);
        imageMetadata.dic.size = size(handles.rawDIC);
        imageMetadata.dic.processing_applied = 'None (raw DIC image)';
        
        fprintf('✓ DIC image exported: %s\n', dicFilename);
        images_exported = images_exported + 1;
        
    catch ME
        fprintf('Warning: Could not export DIC image: %s\n', ME.message);
        imageMetadata.dic = struct('error', ME.message);
    end
end

% Save metadata file if any images were exported
if images_exported > 0
    metadataFilename = fullfile(processedImagesDir, sprintf('image_metadata_%s.mat', timestamp));
    save(metadataFilename, 'imageMetadata');
    
    % Create human-readable summary
    summaryFilename = fullfile(processedImagesDir, sprintf('processing_summary_%s.txt', timestamp));
    createProcessedImageSummary(imageMetadata, summaryFilename);
    
    fprintf('\n✅ Export complete! %d images exported to: %s\n', images_exported, processedImagesDir);
else
    fprintf('\n❌ No images exported. Please run "Update Previews" first to process your images.\n');
    % Remove empty directory
    if exist(processedImagesDir, 'dir')
        rmdir(processedImagesDir);
    end
end
end

function saveProcessedImage(imageData, filename, format)
% Save processed image data in the specified format
% Handles both 2D and 3D data appropriately

% Ensure image data is in appropriate format for saving
if isa(imageData, 'double')
    % Convert double to appropriate integer type for saving
    % Preserve the full dynamic range by scaling to uint16
    minVal = min(imageData(:));
    maxVal = max(imageData(:));
    
    if maxVal > minVal
        % Scale to uint16 range (0-65535) to preserve precision
        imageData_scaled = uint16(65535 * (imageData - minVal) / (maxVal - minVal));
    else
        % Constant image
        imageData_scaled = uint16(zeros(size(imageData)));
    end
else
    imageData_scaled = imageData;
end

switch lower(format)
    case {'tiff', 'tif'}
        % Use imwrite for 2D, custom function for 3D
        if ndims(imageData_scaled) == 2
            imwrite(imageData_scaled, filename, 'TIFF', 'Compression', 'none');
        else
            % 3D TIFF stack
            for z = 1:size(imageData_scaled, 3)
                if z == 1
                    imwrite(imageData_scaled(:,:,z), filename, 'TIFF', 'Compression', 'none');
                else
                    imwrite(imageData_scaled(:,:,z), filename, 'TIFF', 'Compression', 'none', 'WriteMode', 'append');
                end
            end
        end
        
    case 'png'
        % PNG only supports 2D, so use max projection for 3D
        if ndims(imageData_scaled) == 3
            imageData_scaled = max(imageData_scaled, [], 3);
        end
        imwrite(imageData_scaled, filename, 'PNG');
        
    case {'jpg', 'jpeg'}
        % JPEG only supports 2D and 8-bit, so convert and project
        if ndims(imageData_scaled) == 3
            imageData_scaled = max(imageData_scaled, [], 3);
        end
        % Convert to 8-bit
        imageData_8bit = uint8(double(imageData_scaled) / 65535 * 255);
        imwrite(imageData_8bit, filename, 'JPEG', 'Quality', 95);
        
    otherwise
        % Default to TIFF
        if ndims(imageData_scaled) == 2
            imwrite(imageData_scaled, filename, 'TIFF', 'Compression', 'none');
        else
            for z = 1:size(imageData_scaled, 3)
                if z == 1
                    imwrite(imageData_scaled(:,:,z), filename, 'TIFF', 'Compression', 'none');
                else
                    imwrite(imageData_scaled(:,:,z), filename, 'TIFF', 'Compression', 'none', 'WriteMode', 'append');
                end
            end
        end
end
end

function description = getNucleiProcessingDescription(handles)
% Generate description of nuclei processing steps applied
description = {};

% Check preprocessing
if isfield(handles, 'nucPreprocEnabledCheck') && isvalid(handles.nucPreprocEnabledCheck) && handles.nucPreprocEnabledCheck.Value
    if isfield(handles, 'nucPreprocMethodDrop') && isvalid(handles.nucPreprocMethodDrop)
        method = handles.nucPreprocMethodDrop.Value;
        if ~strcmp(method, 'None')
            description{end+1} = sprintf('Preprocessing: %s', method);
            
            % Add method-specific parameters
            switch method
                case 'Gaussian'
                    if isfield(handles, 'nucGaussInput') && isvalid(handles.nucGaussInput)
                        description{end} = sprintf('%s (sigma=%.2f)', description{end}, handles.nucGaussInput.Value);
                    end
                case 'Median'
                    if isfield(handles, 'nucMedianInput') && isvalid(handles.nucMedianInput)
                        description{end} = sprintf('%s (size=%.1f)', description{end}, handles.nucMedianInput.Value);
                    end
                case 'Non-Local Means'
                    if isfield(handles, 'nucNlmFilterStrengthInput') && isvalid(handles.nucNlmFilterStrengthInput)
                        description{end} = sprintf('%s (strength=%.2f)', description{end}, handles.nucNlmFilterStrengthInput.Value);
                    end
            end
        end
    end
end

% Check background correction
if isfield(handles, 'nucBgCorrEnabledCheck') && isvalid(handles.nucBgCorrEnabledCheck) && handles.nucBgCorrEnabledCheck.Value
    if isfield(handles, 'nucBgMethodDrop') && isvalid(handles.nucBgMethodDrop)
        method = handles.nucBgMethodDrop.Value;
        if ~strcmp(method, 'None')
            bgDesc = sprintf('Background Correction: %s', method);
            if isfield(handles, 'nucBgParamInput') && isvalid(handles.nucBgParamInput)
                bgDesc = sprintf('%s (param=%.1f)', bgDesc, handles.nucBgParamInput.Value);
            end
            description{end+1} = bgDesc;
        end
    end
end

% Check clipping
if isfield(handles, 'nucBgCorrClipChecks') && isvalid(handles.nucBgCorrClipChecks) && handles.nucBgCorrClipChecks.Value
    description{end+1} = 'Clipping: Negative values set to zero';
end

if isempty(description)
    description = {'No processing applied (raw image)'};
end
end

function description = getChannelProcessingDescription(handles, ch)
% Generate description of channel processing steps applied
description = {};

% Check preprocessing
if isfield(handles, 'preprocEnabledChecks') && length(handles.preprocEnabledChecks) >= ch && handles.preprocEnabledChecks(ch).Value
    if isfield(handles, 'preprocMethodDrops') && length(handles.preprocMethodDrops) >= ch
        method = handles.preprocMethodDrops(ch).Value;
        if ~strcmp(method, 'None')
            description{end+1} = sprintf('Preprocessing: %s', method);
            
            % Add method-specific parameters
            switch method
                case 'Gaussian'
                    if isfield(handles, 'gaussInputs') && length(handles.gaussInputs) >= ch
                        description{end} = sprintf('%s (sigma=%.2f)', description{end}, handles.gaussInputs(ch).Value);
                    end
                case 'Median'
                    if isfield(handles, 'medianInputs') && length(handles.medianInputs) >= ch
                        description{end} = sprintf('%s (size=%.1f)', description{end}, handles.medianInputs(ch).Value);
                    end
                case 'Non-Local Means'
                    if isfield(handles, 'nlmFilterStrengthInputs') && length(handles.nlmFilterStrengthInputs) >= ch
                        description{end} = sprintf('%s (strength=%.2f)', description{end}, handles.nlmFilterStrengthInputs(ch).Value);
                    end
            end
        end
    end
end

% Check background correction
if isfield(handles, 'bgCorrEnabledChecks') && length(handles.bgCorrEnabledChecks) >= ch && handles.bgCorrEnabledChecks(ch).Value
    if isfield(handles, 'bgMethodDrops') && length(handles.bgMethodDrops) >= ch
        method = handles.bgMethodDrops(ch).Value;
        if ~strcmp(method, 'None')
            bgDesc = sprintf('Background Correction: %s', method);
            if isfield(handles, 'bgParamInputs') && length(handles.bgParamInputs) >= ch
                bgDesc = sprintf('%s (param=%.1f)', bgDesc, handles.bgParamInputs(ch).Value);
            end
            description{end+1} = bgDesc;
        end
    end
end

% Check clipping
if isfield(handles, 'bgCorrClipChecks') && length(handles.bgCorrClipChecks) >= ch && handles.bgCorrClipChecks(ch).Value
    description{end+1} = 'Clipping: Negative values set to zero';
end

if isempty(description)
    description = {'No processing applied (raw image)'};
end
end

function createProcessedImageSummary(metadata, filename)
% Create human-readable summary of exported processed images
fid = fopen(filename, 'w');
if fid == -1
    warning('Could not create processed image summary file: %s', filename);
    return;
end

fprintf(fid, 'Processed Images Export Summary\n');
fprintf(fid, '==============================\n\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Export Timestamp: %s\n\n', metadata.exportTimestamp);

fprintf(fid, 'Description:\n');
fprintf(fid, '%s\n\n', metadata.description);

% Nuclei information
if isfield(metadata, 'nuclei')
    fprintf(fid, 'NUCLEI IMAGE:\n');
    if isfield(metadata.nuclei, 'filename')
        fprintf(fid, '  File: %s\n', metadata.nuclei.filename);
        fprintf(fid, '  Original Size: %s\n', mat2str(metadata.nuclei.original_size));
        fprintf(fid, '  Processed Size: %s\n', mat2str(metadata.nuclei.processed_size));
        
        if isfield(metadata.nuclei, 'xy_spacing')
            fprintf(fid, '  XY Spacing: %.3f µm/pixel\n', metadata.nuclei.xy_spacing);
        end
        if isfield(metadata.nuclei, 'z_spacing')
            fprintf(fid, '  Z Spacing: %.3f µm/slice\n', metadata.nuclei.z_spacing);
        end
        
        fprintf(fid, '  Processing Applied:\n');
        for i = 1:length(metadata.nuclei.processing_applied)
            fprintf(fid, '    - %s\n', metadata.nuclei.processing_applied{i});
        end
    elseif isfield(metadata.nuclei, 'error')
        fprintf(fid, '  ERROR: %s\n', metadata.nuclei.error);
    end
    fprintf(fid, '\n');
end

% Channel information
if isfield(metadata, 'channels')
    channelFields = fieldnames(metadata.channels);
    fprintf(fid, 'FLUORESCENT CHANNELS:\n');
    
    for i = 1:length(channelFields)
        channelName = channelFields{i};
        channelData = metadata.channels.(channelName);
        
        fprintf(fid, '  %s:\n', upper(strrep(channelName, '_', ' ')));
        
        if isfield(channelData, 'filename')
            fprintf(fid, '    File: %s\n', channelData.filename);
            fprintf(fid, '    Original Size: %s\n', mat2str(channelData.original_size));
            fprintf(fid, '    Processed Size: %s\n', mat2str(channelData.processed_size));
            
            if isfield(channelData, 'xy_spacing')
                fprintf(fid, '    XY Spacing: %.3f µm/pixel\n', channelData.xy_spacing);
            end
            if isfield(channelData, 'z_spacing')
                fprintf(fid, '    Z Spacing: %.3f µm/slice\n', channelData.z_spacing);
            end
            
            fprintf(fid, '    Processing Applied:\n');
            for j = 1:length(channelData.processing_applied)
                fprintf(fid, '      - %s\n', channelData.processing_applied{j});
            end
        elseif isfield(channelData, 'error')
            fprintf(fid, '    ERROR: %s\n', channelData.error);
        elseif isfield(channelData, 'status')
            fprintf(fid, '    STATUS: %s\n', channelData.status);
        end
        fprintf(fid, '\n');
    end
end

% DIC information
if isfield(metadata, 'dic')
    fprintf(fid, 'DIC IMAGE:\n');
    if isfield(metadata.dic, 'filename')
        fprintf(fid, '  File: %s\n', metadata.dic.filename);
        fprintf(fid, '  Size: %s\n', mat2str(metadata.dic.size));
        fprintf(fid, '  Processing: %s\n', metadata.dic.processing_applied);
    elseif isfield(metadata.dic, 'error')
        fprintf(fid, '  ERROR: %s\n', metadata.dic.error);
    end
    fprintf(fid, '\n');
end

fprintf(fid, 'USAGE NOTES:\n');
fprintf(fid, '- All images are saved in the same format and bit depth\n');
fprintf(fid, '- 3D images are saved as multi-page TIFF stacks\n');
fprintf(fid, '- Processing parameters are preserved in the metadata file\n');
fprintf(fid, '- These images match exactly what is displayed in the preview panels\n');
fprintf(fid, '- Intensity values are scaled to preserve the full dynamic range\n');

fclose(fid);
fprintf('Processed image summary exported to: %s\n', filename);
end

function exportNucleiImage(handles, exportDir, timestamp)
% Export processed nuclei image (uses cached data from Update Previews)
fprintf('Exporting nuclei image...\n');

img = [];

% PRIMARY SOURCE: Get from handles.processedNuclei (populated by Update Previews)
if isfield(handles, 'processedNuclei') && ~isempty(handles.processedNuclei)
    img = handles.processedNuclei;
% FALLBACK 1: Try preview cache
elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ...
       isfield(handles.previewCache.nuclei, 'processed') && ~isempty(handles.previewCache.nuclei.processed)
    img = handles.previewCache.nuclei.processed;
% FALLBACK 2: Recompute only if necessary
else
    fprintf('  Warning: No cached nuclei data, recomputing...\n');
    try
        img = snap_helpers.preprocessNucleiWithBgCorr(handles);
    catch ME
        warning('Could not process nuclei image: %s', ME.message);
        return;
    end
end

if isempty(img)
    warning('No nuclei image available.');
    return;
end

exportFormat = lower(handles.exportImageFormatDrop.Value);
filename = fullfile(exportDir, sprintf('export_nuclei_image.%s', exportFormat));
saveProcessedImage(img, filename, exportFormat);
fprintf('Nuclei image exported to: %s\n', filename);
end

function exportChannelImage(handles, exportDir, timestamp, ch)
% Export processed channel image (uses cached data from Update Previews)
fprintf('Exporting Channel %d image...\n', ch);

img = [];

% PRIMARY SOURCE: Get from handles.processedChannel (populated by Update Previews)
if isfield(handles, 'processedChannel') && ch <= length(handles.processedChannel) && ~isempty(handles.processedChannel{ch})
    img = handles.processedChannel{ch};
% FALLBACK 1: Try preview cache
elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
       length(handles.previewCache.channels) >= ch && ...
       isfield(handles.previewCache.channels{ch}, 'processed') && ~isempty(handles.previewCache.channels{ch}.processed)
    img = handles.previewCache.channels{ch}.processed;
% FALLBACK 2: Recompute only if necessary
else
    fprintf('  Warning: No cached channel %d data, recomputing...\n', ch);
    try
        img = snap_helpers.processImage(handles, ch);
    catch ME
        warning('Could not process Channel %d image: %s', ch, ME.message);
        return;
    end
end

if isempty(img)
    warning('No Channel %d image available.', ch);
    return;
end

exportFormat = lower(handles.exportImageFormatDrop.Value);
filename = fullfile(exportDir, sprintf('export_ch%d_image.%s', ch, exportFormat));
saveProcessedImage(img, filename, exportFormat);
fprintf('Channel %d image exported to: %s\n', ch, filename);
end

function lastUsed = collectCurrentParameters(handles)
% Collect ALL current parameter values directly from UI elements
% This reads the CURRENT state of all UI controls, not cached values
% Returns a lastUsed structure compatible with SNAP_batch loading

% Start with defaults to ensure all required fields exist
lastUsed = snap_helpers.initializeParameters(handles.Nmax);

try
    % === BASIC SETTINGS ===
    if isfield(handles, 'numChanDrop')
        lastUsed.numChannels = str2double(handles.numChanDrop.Value);
    end
    
    % === NUCLEI SPACING ===
    if isfield(handles, 'nucXYSpacingInput')
        lastUsed.nucXYSpacing = handles.nucXYSpacingInput.Value;
    end
    if isfield(handles, 'nucZSpacingInput')
        lastUsed.nucZSpacing = handles.nucZSpacingInput.Value;
    end
    
    % === NUCLEI DECONVOLUTION ===
    if isfield(handles, 'nucDeconvEnabledCheck')
        lastUsed.nucDeconvEnabled = handles.nucDeconvEnabledCheck.Value;
    end
    if isfield(handles, 'nucDeconvMethodDrop')
        lastUsed.nucDeconvMethod = handles.nucDeconvMethodDrop.Value;
    end
    if isfield(handles, 'nucDeconvLRIterationsInput')
        lastUsed.nucDeconvLRIterations = handles.nucDeconvLRIterationsInput.Value;
    end
    if isfield(handles, 'nucDeconvLRDampingInput')
        lastUsed.nucDeconvLRDamping = handles.nucDeconvLRDampingInput.Value;
    end
    if isfield(handles, 'nucDeconvWienerNSRInput')
        lastUsed.nucDeconvWienerNSR = handles.nucDeconvWienerNSRInput.Value;
    end
    if isfield(handles, 'nucDeconvBlindIterationsInput')
        lastUsed.nucDeconvBlindIterations = handles.nucDeconvBlindIterationsInput.Value;
    end
    if isfield(handles, 'nucDeconvBlindUnderRelaxInput')
        lastUsed.nucDeconvBlindUnderRelax = handles.nucDeconvBlindUnderRelaxInput.Value;
    end
    if isfield(handles, 'nucDeconvPSFSourceDrop')
        lastUsed.nucDeconvPSFSource = handles.nucDeconvPSFSourceDrop.Value;
    end
    if isfield(handles, 'nucDeconvPSFPathText')
        lastUsed.nucDeconvPSFFilePath = handles.nucDeconvPSFPathText.Value;
    end
    if isfield(handles, 'nucDeconvPSFSigmaXYInput')
        lastUsed.nucDeconvPSFSigmaXY = handles.nucDeconvPSFSigmaXYInput.Value;
    end
    if isfield(handles, 'nucDeconvPSFSigmaZInput')
        lastUsed.nucDeconvPSFSigmaZ = handles.nucDeconvPSFSigmaZInput.Value;
    end
    if isfield(handles, 'nucDeconvPSFSizeXYInput')
        lastUsed.nucDeconvPSFSizeXY = handles.nucDeconvPSFSizeXYInput.Value;
    end
    if isfield(handles, 'nucDeconvPSFSizeZInput')
        lastUsed.nucDeconvPSFSizeZ = handles.nucDeconvPSFSizeZInput.Value;
    end
    
    % === NUCLEI PREPROCESSING ===
    if isfield(handles, 'nucPreprocEnabledCheck')
        lastUsed.nucPreprocEnabled = handles.nucPreprocEnabledCheck.Value;
    end
    if isfield(handles, 'nucPreprocessModeDrop')
        lastUsed.nucPreProcMode = handles.nucPreprocessModeDrop.Value;
    end
    if isfield(handles, 'nucPreprocessScaleCheck')
        lastUsed.nucPreProcScale = handles.nucPreprocessScaleCheck.Value;
    end
    if isfield(handles, 'nucPreprocessProjectionDrop')
        lastUsed.nucPreProcProjection = handles.nucPreprocessProjectionDrop.Value;
    end
    if isfield(handles, 'nucPreprocMethodDrop')
        lastUsed.nucPreProcMethod = handles.nucPreprocMethodDrop.Value;
    end
    if isfield(handles, 'nucPreprocClipChecks')
        lastUsed.nucPreprocClipAtZero = handles.nucPreprocClipChecks.Value;
    end
    
    % Nuclei preprocessing method-specific parameters
    if isfield(handles, 'nucPreprocParam1Inputs')
        lastUsed.nucSmoothGaussianValue = handles.nucPreprocParam1Inputs.Value;
    end
    if isfield(handles, 'nucPreprocParam2Inputs')
        lastUsed.nucSmoothMedianValue = handles.nucPreprocParam2Inputs.Value;
    end
    if isfield(handles, 'nucNlmFilterStrengthInput')
        lastUsed.nucNlmFilterStrength = handles.nucNlmFilterStrengthInput.Value;
    end
    if isfield(handles, 'nucNlmSearchWindowInput')
        lastUsed.nucNlmSearchWindow = handles.nucNlmSearchWindowInput.Value;
    end
    if isfield(handles, 'nucNlmComparisonWindowInput')
        lastUsed.nucNlmComparisonWindow = handles.nucNlmComparisonWindowInput.Value;
    end
    if isfield(handles, 'nucWaveletNameDrop')
        lastUsed.nucWaveletName = handles.nucWaveletNameDrop.Value;
    end
    if isfield(handles, 'nucWaveletLevelInput')
        lastUsed.nucWaveletLevel = handles.nucWaveletLevelInput.Value;
    end
    if isfield(handles, 'nucWaveletThresholdRuleDrop')
        lastUsed.nucWaveletThresholdRule = handles.nucWaveletThresholdRuleDrop.Value;
    end
    if isfield(handles, 'nucWaveletThresholdMethodDrop')
        lastUsed.nucWaveletThresholdMethod = handles.nucWaveletThresholdMethodDrop.Value;
    end
    
    % === NUCLEI BACKGROUND CORRECTION ===
    if isfield(handles, 'nucBgCorrEnabledCheck')
        lastUsed.nucBgCorrEnabled = handles.nucBgCorrEnabledCheck.Value;
    end
    if isfield(handles, 'nucBgCorrModeDrop')
        lastUsed.nucBgCorrMode = handles.nucBgCorrModeDrop.Value;
    end
    if isfield(handles, 'nucBgMethodDrop')
        lastUsed.nucBgMethod = handles.nucBgMethodDrop.Value;
    end
    if isfield(handles, 'nucBgParamInput')
        lastUsed.nucBgParam = handles.nucBgParamInput.Value;
    end
    if isfield(handles, 'nucBgCorrScaleCheck')
        lastUsed.nucBgCorrScale = handles.nucBgCorrScaleCheck.Value;
    end
    if isfield(handles, 'nucBgCorrProjectionDrop')
        lastUsed.nucBgCorrProjection = handles.nucBgCorrProjectionDrop.Value;
    end
    if isfield(handles, 'nucBgCorrClipCheck')
        lastUsed.nucBgCorrClipAtZero = handles.nucBgCorrClipCheck.Value;
    end
    
    % === NUCLEI SEGMENTATION ===
    if isfield(handles, 'nucSegEnabledCheck')
        lastUsed.nucSegEnabled = handles.nucSegEnabledCheck.Value;
    end
    if isfield(handles, 'nucSegModeDrop')
        lastUsed.nucSegMode = handles.nucSegModeDrop.Value;
    end
    if isfield(handles, 'nucSegProjectionDrop')
        lastUsed.nucSegProjection = handles.nucSegProjectionDrop.Value;
    end
    if isfield(handles, 'nucSegMainMethodDrop')
        lastUsed.nucSegMainMethod = handles.nucSegMainMethodDrop.Value;
    end
    if isfield(handles, 'nucSegSubMethodDrop')
        lastUsed.nucSegSubMethod = handles.nucSegSubMethodDrop.Value;
    end
    if isfield(handles, 'nucSegAbsoluteThresholdInput')
        lastUsed.nucSegAbsoluteThreshold = handles.nucSegAbsoluteThresholdInput.Value;
    end
    if isfield(handles, 'nucSegStdMultiplierInput')
        lastUsed.nucSegStdMultiplier = handles.nucSegStdMultiplierInput.Value;
    end
    if isfield(handles, 'nucSegOffsetInput')
        lastUsed.nucSegOffset = handles.nucSegOffsetInput.Value;
    end
    if isfield(handles, 'nucSegLocalAlgorithmDrop')
        lastUsed.nucSegLocalAlgorithm = handles.nucSegLocalAlgorithmDrop.Value;
    end
    if isfield(handles, 'nucSegAlgParamInput')
        lastUsed.nucSegAlgParam1 = handles.nucSegAlgParamInput.Value;
    end
    if isfield(handles, 'nucSegAlgParam2Input')
        lastUsed.nucSegAlgParam2 = handles.nucSegAlgParam2Input.Value;
    end
    if isfield(handles, 'nucShowSegCheck')
        lastUsed.nucShowSeg = handles.nucShowSegCheck.Value;
    end
    
    % Algorithm-specific parameters
    % NOTE: These are mapped from nucSegAlgParamInput/nucSegAlgParam2Input
    % based on the selected algorithm in nucSegLocalAlgorithmDrop
    if isfield(handles, 'nucSegLocalAlgorithmDrop')
        algorithm = handles.nucSegLocalAlgorithmDrop.Value;
        if isfield(handles, 'nucSegAlgParamInput')
            algParam1 = handles.nucSegAlgParamInput.Value;
        else
            algParam1 = 0;
        end
        if isfield(handles, 'nucSegAlgParam2Input')
            algParam2 = handles.nucSegAlgParam2Input.Value;
        else
            algParam2 = 0;
        end
        
        % Map algorithm-specific parameters to their named fields
        switch algorithm
            case 'Bernsen'
                lastUsed.nucSegBernsenContrast = algParam1;
            case 'Mean'
                lastUsed.nucSegMeanC = algParam1;
            case 'Median'
                lastUsed.nucSegMedianC = algParam1;
            case 'MidGrey'
                lastUsed.nucSegMidGreyC = algParam1;
            case 'Niblack'
                lastUsed.nucSegNiblackK = algParam1;
                lastUsed.nucSegNiblackC = algParam2;
            case 'Phansalkar'
                lastUsed.nucSegPhansalkarK = algParam1;
                lastUsed.nucSegPhansalkarR = algParam2;
            case 'Sauvola'
                lastUsed.nucSegSauvolaK = algParam1;
                lastUsed.nucSegSauvolaR = algParam2;
            % Otsu has no parameters
        end
    end
    
    % Default checkbox states
    % NOTE: These are mapped from nucSegAlgParamDefaultCheck/nucSegAlgParam2DefaultCheck
    % based on the selected algorithm in nucSegLocalAlgorithmDrop
    if isfield(handles, 'nucSegLocalAlgorithmDrop')
        algorithm = handles.nucSegLocalAlgorithmDrop.Value;
        if isfield(handles, 'nucSegAlgParamDefaultCheck')
            algParam1Default = handles.nucSegAlgParamDefaultCheck.Value;
        else
            algParam1Default = false;
        end
        if isfield(handles, 'nucSegAlgParam2DefaultCheck')
            algParam2Default = handles.nucSegAlgParam2DefaultCheck.Value;
        else
            algParam2Default = false;
        end
        
        % Map algorithm-specific default checkbox states to their named fields
        switch algorithm
            case 'Bernsen'
                lastUsed.nucSegBernsenContrastDefault = algParam1Default;
            case 'Mean'
                lastUsed.nucSegMeanCDefault = algParam1Default;
            case 'Median'
                lastUsed.nucSegMedianCDefault = algParam1Default;
            case 'MidGrey'
                lastUsed.nucSegMidGreyCDefault = algParam1Default;
            case 'Niblack'
                lastUsed.nucSegNiblackKDefault = algParam1Default;
                lastUsed.nucSegNiblackCDefault = algParam2Default;
            case 'Phansalkar'
                lastUsed.nucSegPhansalkarKDefault = algParam1Default;
                lastUsed.nucSegPhansalkarRDefault = algParam2Default;
            case 'Sauvola'
                lastUsed.nucSegSauvolaKDefault = algParam1Default;
                lastUsed.nucSegSauvolaRDefault = algParam2Default;
            % Otsu has no parameters
        end
    end
    
    % === NUCLEI FILTERING ===
    if isfield(handles, 'nucFilterEnabledCheck')
        lastUsed.nucFilterEnabled = handles.nucFilterEnabledCheck.Value;
    end
    if isfield(handles, 'nucFilterSizeEnabledCheck')
        lastUsed.nucFilterSizeEnabled = handles.nucFilterSizeEnabledCheck.Value;
    end
    if isfield(handles, 'nucFilterMinSizeInput')
        lastUsed.nucFilterMinSize = handles.nucFilterMinSizeInput.Value;
    end
    if isfield(handles, 'nucFilterSizeUnitDrop')
        lastUsed.nucFilterSizeUnit = handles.nucFilterSizeUnitDrop.Value;
    end
    if isfield(handles, 'nucFilterCircularityEnabledCheck')
        lastUsed.nucFilterCircularityEnabled = handles.nucFilterCircularityEnabledCheck.Value;
    end
    if isfield(handles, 'nucFilterMinCircularityInput')
        lastUsed.nucFilterMinCircularity = handles.nucFilterMinCircularityInput.Value;
    end
    if isfield(handles, 'nucFilterSolidityEnabledCheck')
        lastUsed.nucFilterSolidityEnabled = handles.nucFilterSolidityEnabledCheck.Value;
    end
    if isfield(handles, 'nucFilterMinSolidityInput')
        lastUsed.nucFilterMinSolidity = handles.nucFilterMinSolidityInput.Value;
    end
    if isfield(handles, 'nucExcludeEdgesCheck')
        lastUsed.nucExcludeEdges = handles.nucExcludeEdgesCheck.Value;
    end
    
    % === NUCLEI INCLUSION/EXCLUSION ===
    if isfield(handles, 'nucInclusionExclusionEnabledCheck')
        lastUsed.nucInclusionExclusionEnabled = handles.nucInclusionExclusionEnabledCheck.Value;
    end
    if isfield(handles, 'nucInclusionExclusionModeDrop')
        lastUsed.nucInclusionExclusionMode = handles.nucInclusionExclusionModeDrop.Value;
    end
    if isfield(handles, 'nucInclusionExclusionApplyDrop')
        lastUsed.nucInclusionExclusionApplyTo = handles.nucInclusionExclusionApplyDrop.Value;
    end
    
    % === CHANNEL PARAMETERS ===
    numChannels = lastUsed.numChannels;
    % Export ALL channels up to Nmax (typically 5), not just active ones
    for ch = 1:handles.Nmax
        try
        % Spacing
        if isfield(handles, 'xySpacingInputs') && length(handles.xySpacingInputs) >= ch
            lastUsed.xySpacing{ch} = handles.xySpacingInputs(ch).Value;
        end
        if isfield(handles, 'zSpacingInputs') && length(handles.zSpacingInputs) >= ch
            lastUsed.zSpacing{ch} = handles.zSpacingInputs(ch).Value;
        end
        
        % Deconvolution
        if isfield(handles, 'deconvEnabledChecks') && length(handles.deconvEnabledChecks) >= ch
            lastUsed.deconvEnabled{ch} = handles.deconvEnabledChecks(ch).Value;
        end
        if isfield(handles, 'deconvMethodDrops') && length(handles.deconvMethodDrops) >= ch
            lastUsed.deconvMethod{ch} = handles.deconvMethodDrops(ch).Value;
        end
        if isfield(handles, 'deconvLRIterationsInputs') && length(handles.deconvLRIterationsInputs) >= ch
            lastUsed.deconvLRIterations{ch} = handles.deconvLRIterationsInputs(ch).Value;
        end
        if isfield(handles, 'deconvLRDampingInputs') && length(handles.deconvLRDampingInputs) >= ch
            lastUsed.deconvLRDamping{ch} = handles.deconvLRDampingInputs(ch).Value;
        end
        if isfield(handles, 'deconvWienerNSRInputs') && length(handles.deconvWienerNSRInputs) >= ch
            lastUsed.deconvWienerNSR{ch} = handles.deconvWienerNSRInputs(ch).Value;
        end
        if isfield(handles, 'deconvBlindIterationsInputs') && length(handles.deconvBlindIterationsInputs) >= ch
            lastUsed.deconvBlindIterations{ch} = handles.deconvBlindIterationsInputs(ch).Value;
        end
        if isfield(handles, 'deconvBlindUnderRelaxInputs') && length(handles.deconvBlindUnderRelaxInputs) >= ch
            lastUsed.deconvBlindUnderRelax{ch} = handles.deconvBlindUnderRelaxInputs(ch).Value;
        end
        if isfield(handles, 'deconvPSFSourceDrops') && length(handles.deconvPSFSourceDrops) >= ch
            lastUsed.deconvPSFSource{ch} = handles.deconvPSFSourceDrops(ch).Value;
        end
        if isfield(handles, 'deconvPSFPathTexts') && length(handles.deconvPSFPathTexts) >= ch
            lastUsed.deconvPSFFilePath{ch} = handles.deconvPSFPathTexts(ch).Value;
        end
        if isfield(handles, 'deconvPSFSigmaXYInputs') && length(handles.deconvPSFSigmaXYInputs) >= ch
            lastUsed.deconvPSFSigmaXY{ch} = handles.deconvPSFSigmaXYInputs(ch).Value;
        end
        if isfield(handles, 'deconvPSFSigmaZInputs') && length(handles.deconvPSFSigmaZInputs) >= ch
            lastUsed.deconvPSFSigmaZ{ch} = handles.deconvPSFSigmaZInputs(ch).Value;
        end
        if isfield(handles, 'deconvPSFSizeXYInputs') && length(handles.deconvPSFSizeXYInputs) >= ch
            lastUsed.deconvPSFSizeXY{ch} = handles.deconvPSFSizeXYInputs(ch).Value;
        end
        if isfield(handles, 'deconvPSFSizeZInputs') && length(handles.deconvPSFSizeZInputs) >= ch
            lastUsed.deconvPSFSizeZ{ch} = handles.deconvPSFSizeZInputs(ch).Value;
        end
        
        % Preprocessing
        if isfield(handles, 'preprocEnabledChecks') && length(handles.preprocEnabledChecks) >= ch
            lastUsed.preprocEnabled{ch} = handles.preprocEnabledChecks(ch).Value;
        end
        if isfield(handles, 'preprocessModeDrops') && length(handles.preprocessModeDrops) >= ch
            lastUsed.preProcMode{ch} = handles.preprocessModeDrops(ch).Value;
        end
        if isfield(handles, 'preprocessProjectionDrops') && length(handles.preprocessProjectionDrops) >= ch
            lastUsed.preProcProjection{ch} = handles.preprocessProjectionDrops(ch).Value;
        end
        if isfield(handles, 'preprocessScaleChecks') && length(handles.preprocessScaleChecks) >= ch
            lastUsed.preProcScale{ch} = handles.preprocessScaleChecks(ch).Value;
        end
        if isfield(handles, 'preprocMethodDrops') && length(handles.preprocMethodDrops) >= ch
            lastUsed.preProcMethod{ch} = handles.preprocMethodDrops(ch).Value;
        end
        if isfield(handles, 'preprocClipChecks') && length(handles.preprocClipChecks) >= ch
            lastUsed.preprocClipAtZero{ch} = handles.preprocClipChecks(ch).Value;
        end
        
        % Method-specific preprocessing parameters
        if isfield(handles, 'gaussInputs') && length(handles.gaussInputs) >= ch
            lastUsed.smoothGaussianValues{ch} = handles.gaussInputs(ch).Value;
        end
        if isfield(handles, 'medianInputs') && length(handles.medianInputs) >= ch
            lastUsed.smoothMedianValues{ch} = handles.medianInputs(ch).Value;
        end
        if isfield(handles, 'nlmFilterStrengthInputs') && length(handles.nlmFilterStrengthInputs) >= ch
            lastUsed.nlmFilterStrength{ch} = handles.nlmFilterStrengthInputs(ch).Value;
        end
        if isfield(handles, 'nlmSearchWindowInputs') && length(handles.nlmSearchWindowInputs) >= ch
            lastUsed.nlmSearchWindow{ch} = handles.nlmSearchWindowInputs(ch).Value;
        end
        if isfield(handles, 'nlmComparisonWindowInputs') && length(handles.nlmComparisonWindowInputs) >= ch
            lastUsed.nlmComparisonWindow{ch} = handles.nlmComparisonWindowInputs(ch).Value;
        end
        if isfield(handles, 'waveletNameDrops') && length(handles.waveletNameDrops) >= ch
            lastUsed.waveletName{ch} = handles.waveletNameDrops(ch).Value;
        end
        if isfield(handles, 'waveletLevelInputs') && length(handles.waveletLevelInputs) >= ch
            lastUsed.waveletLevel{ch} = handles.waveletLevelInputs(ch).Value;
        end
        if isfield(handles, 'waveletThresholdRuleDrops') && length(handles.waveletThresholdRuleDrops) >= ch
            lastUsed.waveletThresholdRule{ch} = handles.waveletThresholdRuleDrops(ch).Value;
        end
        if isfield(handles, 'waveletThresholdMethodDrops') && length(handles.waveletThresholdMethodDrops) >= ch
            lastUsed.waveletThresholdMethod{ch} = handles.waveletThresholdMethodDrops(ch).Value;
        end
        
        % Background Correction
        if isfield(handles, 'bgCorrEnabledChecks') && length(handles.bgCorrEnabledChecks) >= ch
            lastUsed.bgCorrEnabled{ch} = handles.bgCorrEnabledChecks(ch).Value;
        end
        if isfield(handles, 'bgCorrModeDrops') && length(handles.bgCorrModeDrops) >= ch
            lastUsed.bgCorrMode{ch} = handles.bgCorrModeDrops(ch).Value;
        end
        if isfield(handles, 'bgCorrScaleChecks') && length(handles.bgCorrScaleChecks) >= ch
            lastUsed.bgCorrScale{ch} = handles.bgCorrScaleChecks(ch).Value;
        end
        if isfield(handles, 'bgCorrProjectionDrops') && length(handles.bgCorrProjectionDrops) >= ch
            lastUsed.bgCorrProjection{ch} = handles.bgCorrProjectionDrops(ch).Value;
        end
        if isfield(handles, 'bgMethodDrops') && length(handles.bgMethodDrops) >= ch
            lastUsed.bgMethod{ch} = handles.bgMethodDrops(ch).Value;
        end
        if isfield(handles, 'bgParamInputs') && length(handles.bgParamInputs) >= ch
            lastUsed.bgParam{ch} = handles.bgParamInputs(ch).Value;
        end
        if isfield(handles, 'bgCorrClipChecks') && length(handles.bgCorrClipChecks) >= ch
            lastUsed.bgCorrClipAtZero{ch} = handles.bgCorrClipChecks(ch).Value;
        end
        
        % Maxima Detection
        if isfield(handles, 'maximaEnabledChecks') && length(handles.maximaEnabledChecks) >= ch
            lastUsed.maximaEnabled{ch} = handles.maximaEnabledChecks(ch).Value;
        end
        if isfield(handles, 'maximaModeDrops') && length(handles.maximaModeDrops) >= ch
            lastUsed.maximaMode{ch} = handles.maximaModeDrops(ch).Value;
        end
        if isfield(handles, 'maximaScaleChecks') && length(handles.maximaScaleChecks) >= ch
            lastUsed.maximaScale{ch} = handles.maximaScaleChecks(ch).Value;
        end
        if isfield(handles, 'maximaProjectionDrops') && length(handles.maximaProjectionDrops) >= ch
            lastUsed.maximaProjection{ch} = handles.maximaProjectionDrops(ch).Value;
        end
        if isfield(handles, 'maximaMethodDrops') && length(handles.maximaMethodDrops) >= ch
            lastUsed.maximaMethod{ch} = handles.maximaMethodDrops(ch).Value;
        end
        if isfield(handles, 'maximaNeighborhoodInputs') && length(handles.maximaNeighborhoodInputs) >= ch
            lastUsed.maximaNeighborhoodSize{ch} = handles.maximaNeighborhoodInputs(ch).Value;
        end
        if isfield(handles, 'hMaxInputs') && length(handles.hMaxInputs) >= ch
            lastUsed.hMaxValue{ch} = handles.hMaxInputs(ch).Value;
        end
        if isfield(handles, 'logSigmaInputs') && length(handles.logSigmaInputs) >= ch
            lastUsed.sigmaValue{ch} = handles.logSigmaInputs(ch).Value;
        end
        if isfield(handles, 'logThresholdInputs') && length(handles.logThresholdInputs) >= ch
            lastUsed.peakThresholdValue{ch} = handles.logThresholdInputs(ch).Value;
        end
        if isfield(handles, 'showMaximaChecks') && length(handles.showMaximaChecks) >= ch && isvalid(handles.showMaximaChecks(ch))
            lastUsed.showMaxima{ch} = handles.showMaximaChecks(ch).Value;
        end
        if isfield(handles, 'maximaColorDrops') && length(handles.maximaColorDrops) >= ch && isvalid(handles.maximaColorDrops(ch))
            lastUsed.maximaColor{ch} = handles.maximaColorDrops(ch).Value;
        end
        if isfield(handles, 'displayOnAllPreviewsChecks') && length(handles.displayOnAllPreviewsChecks) >= ch && isvalid(handles.displayOnAllPreviewsChecks(ch))
            lastUsed.displayOnAllPreviews{ch} = handles.displayOnAllPreviewsChecks(ch).Value;
        end
        
        % Gaussian Fitting
        if isfield(handles, 'gaussFitEnabledChecks') && length(handles.gaussFitEnabledChecks) >= ch
            lastUsed.gaussFitEnabled{ch} = handles.gaussFitEnabledChecks(ch).Value;
        end
        if isfield(handles, 'gaussFitMethodDrop') && length(handles.gaussFitMethodDrop) >= ch
            lastUsed.gaussFitMethod{ch} = handles.gaussFitMethodDrop(ch).Value;
        end
        if isfield(handles, 'gaussFitVoxelWindowSlider') && length(handles.gaussFitVoxelWindowSlider) >= ch
            lastUsed.gaussFitVoxelWindowSize{ch} = handles.gaussFitVoxelWindowSlider(ch).Value;
        end
        if isfield(handles, 'gaussFitBgCorrMethodDrop') && length(handles.gaussFitBgCorrMethodDrop) >= ch
            lastUsed.gaussFitBgCorrMethod{ch} = handles.gaussFitBgCorrMethodDrop(ch).Value;
        end
        if isfield(handles, 'gaussFitBgCorrWidthEdit') && length(handles.gaussFitBgCorrWidthEdit) >= ch
            lastUsed.gaussFitBgCorrWidth{ch} = handles.gaussFitBgCorrWidthEdit(ch).Value;
        end
        if isfield(handles, 'gaussFitPolyDegreeEdit') && length(handles.gaussFitPolyDegreeEdit) >= ch
            lastUsed.gaussFitPolyDegree{ch} = handles.gaussFitPolyDegreeEdit(ch).Value;
        end
        if isfield(handles, 'gaussFitMaxIterationsEdit') && length(handles.gaussFitMaxIterationsEdit) >= ch
            lastUsed.gaussFitMaxIterations{ch} = handles.gaussFitMaxIterationsEdit(ch).Value;
        end
        if isfield(handles, 'gaussFitToleranceEdit') && length(handles.gaussFitToleranceEdit) >= ch
            lastUsed.gaussFitTolerance{ch} = handles.gaussFitToleranceEdit(ch).Value;
        end
        if isfield(handles, 'gaussFitRadialRadiusEdit') && length(handles.gaussFitRadialRadiusEdit) >= ch
            lastUsed.gaussFitRadialRadius{ch} = handles.gaussFitRadialRadiusEdit(ch).Value;
        end
        if isfield(handles, 'gaussFitPlotCheck') && length(handles.gaussFitPlotCheck) >= ch && isvalid(handles.gaussFitPlotCheck(ch))
            lastUsed.gaussFitPlotCheck{ch} = handles.gaussFitPlotCheck(ch).Value;
        end
        
        % Fit Filtering
        if isfield(handles, 'fitFilterEnabledChecks') && length(handles.fitFilterEnabledChecks) >= ch
            lastUsed.fitFilterEnabled{ch} = handles.fitFilterEnabledChecks(ch).Value;
        end
        if isfield(handles, 'fitFilterRSquaredEnabledChecks') && length(handles.fitFilterRSquaredEnabledChecks) >= ch
            lastUsed.fitFilterRSquaredEnabled{ch} = handles.fitFilterRSquaredEnabledChecks(ch).Value;
        end
        if isfield(handles, 'fitFilterRSquaredMinInputs') && length(handles.fitFilterRSquaredMinInputs) >= ch
            lastUsed.fitFilterRSquaredMin{ch} = handles.fitFilterRSquaredMinInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterRSquaredMaxInputs') && length(handles.fitFilterRSquaredMaxInputs) >= ch
            lastUsed.fitFilterRSquaredMax{ch} = handles.fitFilterRSquaredMaxInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterSigmaSumEnabledChecks') && length(handles.fitFilterSigmaSumEnabledChecks) >= ch
            lastUsed.fitFilterSigmaSumEnabled{ch} = handles.fitFilterSigmaSumEnabledChecks(ch).Value;
        end
        if isfield(handles, 'fitFilterSigmaSumMinInputs') && length(handles.fitFilterSigmaSumMinInputs) >= ch
            lastUsed.fitFilterSigmaSumMin{ch} = handles.fitFilterSigmaSumMinInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterSigmaSumMaxInputs') && length(handles.fitFilterSigmaSumMaxInputs) >= ch
            lastUsed.fitFilterSigmaSumMax{ch} = handles.fitFilterSigmaSumMaxInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterAmplitudeEnabledChecks') && length(handles.fitFilterAmplitudeEnabledChecks) >= ch
            lastUsed.fitFilterAmplitudeEnabled{ch} = handles.fitFilterAmplitudeEnabledChecks(ch).Value;
        end
        if isfield(handles, 'fitFilterAmplitudeMinInputs') && length(handles.fitFilterAmplitudeMinInputs) >= ch
            lastUsed.fitFilterAmplitudeMin{ch} = handles.fitFilterAmplitudeMinInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterAmplitudeMaxInputs') && length(handles.fitFilterAmplitudeMaxInputs) >= ch
            lastUsed.fitFilterAmplitudeMax{ch} = handles.fitFilterAmplitudeMaxInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterIntensityEnabledChecks') && length(handles.fitFilterIntensityEnabledChecks) >= ch
            lastUsed.fitFilterIntensityEnabled{ch} = handles.fitFilterIntensityEnabledChecks(ch).Value;
        end
        if isfield(handles, 'fitFilterIntensityMinInputs') && length(handles.fitFilterIntensityMinInputs) >= ch
            lastUsed.fitFilterIntensityMin{ch} = handles.fitFilterIntensityMinInputs(ch).Value;
        end
        if isfield(handles, 'fitFilterIntensityMaxInputs') && length(handles.fitFilterIntensityMaxInputs) >= ch
            lastUsed.fitFilterIntensityMax{ch} = handles.fitFilterIntensityMaxInputs(ch).Value;
        end
        catch ME
            % Skip this channel if UI elements are invalid (e.g., GraphicsPlaceholder)
            % This can happen for channels that weren't fully initialized
            % fprintf('Warning: Could not export parameters for channel %d: %s\n', ch, ME.message);
        end
    end
    
    % Preview settings
    if isfield(handles, 'previewContentDrops') && length(handles.previewContentDrops) >= 5
        for i = 1:5
            lastUsed.previewContents{i} = handles.previewContentDrops(i).Value;
        end
    end
    if isfield(handles, 'previewModeDrops') && length(handles.previewModeDrops) >= 5
        for i = 1:5
            lastUsed.previewModes{i} = handles.previewModeDrops(i).Value;
        end
    end
    if isfield(handles, 'previewProjectionDrops') && length(handles.previewProjectionDrops) >= 5
        for i = 1:5
            lastUsed.previewProjections{i} = handles.previewProjectionDrops(i).Value;
        end
    end
    
    % Export settings
    if isfield(handles, 'exportImageFormatDrop')
        lastUsed.exportImageFormat = handles.exportImageFormatDrop.Value;
    end
    
catch ME
    warning('Error collecting parameters: %s', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).file, ME.stack(i).line);
    end
    % Return a basic structure instead of failing completely
end

% Note: We deliberately exclude some file paths and panel states from export
% These should not be exported with parameters
% - dicFilePath, nucFilePath, channelFilePaths are excluded (except PSF file paths)
% - nucNavPanelIndex, navPanelIndex are excluded
% Preview settings (previewContents, previewModes, previewProjections) ARE exported

end
