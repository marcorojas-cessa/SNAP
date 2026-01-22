function exportNucleiSignalDataStandardized(outputPath, imageName, nucleusLabels, nucleiMask, signalComposition, options, precomputedMeasurements)
% exportNucleiSignalDataStandardized - Combined nuclei-signal export
%
% ============================================================================
% EXPORT DESIGN PHILOSOPHY (Critical for Understanding)
% ============================================================================
%
% ZERO REDUNDANCY PRINCIPLE:
%   This function NEVER recomputes data that was already calculated elsewhere.
%   Instead, it assembles pre-computed data from shared sources:
%
%   1. NUCLEUS DATA: From computeNucleusMeasurements()
%      - If precomputedMeasurements provided → use it (CONSISTENCY MODE)
%      - If not provided → compute now (BACKWARD COMPATIBILITY)
%      - Result: Nucleus data IDENTICAL to exportNucleiDataStandardized()
%
%   2. SIGNAL DATA: From signalComposition.nuclei_data
%      - Already contains FULL signal data (via computeSignalMeasurements())
%      - This was computed during analyzeNucleiSignalComposition()
%      - Result: Signal data IDENTICAL to exportChannelDataStandardized()
%
% TWO CSV FORMATS:
%
%   SIMPLE CSV (one row per nucleus):
%   - Purpose: Quick overview of nuclei properties + signal counts
%   - Contains: ALL nucleus measurements + per-channel signal counts
%   - Same nucleus fields as exportNucleiDataStandardized()
%
%   EXPANDED CSV (one row per signal):
%   - Purpose: Detailed analysis linking each signal to its parent nucleus
%   - Contains: COMPLETE nucleus data + COMPLETE signal data for each signal
%   - Each signal row includes ALL parent nucleus measurements
%   - Allows filtering/grouping by nucleus properties in downstream analysis
%
% DATA CONSISTENCY:
%   - computation_id links this export to nuclei-only export
%   - Same ID = guaranteed identical nucleus measurements
%   - ID stored in MAT metadata and receipt files
%
% ============================================================================
%
% INPUTS:
%   outputPath        - Base path (without extension)
%   imageName         - Image identifier
%   nucleusLabels     - From segmentation
%   nucleiMask        - Binary/labeled mask
%   signalComposition - From analyzeNucleiSignalComposition (contains full signal data)
%   options           - Struct with is_3d, xy_spacing, z_spacing
%   precomputedMeasurements - (OPTIONAL) Pre-computed nuclei measurements from computeNucleusMeasurements()
%                             If provided, uses this instead of recomputing (ensures consistency)
%                             This is the KEY parameter for consistency!
%
% OUTPUTS:
%   {outputPath}.mat - Complete nested structure with metadata.computation_id
%   {outputPath}_simple.csv - One row per nucleus with COMPLETE nucleus data + signal counts
%   {outputPath}_expanded.csv - One row per signal with COMPLETE nucleus + signal data
%   {outputPath}_receipt.txt - Statistics receipt with computation_id for verification

% Silent export (no console output)

% Check if we have data
if isempty(signalComposition) || ~isfield(signalComposition, 'nuclei_data') || isempty(signalComposition.nuclei_data)
    warning('No signal composition data to export');
    return;
end

% Parse options
if ~isfield(options, 'is_3d'), options.is_3d = false; end
if ~isfield(options, 'xy_spacing'), options.xy_spacing = []; end
if ~isfield(options, 'z_spacing'), options.z_spacing = []; end

has_spacing = ~isempty(options.xy_spacing) && options.xy_spacing > 0;
if options.is_3d && (~isempty(options.z_spacing) && options.z_spacing > 0)
    has_spacing = has_spacing;
else
    has_spacing = false;
end

%% GET NUCLEUS MEASUREMENTS (use pre-computed if available for consistency)
%  ========================================================================
%  DECISION POINT: Use shared measurements or compute now?
%  ========================================================================
%
%  IF precomputedMeasurements provided (consistency mode):
%    → Use it directly (NEVER recompute)
%    → Contains same computation_id as nuclei-only export
%    → GUARANTEES identical nucleus data between exports
%
%  IF NOT provided (standalone mode):
%    → Call computeNucleusMeasurements() now
%    → Generates new computation_id
%    → Works but can't guarantee consistency with nuclei-only export
%
%  WHY THIS MATTERS:
%    If both nuclei and nuclei+signal exports are generated, users expect
%    nucleus measurements (volume, area, etc.) to be IDENTICAL in both.
%    This decision point ensures that happens.
%
%  ========================================================================

if nargin >= 7 && ~isempty(precomputedMeasurements)
    % CONSISTENCY MODE: Use pre-computed measurements
    nucleiData = precomputedMeasurements;
    fprintf('    Using pre-computed measurements (ensures consistency)\n');
else
    % STANDALONE MODE: Compute now (backward compatibility)
nucleiData = snap_helpers.computeNucleusMeasurements(nucleusLabels, nucleiMask, options);
end

%% Add data consistency metadata
% Store computation ID from nucleiData to verify consistency
if isfield(nucleiData, 'computation_id')
    computation_id = nucleiData.computation_id;
else
    % Generate new ID if not present
    computation_id = sprintf('nucleisignal_%s', datestr(now, 'yyyymmdd_HHMMSS_FFF'));
    nucleiData.computation_id = computation_id;
end

%% Build combined structure for MATLAB export
matData = struct();
matData.metadata.description = 'Combined nuclei and signal data with complete associations';
matData.metadata.coordinate_convention = 'ARRAY_CONVENTION';
matData.metadata.coordinate_format = '[row, col, slice]';
matData.metadata.image_name = imageName;
matData.metadata.num_nuclei = nucleiData.num_nuclei;
matData.metadata.is_3d = options.is_3d;
matData.metadata.computation_id = computation_id;  % For consistency verification

% Get number of channels
num_channels = length(signalComposition.active_channels);
matData.metadata.num_channels = num_channels;
matData.metadata.active_channels = signalComposition.active_channels;

% Build nuclei array with embedded signals
% Pre-allocate struct array to avoid dissimilar structures error
total_signals = 0;

% Create first nucleus to establish structure
nuclei_array = [];

for i = 1:nucleiData.num_nuclei
    % CRITICAL: Create ALL fields in the SAME ORDER for every nucleus
    % This prevents "dissimilar structures" error
    nuc = struct();
    
    % Nucleus ID and centroid
    nuc.nucleus_id = i;
    nuc.centroid_x = nucleiData.centroids(i, 1);  % col
    nuc.centroid_y = nucleiData.centroids(i, 2);  % row
    nuc.centroid_z = nucleiData.centroids(i, 3);  % slice
    
    % ALL measurement fields - ALWAYS in this exact order
    % Populate with actual data or NaN based on what's available
    if nucleiData.is_3d && isfield(nucleiData, 'volume_voxels')
        nuc.volume_voxels = nucleiData.volume_voxels(i);
        nuc.area_pixels = NaN;
    else
        nuc.volume_voxels = NaN;
        nuc.area_pixels = nucleiData.area_pixels(i);
    end
    
    if nucleiData.is_3d && isfield(nucleiData, 'surface_area_voxels')
        nuc.surface_area_voxels = nucleiData.surface_area_voxels(i);
        nuc.perimeter_pixels = NaN;
    else
        nuc.surface_area_voxels = NaN;
        nuc.perimeter_pixels = nucleiData.perimeter_pixels(i);
    end
    
    nuc.equivalent_diameter_pixels = nucleiData.equivalent_diameter_pixels(i);
    nuc.solidity = nucleiData.solidity(i);
    
    if nucleiData.is_3d && isfield(nucleiData, 'sphericity')
        nuc.sphericity = nucleiData.sphericity(i);
        nuc.circularity = NaN;
    else
        nuc.sphericity = NaN;
        nuc.circularity = nucleiData.circularity(i);
    end
    
    % Spacing fields - ALWAYS in this exact order
    if has_spacing
        nuc.equivalent_diameter_microns = nucleiData.equivalent_diameter_microns(i);
        if nucleiData.is_3d
            nuc.volume_microns3 = nucleiData.volume_microns3(i);
            nuc.area_microns2 = NaN;
            nuc.surface_area_microns2 = nucleiData.surface_area_microns2(i);
            nuc.perimeter_microns = NaN;
        else
            nuc.volume_microns3 = NaN;
            nuc.area_microns2 = nucleiData.area_microns2(i);
            nuc.surface_area_microns2 = NaN;
            nuc.perimeter_microns = nucleiData.perimeter_microns(i);
        end
    else
        nuc.equivalent_diameter_microns = NaN;
        nuc.volume_microns3 = NaN;
        nuc.area_microns2 = NaN;
        nuc.surface_area_microns2 = NaN;
        nuc.perimeter_microns = NaN;
    end
    
    % Extract signal data from signalComposition (already has full data!)
    nuc.signals = cell(1, num_channels);
    nuc.signal_counts = zeros(1, num_channels);
    
    if i <= length(signalComposition.nuclei_data)
        comp_nuc = signalComposition.nuclei_data(i);
        
        for ch_idx = 1:num_channels
            ch = signalComposition.active_channels(ch_idx);
            ch_field = sprintf('channel_%d', ch);
            
            if isfield(comp_nuc.channels, ch_field)
                ch_data = comp_nuc.channels.(ch_field);
                
                if isfield(ch_data, 'signals') && ~isempty(ch_data.signals)
                    % Signals already have full data from getSignalsInNucleus!
                    nuc.signals{ch_idx} = ch_data.signals;
                    nuc.signal_counts(ch_idx) = length(ch_data.signals);
                    total_signals = total_signals + length(ch_data.signals);
                else
                    nuc.signals{ch_idx} = [];
                    nuc.signal_counts(ch_idx) = 0;
                end
            else
                nuc.signals{ch_idx} = [];
                nuc.signal_counts(ch_idx) = 0;
            end
        end
    end
    
    nuc.total_signal_count = sum(nuc.signal_counts);
    
    % Append to array
    if i == 1
        nuclei_array = nuc;  % First element establishes structure
    else
        nuclei_array(i) = nuc;  % Subsequent elements must match
    end
end

% Store in matData
matData.nuclei = nuclei_array;

%% Save MATLAB file
matFilename = [outputPath '.mat'];
save(matFilename, '-struct', 'matData');
% MAT file saved silently

%% Export SIMPLE CSV (one row per nucleus) - COMPLETE NUCLEUS DATA
csvSimpleFilename = [outputPath '_simple.csv'];
fid = fopen(csvSimpleFilename, 'w');

% Use SHARED helper to build nucleus column list (ZERO REDUNDANCY!)
% This ensures columns EXACTLY match exportNucleiDataStandardized()
[nucleus_cols, include_flags] = snap_helpers.buildNucleusColumnList(nucleiData, has_spacing, '');

% Add image_name column ONLY if image name provided (for batch exports)
include_image_name = ~isempty(imageName);
if include_image_name
    header_cols = {'image_name'};
    header_cols = [header_cols, nucleus_cols];
else
    header_cols = nucleus_cols;
end

% Add per-channel signal counts
for ch_idx = 1:num_channels
    ch = signalComposition.active_channels(ch_idx);
    header_cols{end+1} = sprintf('ch%d_signal_count', ch);
end

fprintf(fid, '%s\n', strjoin(header_cols, ','));

% Data rows - Use SHARED logic for writing nucleus data
for i = 1:nucleiData.num_nuclei
    nuc = nuclei_array(i);
    
    % Use SHARED helper to write nucleus data (ZERO REDUNDANCY!)
    % This ensures data order EXACTLY matches the header columns
    row_data = snap_helpers.writeNucleusDataRow(nuc, nucleiData, has_spacing, include_flags);
    
    % Prepend image_name ONLY if included in header
    if include_image_name
        row_data = [{imageName}, row_data];
    end
    
    % Append signal counts per channel
    for ch_idx = 1:num_channels
        row_data{end+1} = sprintf('%d', nuc.signal_counts(ch_idx));
    end
    
    fprintf(fid, '%s\n', strjoin(row_data, ','));
end

fclose(fid);
% Simple CSV saved silently

%% Export EXPANDED CSV (one row per signal) - COMPLETE NUCLEUS + SIGNAL DATA
csvExpandedFilename = [outputPath '_expanded.csv'];
fid = fopen(csvExpandedFilename, 'w');

% Build header using SHARED helpers (ZERO REDUNDANCY!)
% This ensures columns EXACTLY match individual exports

% Nucleus columns with 'nucleus_' prefix (uses SAME logic as exportNucleiDataStandardized)
[nucleus_cols_prefixed, nuc_include_flags] = snap_helpers.buildNucleusColumnList(nucleiData, has_spacing, 'nucleus_');

% Determine fit method from signal data to build signal columns correctly
fitMethod = '';
hasFitting = false;
if ~isempty(signalComposition.nuclei_data)
    % Try to get fit method from first signal found
    for nuc_idx = 1:length(signalComposition.nuclei_data)
        for ch_idx = 1:num_channels
            ch = signalComposition.active_channels(ch_idx);
            ch_field = sprintf('channel_%d', ch);
            if isfield(signalComposition.nuclei_data(nuc_idx).channels, ch_field)
                ch_data = signalComposition.nuclei_data(nuc_idx).channels.(ch_field);
                if isfield(ch_data, 'signals') && ~isempty(ch_data.signals) && isstruct(ch_data.signals)
                    % Get fit method from first signal
                    sig = ch_data.signals(1);
                    if isfield(sig, 'fitted_coords') && any(~isnan(sig.fitted_coords))
                        hasFitting = true;
                        % Infer method from which amplitude fields are present
                        if ~isnan(sig.amplitude_x)
                            fitMethod = '1D (X,Y,Z) Gaussian';
                        elseif ~isnan(sig.amplitude_xy)
                            fitMethod = '2D (XY) + 1D (Z) Gaussian';
                        elseif ~isnan(sig.radial_symmetry_score)
                            fitMethod = 'Radial Symmetry';
                        else
                            fitMethod = 'Gaussian';
                        end
                        break;
                    end
                end
            end
        end
        if hasFitting, break; end
    end
end

% Signal columns using SHARED helper (uses SAME logic as exportChannelDataStandardized)
% NOTE: Nuclei+signal export NEEDS channel_id because it combines multiple channels
% Pass true to include channel_id column
[signal_cols, sig_include_flags] = snap_helpers.buildSignalColumnList(hasFitting, fitMethod, nucleiData.is_3d, true);

% Combine: image_name (if provided) + nucleus columns + signal columns
if include_image_name
    header_cols = {'image_name'};
    header_cols = [header_cols, nucleus_cols_prefixed, signal_cols];
else
    header_cols = [nucleus_cols_prefixed, signal_cols];
end

fprintf(fid, '%s\n', strjoin(header_cols, ','));

% Data rows (one per signal) - Use SHARED logic for both nucleus and signal data
for i = 1:nucleiData.num_nuclei
    nuc = nuclei_array(i);
    
    % Use SHARED helper to write nucleus data (ZERO REDUNDANCY!)
    nuc_data = snap_helpers.writeNucleusDataRow(nuc, nucleiData, has_spacing, nuc_include_flags);
    
    % Prepend image_name ONLY if included in header
    if include_image_name
        nuc_data = [{imageName}, nuc_data];
    end
    
    % Write a row for each signal in this nucleus
    for ch_idx = 1:num_channels
        signals_in_channel = nuc.signals{ch_idx};
        if ~isempty(signals_in_channel) && isstruct(signals_in_channel)
            ch = signalComposition.active_channels(ch_idx);
            
            num_signals = length(signals_in_channel);
            for sig_idx = 1:num_signals
                sig = signals_in_channel(sig_idx);
                
                % Start with nucleus data
                row_data = nuc_data;
                
                % Add signal data using SHARED helper (ZERO REDUNDANCY!)
                sig.channel_id = ch;  % Add channel ID to signal struct for helper
                signal_data = snap_helpers.writeSignalDataRow(sig, sig_include_flags);
                row_data = [row_data, signal_data];
                
                fprintf(fid, '%s\n', strjoin(row_data, ','));
            end
        end
    end
end

fclose(fid);
% Expanded CSV saved silently

%% Export TXT receipt
txtFilename = [outputPath '_receipt.txt'];
fid = fopen(txtFilename, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'SNAP Nuclei-Signal Composition Analysis\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Image: %s\n\n', imageName);

fprintf(fid, '--- Summary Statistics ---\n');
fprintf(fid, 'Total Nuclei: %d\n', nucleiData.num_nuclei);
fprintf(fid, 'Total Signals: %d\n', total_signals);
fprintf(fid, 'Channels: %d\n', num_channels);
fprintf(fid, 'Active Channels: %s\n\n', mat2str(signalComposition.active_channels));

fprintf(fid, '--- Export Format ---\n');
fprintf(fid, 'Simple CSV: One row per nucleus with COMPLETE nucleus data + signal counts\n');
fprintf(fid, 'Expanded CSV: One row per signal with COMPLETE nucleus + signal data\n\n');

fprintf(fid, '--- Nucleus Fields Exported ---\n');
fprintf(fid, 'Centroids: centroid_x, centroid_y, centroid_z\n');
if nucleiData.is_3d
    fprintf(fid, 'Size: volume_voxels');
    if has_spacing, fprintf(fid, ', volume_microns3'); end
    fprintf(fid, '\n');
    fprintf(fid, 'Surface: surface_area_voxels');
    if has_spacing, fprintf(fid, ', surface_area_microns2'); end
    fprintf(fid, '\n');
    fprintf(fid, 'Shape: sphericity, solidity\n');
else
    fprintf(fid, 'Size: area_pixels');
    if has_spacing, fprintf(fid, ', area_microns2'); end
    fprintf(fid, '\n');
    fprintf(fid, 'Perimeter: perimeter_pixels');
    if has_spacing, fprintf(fid, ', perimeter_microns'); end
    fprintf(fid, '\n');
    fprintf(fid, 'Shape: circularity, solidity\n');
end
fprintf(fid, 'Diameter: equivalent_diameter_pixels');
if has_spacing, fprintf(fid, ', equivalent_diameter_microns'); end
fprintf(fid, '\n\n');

fprintf(fid, '--- Signal Fields Exported (Expanded CSV) ---\n');
fprintf(fid, 'Coordinates: maxima_x/y/z, fitted_x/y/z\n');
fprintf(fid, 'Intensity: amplitude(s), integrated_intensity, background\n');
fprintf(fid, 'Quality: r_squared or radial_symmetry_score\n');
fprintf(fid, 'Shape: sigma_x/y/z, rho_xy/xz/yz (if distorted), alpha_x/y/z (if skewed)\n\n');

fprintf(fid, '--- Coordinate Convention ---\n');
fprintf(fid, 'Convention: ARRAY (row, col, slice)\n');
fprintf(fid, 'Access: imageData(row, col, slice)\n\n');

fprintf(fid, '--- Data Consistency ---\n');
fprintf(fid, 'Computation ID: %s\n', computation_id);
fprintf(fid, 'Note: All exports with the same Computation ID are guaranteed to have\n');
fprintf(fid, '      identical nucleus measurements (computed once and shared).\n\n');

fprintf(fid, '--- Export Date & Time ---\n');
fprintf(fid, '%s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fclose(fid);
    % Receipt saved silently
end

function str = formatValue(val)
    if isnan(val) || isempty(val)
        str = 'NaN';
    else
        str = sprintf('%.4f', val);
    end
end

