function exportNucleiDataStandardized(nucleusLabels, nucleiMask, outputPath, options, precomputedMeasurements)
% Export nuclei segmentation data in standardized format
% Used by both SNAP and SNAP_batch to ensure identical export format
%
% ============================================================================
% CONSISTENCY GUARANTEE (Critical for Understanding)
% ============================================================================
%
% This function is part of a SHARED MEASUREMENT SYSTEM that ensures nucleus
% data is IDENTICAL across all export types (nuclei, channel, nuclei+signal).
%
% HOW IT WORKS:
%   - When called WITH precomputedMeasurements:
%     → Uses the provided data (CONSISTENCY MODE)
%     → NEVER recomputes (prevents any possibility of discrepancy)
%     → Includes computation_id for verification
%
%   - When called WITHOUT precomputedMeasurements:
%     → Calls computeNucleusMeasurements() itself (STANDALONE MODE)
%     → Generates new computation_id
%     → Backward compatible with old code
%
% CALLING PATTERN:
%   From exportData() or SNAP_batch:
%     1. Compute measurements ONCE before export loop
%     2. Pass to ALL nuclei-related exports
%     3. Result: Perfect consistency across all files
%
% ============================================================================
%
% INPUTS:
%   nucleusLabels - Struct with fields:
%       .num_nuclei
%       .centroids_3d - [N x 3] array (Note: in Cartesian [col, row, slice])
%       .labeled_mask - Indexed mask where each nucleus has unique ID
%       .areas (2D) or .volumes (3D) - Optional
%       .circularities (2D) or .sphericities (3D) - Optional
%   nucleiMask    - Binary or labeled mask
%   outputPath    - Base path for output (without extension)
%   options       - Struct with fields:
%       .is_3d          - Boolean (3D vs 2D segmentation)
%       .xy_spacing     - Physical XY spacing (optional, microns/pixel)
%       .z_spacing      - Physical Z spacing (optional, microns/slice)
%       .image_name     - Image identifier (for batch, optional)
%   precomputedMeasurements - (OPTIONAL) Pre-computed nuclei measurements from computeNucleusMeasurements()
%                             If provided, uses this instead of recomputing (ensures consistency)
%                             This is the KEY to the consistency system!
%
% OUTPUTS:
%   Creates three files:
%   - {outputPath}.csv  - CSV with all measurements
%   - {outputPath}.mat  - MAT with data + labeled_mask + computation_id
%   - {outputPath}_receipt.txt - Parameter receipt with computation_id

    % Silent export (no console output)
    
    % Validate inputs
    if isempty(nucleusLabels) || nucleusLabels.num_nuclei == 0
        warning('No nuclei to export');
        return;
    end
    
    % Check options
    if ~isfield(options, 'is_3d'), options.is_3d = false; end
    if ~isfield(options, 'xy_spacing'), options.xy_spacing = []; end
    if ~isfield(options, 'z_spacing'), options.z_spacing = []; end
    if ~isfield(options, 'image_name'), options.image_name = ''; end
    
    has_spacing = ~isempty(options.xy_spacing) && options.xy_spacing > 0;
    if options.is_3d
        has_spacing = has_spacing && ~isempty(options.z_spacing) && options.z_spacing > 0;
    end
    
    %% GET NUCLEUS MEASUREMENTS (use pre-computed if available for consistency)
    %  ========================================================================
    %  DECISION POINT: Use shared measurements or compute now?
    %  ========================================================================
    %
    %  IF precomputedMeasurements provided (consistency mode):
    %    → Use it directly (NEVER recompute)
    %    → Contains computation_id from original computation
    %    → GUARANTEES identical data with other exports
    %
    %  IF NOT provided (standalone mode):
    %    → Call computeNucleusMeasurements() now
    %    → Generates new computation_id
    %    → Works but can't guarantee consistency with other exports
    %
    %  ========================================================================
    
    if nargin >= 5 && ~isempty(precomputedMeasurements)
        % CONSISTENCY MODE: Use pre-computed measurements
        nucleiData = precomputedMeasurements;
        fprintf('    Using pre-computed measurements (ensures consistency)\n');
    else
        % STANDALONE MODE: Compute now (backward compatibility)
    nucleiData = snap_helpers.computeNucleusMeasurements(nucleusLabels, nucleiMask, options);
    end
    
    % Extract variables for CSV export
    num_nuclei = nucleiData.num_nuclei;
    centroids = nucleiData.centroids;  % [col, row, slice] Cartesian
    labeled_mask = nucleiData.labeled_mask;
    
    % Arrays for CSV export
    if ~isempty(options.image_name)
        image_names = repmat({options.image_name}, num_nuclei, 1);
    else
        image_names = [];
    end
    
    nucleus_ids = (1:num_nuclei)';
    centroid_x = centroids(:, 1);  % col
    centroid_y = centroids(:, 2);  % row
    centroid_z = centroids(:, 3);  % slice
    
    %% Add data consistency metadata
    % Generate unique computation ID for this measurement set
    % This allows verification that all exports came from the same computation
    if ~isfield(nucleiData, 'computation_id') || isempty(nucleiData.computation_id)
        nucleiData.computation_id = sprintf('nuclei_%s', datestr(now, 'yyyymmdd_HHMMSS_FFF'));
    end
    
    %% Export MAT and CSV files
    matFilename = [outputPath '.mat'];
    csvFilename = [outputPath '.csv'];
    
    % Build table with controlled column order using SHARED COLUMN LOGIC
    % This ensures columns match exactly what other exports would generate
    exportTable = table();
    
    % Add image name ONLY if provided (for batch exports)
    % SNAP exports don't need image_name column (single image context)
    include_image_name = ~isempty(options.image_name);
    if include_image_name
        exportTable.image_name = image_names;
    end
    
    % Use SHARED helper to determine which columns to include (ZERO REDUNDANCY!)
    [~, include_flags] = snap_helpers.buildNucleusColumnList(nucleiData, has_spacing, '');
    
    % Core fields (always)
    exportTable.nucleus_id = nucleus_ids;
    exportTable.centroid_x = centroid_x;
    exportTable.centroid_y = centroid_y;
    
    % Centroid Z (conditional based on mask dimensionality)
    if include_flags.centroid_z
        exportTable.centroid_z = centroid_z;
    end
    
    % Size metrics (conditional based on dimensionality)
    if include_flags.is_3d
        % 3D fields
        exportTable.volume_voxels = nucleiData.volume_voxels;
        if has_spacing
            exportTable.volume_microns3 = nucleiData.volume_microns3;
        end
        exportTable.surface_area_voxels = nucleiData.surface_area_voxels;
        if has_spacing
            exportTable.surface_area_microns2 = nucleiData.surface_area_microns2;
        end
    else
        % 2D fields
        exportTable.area_pixels = nucleiData.area_pixels;
        if has_spacing
            exportTable.area_microns2 = nucleiData.area_microns2;
        end
        exportTable.perimeter_pixels = nucleiData.perimeter_pixels;
        if has_spacing
            exportTable.perimeter_microns = nucleiData.perimeter_microns;
        end
        end
        
    % Common measurements (always)
        exportTable.equivalent_diameter_pixels = nucleiData.equivalent_diameter_pixels;
        if has_spacing
            exportTable.equivalent_diameter_microns = nucleiData.equivalent_diameter_microns;
        end
        
    % Shape metrics (dimension-specific)
    if include_flags.is_3d
        exportTable.sphericity = nucleiData.sphericity;
        exportTable.solidity = nucleiData.solidity;
    else
        exportTable.circularity = nucleiData.circularity;
        exportTable.solidity = nucleiData.solidity;
    end
    
    % Save files
    save(matFilename, 'nucleiData');
    writetable(exportTable, csvFilename);
    % Files saved silently
    
    %% Create parameter receipt (TXT)
    txtFilename = [outputPath '_receipt.txt'];
    fid = fopen(txtFilename, 'w');
    
    fprintf(fid, '========================================\n');
    fprintf(fid, 'SNAP Nuclei Segmentation Parameters\n');
    fprintf(fid, '========================================\n\n');
    
    if ~isempty(options.image_name)
        fprintf(fid, 'Image: %s\n\n', options.image_name);
    end
    
    fprintf(fid, '--- Segmentation Type ---\n');
    if options.is_3d
        fprintf(fid, 'Type: 3D Volume Segmentation\n');
        fprintf(fid, 'Metrics: Volume, Surface Area, Sphericity\n\n');
    else
        fprintf(fid, 'Type: 2D Segmentation\n');
        fprintf(fid, 'Metrics: Area, Perimeter, Circularity\n\n');
    end
    
    fprintf(fid, '--- Physical Spacing ---\n');
    if has_spacing
        fprintf(fid, 'XY Spacing: %.4f microns/pixel\n', options.xy_spacing);
        if options.is_3d && ~isempty(options.z_spacing)
            fprintf(fid, 'Z Spacing: %.4f microns/slice\n', options.z_spacing);
        end
    else
        fprintf(fid, 'Physical spacing not provided\n');
        fprintf(fid, 'All measurements in pixels/voxels only\n');
    end
    fprintf(fid, '\n');
    
    fprintf(fid, '--- Exported Fields ---\n');
    fprintf(fid, 'Mandatory: image_name, nucleus_id, centroid_x, centroid_y, centroid_z\n');
    if options.is_3d
        fprintf(fid, 'Size: volume_voxels');
        if has_spacing
            fprintf(fid, ', volume_microns3');
        end
        fprintf(fid, '\n');
        fprintf(fid, 'Surface: surface_area_voxels');
        if has_spacing
            fprintf(fid, ', surface_area_microns2');
        end
        fprintf(fid, '\n');
    else
        fprintf(fid, 'Size: area_pixels');
        if has_spacing
            fprintf(fid, ', area_microns2');
        end
        fprintf(fid, '\n');
        fprintf(fid, 'Perimeter: perimeter_pixels');
        if has_spacing
            fprintf(fid, ', perimeter_microns');
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, 'Equivalent Diameter: equivalent_diameter_pixels');
    if has_spacing
        fprintf(fid, ', equivalent_diameter_microns');
    end
    fprintf(fid, '\n');
    
    if options.is_3d
        fprintf(fid, 'Shape: sphericity, solidity\n\n');
    else
        fprintf(fid, 'Shape: circularity, solidity\n\n');
    end
    
    fprintf(fid, '--- Coordinate Convention ---\n');
    fprintf(fid, 'Convention: ARRAY (row, col, slice)\n');
    fprintf(fid, 'Centroid coordinates: [x=col, y=row, z=slice] Cartesian\n');
    fprintf(fid, 'Indexed mask: imageData(row, col, slice)\n\n');
    
    fprintf(fid, '--- Export Statistics ---\n');
    fprintf(fid, 'Total Nuclei: %d\n', num_nuclei);
    if options.is_3d
        fprintf(fid, 'Mean Volume: %.2f voxels\n', mean(nucleiData.volume_voxels));
        if has_spacing
            fprintf(fid, 'Mean Volume: %.2f microns³\n', mean(nucleiData.volume_microns3));
        end
    else
        fprintf(fid, 'Mean Area: %.2f pixels\n', mean(nucleiData.area_pixels));
        if has_spacing
            fprintf(fid, 'Mean Area: %.2f microns²\n', mean(nucleiData.area_microns2));
        end
    end
    fprintf(fid, '\n');
    
    fprintf(fid, '--- Data Consistency ---\n');
    fprintf(fid, 'Computation ID: %s\n', nucleiData.computation_id);
    fprintf(fid, 'Note: All exports with the same Computation ID are guaranteed to have\n');
    fprintf(fid, '      identical measurements (computed once and shared).\n\n');
    
    fprintf(fid, '--- Export Date & Time ---\n');
    fprintf(fid, '%s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    
    fclose(fid);
    % Receipt saved silently
end

