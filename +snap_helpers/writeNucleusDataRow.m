function row_data = writeNucleusDataRow(nuc, nucleiData, has_spacing, include_flags)
% writeNucleusDataRow - Build nucleus data row for CSV using shared logic
%
% ============================================================================
% SINGLE SOURCE OF TRUTH for nucleus data row writing
% ============================================================================
%
% This function writes nucleus measurements to CSV in the EXACT order
% and with the EXACT conditional logic as buildNucleusColumnList().
%
% CRITICAL: This function MUST write fields in the same order and with
%           the same conditionals as buildNucleusColumnList() builds headers.
%
% Used by:
%   1. exportNucleiSignalDataStandardized (simple CSV)
%   2. exportNucleiSignalDataStandardized (expanded CSV)
%
% CONSISTENCY GUARANTEE:
%   - include_flags comes from buildNucleusColumnList()
%   - This function follows the EXACT same conditionals
%   - Result: Data columns perfectly match header columns
%
%   If this function and buildNucleusColumnList() ever get out of sync,
%   you'll get CSV parsing errors or misaligned columns. That's why they
%   must be updated together.
%
% ============================================================================
%
% INPUTS:
%   nuc          - Nucleus struct with all measurement fields
%   nucleiData   - Full nucleiData struct (for accessing labeled_mask)
%   has_spacing  - Boolean indicating if physical spacing available
%   include_flags - Flags from buildNucleusColumnList() indicating which fields to include
%
% OUTPUTS:
%   row_data - Cell array of formatted values matching header order

row_data = {};

% Core fields (always included)
row_data{end+1} = sprintf('%d', nuc.nucleus_id);
row_data{end+1} = sprintf('%.4f', nuc.centroid_x);
row_data{end+1} = sprintf('%.4f', nuc.centroid_y);

% Centroid Z (conditional)
if include_flags.centroid_z
    row_data{end+1} = sprintf('%.4f', nuc.centroid_z);
end

% Size measurements (2D vs 3D) - ONLY write applicable fields
if include_flags.is_3d
    % 3D fields
    row_data{end+1} = sprintf('%.2f', nuc.volume_voxels);
    if has_spacing
        row_data{end+1} = sprintf('%.4f', nuc.volume_microns3);
    end
    row_data{end+1} = sprintf('%.2f', nuc.surface_area_voxels);
    if has_spacing
        row_data{end+1} = sprintf('%.4f', nuc.surface_area_microns2);
    end
else
    % 2D fields
    row_data{end+1} = sprintf('%.2f', nuc.area_pixels);
    if has_spacing
        row_data{end+1} = sprintf('%.4f', nuc.area_microns2);
    end
    row_data{end+1} = sprintf('%.2f', nuc.perimeter_pixels);
    if has_spacing
        row_data{end+1} = sprintf('%.4f', nuc.perimeter_microns);
    end
end

% Common measurements (always included)
row_data{end+1} = sprintf('%.2f', nuc.equivalent_diameter_pixels);
if has_spacing
    row_data{end+1} = sprintf('%.4f', nuc.equivalent_diameter_microns);
end

% Shape metrics (dimension-specific)
if include_flags.is_3d
    row_data{end+1} = sprintf('%.4f', nuc.sphericity);
    row_data{end+1} = sprintf('%.4f', nuc.solidity);
else
    row_data{end+1} = sprintf('%.4f', nuc.circularity);
    row_data{end+1} = sprintf('%.4f', nuc.solidity);
end

end

