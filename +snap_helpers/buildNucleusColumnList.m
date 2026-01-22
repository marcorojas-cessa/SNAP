function [header_cols, include_flags] = buildNucleusColumnList(nucleiData, has_spacing, prefix)
% buildNucleusColumnList - Determine which nucleus columns to include in CSV
%
% ============================================================================
% SINGLE SOURCE OF TRUTH for nucleus column selection
% ============================================================================
%
% This function encapsulates the EXACT logic for deciding which nucleus
% columns to include in CSV exports. Used by:
%   1. exportNucleiDataStandardized (no prefix)
%   2. exportNucleiSignalDataStandardized simple CSV (no prefix)
%   3. exportNucleiSignalDataStandardized expanded CSV (with 'nucleus_' prefix)
%
% WHY THIS EXISTS:
%   Ensures all three exports have IDENTICAL column selection logic.
%   Without this, we'd have to maintain the same conditional logic in
%   three places, risking inconsistencies.
%
%   PROBLEM THIS SOLVES:
%     Before: Each export had its own column selection code
%     - exportNucleiData: volume_voxels column present for 3D
%     - exportNucleiSignal: volume_voxels AND area_pixels columns present (with NaNs)
%     - User sees different column sets (CONFUSING!)
%
%     Now: All exports call this function
%     - Same input → same columns
%     - exportNucleiData: only applicable columns
%     - exportNucleiSignal: only applicable columns (MATCHES!)
%     - User sees consistent column sets (CLEAR!)
%
% CONDITIONAL LOGIC:
%   - Always include: nucleus_id, centroid_x, centroid_y
%   - Include centroid_z ONLY IF mask is 3D
%   - Include volume fields ONLY IF is_3d=true
%   - Include area fields ONLY IF is_3d=false
%   - Include micron fields ONLY IF has_spacing=true
%   - Include sphericity ONLY IF is_3d=true
%   - Include circularity ONLY IF is_3d=false
%
% ============================================================================
%
% INPUTS:
%   nucleiData  - Struct from computeNucleusMeasurements with .is_3d, .labeled_mask
%   has_spacing - Boolean indicating if physical spacing is available
%   prefix      - String prefix for column names (e.g., '' or 'nucleus_')
%
% OUTPUTS:
%   header_cols   - Cell array of column names to include
%   include_flags - Struct with boolean flags for optional columns
%                   (useful for data writing logic)

% Add prefix to column name helper
addPrefix = @(name) [prefix name];

% Initialize output
header_cols = {};
include_flags = struct();

% Always include these core fields
header_cols = {addPrefix('nucleus_id'), addPrefix('centroid_x'), addPrefix('centroid_y')};

% Include centroid_z if labeled_mask is 3D (covers 3D volume and 2D slice-by-slice modes)
if ndims(nucleiData.labeled_mask) == 3
    header_cols = [header_cols, {addPrefix('centroid_z')}];
    include_flags.centroid_z = true;
else
    include_flags.centroid_z = false;
end

% Size measurements (2D vs 3D) - ONLY include applicable fields
if nucleiData.is_3d
    % 3D fields only
    header_cols = [header_cols, {addPrefix('volume_voxels')}];
    if has_spacing
        header_cols = [header_cols, {addPrefix('volume_microns3')}];
    end
    header_cols = [header_cols, {addPrefix('surface_area_voxels')}];
    if has_spacing
        header_cols = [header_cols, {addPrefix('surface_area_microns2')}];
    end
    include_flags.is_3d = true;
else
    % 2D fields only
    header_cols = [header_cols, {addPrefix('area_pixels')}];
    if has_spacing
        header_cols = [header_cols, {addPrefix('area_microns2')}];
    end
    header_cols = [header_cols, {addPrefix('perimeter_pixels')}];
    if has_spacing
        header_cols = [header_cols, {addPrefix('perimeter_microns')}];
    end
    include_flags.is_3d = false;
end

% Common measurements (always include)
header_cols = [header_cols, {addPrefix('equivalent_diameter_pixels')}];
if has_spacing
    header_cols = [header_cols, {addPrefix('equivalent_diameter_microns')}];
end

% Shape metrics (dimension-specific)
if nucleiData.is_3d
    header_cols = [header_cols, {addPrefix('sphericity'), addPrefix('solidity')}];
else
    header_cols = [header_cols, {addPrefix('circularity'), addPrefix('solidity')}];
end

% Store flags for data writing
include_flags.has_spacing = has_spacing;
end

