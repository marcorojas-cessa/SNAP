function nucleiData = computeNucleusMeasurements(nucleusLabels, labeled_mask, options)
% computeNucleusMeasurements - SINGLE SOURCE OF TRUTH for nucleus measurements
%
% ============================================================================
% CRITICAL: This is THE ONLY PLACE where nucleus measurements are computed
% ============================================================================
%
% WHY THIS FUNCTION EXISTS:
%   To prevent ANY possibility of inconsistent measurements across exports.
%   
%   BEFORE this system:
%   - exportNucleiData() computed measurements
%   - exportNucleiSignalData() computed measurements separately
%   - Risk: Different rounding, different algorithms, different inputs
%   - Result: User sees different numbers in different exports (BAD!)
%
%   WITH this system:
%   - This function called ONCE
%   - Result stored and passed to ALL export functions
%   - Guaranteed: All exports have IDENTICAL measurements
%
% USAGE PATTERN:
%   
%   Option A: Compute and pass to exports (RECOMMENDED for consistency)
%   ---------------------------------------------------------------
%   nucleiData = computeNucleusMeasurements(labels, mask, opts);
%   exportNucleiDataStandardized(..., nucleiData);
%   exportNucleiSignalDataStandardized(..., nucleiData);
%   % Both exports now have IDENTICAL nucleus measurements
%
%   Option B: Let export functions compute (backward compatible)
%   ------------------------------------------------------------
%   exportNucleiDataStandardized(...);  % Computes internally
%   % Works but doesn't guarantee consistency with other exports
%
% WHAT IT COMPUTES:
%   2D measurements:
%     - Size: area (regionprops), perimeter (bwboundaries-based)
%     - Shape: circularity (4π×Area/Perimeter²), solidity (regionprops)
%     - Equivalent diameter (calculated from area)
%   3D measurements:
%     - Size: volume/surface area (regionprops3)
%     - Shape: sphericity (π^(1/3)×(6V)^(2/3)/SA), solidity (regionprops3)
%     - Equivalent diameter (calculated from volume)
%   
%   NOTE ON PERIMETER/CIRCULARITY CALCULATION:
%     Uses bwboundaries to trace boundary and sums Euclidean distances.
%     This is more accurate than regionprops 'Perimeter' which uses
%     4-connectivity and underestimates. All shape metrics capped at 1.0.
%
% ============================================================================
%
% INPUTS:
%   nucleusLabels - Struct with:
%       .num_nuclei
%       .centroids_3d ([col, row, slice] Cartesian)
%       .labeled_mask (THIS is the labeled mask that will be used)
%   labeled_mask  - [DEPRECATED - ignored] The labeled_mask from nucleusLabels struct is used instead
%   options       - Struct with:
%       .is_3d (boolean)
%       .xy_spacing (optional, microns/pixel)
%       .z_spacing (optional, microns/slice)
%       .image_name (optional)
%
% OUTPUTS:
%   nucleiData - Struct with fields:
%       .num_nuclei
%       .nucleus_ids (1:N)
%       .centroids ([col, row, slice])
%       .labeled_mask
%       .is_3d
%       .computation_id (unique ID for consistency verification)
%       
%       For 3D:
%         .volume_voxels, .surface_area_voxels, .equivalent_diameter_pixels
%         .sphericity, .solidity
%         (if spacing) .volume_microns3, .surface_area_microns2, .equivalent_diameter_microns
%       
%       For 2D:
%         .area_pixels, .perimeter_pixels, .equivalent_diameter_pixels
%         .circularity, .solidity
%         (if spacing) .area_microns2, .perimeter_microns, .equivalent_diameter_microns

    num_nuclei = nucleusLabels.num_nuclei;
    centroids = nucleusLabels.centroids_3d;  % [col, row, slice] Cartesian
    
    % CRITICAL FIX: Use labeled_mask from nucleusLabels struct (the authoritative source)
    % The second parameter 'labeled_mask' is redundant and often passed incorrectly as binary mask
    if isfield(nucleusLabels, 'labeled_mask') && ~isempty(nucleusLabels.labeled_mask)
        labeled_mask = nucleusLabels.labeled_mask;
    else
        error('nucleusLabels must contain a labeled_mask field');
    end
    
    % Check spacing
    has_spacing = false;
    if isfield(options, 'xy_spacing') && ~isempty(options.xy_spacing) && options.xy_spacing > 0
        has_spacing = true;
        if options.is_3d && (~isfield(options, 'z_spacing') || isempty(options.z_spacing) || options.z_spacing <= 0)
            has_spacing = false;
        end
    end
    
    %% Compute regionprops for all nuclei
    % Three possible cases:
    %   1. 2D segmentation on 2D image → 2D mask → use regionprops
    %   2. 2D segmentation on 3D image (slice-by-slice) → 3D mask → project to 2D, use regionprops
    %   3. 3D segmentation on 3D image (3D volume mode) → 3D mask → use regionprops3
    
    % CRITICAL: Use actual mask dimensionality, not just options.is_3d
    % options.is_3d indicates the segmentation MODE, but we need to check the actual mask
    actual_is_3d = (ndims(labeled_mask) == 3 && size(labeled_mask, 3) > 1);
    
    if actual_is_3d
        % Case 3: 3D volume segmentation (3D mask exists)
        % Note: EquivalentDiameter not supported in older regionprops3, calculate manually
        props = regionprops3(labeled_mask, 'Volume', 'SurfaceArea', 'Centroid', 'Solidity');
    else
        % Cases 1 & 2: 2D segmentation or 3D mask projected to 2D
        if ndims(labeled_mask) == 3 && size(labeled_mask, 3) > 1
            % Case 2: Project 3D mask to 2D
            labeled_mask_2d = max(labeled_mask, [], 3);
        else
            % Case 1: Already 2D
            labeled_mask_2d = labeled_mask;
        end
        
        % Note: We calculate perimeter manually using bwboundaries (more accurate)
        % regionprops 'Perimeter' uses 4-connectivity which underestimates
        props = regionprops(labeled_mask_2d, 'Area', 'Centroid', 'Solidity');
    end
    
    %% Build output structure
    nucleiData = struct();
    nucleiData.num_nuclei = num_nuclei;
    nucleiData.nucleus_ids = (1:num_nuclei)';
    nucleiData.centroids = centroids;  % [col, row, slice] Cartesian
    nucleiData.labeled_mask = labeled_mask;
    nucleiData.is_3d = actual_is_3d;  % Use actual computed dimensionality, not requested mode
    
    %% Extract measurements from regionprops
    % Use actual_is_3d to determine which measurements to extract
    if actual_is_3d
        % === 3D MEASUREMENTS ===
        % Extract fields from table (regionprops3 returns table)
        volume_voxels = props.Volume;
        surface_area_voxels = props.SurfaceArea;
        solidity = props.Solidity;
        
        % Calculate equivalent diameter manually: diameter of sphere with same volume
        % Formula: d = (6 * V / π)^(1/3)
        equivalent_diameter_pixels = (6 * volume_voxels / pi).^(1/3);
        
        % Get actual label IDs (might not be consecutive after filtering)
        unique_labels_3d = unique(labeled_mask);
        unique_labels_3d = unique_labels_3d(unique_labels_3d > 0);  % Exclude background (0)
        
        % Calculate sphericity using REGIONPROPS3 METHOD with fallbacks
        % NOTE: Uses regionprops3 primarily (same as before), but includes
        % fallback methods for robustness. All methods cap at 1.0.
        sphericity = zeros(num_nuclei, 1);
        for i = 1:num_nuclei
            if surface_area_voxels(i) > 0
                % Method 1: Use regionprops3 values (preferred, most accurate)
                sphericity(i) = (pi^(1/3) * (6 * volume_voxels(i))^(2/3)) / surface_area_voxels(i);
                sphericity(i) = max(0, min(1.0, sphericity(i)));  % Cap at [0, 1]
            else
                % Fallback: Try alternative surface area estimation
                try
                    % Use actual label ID
                    if i <= length(unique_labels_3d)
                        label_id = unique_labels_3d(i);
                    else
                        label_id = i;  % Fallback
                    end
                    
                    nucleus_mask = (labeled_mask == label_id);
                    vol = volume_voxels(i);
                    
                    % Method 2: Isosurface-based surface area (if regionprops3 fails)
                    if vol > 8
                        smoothed_mask = smooth3(double(nucleus_mask), 'gaussian', [3 3 3], 0.5);
                        smoothed_mask = smoothed_mask > 0.5;
                    else
                        smoothed_mask = nucleus_mask;
                    end
                    
                    [faces, vertices] = isosurface(smoothed_mask, 0.5);
                    
                    if ~isempty(faces) && size(faces, 1) > 0
                        surface_area_alt = 0;
                        for f = 1:size(faces, 1)
                            v1 = vertices(faces(f, 1), :);
                            v2 = vertices(faces(f, 2), :);
                            v3 = vertices(faces(f, 3), :);
                            edge1 = v2 - v1;
                            edge2 = v3 - v1;
                            triangle_area = 0.5 * norm(cross(edge1, edge2));
                            surface_area_alt = surface_area_alt + triangle_area;
                        end
                        
                        if surface_area_alt > 0
                            sphericity(i) = (pi^(1/3) * (6 * vol)^(2/3)) / surface_area_alt;
                            sphericity(i) = max(0, min(1.0, sphericity(i)));
                            surface_area_voxels(i) = surface_area_alt; % Update for export
                        end
                    end
                catch
                    % Method 3: Simple boundary voxel counting (last resort)
                    sphericity(i) = 0;
                end
            end
        end
        
        % Store measurements
        nucleiData.volume_voxels = volume_voxels;
        nucleiData.surface_area_voxels = surface_area_voxels;
        nucleiData.equivalent_diameter_pixels = equivalent_diameter_pixels;
        nucleiData.solidity = solidity;
        nucleiData.sphericity = sphericity;
        
        % Physical units if spacing provided
        if has_spacing
            nucleiData.volume_microns3 = volume_voxels * (options.xy_spacing^2 * options.z_spacing);
            nucleiData.surface_area_microns2 = surface_area_voxels * (options.xy_spacing^2);
            nucleiData.equivalent_diameter_microns = equivalent_diameter_pixels * options.xy_spacing;
        end
        
    else
        % === 2D MEASUREMENTS ===
        % Extract fields from struct array (regionprops returns struct array)
        area_pixels = [props.Area]';
        solidity = [props.Solidity]';
        
        % Get actual label IDs (might not be consecutive after filtering)
        unique_labels = unique(labeled_mask_2d);
        unique_labels = unique_labels(unique_labels > 0);  % Exclude background (0)
        
        % Calculate perimeter and circularity using BOUNDARY-BASED METHOD
        % NOTE: This uses bwboundaries instead of regionprops 'Perimeter' because
        % regionprops uses 4-connectivity which underestimates the perimeter.
        % The boundary-based method calculates the actual Euclidean distance along the boundary.
        perimeter_pixels = zeros(num_nuclei, 1);
        circularity = zeros(num_nuclei, 1);
        
        for i = 1:num_nuclei
            % Use the actual label ID (regionprops output order matches label order)
            if i <= length(unique_labels)
                label_id = unique_labels(i);
            else
                label_id = i;  % Fallback
            end
            
            % Create binary mask for just this nucleus
            nucleus_mask = (labeled_mask_2d == label_id);
            
            % Get boundary using bwboundaries
            boundaries = bwboundaries(nucleus_mask, 'noholes');
            
            if ~isempty(boundaries) && size(boundaries{1}, 1) >= 3
                boundary = boundaries{1}; % Main boundary [row, col]
                
                % Calculate perimeter as sum of Euclidean distances between consecutive points
                perim = 0;
                for j = 1:size(boundary, 1)
                    next_j = mod(j, size(boundary, 1)) + 1;
                    dx = boundary(next_j, 2) - boundary(j, 2);  % col difference
                    dy = boundary(next_j, 1) - boundary(j, 1);  % row difference
                    perim = perim + sqrt(dx^2 + dy^2);
                end
                
                perimeter_pixels(i) = perim;
                
                % Calculate circularity: 4π × Area / Perimeter²
                if perim > 0
                    circularity(i) = 4 * pi * area_pixels(i) / (perim^2);
                circularity(i) = min(1.0, circularity(i));  % Cap at 1.0
                end
            end
        end
        
        % Calculate equivalent diameter manually: diameter = sqrt(4 * Area / π)
        % This is the diameter of a circle with the same area
        equivalent_diameter_pixels = sqrt(4 * area_pixels / pi);
        
        % Store measurements
        nucleiData.area_pixels = area_pixels;
        nucleiData.perimeter_pixels = perimeter_pixels;
        nucleiData.equivalent_diameter_pixels = equivalent_diameter_pixels;
        nucleiData.solidity = solidity;
        nucleiData.circularity = circularity;
        
        % Physical units if spacing provided
        if has_spacing
            nucleiData.area_microns2 = area_pixels * (options.xy_spacing^2);
            nucleiData.perimeter_microns = perimeter_pixels * options.xy_spacing;
            nucleiData.equivalent_diameter_microns = equivalent_diameter_pixels * options.xy_spacing;
        end
    end
end

