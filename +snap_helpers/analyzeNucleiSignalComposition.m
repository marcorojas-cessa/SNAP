function analysis_results = analyzeNucleiSignalComposition(handles)
    % Analyze signal composition within each segmented nucleus
    % Returns a structure with quantitative metrics for all nuclei and channels
%
% ============================================================================
% SIGNAL COMPOSITION ANALYSIS (Part of Export Consistency System)
% ============================================================================
%
% This function creates the DATA SOURCE for nuclei+signal export.
% It links signals to their parent nuclei and embeds FULL signal data.
%
% KEY PRINCIPLE - ZERO REDUNDANCY:
%   Instead of storing just signal counts or IDs, this function embeds
%   the COMPLETE signal data (from computeSignalMeasurements()) within
%   each nucleus's structure.
%
% DATA FLOW:
%   1. Get nucleus labels (from cache or compute)
%   2. Get signal data per channel (from cache: handles.maximaCoords, handles.gaussFitResults)
%   3. For each nucleus:
%      → Find signals within its boundaries
%      → Call computeSignalMeasurements() to get FULL signal data
%      → Embed complete signal data in nucleus structure
%   4. Return analysis_results with:
%      - .nuclei_data(i).channels.channel_X.signals (FULL data per signal)
%
% USAGE IN EXPORTS:
%   exportNucleiSignalDataStandardized() receives this structure and:
%   - Extracts signal data WITHOUT recomputing anything
%   - Combines with nucleus measurements (from computeNucleusMeasurements)
%   - Result: Complete data with zero redundant computation
%
% CONSISTENCY WITH CHANNEL EXPORT:
%   The signal data embedded here came from computeSignalMeasurements(),
%   the SAME function that exportChannelDataStandardized() uses.
%   Therefore, signal measurements are IDENTICAL in both exports.
%
% ============================================================================
    
    analysis_results = [];
    
    % Check if we have actual segmented nuclei objects to analyze
    % This is the definitive test - if nuclei were segmented, we'll have labeled objects
    nucleus_labels = [];
    nuc_mask = [];
    
    % Try to get nuclei labels from cache first
    if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ...
       ~isempty(handles.previewCache.nuclei) && isfield(handles.previewCache.nuclei, 'labels')
        nucleus_labels = handles.previewCache.nuclei.labels;
        nuc_mask = handles.previewCache.nuclei.mask;
    else
    end
    
    % If no cached labels, check if we should generate them
    if isempty(nucleus_labels) && isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
        try
            processed_nuc = snap_helpers.preprocessNucleiWithBgCorr(handles);
            [nuc_mask, nucleus_labels] = snap_helpers.segmentNuclei(processed_nuc, handles);
        catch
            % If segmentation fails, no analysis possible
            return;
        end
    end
    
    % If we still don't have nuclei labels or they're empty, no analysis possible
    if isempty(nucleus_labels) || ~isfield(nucleus_labels, 'centroids_3d') || ...
       isempty(nucleus_labels.centroids_3d) || nucleus_labels.num_nuclei == 0
        return;
    end
    
    try
        % Use the segmented nuclei objects for analysis
        centroids_3d = nucleus_labels.centroids_3d;
        cc = nucleus_labels.connected_components;
        
        % Double-check we have actual nuclei objects
        if isempty(centroids_3d) || isempty(cc) || cc.NumObjects == 0
            return;
        end
        
        % Get active channels - check BOTH handles.rawChannel AND previewCache
        numActiveChannels = str2double(handles.numChanDrop.Value);
        active_channels = [];
        for k = 1:numActiveChannels
            % Check if channel has data in rawChannel OR in preview cache
            has_raw = (k <= numel(handles.rawChannel) && ~isempty(handles.rawChannel{k}));
            has_cached = (isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
                         k <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{k}));
            
            if has_raw || has_cached
                active_channels(end+1) = k;
            end
        end
        
        if isempty(active_channels)
            return;
        end
        
        % COHESION VALIDATION: Check data consistency between preview and analysis
        
        % Initialize results structure
        analysis_results = struct();
        analysis_results.num_nuclei = size(centroids_3d, 1);
        analysis_results.active_channels = active_channels;
        analysis_results.nuclei_data = [];
        analysis_results.centroids_3d = centroids_3d; % Store centroids for consistency
        
        % Get voxel spacing for volume calculations
        xy_spacing = handles.nucXYSpacingInput.Value;
        z_spacing = handles.nucZSpacingInput.Value;
        if xy_spacing <= 0, xy_spacing = 1; end
        if z_spacing <= 0, z_spacing = 1; end
        
        % Analyze each nucleus using consistent numbering
        nuclei_data = [];
        for nuc_idx = 1:size(centroids_3d, 1)
            nucleus_data = struct();
            nucleus_data.nucleus_id = nuc_idx;  % This now matches the preview labels
            
            % Use centroid from the consistent labeling system
            nucleus_data.centroid = centroids_3d(nuc_idx, :);
            
            % Get pixel indices for this nucleus by finding the connected component
            % that contains the centroid coordinates
            pixel_indices = [];
            centroid_coords = round(nucleus_data.centroid);
            
            % Find which connected component contains this centroid
            for cc_idx = 1:cc.NumObjects
                test_indices = cc.PixelIdxList{cc_idx};
                if ndims(nuc_mask) == 3
                    [y_coords, x_coords, z_coords] = ind2sub(size(nuc_mask), test_indices);
                    if any(x_coords == centroid_coords(1) & y_coords == centroid_coords(2) & z_coords == centroid_coords(3))
                        pixel_indices = test_indices;
                        break;
                    end
                else
                    [y_coords, x_coords] = ind2sub(size(nuc_mask), test_indices);
                    if any(x_coords == centroid_coords(1) & y_coords == centroid_coords(2))
                        pixel_indices = test_indices;
                        break;
                    end
                end
            end
            
            % If we couldn't match by exact centroid, use the cc_idx as fallback
            if isempty(pixel_indices) && nuc_idx <= cc.NumObjects
                pixel_indices = cc.PixelIdxList{nuc_idx};
            end
            
            if isempty(pixel_indices)
                continue; % Skip this nucleus if we can't find its pixels
            end
            
            
            % Calculate nucleus volume/area (ensure integers)
            if ndims(nuc_mask) == 3
                nucleus_data.volume_voxels = round(length(pixel_indices));  % Ensure integer
                nucleus_data.volume_microns3 = nucleus_data.volume_voxels * (xy_spacing^2 * z_spacing);
            else
                nucleus_data.area_pixels = round(length(pixel_indices));  % Ensure integer
                nucleus_data.area_microns2 = nucleus_data.area_pixels * (xy_spacing^2);
            end
            
            % Calculate circularity/sphericity based on segmentation mode
            % CRITICAL FIX: Check segmentation mode, not just mask dimensions
            seg_mode = handles.nucSegModeDrop.Value;
            
            % DEEP DEBUG: Investigate the actual data being passed
            
            % ================================================================
            % CRITICAL CONSISTENCY FIX (Do NOT change without understanding!)
            % ================================================================
            %
            % SHAPE METRICS ARE NOT COMPUTED HERE
            %
            % WHY:
            %   This function previously had its own calculateCircularity2D()
            %   and calculateSphericity3D() functions that used custom perimeter
            %   calculations (via bwboundaries). This created DISCREPANCIES
            %   with exports which use computeNucleusMeasurements() (uses
            %   regionprops perimeter).
            %
            %   Example discrepancy:
            %   - Analysis table: 0.987 (from custom calculation)
            %   - Export CSV: 0.992 (from regionprops)
            %   - User sees different values (BAD!)
            %
            % SOLUTION:
            %   analyzeNucleiSignalComposition() focuses ONLY on signal analysis
            %   (signal counts, intensities, distributions). Shape metrics are
            %   computed by computeNucleusMeasurements() which is the single
            %   source of truth for ALL exports and the analysis table.
            %
            % RESULT:
            %   - Analysis table: gets shape metrics from computeNucleusMeasurements()
            %   - All exports: use computeNucleusMeasurements()
            %   - Perfect consistency: SAME values everywhere!
            %
            % ================================================================
            
            % Store placeholders (analysis table gets real values from shared source)
                        nucleus_data.circularity = NaN;
                        nucleus_data.sphericity = NaN;
            nucleus_data.solidity = NaN;
            
            % Analyze signal in each channel
            channel_data = struct();
            
            for ch = active_channels
                % Get channel image - check both raw data and cache
                channel_img = [];
                if ch <= numel(handles.rawChannel) && ~isempty(handles.rawChannel{ch})
                    channel_img = double(handles.rawChannel{ch});
                elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
                       ch <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{ch}) && ...
                       isfield(handles.previewCache.channels{ch}, 'processed') && ~isempty(handles.previewCache.channels{ch}.processed)
                    channel_img = double(handles.previewCache.channels{ch}.processed);
                end
                
                if isempty(channel_img)
                    continue;  % Skip this channel if no data available
                end
                
                % Extract intensities within this nucleus
                % CRITICAL FIX: Handle dimension mismatch properly
                
                if ndims(channel_img) == ndims(nuc_mask)
                    % Same dimensions - direct indexing
                    intensities = channel_img(pixel_indices);
                    
                elseif ndims(channel_img) == 2 && ndims(nuc_mask) == 3
                    % 2D channel, 3D nucleus mask - project nucleus mask to 2D
                    nuc_mask_2d = max(nuc_mask, [], 3) > 0;
                    [y_2d, x_2d] = find(nuc_mask_2d);
                    nucleus_pixels_2d = sub2ind(size(channel_img), y_2d, x_2d);
                    intensities = channel_img(nucleus_pixels_2d);
                    
                elseif ndims(channel_img) == 3 && ndims(nuc_mask) == 2
                    % 3D channel, 2D nucleus mask - this is your case!
                    % We need to extract intensities from all Z-slices where the nucleus exists
                    
                    % Create a 2D mask for just this specific nucleus
                    this_nucleus_mask_2d = false(size(nuc_mask));
                    this_nucleus_mask_2d(pixel_indices) = true;
                    
                    % Get the 2D coordinates of this nucleus
                    [y_coords, x_coords] = find(this_nucleus_mask_2d);
                    
                    % Extract intensities from all Z-slices at these X,Y coordinates
                    intensities = [];
                    for z = 1:size(channel_img, 3)
                        for i = 1:length(y_coords)
                            intensities(end+1) = channel_img(y_coords(i), x_coords(i), z);
                        end
                    end
                    
                    % Keep original pixel_indices for signal counting (they're correct for 2D)
                    
                else
                    continue; % Skip this channel if dimensions don't match
                end
                
                % Calculate statistics
                ch_stats = struct();
                ch_stats.mean_intensity = mean(intensities);
                ch_stats.std_intensity = std(intensities);
                ch_stats.max_intensity = max(intensities);
                ch_stats.min_intensity = min(intensities);
                ch_stats.total_intensity = sum(intensities);
                ch_stats.median_intensity = median(intensities);
                
                % Get signals (local maxima) within this nucleus with FULL fit data
                % USE SHARED HELPER FOR ZERO REDUNDANCY!
                [signal_indices, signal_data_array] = getSignalsInNucleus(handles, ch, pixel_indices, nuc_mask);
                
                % Store signal data (ensure consistent structure)
                if ~isempty(signal_data_array) && isstruct(signal_data_array)
                    ch_stats.signal_count = length(signal_data_array);
                    ch_stats.signals = signal_data_array;
                else
                    ch_stats.signal_count = 0;
                    ch_stats.signals = struct([]);  % Empty struct array
                end
                
                % Calculate background-corrected metrics if possible
                % Use pixels just outside the nucleus as background estimate
                try
                    if ndims(nuc_mask) == 3
                        dilated_mask = imdilate(nuc_mask, ones(3,3,3));
                    else
                        dilated_mask = imdilate(nuc_mask, ones(3,3));
                    end
                    background_mask = dilated_mask & ~nuc_mask;
                    
                    if sum(background_mask(:)) > 0
                        if ndims(channel_img) == ndims(background_mask)
                            bg_intensities = channel_img(background_mask);
                            background_mean = mean(bg_intensities);
                            ch_stats.background_mean = background_mean;
                            ch_stats.signal_to_background = ch_stats.mean_intensity / max(background_mean, 1);
                            ch_stats.background_corrected_total = ch_stats.total_intensity - (background_mean * length(intensities));
                        end
                    end
                catch
                    % If background calculation fails, continue without it
                end
                
                channel_data.(sprintf('channel_%d', ch)) = ch_stats;
            end
            
            nucleus_data.channels = channel_data;
            
            % Calculate channel ratios if multiple channels available
            if length(active_channels) > 1
                ratios = struct();
                for i = 1:length(active_channels)
                    for j = i+1:length(active_channels)
                        ch1 = active_channels(i);
                        ch2 = active_channels(j);
                        
                        ch1_field = sprintf('channel_%d', ch1);
                        ch2_field = sprintf('channel_%d', ch2);
                        
                        if isfield(channel_data, ch1_field) && isfield(channel_data, ch2_field)
                            mean1 = channel_data.(ch1_field).mean_intensity;
                            mean2 = channel_data.(ch2_field).mean_intensity;
                            
                            if mean2 > 0
                                ratio_name = sprintf('ch%d_ch%d_ratio', ch1, ch2);
                                ratios.(ratio_name) = mean1 / mean2;
                            end
                        end
                    end
                end
                nucleus_data.ratios = ratios;
            end
            
            nuclei_data = [nuclei_data; nucleus_data];
        end
        
        analysis_results.nuclei_data = nuclei_data;
        analysis_results.timestamp = now();
        
    catch ME
        warning('Failed to analyze nuclei signal composition: %s', ME.message);
        analysis_results = [];
    end
end

function circularity = calculateCircularity2D(pixel_indices, mask)
    % Calculate circularity for 2D nucleus: 4*pi*Area/Perimeter^2
    % Perfect circle = 1.0, less circular shapes < 1.0
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        nucleus_mask(pixel_indices) = true;
        
        % Use regionprops for accurate measurements
        props = regionprops(nucleus_mask, 'Area', 'Perimeter');
        
        if isempty(props) || props.Area == 0 || props.Perimeter == 0
            circularity = 0;
            return;
        end
        
        % CRITICAL FIX: regionprops 'Perimeter' uses 4-connectivity which underestimates
        % Use boundary-based perimeter calculation for more accurate results
        boundaries = bwboundaries(nucleus_mask, 'noholes');
        if isempty(boundaries)
            circularity = 0;
            return;
        end
        
        % Calculate actual perimeter from boundary coordinates
        boundary = boundaries{1}; % Main boundary
        if size(boundary, 1) < 3
            circularity = 0;
            return;
        end
        
        % Calculate perimeter as sum of distances between consecutive boundary points
        perimeter = 0;
        for i = 1:size(boundary, 1)
            next_i = mod(i, size(boundary, 1)) + 1;
            dx = boundary(next_i, 2) - boundary(i, 2);
            dy = boundary(next_i, 1) - boundary(i, 1);
            perimeter = perimeter + sqrt(dx^2 + dy^2);
        end
        
        area = props.Area;
        
        % Circularity formula: 4*pi*Area/Perimeter^2
        circularity = 4 * pi * area / (perimeter^2);
        circularity = min(1.0, circularity); % Cap at 1.0
        
    catch ME
        circularity = 0;
    end
end

function circularity = calculateCircularity2D_debug(pixel_indices, mask, nuc_id)
    % DEBUG VERSION: Calculate circularity with extensive debugging
    fprintf('  [Nucleus %d] Starting circularity calculation...\n', nuc_id);
    
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        fprintf('  [Nucleus %d] Created empty mask of size %s\n', nuc_id, mat2str(size(nucleus_mask)));
        
        % Check pixel indices validity
        if any(pixel_indices > numel(nucleus_mask)) || any(pixel_indices < 1)
            fprintf('  [Nucleus %d] ERROR: Invalid pixel indices!\n', nuc_id);
            circularity = 0;
            return;
        end
        
        nucleus_mask(pixel_indices) = true;
        fprintf('  [Nucleus %d] Set %d pixels to true\n', nuc_id, length(pixel_indices));
        
        % Check if we actually have a connected shape
        cc = bwconncomp(nucleus_mask);
        fprintf('  [Nucleus %d] Found %d connected components\n', nuc_id, cc.NumObjects);
        
        if cc.NumObjects == 0
            fprintf('  [Nucleus %d] ERROR: No connected components found!\n', nuc_id);
            circularity = 0;
            return;
        end
        
        % Use regionprops for accurate measurements
        props = regionprops(nucleus_mask, 'Area', 'Perimeter', 'BoundingBox', 'Centroid');
        fprintf('  [Nucleus %d] regionprops returned %d objects\n', nuc_id, length(props));
        
        if isempty(props)
            fprintf('  [Nucleus %d] ERROR: regionprops returned empty!\n', nuc_id);
            circularity = 0;
            return;
        end
        
        % Use the largest object if multiple exist
        if length(props) > 1
            areas = [props.Area];
            [~, max_idx] = max(areas);
            props = props(max_idx);
            fprintf('  [Nucleus %d] Using largest object (index %d) out of %d\n', nuc_id, max_idx, length(areas));
        end
        
        % CRITICAL FIX: regionprops 'Perimeter' uses 4-connectivity which underestimates
        % Use boundary-based perimeter calculation for more accurate results
        boundaries = bwboundaries(nucleus_mask, 'noholes');
        if isempty(boundaries)
            fprintf('  [Nucleus %d] ERROR: No boundaries found!\n', nuc_id);
            circularity = 0;
            return;
        end
        
        % Calculate actual perimeter from boundary coordinates
        boundary = boundaries{1}; % Main boundary
        if size(boundary, 1) < 3
            fprintf('  [Nucleus %d] ERROR: Boundary too small (%d points)!\n', nuc_id, size(boundary, 1));
            circularity = 0;
            return;
        end
        
        % Calculate perimeter using more accurate method
        % Method 1: Try bwperim + sum for pixel-accurate perimeter
        perim_mask = bwperim(nucleus_mask);
        pixel_perimeter = sum(perim_mask(:));
        
        % Method 2: Boundary-based calculation with proper edge weighting
        boundary_perimeter = 0;
        for i = 1:size(boundary, 1)
            next_i = mod(i, size(boundary, 1)) + 1;
            dx = boundary(next_i, 2) - boundary(i, 2);
            dy = boundary(next_i, 1) - boundary(i, 1);
            dist = sqrt(dx^2 + dy^2);
            boundary_perimeter = boundary_perimeter + dist;
        end
        
        % Use the larger of the two estimates (more conservative)
        perimeter = max(pixel_perimeter, boundary_perimeter);
        
        area = props.Area;
        regionprops_perimeter = props.Perimeter;
        
        fprintf('  [Nucleus %d] Raw measurements - Area: %.2f, RegionProps Perimeter: %.2f, Pixel Perimeter: %.2f, Boundary Perimeter: %.2f, Final Perimeter: %.2f\n', ...
            nuc_id, area, regionprops_perimeter, pixel_perimeter, boundary_perimeter, perimeter);
        
        if area == 0 || perimeter == 0
            fprintf('  [Nucleus %d] ERROR: Zero area (%.2f) or perimeter (%.2f)!\n', nuc_id, area, perimeter);
            circularity = 0;
            return;
        end
        
        % Check for suspicious values
        if area == 1 && perimeter < 4
            fprintf('  [Nucleus %d] WARNING: Single pixel detected (Area=1, Perimeter=%.2f)\n', nuc_id, perimeter);
        end
        
        % Circularity formula: 4*pi*Area/Perimeter^2
        raw_circularity = 4 * pi * area / (perimeter^2);
        raw_circularity = min(1.0, raw_circularity); % Cap at 1.0
        fprintf('  [Nucleus %d] Circularity: %.6f\n', nuc_id, raw_circularity);
        
        % Check for perfect circle (suspicious)
        if raw_circularity > 0.99
            fprintf('  [Nucleus %d] WARNING: Near-perfect circularity detected! This may indicate an error.\n', nuc_id);
            fprintf('  [Nucleus %d] Centroid: [%.2f, %.2f], BoundingBox: [%.1f %.1f %.1f %.1f]\n', ...
                nuc_id, props.Centroid(1), props.Centroid(2), props.BoundingBox(1), props.BoundingBox(2), props.BoundingBox(3), props.BoundingBox(4));
        end
        
        % Clamp to maximum of 1.0 (perfect circle) since values > 1.0 are mathematically impossible
        circularity = min(1.0, raw_circularity);
        
        if raw_circularity > 1.0
            fprintf('  [Nucleus %d] NOTE: Clamped circularity from %.6f to 1.0000 (measurement precision limit)\n', nuc_id, raw_circularity);
        end
        
        fprintf('  [Nucleus %d] Final circularity: %.3f\n', nuc_id, circularity);
        
    catch ME
        fprintf('  [Nucleus %d] ERROR in circularity calculation: %s\n', nuc_id, ME.message);
        circularity = 0;
    end
end

function sphericity = calculateSphericity3D(pixel_indices, mask)
    % Calculate sphericity for 3D nucleus: pi^(1/3) * (6*Volume)^(2/3) / SurfaceArea
    % Perfect sphere = 1.0, less spherical shapes < 1.0
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        nucleus_mask(pixel_indices) = true;
        
        % Calculate volume (number of voxels)
        volume = length(pixel_indices);
        
        if volume == 0
            sphericity = 0;
            return;
        end
        
        % Method 1: Use regionprops3 if available (MATLAB R2017b+)
        if exist('regionprops3', 'file')
            try
                props = regionprops3(nucleus_mask, 'Volume', 'SurfaceArea');
                if ~isempty(props) && props.SurfaceArea > 0
                    volume_measured = props.Volume;
                    surface_area = props.SurfaceArea;
                    
                    % Calculate sphericity
                    sphericity = (pi^(1/3)) * ((6 * volume_measured)^(2/3)) / surface_area;
                    
                    fprintf('DEBUG Sphericity (regionprops3): Volume=%.1f, SurfaceArea=%.1f, Sphericity=%.3f\n', ...
                        volume_measured, surface_area, sphericity);
                    
                    % Clamp to [0, 1] range
                    sphericity = max(0, min(1, sphericity));
                    return;
                end
            catch
                % Fall through to alternative method
            end
        end
        
        % Method 2: Improved surface area estimation using boundary voxels
        % This is more accurate than the face-counting method
        try
            % Smooth the mask slightly to get better surface estimation
            if volume > 8  % Only smooth if we have enough voxels
                smoothed_mask = smooth3(double(nucleus_mask), 'gaussian', [3 3 3], 0.5);
                smoothed_mask = smoothed_mask > 0.5;
            else
                smoothed_mask = nucleus_mask;
            end
            
            % Use isosurface to estimate surface area
            [faces, vertices] = isosurface(smoothed_mask, 0.5);
            
            if ~isempty(faces) && size(faces, 1) > 0
                % Calculate surface area from triangular faces
                surface_area = 0;
                for i = 1:size(faces, 1)
                    v1 = vertices(faces(i, 1), :);
                    v2 = vertices(faces(i, 2), :);
                    v3 = vertices(faces(i, 3), :);
                    
                    % Calculate triangle area using cross product
                    edge1 = v2 - v1;
                    edge2 = v3 - v1;
                    triangle_area = 0.5 * norm(cross(edge1, edge2));
                    surface_area = surface_area + triangle_area;
                end
                
                if surface_area > 0
                    sphericity = (pi^(1/3)) * ((6 * volume)^(2/3)) / surface_area;
                    
                    fprintf('DEBUG Sphericity (isosurface): Volume=%.1f, SurfaceArea=%.1f, Sphericity=%.3f\n', ...
                        volume, surface_area, sphericity);
                    
                    % Clamp to [0, 1] range
                    sphericity = max(0, min(1, sphericity));
                    return;
                end
            end
        catch
            % Fall through to simple method
        end
        
        % Method 3: Simple boundary-based estimation (fallback)
        % Count boundary voxels (voxels with at least one non-nucleus neighbor)
        boundary_voxels = 0;
        [y_coords, x_coords, z_coords] = ind2sub(size(mask), pixel_indices);
        
        for i = 1:length(pixel_indices)
            x = x_coords(i);
            y = y_coords(i);
            z = z_coords(i);
            
            % Check if this voxel is on the boundary (has non-nucleus neighbors)
            is_boundary = false;
            
            % Check 6-connected neighbors
            neighbors = [x+1,y,z; x-1,y,z; x,y+1,z; x,y-1,z; x,y,z+1; x,y,z-1];
            
            for j = 1:size(neighbors, 1)
                nx = neighbors(j, 1);
                ny = neighbors(j, 2);
                nz = neighbors(j, 3);
                
                % Check if neighbor is outside bounds or not part of nucleus
                if nx < 1 || nx > size(mask, 2) || ...
                   ny < 1 || ny > size(mask, 1) || ...
                   nz < 1 || nz > size(mask, 3) || ...
                   ~nucleus_mask(ny, nx, nz)
                    is_boundary = true;
                    break;
                end
            end
            
            if is_boundary
                boundary_voxels = boundary_voxels + 1;
            end
        end
        
        % Estimate surface area (each boundary voxel contributes ~1 unit of surface)
        surface_area = boundary_voxels;
        
        if surface_area > 0
            sphericity = (pi^(1/3)) * ((6 * volume)^(2/3)) / surface_area;
            
            fprintf('DEBUG Sphericity (boundary): Volume=%.1f, BoundaryVoxels=%.1f, Sphericity=%.3f\n', ...
                volume, boundary_voxels, sphericity);
        else
            sphericity = 0;
        end
        
        % Clamp to [0, 1] range
        sphericity = max(0, min(1, sphericity));
        
    catch ME
        fprintf('DEBUG: Error calculating sphericity: %s\n', ME.message);
        sphericity = 0;
    end
end

function sphericity = calculateSphericity3D_debug(pixel_indices, mask, nuc_id)
    % DEBUG VERSION: Calculate sphericity with extensive debugging
    fprintf('  [Nucleus %d] Starting sphericity calculation...\n', nuc_id);
    
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        fprintf('  [Nucleus %d] Created empty 3D mask of size %s\n', nuc_id, mat2str(size(nucleus_mask)));
        
        % Check pixel indices validity
        if any(pixel_indices > numel(nucleus_mask)) || any(pixel_indices < 1)
            fprintf('  [Nucleus %d] ERROR: Invalid pixel indices!\n', nuc_id);
            sphericity = 0;
            return;
        end
        
        nucleus_mask(pixel_indices) = true;
        volume = length(pixel_indices);
        fprintf('  [Nucleus %d] Set %d voxels to true (volume=%d)\n', nuc_id, length(pixel_indices), volume);
        
        if volume == 0
            fprintf('  [Nucleus %d] ERROR: Zero volume!\n', nuc_id);
            sphericity = 0;
            return;
        end
        
        % Check if we actually have a connected shape
        cc = bwconncomp(nucleus_mask, 26);
        fprintf('  [Nucleus %d] Found %d connected components in 3D\n', nuc_id, cc.NumObjects);
        
        if cc.NumObjects == 0
            fprintf('  [Nucleus %d] ERROR: No connected components found!\n', nuc_id);
            sphericity = 0;
            return;
        end
        
        % Method 1: Try regionprops3 if available
        if exist('regionprops3', 'file')
            try
                props = regionprops3(nucleus_mask, 'Volume', 'SurfaceArea');
                fprintf('  [Nucleus %d] regionprops3 returned %d objects\n', nuc_id, height(props));
                
                if ~isempty(props) && props.SurfaceArea > 0
                    volume_measured = props.Volume;
                    surface_area = props.SurfaceArea;
                    
                    fprintf('  [Nucleus %d] regionprops3 measurements - Volume: %.2f, SurfaceArea: %.2f\n', ...
                        nuc_id, volume_measured, surface_area);
                    
                    % Calculate sphericity
                    raw_sphericity = (pi^(1/3)) * ((6 * volume_measured)^(2/3)) / surface_area;
                    fprintf('  [Nucleus %d] Raw sphericity (regionprops3): %.6f\n', nuc_id, raw_sphericity);
                    
                    % Check for perfect sphere (suspicious)
                    if raw_sphericity > 0.99
                        fprintf('  [Nucleus %d] WARNING: Near-perfect sphericity detected! This may indicate an error.\n', nuc_id);
                    end
                    
                    % Clamp to maximum of 1.0 (perfect sphere) since values > 1.0 are mathematically impossible
                    sphericity = min(1.0, raw_sphericity);
                    
                    if raw_sphericity > 1.0
                        fprintf('  [Nucleus %d] NOTE: Clamped sphericity from %.6f to 1.0000 (measurement precision limit)\n', nuc_id, raw_sphericity);
                    end
                    fprintf('  [Nucleus %d] Final sphericity (regionprops3): %.3f\n', nuc_id, sphericity);
                    return;
                end
            catch ME
                fprintf('  [Nucleus %d] regionprops3 failed: %s\n', nuc_id, ME.message);
            end
        else
            fprintf('  [Nucleus %d] regionprops3 not available, using fallback method\n', nuc_id);
        end
        
        % Method 2: Boundary voxel counting (fallback)
        fprintf('  [Nucleus %d] Using boundary voxel counting method...\n', nuc_id);
        
        boundary_voxels = 0;
        [y_coords, x_coords, z_coords] = ind2sub(size(mask), pixel_indices);
        
        for i = 1:length(pixel_indices)
            x = x_coords(i);
            y = y_coords(i);
            z = z_coords(i);
            
            % Check 6-connected neighbors
            neighbors = [x+1,y,z; x-1,y,z; x,y+1,z; x,y-1,z; x,y,z+1; x,y,z-1];
            is_boundary = false;
            
            for j = 1:size(neighbors, 1)
                nx = neighbors(j, 1);
                ny = neighbors(j, 2);
                nz = neighbors(j, 3);
                
                % Check if neighbor is outside bounds or not part of nucleus
                if nx < 1 || nx > size(mask, 2) || ...
                   ny < 1 || ny > size(mask, 1) || ...
                   nz < 1 || nz > size(mask, 3) || ...
                   ~nucleus_mask(ny, nx, nz)
                    is_boundary = true;
                    break;
                end
            end
            
            if is_boundary
                boundary_voxels = boundary_voxels + 1;
            end
        end
        
        fprintf('  [Nucleus %d] Boundary voxel count: %d (out of %d total voxels)\n', ...
            nuc_id, boundary_voxels, volume);
        
        if boundary_voxels == 0
            fprintf('  [Nucleus %d] ERROR: No boundary voxels found!\n', nuc_id);
            sphericity = 0;
            return;
        end
        
        % Calculate sphericity using boundary voxels as surface area estimate
        surface_area = boundary_voxels;
        raw_sphericity = (pi^(1/3)) * ((6 * volume)^(2/3)) / surface_area;
        
        fprintf('  [Nucleus %d] Raw sphericity (boundary): %.6f\n', nuc_id, raw_sphericity);
        
        % Check for perfect sphere (suspicious)
        if raw_sphericity > 0.99
            fprintf('  [Nucleus %d] WARNING: Near-perfect sphericity detected! This may indicate an error.\n', nuc_id);
        end
        
        sphericity = max(0, min(1, raw_sphericity));
        fprintf('  [Nucleus %d] Final sphericity: %.3f\n', nuc_id, sphericity);
        
    catch ME
        fprintf('  [Nucleus %d] ERROR in sphericity calculation: %s\n', nuc_id, ME.message);
        sphericity = 0;
    end
end

function signal_count = countCachedSignalsInNucleus(handles, channel_idx, pixel_indices, nuc_mask)
    % Count cached local maxima (filtered signals) within the nucleus region
    % UNIFIED DATA ACCESS: Uses the same data source as preview system
    signal_count = 0;
    
    fprintf('DEBUG: *** ENTERING countCachedSignalsInNucleus for channel %d ***\n', channel_idx);
    fprintf('DEBUG: Nucleus has %d pixels, mask dimensions: %s\n', length(pixel_indices), mat2str(size(nuc_mask)));
    
    try
        % CRITICAL FIX: Try multiple data sources for maxima coordinates
        % This ensures consistency between preview and signal analysis
        cached_maxima = [];
        
        % First priority: Use handles.maximaCoords (contains filtered coordinates)
        if isfield(handles, 'maximaCoords') && channel_idx <= length(handles.maximaCoords) && ...
           ~isempty(handles.maximaCoords{channel_idx})
            cached_maxima = handles.maximaCoords{channel_idx};
            fprintf('DEBUG: Using maximaCoords for channel %d signal counting (%d maxima)\n', channel_idx, size(cached_maxima, 1));
            fprintf('DEBUG: First few maxima coordinates: %s\n', mat2str(cached_maxima(1:min(3,size(cached_maxima,1)), :)));
        
        % Second priority: Use preview cache (fallback)
        elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
               channel_idx <= length(handles.previewCache.channels) && ...
               ~isempty(handles.previewCache.channels{channel_idx}) && ...
               isfield(handles.previewCache.channels{channel_idx}, 'allMaxima')
            cached_maxima = handles.previewCache.channels{channel_idx}.allMaxima;
            fprintf('DEBUG: Using previewCache for channel %d signal counting (%d maxima)\n', channel_idx, size(cached_maxima, 1));
            fprintf('DEBUG: First few maxima coordinates: %s\n', mat2str(cached_maxima(1:min(3,size(cached_maxima,1)), :)));
        
        else
            fprintf('DEBUG: No maxima data found for channel %d\n', channel_idx);
            fprintf('DEBUG: handles.maximaCoords exists: %d\n', isfield(handles, 'maximaCoords'));
            if isfield(handles, 'maximaCoords')
                fprintf('DEBUG: maximaCoords length: %d, channel_idx: %d\n', length(handles.maximaCoords), channel_idx);
            end
            return;
        end
        
        if isempty(cached_maxima)
            return;
        end
        
        % Create a mask for just this nucleus
        % CRITICAL FIX: Handle dimension mismatch between nucleus mask and maxima coordinates
        if ndims(nuc_mask) == 2 && size(cached_maxima, 2) >= 3
            % 2D nucleus mask but 3D maxima coordinates
            % Create a 3D nucleus mask by replicating across Z-slices
            fprintf('DEBUG: Creating 3D nucleus mask from 2D mask for 3D maxima\n');
            
            % Determine Z-dimension from maxima coordinates
            max_z = max(cached_maxima(:, 3));
            min_z = min(cached_maxima(:, 3));
            z_slices = ceil(max_z);
            
            fprintf('DEBUG: Maxima Z range: %.1f to %.1f, creating %d Z-slices\n', min_z, max_z, z_slices);
            
            % Create 3D mask
            nucleus_mask_3d = repmat(nuc_mask, [1, 1, z_slices]);
            nucleus_mask = false(size(nucleus_mask_3d));
            nucleus_mask(pixel_indices) = true;
            
            fprintf('DEBUG: Created 3D nucleus mask, size: %s\n', mat2str(size(nucleus_mask)));
        else
            % Standard case - same dimensions
            nucleus_mask = false(size(nuc_mask));
            nucleus_mask(pixel_indices) = true;
            fprintf('DEBUG: Using standard nucleus mask, size: %s\n', mat2str(size(nucleus_mask)));
        end
        
        % DEBUG: Add detailed information about this nucleus
        fprintf('DEBUG: Counting signals in nucleus with %d pixels, mask size: %s\n', ...
            length(pixel_indices), mat2str(size(nucleus_mask)));
        
        % Count how many cached maxima fall within this nucleus
        % cached_maxima is in ARRAY CONVENTION [row, col, slice]
        signals_found = [];
        for i = 1:size(cached_maxima, 1)
            row = round(cached_maxima(i, 1));
            col = round(cached_maxima(i, 2));
            
            if size(cached_maxima, 2) >= 3
                slice = round(cached_maxima(i, 3));
                % Check bounds for 3D (ARRAY CONVENTION)
                if row >= 1 && row <= size(nucleus_mask, 1) && ...
                   col >= 1 && col <= size(nucleus_mask, 2) && ...
                   slice >= 1 && slice <= size(nucleus_mask, 3)
                    
                    % DIRECT array access with ARRAY CONVENTION
                    is_inside = nucleus_mask(row, col, slice);
                    if is_inside
                        signal_count = signal_count + 1;
                        signals_found(end+1) = i;
                        fprintf('DEBUG: Signal %d at [%d,%d,%d] is INSIDE nucleus\n', i, row, col, slice);
                    else
                        fprintf('DEBUG: Signal %d at [%d,%d,%d] is outside nucleus\n', i, row, col, slice);
                    end
                else
                    fprintf('DEBUG: Signal %d at [%d,%d,%d] is OUT OF BOUNDS (mask size: %s)\n', ...
                        i, row, col, slice, mat2str(size(nucleus_mask)));
                end
            else
                % 2D coordinates - check if nucleus mask is 3D and project it
                if ndims(nucleus_mask) == 3
                    nucleus_mask_2d = max(nucleus_mask, [], 3) > 0;
                else
                    nucleus_mask_2d = nucleus_mask;
                end
                
                % Check bounds for 2D (ARRAY CONVENTION)
                if row >= 1 && row <= size(nucleus_mask_2d, 1) && ...
                   col >= 1 && col <= size(nucleus_mask_2d, 2)
                    
                    % DIRECT array access with ARRAY CONVENTION
                    is_inside = nucleus_mask_2d(row, col);
                    if is_inside
                        signal_count = signal_count + 1;
                        signals_found(end+1) = i;
                        fprintf('DEBUG: Signal %d at [%d,%d] is INSIDE nucleus (2D mode)\n', i, row, col);
                    else
                        fprintf('DEBUG: Signal %d at [%d,%d] is outside nucleus (2D mode)\n', i, row, col);
                    end
                else
                    fprintf('DEBUG: Signal %d at [%d,%d] is OUT OF BOUNDS (2D mask size: %s)\n', ...
                        i, row, col, mat2str(size(nucleus_mask_2d)));
                end
            end
        end
        
        fprintf('DEBUG: Found %d signals in this nucleus (indices: %s)\n', ...
            signal_count, mat2str(signals_found));
        
    catch ME
        fprintf('DEBUG: Error counting cached signals: %s\n', ME.message);
        signal_count = 0;
    end
end

function validateDataConsistency(handles, active_channels)
    % Validate data consistency between preview system and signal analysis
    % Provides debugging information about data sources
    %
    % TESTING GUIDE FOR COHESIVE SYSTEM:
    % 1. Load images and run "Update Live Preview"
    % 2. Enable nuclei filtering (Include/Exclude Inside Nuclei)
    % 3. Check that preview shows filtered maxima
    % 4. Run nuclei signal analysis
    % 5. Verify that signal counts in analysis match what's shown in preview
    % 6. Look for "✓ CONSISTENT" messages in console output
    
    fprintf('\n=== DATA CONSISTENCY VALIDATION ===\n');
    
    for ch = active_channels
        maxima_coords_count = 0;
        cache_count = 0;
        
        % Check handles.maximaCoords
        if isfield(handles, 'maximaCoords') && ch <= length(handles.maximaCoords) && ...
           ~isempty(handles.maximaCoords{ch})
            maxima_coords_count = size(handles.maximaCoords{ch}, 1);
        end
        
        % Check preview cache
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
           ch <= length(handles.previewCache.channels) && ...
           ~isempty(handles.previewCache.channels{ch}) && ...
           isfield(handles.previewCache.channels{ch}, 'allMaxima')
            cache_count = size(handles.previewCache.channels{ch}.allMaxima, 1);
        end
        
        fprintf('Channel %d: maximaCoords=%d, previewCache=%d', ch, maxima_coords_count, cache_count);
        
        if maxima_coords_count == cache_count
            fprintf(' ✓ CONSISTENT\n');
        else
            fprintf(' ⚠ INCONSISTENT\n');
        end
    end
    
    fprintf('=====================================\n\n');
end

function signal_count = countCachedSignalsInNucleus_fixed(handles, channel_idx, pixel_indices, nuc_mask)
    % Count cached local maxima (filtered signals) within the nucleus region
    % UNIFIED DATA ACCESS: Uses the same data source as preview system
    % FIXED VERSION: Properly handles Z-projection nuclei with 3D maxima
    signal_count = 0;
    
    try
        % CRITICAL FIX: Try multiple data sources for maxima coordinates
        % This ensures consistency between preview and signal analysis
        cached_maxima = [];
        
        % First priority: Use handles.maximaCoords (contains filtered coordinates)
        if isfield(handles, 'maximaCoords') && channel_idx <= length(handles.maximaCoords) && ...
           ~isempty(handles.maximaCoords{channel_idx})
            cached_maxima = handles.maximaCoords{channel_idx};
        
        % Fallback: Use preview cache
        elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
               channel_idx <= length(handles.previewCache.channels) && ...
               ~isempty(handles.previewCache.channels{channel_idx}) && ...
               isfield(handles.previewCache.channels{channel_idx}, 'allMaxima')
            cached_maxima = handles.previewCache.channels{channel_idx}.allMaxima;
        else
            return;
        end
        
        if isempty(cached_maxima)
            return;
        end
        
        % CRITICAL FIX: For Z-projection nuclei segmentation
        % The nucleus has the same X,Y shape across ALL Z-slices
        % So we only need to check if the maxima's X,Y coordinates fall within the nucleus shape
        
        % SANITY CHECK: If pixel_indices is too large, something is wrong
        max_expected_pixels = size(nuc_mask, 1) * size(nuc_mask, 2);
        if length(pixel_indices) > max_expected_pixels * 0.1  % More than 10% of image
            fprintf('DEBUG: WARNING - pixel_indices seems too large (%d pixels), skipping this nucleus\n', length(pixel_indices));
            return;
        end
        
        % Create a 2D nucleus mask from pixel indices
        % IMPORTANT: pixel_indices are linear indices, need to convert to 2D coordinates
        nucleus_mask_2d = false(size(nuc_mask, 1), size(nuc_mask, 2));
        
        try
            if ndims(nuc_mask) == 2
                % For 2D nucleus mask, pixel_indices should be correct
                % But validate they're within bounds
                valid_indices = pixel_indices(pixel_indices >= 1 & pixel_indices <= numel(nucleus_mask_2d));
                nucleus_mask_2d(valid_indices) = true;
            else
                % For 3D nucleus mask, convert linear indices to subscripts and project to 2D
                [y_coords, x_coords, ~] = ind2sub(size(nuc_mask), pixel_indices);
                % Remove duplicates and ensure we stay within 2D bounds
                unique_coords = unique([y_coords(:), x_coords(:)], 'rows');
                valid_idx = unique_coords(:,1) >= 1 & unique_coords(:,1) <= size(nucleus_mask_2d, 1) & ...
                           unique_coords(:,2) >= 1 & unique_coords(:,2) <= size(nucleus_mask_2d, 2);
                valid_coords = unique_coords(valid_idx, :);
                
                % Convert back to linear indices for 2D mask
                if ~isempty(valid_coords)
                    linear_indices_2d = sub2ind(size(nucleus_mask_2d), valid_coords(:,1), valid_coords(:,2));
                    nucleus_mask_2d(linear_indices_2d) = true;
                end
            end
        catch ME
            return;
        end
        
        % Count how many cached maxima fall within this nucleus (row, col coordinates only)
        % cached_maxima is in ARRAY CONVENTION [row, col, slice]
        signals_found = [];
        for i = 1:size(cached_maxima, 1)
            row = round(cached_maxima(i, 1));
            col = round(cached_maxima(i, 2));
            
            % For Z-projection nuclei, we ONLY check row, col coordinates
            % The nucleus is assumed to exist at ALL Z-levels with the same row, col shape
            if row >= 1 && row <= size(nucleus_mask_2d, 1) && ...
               col >= 1 && col <= size(nucleus_mask_2d, 2)
                
                % DIRECT array access with ARRAY CONVENTION
                if nucleus_mask_2d(row, col)
                    signal_count = signal_count + 1;
                    if isempty(signals_found)
                        signals_found = i;
                    else
                        signals_found(end+1) = i;
                    end
                end
            end
        end
        
        % Signal counting complete
        
    catch ME
        signal_count = 0;
    end
end

function [signal_indices, signal_data_array] = getSignalsInNucleus(handles, channel_idx, pixel_indices, nuc_mask)
    % Get signals within nucleus with FULL fit data using shared helpers
    % Returns signal indices and complete signal data array
    % USES computeSignalMeasurements for ZERO REDUNDANCY
    
    signal_indices = [];
    signal_data_array = struct([]);  % Initialize as empty struct array
    
    try
        % Get cached maxima coordinates - PRIORITY: previewCache (more reliable during cache creation)
        cached_maxima = [];
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
           channel_idx <= length(handles.previewCache.channels) && ...
           ~isempty(handles.previewCache.channels{channel_idx}) && ...
           isfield(handles.previewCache.channels{channel_idx}, 'allMaxima')
            cached_maxima = handles.previewCache.channels{channel_idx}.allMaxima;
        elseif isfield(handles, 'maximaCoords') && channel_idx <= length(handles.maximaCoords) && ...
           ~isempty(handles.maximaCoords{channel_idx})
            cached_maxima = handles.maximaCoords{channel_idx};
        end
        
        if isempty(cached_maxima)
            return;  % Returns [], struct([])
        end
        
        % Get fit results - PRIORITY: previewCache
        fit_results = [];
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
           channel_idx <= length(handles.previewCache.channels) && ...
           ~isempty(handles.previewCache.channels{channel_idx}) && ...
           isfield(handles.previewCache.channels{channel_idx}, 'fit_results')
            fit_results = handles.previewCache.channels{channel_idx}.fit_results;
        elseif isfield(handles, 'gaussFitResults') && channel_idx <= length(handles.gaussFitResults)
            fit_results = handles.gaussFitResults{channel_idx};
        end
        
        % Get fitting method
        fit_method = 'None';
        if isfield(handles, 'gaussFitMethodDrop') && channel_idx <= length(handles.gaussFitMethodDrop) && ...
           isvalid(handles.gaussFitMethodDrop(channel_idx))
            fit_method = handles.gaussFitMethodDrop(channel_idx).Value;
        end
        
        % Create 2D nucleus mask from pixel indices
        nucleus_mask_2d = false(size(nuc_mask, 1), size(nuc_mask, 2));
        
        if ndims(nuc_mask) == 2
            valid_indices = pixel_indices(pixel_indices >= 1 & pixel_indices <= numel(nucleus_mask_2d));
            nucleus_mask_2d(valid_indices) = true;
        else
            [y_coords, x_coords, ~] = ind2sub(size(nuc_mask), pixel_indices);
            unique_coords = unique([y_coords(:), x_coords(:)], 'rows');
            valid_idx = unique_coords(:,1) >= 1 & unique_coords(:,1) <= size(nucleus_mask_2d, 1) & ...
                       unique_coords(:,2) >= 1 & unique_coords(:,2) <= size(nucleus_mask_2d, 2);
            valid_coords = unique_coords(valid_idx, :);
            
            if ~isempty(valid_coords)
                linear_indices_2d = sub2ind(size(nucleus_mask_2d), valid_coords(:,1), valid_coords(:,2));
                nucleus_mask_2d(linear_indices_2d) = true;
            end
        end
        
        % Find which signals fall within this nucleus
        signals_found = [];
        for i = 1:size(cached_maxima, 1)
            row = round(cached_maxima(i, 1));
            col = round(cached_maxima(i, 2));
            
            if row >= 1 && row <= size(nucleus_mask_2d, 1) && ...
               col >= 1 && col <= size(nucleus_mask_2d, 2)
                if nucleus_mask_2d(row, col)
                    signals_found(end+1) = i;
                end
            end
        end
        
        signal_indices = signals_found;
        
        % If no signals found, return empty consistently
        if isempty(signal_indices)
            signal_data_array = struct([]);
            return;  % Returns [], struct([])
        end
        
        % USE SHARED HELPER to get full signal data (ZERO REDUNDANCY!)
        % This guarantees consistency with channel export
        [all_signal_data, ~] = snap_helpers.computeSignalMeasurements(cached_maxima, fit_results, fit_method);
        
        % Extract only the signals that belong to this nucleus
        % MATLAB struct array indexing - ensure proper extraction
        if ~isempty(all_signal_data)
            signal_data_array = all_signal_data(signal_indices);
        else
            signal_data_array = struct([]);  % Empty struct array
        end
        
    catch ME
        warning('Error getting signals in nucleus: %s', ME.message);
        signal_indices = [];
        signal_data_array = struct([]);  % Return empty struct array, not empty matrix
    end
end

function solidity = calculateSolidity2D(pixel_indices, mask)
    % Calculate 2D solidity: Area / Convex Hull Area
    % Perfect convexity = 1.0, concave shapes < 1.0
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        nucleus_mask(pixel_indices) = true;
        
        % Use regionprops for accurate measurements
        props = regionprops(nucleus_mask, 'Area', 'Solidity');
        
        if isempty(props) || props.Area == 0
            solidity = 0;
            return;
        end
        
        solidity = props.Solidity;
        
        % Ensure solidity is in valid range
        solidity = max(0, min(1.0, solidity));
        
    catch ME
        solidity = 0;
    end
end

function solidity = calculateSolidity3D(pixel_indices, mask)
    % Calculate 3D solidity: Volume / Convex Hull Volume
    % Perfect convexity = 1.0, concave shapes < 1.0
    try
        % Create a binary mask for just this nucleus
        nucleus_mask = false(size(mask));
        nucleus_mask(pixel_indices) = true;
        
        % Use regionprops3 if available (MATLAB R2017b+)
        if exist('regionprops3', 'file')
            try
                props = regionprops3(nucleus_mask, 'Volume', 'Solidity');
                if ~isempty(props) && height(props) > 0
                    solidity = props.Solidity(1);
                    
                    % Ensure solidity is in valid range
                    solidity = max(0, min(1.0, solidity));
                    return;
                end
            catch ME
                % regionprops3 failed, fall through to fallback
            end
        end
        
        % Fallback: Return 1.0 if regionprops3 is not available
        % This is better than returning 0 (which would filter out all nuclei)
        solidity = 1.0;
        
    catch ME
        solidity = 1.0; % Conservative fallback
    end
end
