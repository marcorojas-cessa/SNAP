function maximaCoords = detectMaxima(processedImg, handles, channelIdx)
% DETECTMAXIMA - Detect local maxima using the same logic as SNAP GUI
% This function replicates the exact maxima detection from updateLivePreview
% to ensure SNAP and SNAP_batch produce identical results.
%
% COORDINATE CONVENTION: This function uses ARRAY CONVENTION
%   Coordinates are stored as [row, col, slice]
%   - row: 1st dimension (vertical index, Y-axis in Cartesian display)
%   - col: 2nd dimension (horizontal index, X-axis in Cartesian display)
%   - slice: 3rd dimension (depth index, Z-axis in Cartesian display)
%
% Array access: imageData(row, col, slice) - DIRECT indexing
% Display conversion: [x, y, z] = [col, row, slice] for Cartesian plots
%
% Inputs:
%   processedImg - Processed image (2D or 3D)
%   handles      - Handles structure with processing parameters
%   channelIdx   - Channel index (1-based)
%
% Outputs:
%   maximaCoords - Nx3 array of [row, col, slice] coordinates (ARRAY CONVENTION)

    mode = handles.maximaModeDrops(channelIdx).Value;
    coords = [];
    
    if strcmp(mode, 'On Z-Projection')
        proj_type = handles.maximaProjectionDrops(channelIdx).Value;
        img_2d = projectZ(processedImg, proj_type);
        bw = find2DMaxima(img_2d, handles, channelIdx);
        % MATLAB find() returns [row, col] - keep this order (ARRAY CONVENTION)
        [row, col] = find(bw);
        coords = [row, col, zeros(size(row)) + 0.5];  % [row, col, slice]
    else % 3D or 2D Slice-by-slice
        is_3d_mode = strcmp(mode, '3D');
        bw = find3DMaxima(processedImg, handles, channelIdx, is_3d_mode);
        % MATLAB ind2sub() returns [row, col, slice] - keep this order (ARRAY CONVENTION)
        [row, col, slice] = ind2sub(size(bw), find(bw));
        coords = [row, col, slice];  % [row, col, slice]
    end
    
    maximaCoords = coords;
    
    % Nested helper functions (exact copies from updateLivePreview)
    
    function bw = find2DMaxima(img_2d, h, c)
        % Get neighborhood size (in pixels or microns if scaled)
        neighborhood_size = h.maximaNeighborhoodInputs(c).Value;
        
        % Check if scaling is enabled for maxima detection
        is_scaled = h.maximaScaleChecks(c).Value;
        if is_scaled && h.xySpacingInputs(c).Value > 0
            % Convert microns to pixels
            neighborhood_size = round(neighborhood_size / h.xySpacingInputs(c).Value);
        end
        
        method = h.maximaMethodDrops(c).Value;
        
        switch method
            case 'Simple Regional'
                bw = findLocalMaximaNeighborhood2D(img_2d, neighborhood_size);
            case 'Extended Maxima'
                h_threshold = h.hMaxInputs(c).Value;
                bw = findExtendedMaximaNeighborhood2D(img_2d, neighborhood_size, h_threshold);
            case 'Laplacian of Gaussian'
                bw = false(size(img_2d));
                sigma_val = h.logSigmaInputs(c).Value;
                threshold = h.logThresholdInputs(c).Value;
                if sigma_val > 0
                    % Check if scaling is enabled for LoG sigma
                    if is_scaled && h.xySpacingInputs(c).Value > 0
                        % Convert sigma from microns to pixels
                        sigma_pixels = sigma_val / h.xySpacingInputs(c).Value;
                    else
                        % Use sigma directly in pixels
                        sigma_pixels = sigma_val;
                    end
                    
                    filt_size = 2 * ceil(3*sigma_pixels) + 1;
                    log_filter = fspecial('log', filt_size, sigma_pixels);
                    img_log = -imfilter(img_2d, log_filter, 'replicate');
                    bw = findLocalMaximaNeighborhood2D(img_log, neighborhood_size) & (img_log > threshold);
                end
            otherwise
                bw = false(size(img_2d));
        end
    end

    function bw = find3DMaxima(img_3d, h, c, use_3d_conn)
        % Get neighborhood size (in pixels/voxels or microns if scaled)
        neighborhood_size = h.maximaNeighborhoodInputs(c).Value;
        
        % Check if scaling is enabled for maxima detection
        is_scaled = h.maximaScaleChecks(c).Value;
        if is_scaled && h.xySpacingInputs(c).Value > 0 && h.zSpacingInputs(c).Value > 0
            % Convert microns to pixels/voxels
            neighborhood_size_xy = round(neighborhood_size / h.xySpacingInputs(c).Value);
            neighborhood_size_z = round(neighborhood_size / h.zSpacingInputs(c).Value);
        else
            % Use same size for all dimensions if not scaled
            neighborhood_size_xy = neighborhood_size;
            neighborhood_size_z = neighborhood_size;
        end
        
        method = h.maximaMethodDrops(c).Value;
        bw = false(size(img_3d));
        
        switch method
            case 'Simple Regional'
                if use_3d_conn
                    bw = findLocalMaximaNeighborhood3D(img_3d, neighborhood_size_xy, neighborhood_size_z);
                else % Slice-by-slice 2D
                    for z_slice = 1:size(img_3d, 3)
                        bw(:,:,z_slice) = findLocalMaximaNeighborhood2D(img_3d(:,:,z_slice), neighborhood_size_xy);
                    end
                end
            case 'Extended Maxima'
                h_threshold = h.hMaxInputs(c).Value;
                if use_3d_conn
                    bw = findExtendedMaximaNeighborhood3D(img_3d, neighborhood_size_xy, neighborhood_size_z, h_threshold);
                else % Slice-by-slice 2D
                    for z_slice = 1:size(img_3d, 3)
                        bw(:,:,z_slice) = findExtendedMaximaNeighborhood2D(img_3d(:,:,z_slice), neighborhood_size_xy, h_threshold);
                    end
                end
            case 'Laplacian of Gaussian'
                sigma_val = h.logSigmaInputs(c).Value;
                threshold = h.logThresholdInputs(c).Value;
                if sigma_val > 0
                    % Use the Scale checkbox to determine if we should apply anisotropic scaling
                    is_scaled = h.maximaScaleChecks(c).Value && h.xySpacingInputs(c).Value > 0 && h.zSpacingInputs(c).Value > 0;

                    if is_scaled
                        % Convert sigma from microns to pixels/voxels
                        sigma_pixels_xy = sigma_val / h.xySpacingInputs(c).Value;
                        sigma_pixels_z = sigma_val / h.zSpacingInputs(c).Value;
                        
                        if ~use_3d_conn % Slice-by-slice 2D with scaled sigma
                            for z_slice = 1:size(img_3d, 3)
                                img_2d_slice = img_3d(:,:,z_slice);
                                filt_size = 2 * ceil(3*sigma_pixels_xy) + 1;
                                log_filter = fspecial('log', filt_size, sigma_pixels_xy);
                                img_log = -imfilter(img_2d_slice, log_filter, 'replicate');
                                bw(:,:,z_slice) = findLocalMaximaNeighborhood2D(img_log, neighborhood_size_xy) & (img_log > threshold);
                            end
                        else % True 3D Anisotropic LoG
                            % Create anisotropic 3D LoG filter
                            filt_size_xy = 2 * ceil(3*sigma_pixels_xy) + 1;
                            filt_size_z = 2 * ceil(3*sigma_pixels_z) + 1;
                            
                            % Create 3D LoG filter with different sigma for XY and Z
                            [X, Y, Z] = meshgrid(-(filt_size_xy-1)/2:(filt_size_xy-1)/2, ...
                                                  -(filt_size_xy-1)/2:(filt_size_xy-1)/2, ...
                                                  -(filt_size_z-1)/2:(filt_size_z-1)/2);
                            
                            % Normalize coordinates by sigma
                            X_norm = X / sigma_pixels_xy;
                            Y_norm = Y / sigma_pixels_xy;
                            Z_norm = Z / sigma_pixels_z;
                            
                            % 3D LoG formula: ∇²G = (r² - 2σ²)/(σ⁴) * exp(-r²/(2σ²))
                            % where r² = (x/σx)² + (y/σy)² + (z/σz)²
                            r_squared = X_norm.^2 + Y_norm.^2 + Z_norm.^2;
                            sigma_eff = 1; % Normalized sigma
                            
                            % 3D LoG kernel - use element-wise multiplication
                            log_kernel = (r_squared - 2*sigma_eff^2) / (sigma_eff^4) .* exp(-r_squared/(2*sigma_eff^2));
                            
                            % Normalize kernel to sum to zero
                            log_kernel = log_kernel - mean(log_kernel(:));
                            
                            % Apply 3D LoG filter
                            img_log = -imfilter(img_3d, log_kernel, 'replicate');
                            
                            % Find maxima in the LoG result
                            bw = findLocalMaximaNeighborhood3D(img_log, neighborhood_size_xy, neighborhood_size_z) & (img_log > threshold);
                        end
                    else % Unscaled (isotropic voxels)
                        if ~use_3d_conn % Slice-by-slice 2D
                            for z_slice = 1:size(img_3d, 3)
                                img_2d_slice = img_3d(:,:,z_slice);
                                filt_size = 2 * ceil(3*sigma_val) + 1;
                                log_filter = fspecial('log', filt_size, sigma_val);
                                img_log = -imfilter(img_2d_slice, log_filter, 'replicate');
                                bw(:,:,z_slice) = findLocalMaximaNeighborhood2D(img_log, neighborhood_size_xy) & (img_log > threshold);
                            end
                        else % 3D isotropic LoG
                            % Create isotropic 3D LoG filter
                            filt_size = 2 * ceil(3*sigma_val) + 1;
                            log_filter = fspecial('log', filt_size, sigma_val);
                            
                            % Apply 2D LoG to each slice and combine
                            img_log_3d = zeros(size(img_3d));
                            for z_slice = 1:size(img_3d, 3)
                                img_2d_slice = img_3d(:,:,z_slice);
                                img_log_3d(:,:,z_slice) = -imfilter(img_2d_slice, log_filter, 'replicate');
                            end
                            
                            % Find maxima in the 3D LoG result
                            bw = findLocalMaximaNeighborhood3D(img_log_3d, neighborhood_size_xy, neighborhood_size_z) & (img_log_3d > threshold);
                        end
                    end
                end
        end
    end
    
    % Helper functions for neighborhood-based maxima detection
    
    function bw = findLocalMaximaNeighborhood2D(img, neighborhood_size)
        % Find local maxima using a neighborhood window approach for 2D images
        % neighborhood_size: radius in pixels (e.g., 1 means 3x3 window, 2 means 5x5 window)
        
        % Ensure neighborhood_size is integer to avoid colon operator errors
        neighborhood_size = round(neighborhood_size);
        
        if neighborhood_size <= 0
            bw = false(size(img));
            return;
        end
        
        [rows, cols] = size(img);
        bw = false(rows, cols);
        
        % Create a more robust maxima detection
        for r = 1:rows
            for c = 1:cols
                center_value = img(r, c);
                
                % Skip if center is zero or negative
                if center_value <= 0
                    continue;
                end
                
                % Define neighborhood boundaries
                r_start = max(1, r - neighborhood_size);
                r_end = min(rows, r + neighborhood_size);
                c_start = max(1, c - neighborhood_size);
                c_end = min(cols, c + neighborhood_size);
                
                % Extract neighborhood
                neighborhood = img(r_start:r_end, c_start:c_end);
                
                % Check if center is maximum in neighborhood
                if center_value == max(neighborhood(:))
                    % Additional check: ensure center is not at the edge of a plateau
                    % Count how many pixels have the same value as center
                    same_value_count = sum(neighborhood(:) == center_value);
                    
                    % If only one pixel has this value, it's definitely a peak
                    if same_value_count == 1
                        bw(r, c) = true;
                    else
                        % If multiple pixels have the same value, check if center is the "first" one
                        % This handles plateaus more gracefully
                        center_idx = sub2ind(size(neighborhood), r - r_start + 1, c - c_start + 1);
                        max_positions = find(neighborhood == center_value);
                        if max_positions(1) == center_idx
                            bw(r, c) = true;
                        end
                    end
                end
            end
        end
    end
    
    function bw = findExtendedMaximaNeighborhood2D(img, neighborhood_size, h_threshold)
        % Find extended maxima using a neighborhood window approach for 2D images
        % neighborhood_size: radius in pixels (e.g., 1 means 3x3 window, 2 means 5x5 window)
        % h_threshold: minimum height difference from surrounding (in intensity units)
        
        if neighborhood_size <= 0 || h_threshold <= 0
            bw = false(size(img));
            return;
        end
        
        [rows, cols] = size(img);
        bw = false(rows, cols);
        
        for r = 1:rows
            for c = 1:cols
                center_value = img(r, c);
                
                % Skip if center is zero or negative
                if center_value <= 0
                    continue;
                end
                
                % Define neighborhood boundaries
                r_start = max(1, r - neighborhood_size);
                r_end = min(rows, r + neighborhood_size);
                c_start = max(1, c - neighborhood_size);
                c_end = min(cols, c + neighborhood_size);
                
                % Extract neighborhood
                neighborhood = img(r_start:r_end, c_start:c_end);
                
                % Check if center is maximum in neighborhood
                if center_value == max(neighborhood(:))
                    % Check height threshold: center must be h_threshold above surrounding
                    surrounding_values = neighborhood(:);
                    center_idx = sub2ind(size(neighborhood), r - r_start + 1, c - c_start + 1);
                    surrounding_values(center_idx) = [];
                    
                    % Calculate the height difference from the highest surrounding pixel
                    max_surrounding = max(surrounding_values);
                    height_diff = center_value - max_surrounding;
                    
                    if height_diff >= h_threshold
                        % Additional check: ensure center is not at the edge of a plateau
                        same_value_count = sum(neighborhood(:) == center_value);
                        
                        if same_value_count == 1
                            bw(r, c) = true;
                        else
                            % If multiple pixels have the same value, check if center is the "first" one
                            max_positions = find(neighborhood == center_value);
                            if max_positions(1) == center_idx
                                bw(r, c) = true;
                            end
                        end
                    end
                end
            end
        end
    end
    
    function bw = findLocalMaximaNeighborhood3D(img, neighborhood_size_xy, neighborhood_size_z)
        % Find local maxima using a neighborhood window approach for 3D images
        % neighborhood_size_xy: radius in XY plane (pixels)
        % neighborhood_size_z: radius in Z direction (voxels)
        
        % Ensure neighborhood sizes are integers to avoid colon operator errors
        neighborhood_size_xy = round(neighborhood_size_xy);
        neighborhood_size_z = round(neighborhood_size_z);
        
        if neighborhood_size_xy <= 0 || neighborhood_size_z <= 0
            bw = false(size(img));
            return;
        end
        
        [rows, cols, slices] = size(img);
        bw = false(rows, cols, slices);
        
        for r = 1:rows
            for c = 1:cols
                for z = 1:slices
                    center_value = img(r, c, z);
                    
                    % Skip if center is zero or negative
                    if center_value <= 0
                        continue;
                    end
                    
                    % Define neighborhood boundaries
                    r_start = max(1, r - neighborhood_size_xy);
                    r_end = min(rows, r + neighborhood_size_xy);
                    c_start = max(1, c - neighborhood_size_xy);
                    c_end = min(cols, c + neighborhood_size_xy);
                    z_start = max(1, z - neighborhood_size_z);
                    z_end = min(slices, z + neighborhood_size_z);
                    
                    % Extract neighborhood
                    neighborhood = img(r_start:r_end, c_start:c_end, z_start:z_end);
                    
                    % Check if center is maximum in neighborhood
                    if center_value == max(neighborhood(:))
                        % Additional check: ensure center is not at the edge of a plateau
                        same_value_count = sum(neighborhood(:) == center_value);
                        
                        if same_value_count == 1
                            bw(r, c, z) = true;
                        else
                            % If multiple pixels have the same value, check if center is the "first" one
                            center_idx = sub2ind(size(neighborhood), r - r_start + 1, c - c_start + 1, z - z_start + 1);
                            max_positions = find(neighborhood == center_value);
                            if max_positions(1) == center_idx
                                bw(r, c, z) = true;
                            end
                        end
                    end
                end
            end
        end
    end
    
    function bw = findExtendedMaximaNeighborhood3D(img, neighborhood_size_xy, neighborhood_size_z, h_threshold)
        % Find extended maxima using a neighborhood window approach for 3D images
        % neighborhood_size_xy: radius in XY plane (pixels)
        % neighborhood_size_z: radius in Z direction (voxels)
        % h_threshold: minimum height difference from surrounding (in intensity units)
        
        if neighborhood_size_xy <= 0 || neighborhood_size_z <= 0 || h_threshold <= 0
            bw = false(size(img));
            return;
        end
        
        [rows, cols, slices] = size(img);
        bw = false(rows, cols, slices);
        
        for r = 1:rows
            for c = 1:cols
                for z = 1:slices
                    center_value = img(r, c, z);
                    
                    % Skip if center is zero or negative
                    if center_value <= 0
                        continue;
                    end
                    
                    % Define neighborhood boundaries
                    r_start = max(1, r - neighborhood_size_xy);
                    r_end = min(rows, r + neighborhood_size_xy);
                    c_start = max(1, c - neighborhood_size_xy);
                    c_end = min(cols, c + neighborhood_size_xy);
                    z_start = max(1, z - neighborhood_size_z);
                    z_end = min(slices, z + neighborhood_size_z);
                    
                    % Extract neighborhood
                    neighborhood = img(r_start:r_end, c_start:c_end, z_start:z_end);
                    
                    % Check if center is maximum in neighborhood
                    if center_value == max(neighborhood(:))
                        % Check height threshold: center must be h_threshold above surrounding
                        surrounding_values = neighborhood(:);
                        center_idx = sub2ind(size(neighborhood), r - r_start + 1, c - c_start + 1, z - z_start + 1);
                        surrounding_values(center_idx) = [];
                        
                        % Calculate the height difference from the highest surrounding pixel
                        max_surrounding = max(surrounding_values);
                        height_diff = center_value - max_surrounding;
                        
                        if height_diff >= h_threshold
                            % Additional check: ensure center is not at the edge of a plateau
                            same_value_count = sum(neighborhood(:) == center_value);
                            
                            if same_value_count == 1
                                bw(r, c, z) = true;
                            else
                                % If multiple pixels have the same value, check if center is the "first" one
                                max_positions = find(neighborhood == center_value);
                                if max_positions(1) == center_idx
                                    bw(r, c, z) = true;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    function proj = projectZ(img, type)
        if ndims(img) < 3, proj = img; return; end
        if strcmp(type, 'Max'), proj = max(img, [], 3);
        elseif strcmp(type, 'Min'), proj = min(img, [], 3);
        elseif strcmp(type, 'Median'), proj = median(img, 3);
        elseif strcmp(type, 'Mean'), proj = mean(img, 3);
        else, proj = max(img, [], 3);
        end
    end
end

