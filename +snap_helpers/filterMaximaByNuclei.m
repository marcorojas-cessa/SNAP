function [filteredCoords, keepMask] = filterMaximaByNuclei(maximaCoords, nucleiMask, mode)
% FILTERMAXIMABYNUCLEI - Filter local maxima based on nuclei mask
% This function replicates the exact nuclei filtering logic from updateLivePreview
% to ensure SNAP and SNAP_batch produce identical results.
%
% COORDINATE CONVENTION: This function uses ARRAY CONVENTION
%   Coordinates are [row, col, slice]
%   - row: 1st dimension (vertical index, Y-axis in display)
%   - col: 2nd dimension (horizontal index, X-axis in display)
%   - slice: 3rd dimension (depth index, Z-axis in display)
%
% Array access: nucleiMask(row, col, slice) - DIRECT indexing
%
% Inputs:
%   maximaCoords - Nx3 array of [row, col, slice] coordinates (ARRAY CONVENTION)
%   nucleiMask   - Binary mask of nuclei regions (can be 2D or 3D)
%   mode         - 'Include Inside Nuclei' or 'Exclude Inside Nuclei'
%
% Outputs:
%   filteredCoords - Filtered coordinates array (ARRAY CONVENTION)
%   keepMask       - Boolean mask indicating which maxima were kept

    if isempty(maximaCoords)
        filteredCoords = maximaCoords;
        keepMask = [];
        return;
    end
    
    % Determine if nuclei mask is 2D or 3D
    is_2d_mask = ndims(nucleiMask) == 2;
    
    % maximaCoords are already in ARRAY CONVENTION [row, col, slice]
    % Can access nucleiMask DIRECTLY without swapping
    num_maxima = size(maximaCoords, 1);
    keep_indices = true(num_maxima, 1);
    
    for i = 1:num_maxima
        row = round(maximaCoords(i, 1));
        col = round(maximaCoords(i, 2));
        slice = round(maximaCoords(i, 3));
        
        % Initialize as excluded (safer default)
        is_inside_nuclei = false;
        is_valid_coordinate = false;
        
        if is_2d_mask
            % For 2D nuclei mask, only check row, col coordinates
            % Check bounds for 2D mask (ARRAY CONVENTION)
            if row >= 1 && row <= size(nucleiMask, 1) && ...
               col >= 1 && col <= size(nucleiMask, 2)
                
                is_valid_coordinate = true;
                % DIRECT array access with ARRAY CONVENTION (ignore slice for 2D)
                mask_value = nucleiMask(row, col);
                is_inside_nuclei = (mask_value > 0);  % Explicit: any non-zero = inside nucleus
            end
        else
            % For 3D nuclei mask, check all coordinates
            % Check bounds for 3D mask (ARRAY CONVENTION)
            if row >= 1 && row <= size(nucleiMask, 1) && ...
               col >= 1 && col <= size(nucleiMask, 2) && ...
               slice >= 1 && slice <= size(nucleiMask, 3)
                
                is_valid_coordinate = true;
                % DIRECT array access with ARRAY CONVENTION
                mask_value = nucleiMask(row, col, slice);
                is_inside_nuclei = (mask_value > 0);  % Explicit: any non-zero = inside nucleus
            end
        end
        
        % Apply filtering based on mode, but only for valid coordinates
        if is_valid_coordinate
            if strcmp(mode, 'Include Inside Nuclei')
                % Keep only maxima inside nuclei
                keep_indices(i) = is_inside_nuclei;
            elseif strcmp(mode, 'Exclude Inside Nuclei')
                % Keep only maxima outside nuclei
                keep_indices(i) = ~is_inside_nuclei;
            else
                % Unknown mode - KEEP all by default (safer for backward compatibility)
                warning('Unknown nuclei filtering mode: %s. No filtering applied.', mode);
                keep_indices(i) = true;
            end
        else
            % Invalid coordinates (out of bounds) - exclude them
            keep_indices(i) = false;
        end
    end
    
    % Return filtered coordinates and the mask
    filteredCoords = maximaCoords(keep_indices, :);
    keepMask = keep_indices;
end

