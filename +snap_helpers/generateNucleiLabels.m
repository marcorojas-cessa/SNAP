function [centroids_3d, labels_by_slice] = generateNucleiLabels(nuc_mask)
    % Generate centroids and labels for each segmented nucleus
    % Returns:
    %   centroids_3d: [N x 3] array of [x, y, z] centroids for all nuclei
    %   labels_by_slice: cell array where each cell contains [x, y, label_num] for that Z-slice
    
    centroids_3d = [];
    labels_by_slice = {};
    
    if isempty(nuc_mask)
        return;
    end
    
    % Label connected components
    if ndims(nuc_mask) == 3
        % 3D case
        cc = bwconncomp(nuc_mask, 26);
        num_slices = size(nuc_mask, 3);
        labels_by_slice = cell(num_slices, 1);
        
        % Initialize each slice
        for z = 1:num_slices
            labels_by_slice{z} = [];
        end
        
        % Process each nucleus
        centroids_3d = zeros(cc.NumObjects, 3);
        for i = 1:cc.NumObjects
            % Get pixel indices for this nucleus
            pixel_indices = cc.PixelIdxList{i};
            
            % Convert linear indices to subscripts
            [y_coords, x_coords, z_coords] = ind2sub(size(nuc_mask), pixel_indices);
            
            % Calculate 3D centroid
            centroid_x = mean(x_coords);
            centroid_y = mean(y_coords);
            centroid_z = mean(z_coords);
            
            centroids_3d(i, :) = [centroid_x, centroid_y, centroid_z];
            
            % Add to labels_by_slice for each Z-slice where this nucleus exists
            unique_z = unique(z_coords);
            for z_slice = unique_z'
                if isempty(labels_by_slice{z_slice})
                    labels_by_slice{z_slice} = [centroid_x, centroid_y, i];
                else
                    labels_by_slice{z_slice} = [labels_by_slice{z_slice}; centroid_x, centroid_y, i];
                end
            end
        end
        
    else
        % 2D case
        cc = bwconncomp(nuc_mask, 8);
        labels_by_slice = cell(1, 1);
        
        % Process each nucleus
        centroids_3d = zeros(cc.NumObjects, 3);
        slice_labels = [];
        
        for i = 1:cc.NumObjects
            % Get pixel indices for this nucleus
            pixel_indices = cc.PixelIdxList{i};
            
            % Convert linear indices to subscripts
            [y_coords, x_coords] = ind2sub(size(nuc_mask), pixel_indices);
            
            % Calculate 2D centroid (Z = 1 for 2D)
            centroid_x = mean(x_coords);
            centroid_y = mean(y_coords);
            
            centroids_3d(i, :) = [centroid_x, centroid_y, 1];
            
            % Add to slice labels
            if isempty(slice_labels)
                slice_labels = [centroid_x, centroid_y, i];
            else
                slice_labels = [slice_labels; centroid_x, centroid_y, i];
            end
        end
        
        labels_by_slice{1} = slice_labels;
    end
end
