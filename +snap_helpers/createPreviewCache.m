function createPreviewCache(handles)
% Creates a comprehensive cache of all preview data to enable fast rendering
% This function pre-computes all expensive operations so preview interactions are instant
    
    try
        handles.statusLabel.Text = 'Status: Creating Preview Cache...';
        drawnow;
        
        % Initialize cache structure
        cache = struct();
        
        % --- Cache Nuclei Data ---
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
            try
                cache.nuclei = struct();
                
                % Pre-process nuclei image
                cache.nuclei.processed = snap_helpers.preprocessNucleiWithBgCorr(handles);
                
                % Segment nuclei with consistent labeling
                [cache.nuclei.mask, cache.nuclei.labels] = snap_helpers.segmentNuclei(cache.nuclei.processed, handles);
                
                % Validate nuclei labeling consistency
                if ~isempty(cache.nuclei.labels) && isfield(cache.nuclei.labels, 'centroids_3d') && ...
                   isfield(cache.nuclei.labels, 'labels_by_slice')
                    
                    % Verify that all slice labels reference valid nucleus IDs
                    max_id = size(cache.nuclei.labels.centroids_3d, 1);
                    total_nuclei = cache.nuclei.labels.num_nuclei;
                    
                    % Check for consistency
                    if max_id ~= total_nuclei
                        warning('Nuclei count mismatch: centroids_3d has %d entries but num_nuclei is %d', max_id, total_nuclei);
                    end
                    
                    for slice_idx = 1:length(cache.nuclei.labels.labels_by_slice)
                        if ~isempty(cache.nuclei.labels.labels_by_slice{slice_idx})
                            slice_labels = cache.nuclei.labels.labels_by_slice{slice_idx};
                            if size(slice_labels, 2) >= 3
                                label_ids = slice_labels(:, 3);
                            if any(label_ids > max_id) || any(label_ids < 1)
                                warning('Invalid nucleus IDs detected in slice %d: min=%d, max=%d (should be 1-%d). This indicates filtering did not properly renumber IDs.', ...
                                    slice_idx, min(label_ids), max(label_ids), max_id);
                            end
                            end
                        end
                    end
                end
                
                % Ensure mask is valid
                if isempty(cache.nuclei.mask) || ~islogical(cache.nuclei.mask)
                    warning('Invalid nuclei mask generated, skipping nuclei caching');
                    % Skip this nuclei caching section
                else
                
                % Pre-compute all z-slices and projections
                if ndims(cache.nuclei.mask) == 3
                    z_max = size(cache.nuclei.mask, 3);
                    cache.nuclei.slices = cell(1, z_max);
                    for z = 1:z_max
                        handles.statusLabel.Text = ['Status: Creating Preview Cache - Nuclei Z-frame ' num2str(z) '/' num2str(z_max) '...'];
                        drawnow;
                        cache.nuclei.slices{z} = cache.nuclei.mask(:,:,z);
                    end
                    
                    % Pre-compute projections
                    cache.nuclei.projectionMax = max(cache.nuclei.mask, [], 3);
                    cache.nuclei.projectionMin = min(cache.nuclei.mask, [], 3);
                    cache.nuclei.projectionMean = mean(cache.nuclei.mask, 3);
                    cache.nuclei.projectionMedian = median(cache.nuclei.mask, 3);
                else
                    % Single slice or projection
                    z_max = 1;
                    cache.nuclei.slices = {cache.nuclei.mask};
                    cache.nuclei.projectionMax = cache.nuclei.mask;
                    cache.nuclei.projectionMin = cache.nuclei.mask;
                    cache.nuclei.projectionMean = cache.nuclei.mask;
                    cache.nuclei.projectionMedian = cache.nuclei.mask;
                end
                
                % Pre-compute boundaries for all slices and projections
                cache.nuclei.boundaries = cell(1, z_max);
                for z = 1:z_max
                    handles.statusLabel.Text = ['Status: Creating Preview Cache - Computing boundaries Z-frame ' num2str(z) '/' num2str(z_max) '...'];
                    drawnow;
                    if z <= length(cache.nuclei.slices) && ~isempty(cache.nuclei.slices{z})
                        try
                            boundaries = bwboundaries(cache.nuclei.slices{z});
                            cache.nuclei.boundaries{z} = boundaries;
                        catch ME
                            warning('Failed to compute boundaries for nuclei slice %d: %s', z, ME.message);
                            cache.nuclei.boundaries{z} = {};
                        end
                    else
                        cache.nuclei.boundaries{z} = {};
                    end
                end
                
                % Pre-compute projection boundaries
                try
                    cache.nuclei.boundariesMax = bwboundaries(cache.nuclei.projectionMax);
                    cache.nuclei.boundariesMin = bwboundaries(cache.nuclei.projectionMin);
                    cache.nuclei.boundariesMean = bwboundaries(cache.nuclei.projectionMean);
                    cache.nuclei.boundariesMedian = bwboundaries(cache.nuclei.projectionMedian);
                catch ME
                    warning('Failed to compute nuclei projection boundaries: %s', ME.message);
                    cache.nuclei.boundariesMax = {};
                    cache.nuclei.boundariesMin = {};
                    cache.nuclei.boundariesMean = {};
                    cache.nuclei.boundariesMedian = {};
                end
                
                % Use the consistent nucleus labels from segmentation
                % NOTE: Don't check panel enable states here - labels already generated from segmentation
                if ~isempty(cache.nuclei.labels)
                    try
                        % Use the authoritative labels from segmentation
                        cache.nuclei.centroids3D = cache.nuclei.labels.centroids_3d;
                        cache.nuclei.labelsBySlice = cache.nuclei.labels.labels_by_slice;
                        
                        % For projections, we need to generate labels based on the projected labeled mask
                        if ndims(cache.nuclei.mask) == 3 && ~isempty(cache.nuclei.labels.labeled_mask)
                            labeled_mask = cache.nuclei.labels.labeled_mask;
                            
                            % Generate projection labels by projecting the labeled mask
                            proj_max_labeled = max(labeled_mask, [], 3);
                            proj_min_labeled = min(labeled_mask, [], 3);
                            proj_mean_labeled = round(mean(labeled_mask, 3));
                            proj_median_labeled = round(median(labeled_mask, 3));
                            
                            % Convert projected labeled masks to label coordinates
                            % Use the updated centroids that reflect the filtered nuclei
                            cache.nuclei.labelsProjectionMax = generateLabelsFromLabeledMask(proj_max_labeled, cache.nuclei.labels.centroids_3d);
                            cache.nuclei.labelsProjectionMin = generateLabelsFromLabeledMask(proj_min_labeled, cache.nuclei.labels.centroids_3d);
                            cache.nuclei.labelsProjectionMean = generateLabelsFromLabeledMask(proj_mean_labeled, cache.nuclei.labels.centroids_3d);
                            cache.nuclei.labelsProjectionMedian = generateLabelsFromLabeledMask(proj_median_labeled, cache.nuclei.labels.centroids_3d);
                        else
                            % For 2D data, projection labels are the same as slice labels
                            cache.nuclei.labelsProjectionMax = cache.nuclei.labelsBySlice;
                            cache.nuclei.labelsProjectionMin = cache.nuclei.labelsBySlice;
                            cache.nuclei.labelsProjectionMean = cache.nuclei.labelsBySlice;
                            cache.nuclei.labelsProjectionMedian = cache.nuclei.labelsBySlice;
                        end
                        
                    catch ME
                        warning('Failed to process nuclei labels: %s', ME.message);
                        cache.nuclei.centroids3D = [];
                        cache.nuclei.labelsBySlice = {};
                        cache.nuclei.labelsProjectionMax = {};
                        cache.nuclei.labelsProjectionMin = {};
                        cache.nuclei.labelsProjectionMean = {};
                        cache.nuclei.labelsProjectionMedian = {};
                    end
                else
                    % If no valid labels available, initialize empty label arrays
                    cache.nuclei.centroids3D = [];
                    cache.nuclei.labelsBySlice = {};
                    cache.nuclei.labelsProjectionMax = {};
                    cache.nuclei.labelsProjectionMin = {};
                    cache.nuclei.labelsProjectionMean = {};
                    cache.nuclei.labelsProjectionMedian = {};
                end
                
                end % End of the else block for valid mask
                
            catch ME
                warning('Failed to create nuclei cache: %s', ME.message);
                % Continue without nuclei cache
            end
        end
        
        % --- Cache Channel Data ---
        if isfield(handles, 'processedChannel') && ~isempty(handles.processedChannel)
            cache.channels = cell(1, length(handles.processedChannel));
            
            for k = 1:length(handles.processedChannel)
                if ~isempty(handles.processedChannel{k})
                    try
                        cache.channels{k} = struct();
                        img_3d = handles.processedChannel{k};
                        z_max = size(img_3d, 3);
                        
                        % Pre-compute all z-slices
                        cache.channels{k}.slices = cell(1, z_max);
                        for z = 1:z_max
                            handles.statusLabel.Text = ['Status: Creating Preview Cache - Channel ' num2str(k) ' Z-frame ' num2str(z) '/' num2str(z_max) '...'];
                            drawnow;
                            cache.channels{k}.slices{z} = img_3d(:,:,z);
                        end
                        
                        % Pre-compute all projections
                        cache.channels{k}.projectionMax = max(img_3d, [], 3);
                        cache.channels{k}.projectionMin = min(img_3d, [], 3);
                        cache.channels{k}.projectionMean = mean(img_3d, 3);
                        cache.channels{k}.projectionMedian = median(img_3d, 3);
                        
                        % Cache maxima data for this channel
                        if isfield(handles, 'maximaCoords') && length(handles.maximaCoords) >= k && ~isempty(handles.maximaCoords{k})
                            coords = handles.maximaCoords{k};
                            
                            % Pre-filter maxima by z-slice for fast lookup
                            cache.channels{k}.maximaBySlice = cell(1, z_max);
                            for z = 1:z_max
                                if size(coords, 2) >= 3
                                    slice_maxima = coords(coords(:,3) == z, :);
                                    cache.channels{k}.maximaBySlice{z} = slice_maxima;
                                else
                                    cache.channels{k}.maximaBySlice{z} = coords; % 2D coordinates
                                end
                            end
                            
                            % Store all maxima for projections
                            cache.channels{k}.allMaxima = coords;
                        else
                            % Initialize empty maxima data
                            cache.channels{k}.maximaBySlice = cell(1, z_max);
                            cache.channels{k}.allMaxima = [];
                        end
                        
                    catch ME
                        warning('Failed to cache channel %d: %s', k, ME.message);
                        cache.channels{k} = [];
                    end
                else
                    cache.channels{k} = [];
                end
            end
        end
        
        % --- Cache DIC Data ---
        if isfield(handles, 'rawDIC') && ~isempty(handles.rawDIC)
            try
                cache.dic = struct();
                img_3d = handles.rawDIC;
                z_max = size(img_3d, 3);
                
                % Pre-compute all z-slices
                cache.dic.slices = cell(1, z_max);
                for z = 1:z_max
                    handles.statusLabel.Text = ['Status: Creating Preview Cache - DIC Z-frame ' num2str(z) '/' num2str(z_max) '...'];
                    drawnow;
                    cache.dic.slices{z} = img_3d(:,:,z);
                end
                
                % Pre-compute all projections
                cache.dic.projectionMax = max(img_3d, [], 3);
                cache.dic.projectionMin = min(img_3d, [], 3);
                cache.dic.projectionMean = mean(img_3d, 3);
                cache.dic.projectionMedian = median(img_3d, 3);
                
            catch ME
                warning('Failed to cache DIC data: %s', ME.message);
                cache.dic = [];
            end
        end
        
        % Store cache in handles first
        handles.previewCache = cache;
        guidata(handles.fig, handles);
        
        % --- Cache Analysis Data (after cache is stored) ---
        try
            handles.statusLabel.Text = 'Status: Analyzing Nuclei Composition...';
            drawnow;
            
            analysis_results = snap_helpers.analyzeNucleiSignalComposition(handles);
            handles.previewCache.analysis = analysis_results;
            
            
            % Save updated cache with analysis
            guidata(handles.fig, handles);
            
        catch ME
            warning('Failed to cache analysis data: %s', ME.message);
            handles.previewCache.analysis = [];
            guidata(handles.fig, handles);
        end
        
        handles.statusLabel.Text = 'Status: Preview Cache Created';
        drawnow;
        
    catch ME
        handles.statusLabel.Text = 'Status: Cache Creation Failed';
        warning('Failed to create preview cache: %s', ME.message);
        % Continue without cache - will fall back to real-time computation
    end
end

function labels_cell = generateLabelsFromLabeledMask(labeled_mask, centroids_3d)
    % Generate label coordinates from a labeled mask using consistent centroids
    
    labels_cell = {[]};
    
    if isempty(labeled_mask) || isempty(centroids_3d)
        return;
    end
    
    % Find unique labels (excluding 0 which is background)
    unique_labels = unique(labeled_mask(:));
    unique_labels = unique_labels(unique_labels > 0);
    
    if isempty(unique_labels)
        return;
    end
    
    % Create label coordinates for the projection
    label_coords = [];
    for i = 1:length(unique_labels)
        label_id = unique_labels(i);
        if label_id <= size(centroids_3d, 1)
            centroid = centroids_3d(label_id, :);
            % For projections, use X,Y coordinates and the nucleus ID
            label_coords = [label_coords; centroid(1), centroid(2), label_id];
        end
    end
    
    labels_cell{1} = label_coords;
end
