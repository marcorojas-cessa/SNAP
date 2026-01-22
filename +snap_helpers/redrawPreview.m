function redrawPreview(fig_handle, preview_idx)
% Redraws a single preview panel based on its current settings.
% Uses cached data for fast rendering when available.
%
% COORDINATE CONVENTION:
%   Data stored in ARRAY CONVENTION [row, col, slice]
%   Display uses CARTESIAN CONVENTION for plotting
%   
%   Conversion for display:
%   - XData (horizontal) = col (2nd dimension)
%   - YData (vertical) = row (1st dimension)
%
% NOTE: imagesc() automatically displays correctly (rows→Y, cols→X)
    handles = guidata(fig_handle);
    
    content = handles.previewContentDrops(preview_idx).Value;
    mode = handles.previewModeDrops(preview_idx).Value;
    ax = handles.previewAxes(preview_idx);
    cla(ax, 'reset');
    ax.XTick = []; ax.YTick = [];

    img_3d = [];
    source_channel_idx = 0; % To track which channel's maxima to show
    
    % --- Data Retrieval (with caching) ---
    if strcmp(content, 'DIC')
        if isfield(handles, 'rawDIC')
            img_3d = handles.rawDIC;
        end
    elseif strcmp(content, 'Nuclei')
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
            raw_nuclei_3d = handles.rawNuclei;
            
            % Use cached processed data if available, otherwise preprocess
            if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ...
               isfield(handles.previewCache.nuclei, 'processed') && ~isempty(handles.previewCache.nuclei.processed)
                img_3d = handles.previewCache.nuclei.processed;
            else
                % Fallback to real-time preprocessing if cache not available
                img_3d = snap_helpers.preprocessNucleiWithBgCorr(handles);
            end
            
            % If preprocessing resulted in 2D projection but raw data is 3D, 
            % we need to handle Z-stack mode differently
            if ndims(raw_nuclei_3d) == 3 && size(raw_nuclei_3d, 3) > 1 && ...
               (ndims(img_3d) < 3 || size(img_3d, 3) == 1)
                % Raw data is 3D but processed data is 2D (projection mode)
                % For Z-stack availability, use the raw data dimensions
                img_3d_for_zstack_check = raw_nuclei_3d;
            else
                img_3d_for_zstack_check = img_3d;
            end
        end
    elseif startsWith(content, 'Channel')
        source_channel_idx = str2double(extractAfter(content, 'Channel '));
        if source_channel_idx > 0 && isfield(handles, 'processedChannel') && ...
           source_channel_idx <= numel(handles.processedChannel) && ...
           ~isempty(handles.processedChannel{source_channel_idx})
            img_3d = handles.processedChannel{source_channel_idx};
        end
    end
    
    % --- Display Logic ---
    if isempty(img_3d)
        handles.zSliders(preview_idx).Visible = 'off';
        handles.zLabels(preview_idx).Visible = 'off';
        handles.previewProjectionDrops(preview_idx).Visible = 'off';
        return;
    end

    % Determine Z-stack availability based on appropriate image data
    if strcmp(content, 'Nuclei') && exist('img_3d_for_zstack_check', 'var')
        z_max = size(img_3d_for_zstack_check, 3);
    else
        z_max = size(img_3d, 3);
    end
    img_to_display = [];
    overlay_slice = []; % For nuclei
    maxima_coords_2d = []; % For maxima
    
    % Z-Stack or Z-Projection (with caching)
    if strcmp(mode, 'Z-Stack') && z_max > 1
        handles.zSliders(preview_idx).Visible = 'on';
        handles.zLabels(preview_idx).Visible = 'on';
        handles.previewProjectionDrops(preview_idx).Visible = 'off';

        handles.zSliders(preview_idx).Limits = [1, z_max];
        z_val = round(handles.zSliders(preview_idx).Value);
        if z_val > z_max, z_val = z_max; end
        handles.zSliders(preview_idx).Value = z_val;
        
        handles.zLabels(preview_idx).Text = sprintf('Z: %d / %d', z_val, z_max);
        
        % Use cached slice if available
        if strcmp(content, 'DIC') && isfield(handles, 'previewCache') && isfield(handles.previewCache, 'dic')
            cache = handles.previewCache.dic;
            if z_val <= length(cache.slices)
                img_to_display = cache.slices{z_val};
            else
                img_to_display = img_3d(:,:,z_val); % Fallback
            end
        elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && source_channel_idx > 0 && source_channel_idx <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{source_channel_idx})
            cache = handles.previewCache.channels{source_channel_idx};
            if z_val <= length(cache.slices)
                img_to_display = cache.slices{z_val};
            else
                img_to_display = img_3d(:,:,z_val); % Fallback
            end
        else
            % Special handling for nuclei in Z-stack mode when preprocessing uses projection
            if strcmp(content, 'Nuclei') && exist('raw_nuclei_3d', 'var') && ...
               ndims(raw_nuclei_3d) == 3 && size(raw_nuclei_3d, 3) > 1 && ...
               (ndims(img_3d) < 3 || size(img_3d, 3) == 1)
                % Show raw nuclei slice instead of projected data
                img_to_display = double(raw_nuclei_3d(:,:,z_val));
            else
                img_to_display = img_3d(:,:,z_val);
            end
        end
        
        % Get maxima for this slice if applicable (with caching)
        if source_channel_idx > 0 && handles.showMaximaChecks(source_channel_idx).Value
            if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && source_channel_idx <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{source_channel_idx})
                cache = handles.previewCache.channels{source_channel_idx};
                if z_val <= length(cache.maximaBySlice)
                    maxima_coords_2d = cache.maximaBySlice{z_val};
                else
                    maxima_coords_2d = [];
                end
            else
                % Fallback to real-time computation
                if isfield(handles, 'maximaCoords') && source_channel_idx <= numel(handles.maximaCoords)
                    coords = handles.maximaCoords{source_channel_idx};
                else
                    coords = [];
                end
                if ~isempty(coords)
                    is_projection_data = all(coords(:,3) == 0.5);
                    if is_projection_data
                        maxima_coords_2d = coords; % Show projection maxima on all slices
                    else
                        maxima_coords_2d = coords(coords(:,3) == z_val, :);
                    end
                end
            end
        end
        
    else % Z-Projection or single-slice image
        handles.zSliders(preview_idx).Visible = 'off';
        handles.zLabels(preview_idx).Visible = 'off';
        handles.previewProjectionDrops(preview_idx).Visible = 'on';

        proj_type = handles.previewProjectionDrops(preview_idx).Value;
        
        % Use cached projection if available
        if strcmp(content, 'DIC') && isfield(handles, 'previewCache') && isfield(handles.previewCache, 'dic')
            cache = handles.previewCache.dic;
            switch proj_type
                case 'Max'
                    img_to_display = cache.projectionMax;
                case 'Min'
                    img_to_display = cache.projectionMin;
                case 'Mean'
                    img_to_display = cache.projectionMean;
                case 'Median'
                    img_to_display = cache.projectionMedian;
                otherwise
                    img_to_display = cache.projectionMax;
            end
        elseif isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && source_channel_idx > 0 && source_channel_idx <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{source_channel_idx})
            cache = handles.previewCache.channels{source_channel_idx};
            switch proj_type
                case 'Max'
                    img_to_display = cache.projectionMax;
                case 'Min'
                    img_to_display = cache.projectionMin;
                case 'Mean'
                    img_to_display = cache.projectionMean;
                case 'Median'
                    img_to_display = cache.projectionMedian;
                otherwise
                    img_to_display = cache.projectionMax;
            end
        else
            img_to_display = projectZ(img_3d, proj_type);
        end

        % Get all maxima for projections (with caching)
        if source_channel_idx > 0 && handles.showMaximaChecks(source_channel_idx).Value
            if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && source_channel_idx <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{source_channel_idx})
                cache = handles.previewCache.channels{source_channel_idx};
                maxima_coords_2d = cache.allMaxima;
            elseif isfield(handles, 'maximaCoords') && source_channel_idx <= numel(handles.maximaCoords)
                maxima_coords_2d = handles.maximaCoords{source_channel_idx};
            else
                maxima_coords_2d = [];
            end
        end
    end
    
    % --- Generate Nuclei Overlay for ALL Previews (with caching) ---
    overlay_slice = [];
    boundaries = [];
    
    if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei) && handles.nucShowSegCheck.Value
        % Use cached data if available
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ~isempty(handles.previewCache.nuclei)
            cache = handles.previewCache.nuclei;
            
            try
                if strcmp(mode, 'Z-Stack') && z_max > 1
                    z_val = round(handles.zSliders(preview_idx).Value);
                    if z_val <= length(cache.slices) && z_val <= length(cache.boundaries)
                        overlay_slice = cache.slices{z_val};
                        boundaries = cache.boundaries{z_val};
                    else
                        % Safe fallback to first slice
                        overlay_slice = cache.slices{1};
                        boundaries = cache.boundaries{1};
                    end
                else
                    % Use cached projection based on current projection type
                    proj_type = handles.previewProjectionDrops(preview_idx).Value;
                    switch proj_type
                        case 'Max'
                            if isfield(cache, 'projectionMax')
                                overlay_slice = cache.projectionMax;
                            else
                                overlay_slice = [];
                            end
                            if isfield(cache, 'boundariesMax')
                                boundaries = cache.boundariesMax;
                            else
                                boundaries = {};
                            end
                        case 'Min'
                            if isfield(cache, 'projectionMin')
                                overlay_slice = cache.projectionMin;
                            else
                                overlay_slice = [];
                            end
                            if isfield(cache, 'boundariesMin')
                                boundaries = cache.boundariesMin;
                            else
                                boundaries = {};
                            end
                        case 'Mean'
                            if isfield(cache, 'projectionMean')
                                overlay_slice = cache.projectionMean;
                            else
                                overlay_slice = [];
                            end
                            if isfield(cache, 'boundariesMean')
                                boundaries = cache.boundariesMean;
                            else
                                boundaries = {};
                            end
                        case 'Median'
                            if isfield(cache, 'projectionMedian')
                                overlay_slice = cache.projectionMedian;
                            else
                                overlay_slice = [];
                            end
                            if isfield(cache, 'boundariesMedian')
                                boundaries = cache.boundariesMedian;
                            else
                                boundaries = {};
                            end
                        otherwise
                            if isfield(cache, 'projectionMax')
                                overlay_slice = cache.projectionMax;
                            else
                                overlay_slice = [];
                            end
                            if isfield(cache, 'boundariesMax')
                                boundaries = cache.boundariesMax;
                            else
                                boundaries = {};
                            end
                    end
                end
                
                % Ensure boundaries is a cell array
                if ~iscell(boundaries)
                    boundaries = {};
                end
                
            catch ME
                % If cache access fails, fall back to real-time computation
                warning('Cache access failed for nuclei overlay: %s. Falling back to real-time computation.', ME.message);
                overlay_slice = [];
                boundaries = [];
            end
        end
        
        % Fallback to real-time computation if cache not available or failed
        if isempty(overlay_slice) || isempty(boundaries)
            try
                % Try to use cached segmentation first
                if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ...
                   isfield(handles.previewCache.nuclei, 'mask') && ~isempty(handles.previewCache.nuclei.mask)
                    nuc_mask = handles.previewCache.nuclei.mask;
                else
                    % Last resort: real-time computation
                    processed_nuc = snap_helpers.preprocessNucleiWithBgCorr(handles);
                    [nuc_mask, ~] = snap_helpers.segmentNuclei(processed_nuc, handles);
                end
                
                if strcmp(mode, 'Z-Stack') && z_max > 1
                    z_val = round(handles.zSliders(preview_idx).Value);
                    if ndims(nuc_mask) == 3 && z_val <= size(nuc_mask, 3)
                        overlay_slice = nuc_mask(:,:,z_val);
                    else % Projection mask or out of bounds
                        overlay_slice = nuc_mask;
                    end
                else
                    overlay_slice = projectZ(nuc_mask, 'Max'); % Max project the mask
                end
                
                if ~isempty(overlay_slice)
                    boundaries = bwboundaries(overlay_slice);
                else
                    boundaries = {};
                end
                
            catch ME
                warning('Failed to generate nuclei overlay: %s', ME.message);
                overlay_slice = [];
                boundaries = {};
            end
        end
    end
    
    % --- Render Everything ---
    if isempty(img_to_display)
        % No image to display - show message
        text(ax, 0.5, 0.5, 'No Data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 14);
        return;
    end
    
    imagesc(ax, img_to_display);
    axis(ax, 'image');
    colormap(ax, 'gray');
    hold(ax, 'on');
    
    % Render Nuclei Overlay (Purple)
    if ~isempty(overlay_slice) && ~isempty(boundaries) && iscell(boundaries)
        for k = 1:length(boundaries)
            if ~isempty(boundaries{k}) && size(boundaries{k}, 2) >= 2
                boundary = boundaries{k};
                plot(ax, boundary(:,2), boundary(:,1), 'm', 'LineWidth', 1.5);
            end
        end
    end
    
    % Render Nuclei Labels (ONLY if showing Nuclei content AND "Show Segmentation overlay" is checked)
    nuclei_labels = [];
    if strcmp(content, 'Nuclei') && isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei) && handles.nucShowSegCheck.Value
        % ALWAYS use cached labels for consistency across all preview types
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei') && ~isempty(handles.previewCache.nuclei)
            cache = handles.previewCache.nuclei;
            
            try
                if strcmp(mode, 'Z-Stack') && z_max > 1
                    % Z-Stack mode: show labels only for nuclei present in current slice
                    z_val = round(handles.zSliders(preview_idx).Value);
                    if isfield(cache, 'labelsBySlice') && z_val <= length(cache.labelsBySlice) && ~isempty(cache.labelsBySlice{z_val})
                        nuclei_labels = cache.labelsBySlice{z_val};
                    end
                else
                    % Projection mode: show all nuclei labels
                    proj_type = handles.previewProjectionDrops(preview_idx).Value;
                    switch proj_type
                        case 'Max'
                            if isfield(cache, 'labelsProjectionMax') && ~isempty(cache.labelsProjectionMax) && ~isempty(cache.labelsProjectionMax{1})
                                nuclei_labels = cache.labelsProjectionMax{1};
                            end
                        case 'Min'
                            if isfield(cache, 'labelsProjectionMin') && ~isempty(cache.labelsProjectionMin) && ~isempty(cache.labelsProjectionMin{1})
                                nuclei_labels = cache.labelsProjectionMin{1};
                            end
                        case 'Mean'
                            if isfield(cache, 'labelsProjectionMean') && ~isempty(cache.labelsProjectionMean) && ~isempty(cache.labelsProjectionMean{1})
                                nuclei_labels = cache.labelsProjectionMean{1};
                            end
                        case 'Median'
                            if isfield(cache, 'labelsProjectionMedian') && ~isempty(cache.labelsProjectionMedian) && ~isempty(cache.labelsProjectionMedian{1})
                                nuclei_labels = cache.labelsProjectionMedian{1};
                            end
                        otherwise
                            if isfield(cache, 'labelsProjectionMax') && ~isempty(cache.labelsProjectionMax) && ~isempty(cache.labelsProjectionMax{1})
                                nuclei_labels = cache.labelsProjectionMax{1};
                            end
                    end
                end
            catch ME
                warning('Failed to retrieve nuclei labels from cache: %s', ME.message);
                nuclei_labels = [];
            end
        else
            % No cache available - this means "Update Previews" needs to be run
            warning('No nuclei cache available for preview %d. Run "Update Previews" to see current segmentation results.', preview_idx);
            nuclei_labels = [];
        end
        
        % Display the nuclei labels
        if ~isempty(nuclei_labels) && size(nuclei_labels, 2) >= 3
            for i = 1:size(nuclei_labels, 1)
                x_pos = nuclei_labels(i, 1);
                y_pos = nuclei_labels(i, 2);
                label_num = nuclei_labels(i, 3);
                
                % Display the number as text
                text(ax, x_pos, y_pos, num2str(label_num), ...
                     'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end
    
    % Render Maxima Overlay
    if ~isempty(maxima_coords_2d)
        colorName = handles.maximaColorDrops(source_channel_idx).Value;
        % Use bright neon colors that are highly visible against dark backgrounds
        colorMap = containers.Map({'Red','Green','Blue','Yellow','Magenta','Cyan'}, ...
                                 {[1, 0.2, 0.2], [0.2, 1, 0.2], [0.3, 0.7, 1], [1, 1, 0.2], [1, 0.2, 1], [0.2, 1, 1]});
        markerColor = colorMap(colorName);

        % After filtering, maxima_coords_2d and gaussFitResults are synchronized
        % So we can use the loop index directly
        
        for i = 1:size(maxima_coords_2d, 1)
            point = maxima_coords_2d(i,:);  % point is [row, col, slice] (ARRAY CONVENTION)
            
            % Convert to Cartesian for display: [x, y] = [col, row]
            x_display = point(2);  % col → X
            y_display = point(1);  % row → Y
            
            % The index i corresponds directly to the fit result index
            % (both have been filtered together)
            original_idx = i;

                            if original_idx > 0 && handles.gaussFitEnabledChecks(source_channel_idx).Value
                % Check if the fit result is valid before making it clickable
                is_valid_fit = false;
                if original_idx <= length(handles.gaussFitResults{source_channel_idx})
                    fit_result = handles.gaussFitResults{source_channel_idx}(original_idx);
                    if ~isempty(fit_result) && isfield(fit_result, 'amplitude')
                        % For radial symmetry, be more permissive - just check if it has the required fields
                        if isfield(fit_result, 'fitMethod') && strcmp(fit_result.fitMethod, 'Radial Symmetry')
                            % For radial symmetry, check if we have center coordinates (even if NaN, still allow clicking)
                            if isfield(fit_result, 'center_x') && isfield(fit_result, 'center_y')
                                is_valid_fit = true;
                            end
                        else
                            % For Gaussian fits, check if amplitude is not NaN
                            if ~isnan(fit_result.amplitude)
                                is_valid_fit = true;
                            end
                        end
                    end
                end
                
                if is_valid_fit
                    % Gaussian fitting is enabled and successful for this channel, make it clickable with X marker
                    % Use CARTESIAN for display: XData=col, YData=row
                    h_line = line(ax, 'XData', x_display, 'YData', y_display, ...
                                  'Marker', 'x', 'Color', markerColor, ...
                                  'MarkerSize', 10, 'LineWidth', 1.5);
                    
                    % Set the callback
                    h_line.ButtonDownFcn = @(~,~) snap_helpers.displayStandardizedFitPlot(handles, source_channel_idx, original_idx);
                else
                    % Gaussian fitting is enabled but failed for this maximum, not clickable, use cross
                    % Use CARTESIAN for display: XData=col, YData=row
                    line(ax, 'XData', x_display, 'YData', y_display, ...
                          'Marker', '+', 'Color', markerColor, ...
                          'MarkerSize', 8, 'LineWidth', 1.0);
                end
            else
                % Gaussian fitting is disabled for this channel, use cross marker
                % Use CARTESIAN for display: XData=col, YData=row
                 line(ax, 'XData', x_display, 'YData', y_display, ...
                       'Marker', '+', 'Color', markerColor, ...
                       'MarkerSize', 8, 'LineWidth', 1.0);
            end
        end
    end
    
    % Render Maxima from Other Channels (Display on All Previews)
    numChannels = str2double(handles.numChanDrop.Value);
    for ch = 1:numChannels
        % Skip the source channel (already rendered above)
        if ch == source_channel_idx
            continue;
        end
        
        % Check if this channel has "Display on All Previews" enabled
        if handles.displayOnAllPreviewsChecks(ch).Value && handles.showMaximaChecks(ch).Value
            other_maxima_coords = [];
            
            % Get maxima coordinates for this channel based on current preview mode
            if strcmp(mode, 'Z-Stack') && z_max > 1
                z_val = round(handles.zSliders(preview_idx).Value);
                if isfield(handles, 'maximaCoords') && ch <= numel(handles.maximaCoords)
                    coords = handles.maximaCoords{ch};
                else
                    coords = [];
                end
                if ~isempty(coords)
                    is_projection_data = all(coords(:,3) == 0.5);
                    if is_projection_data
                        other_maxima_coords = coords; % Show projection maxima on all slices
                    else
                        other_maxima_coords = coords(coords(:,3) == z_val, :);
                    end
                end
            else
                % Z-Projection or single-slice image
                if isfield(handles, 'maximaCoords') && ch <= numel(handles.maximaCoords)
                    other_maxima_coords = handles.maximaCoords{ch};
                else
                    other_maxima_coords = [];
                end
            end
            
            % Render maxima from this channel
            if ~isempty(other_maxima_coords)
                colorName = handles.maximaColorDrops(ch).Value;
                % Use the same bright neon color map as source channel for consistency
                colorMap = containers.Map({'Red','Green','Blue','Yellow','Magenta','Cyan'}, ...
                                         {[1, 0.2, 0.2], [0.2, 1, 0.2], [0.3, 0.7, 1], [1, 1, 0.2], [1, 0.2, 1], [0.2, 1, 1]});
                markerColor = colorMap(colorName);
                
                % Always use crosses (+) for maxima from other channels (not clickable)
                for i = 1:size(other_maxima_coords, 1)
                    point = other_maxima_coords(i,:);  % [row, col, slice] (ARRAY CONVENTION)
                    
                    % Convert to Cartesian for display: [x, y] = [col, row]
                    x_disp = point(2);  % col → X
                    y_disp = point(1);  % row → Y
                    
                    % Other channels' maxima are always displayed as crosses and are not clickable
                    line(ax, 'XData', x_disp, 'YData', y_disp, ...
                          'Marker', '+', 'Color', markerColor, ...
                          'MarkerSize', 8, 'LineWidth', 1.0);
                end
            end
        end
    end
    
    hold(ax, 'off');
    
    ax.UserData.OriginalCLim = ax.CLim;
    snap_helpers.updateBrightness(handles.brightnessSliders(preview_idx), handles.fig, preview_idx);
end

function proj = projectZ(img, type)
    if ndims(img) < 3, proj = img; return; end
    switch type
        case 'Max', proj = max(img, [], 3);
        case 'Min', proj = min(img, [], 3);
        case 'Median', proj = median(img, 3);
        case 'Mean', proj = mean(img, 3);
        otherwise, proj = max(img, [], 3);
    end
end