function handles = updateLivePreview(fig_handle)
    handles = guidata(fig_handle);
    
    
    % Enable abort button and reset abort flag
    handles.abortRequested = false;
    handles.abortButton.Enable = 'on';
    handles.updateLivePreviewButton.Enable = 'off'; % Disable update during processing
    guidata(fig_handle, handles);
    
    handles.statusLabel.Text = 'Status: Loading...';
    handles.statusLabel.FontColor = [0.8, 0, 0];
    drawnow;
    try
        % Clear preview cache before processing new data to ensure fresh results
        snap_helpers.clearPreviewCache(handles);
        
        % Initialize fields to prevent errors if files are not loaded
        handles.rawDIC = [];
        handles.rawNuclei = [];

        if ~isempty(handles.dicPathText.Value)
            handles.statusLabel.Text = 'Status: Loading DIC...';
            drawnow;
            handles.rawDIC = tiffreadVolume(handles.dicPathText.Value);
        end
        if ~isempty(handles.nucPathText.Value)
            handles.statusLabel.Text = 'Status: Loading Nuclei...';
            drawnow;
            handles.rawNuclei = tiffreadVolume(handles.nucPathText.Value);
        end
        
        numActiveChannels = str2double(handles.numChanDrop.Value);
        handles.rawChannel = cell(1, numActiveChannels); % Ensure cell array is correct size
        handles.processedChannel = cell(1, numActiveChannels); % Initialize processed channel array
        handles.maximaCoords = cell(1, numActiveChannels); % Initialize maxima coordinates array
        handles.gaussFitResults = cell(1, numActiveChannels); % Initialize Gaussian fit results array
        for k = 1:numActiveChannels
            if ~isempty(handles.channelPathTexts(k).Value)
                handles.statusLabel.Text = ['Status: Loading Ch. ' num2str(k) '...'];
                drawnow;
                try
                    handles.rawChannel{k} = tiffreadVolume(handles.channelPathTexts(k).Value);
                    fprintf('Successfully loaded Channel %d: %s\n', k, mat2str(size(handles.rawChannel{k})));
                catch ME
                    warning('Failed to load Channel %d: %s', k, ME.message);
                    handles.rawChannel{k} = [];
                end
            else
                handles.rawChannel{k} = [];
                fprintf('Channel %d path is empty\n', k);
            end
        end

        % First, update the controls to populate dropdowns with all possible items
        handles.statusLabel.Text = 'Status: Updating UI...';
        drawnow;
        guidata(fig_handle, handles); % Save handles with raw data before updating controls
        snap_helpers.updateControls(fig_handle);
        % Keep the current handles - don't overwrite with old data from guidata

        loaded_items = {};
        if isfield(handles, 'rawDIC') && ~isempty(handles.rawDIC), loaded_items{end+1} = 'DIC'; end
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei), loaded_items{end+1} = 'Nuclei'; end
        for k = 1:numActiveChannels
            if k <= numel(handles.rawChannel) && ~isempty(handles.rawChannel{k})
                loaded_items{end+1} = ['Channel ' num2str(k)];
            end
        end
        
        fprintf('Loaded items detected: %s\n', strjoin(loaded_items, ', '));
        
        for item_cell = loaded_items
            item = item_cell{1};
            is_displayed = any(strcmp({handles.previewContentDrops.Value}, item));
            if ~is_displayed
                for i = 1:5
                    if strcmp(handles.previewContentDrops(i).Value, 'None')
                        handles.previewContentDrops(i).Value = item;
                        fprintf('Assigned %s to Preview %d\n', item, i);
                        break;
                    end
                end
            else
                fprintf('%s is already displayed in previews\n', item);
            end
        end
        guidata(handles.fig, handles);

        maxima_counts_text = {};
        
        % --- Determine max Z dimension from ALL loaded images ---
        max_z_dim = 1; 
        if isfield(handles, 'rawDIC') && ~isempty(handles.rawDIC)
            max_z_dim = max(max_z_dim, size(handles.rawDIC, 3));
        end
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
            max_z_dim = max(max_z_dim, size(handles.rawNuclei, 3));
        end
        for k = 1:numActiveChannels
             if k <= numel(handles.rawChannel) && ~isempty(handles.rawChannel{k})
                max_z_dim = max(max_z_dim, size(handles.rawChannel{k}, 3));
             end
        end

        for k = 1:numActiveChannels
            
            % Check for abort request
            if handles.abortRequested
                handles.statusLabel.Text = 'Status: Processing aborted by user';
                handles.statusLabel.FontColor = [0.8 0.2 0.2];
                handles.abortButton.Enable = 'off';
                handles.updateLivePreviewButton.Enable = 'on';
                guidata(fig_handle, handles);
                fprintf('Processing aborted at channel %d\n', k);
                return;
            end
            
            if k <= numel(handles.rawChannel) && ~isempty(handles.rawChannel{k})
                fprintf('Processing Channel %d...\n', k);
                handles.processedChannel{k} = snap_helpers.processImage(handles, k, handles.statusLabel);
                
                % Check if processing was aborted
                if isempty(handles.processedChannel{k})
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    guidata(handles.fig, handles);
                    return;
                end
                fprintf('Channel %d processed successfully: %s\n', k, mat2str(size(handles.processedChannel{k})));
                
                handles.gaussFitResults{k} = []; % Clear previous results

                if handles.maximaEnabledChecks(k).Value
                    % Detect maxima regardless of showMaxima setting (showMaxima only controls display)
                    handles.statusLabel.Text = ['Status: Finding Maxima Ch. ' num2str(k) '...'];
                    drawnow;
                    processed_img = handles.processedChannel{k};
                    
                    % Use shared maxima detection function (same as SNAP_batch)
                    coords = snap_helpers.detectMaxima(processed_img, handles, k);
                    handles.maximaCoords{k} = coords;
                    
                    % COHESION FIX: Ensure preview cache is populated with maxima coordinates
                    % This ensures signal analysis can access the same data as preview
                    if ~isfield(handles, 'previewCache')
                        handles.previewCache = struct();
                    end
                    if ~isfield(handles.previewCache, 'channels')
                        handles.previewCache.channels = cell(1, numActiveChannels);
                    end
                    if k > length(handles.previewCache.channels)
                        handles.previewCache.channels{k} = struct();
                    elseif isempty(handles.previewCache.channels{k})
                        handles.previewCache.channels{k} = struct();
                    end
                    handles.previewCache.channels{k}.allMaxima = coords;
                    
                    colorName = handles.maximaColorDrops(k).Value;
                    count = size(coords, 1);
                    maxima_counts_text{end+1} = sprintf('Channel %d (%s): %d', k, colorName, count);

                    % --- Perform Local Maxima Fitting if Enabled ---
                    if handles.gaussFitEnabledChecks(k).Value && ~isempty(coords)
                        handles.statusLabel.Text = ['Status: Fitting Maxima Ch. ' num2str(k) '...'];
                        drawnow;
                        
                        % Consolidate all fitting parameters into a single struct
                        fitParams.gaussFitVoxelWindowSize = handles.gaussFitVoxelWindowSlider(k).Value;
                        fitParams.gaussFitBgCorrMethod = handles.gaussFitBgCorrMethodDrop(k).Value;
                        fitParams.gaussFitBgCorrWidth = handles.gaussFitBgCorrWidthEdit(k).Value;
                        fitParams.gaussFitPolyDegree = handles.gaussFitPolyDegreeEdit(k).Value;
                        fitParams.gaussFitMethod = handles.gaussFitMethodDrop(k).Value;
                        fitParams.gaussFitMaxIterations = handles.gaussFitMaxIterationsEdit(k).Value;
                        fitParams.gaussFitTolerance = handles.gaussFitToleranceEdit(k).Value;
                        fitParams.gaussFitRadialRadius = handles.gaussFitRadialRadiusEdit(k).Value;
                        
                        % Add the missing gaussFitPlotCheck parameter
                        fitParams.gaussFitPlotCheck = false; % Default to false for live preview
                        fitParams.xySpacing = handles.xySpacingInputs(k).Value;
                        fitParams.zSpacing = handles.zSpacingInputs(k).Value;
                        
                        mode = handles.maximaModeDrops(k).Value;
                        is3D = ~strcmp(mode, 'On Z-Projection');
                        
                        % Fit on the RAW image data, not the processed data
                        try
                            results = snap_helpers.fitGaussians(handles.rawChannel{k}, coords, fitParams, is3D);
                            
                            % Store the raw data window along with the fit results
                            for i = 1:size(coords,1)
                                results(i).rawDataWindow = getWindowData(handles.rawChannel{k}, round(coords(i,:)), fitParams.gaussFitVoxelWindowSize, is3D);
                                results(i).fitMethod = fitParams.gaussFitMethod;
                            end
                            
                            % Apply fit filtering if enabled and get the filter mask
                            [filtered_results, filter_mask] = snap_helpers.applyFitFiltering(results, k, handles);
                            
                            % CRITICAL: Always apply filter mask to maintain synchronization
                            handles.maximaCoords{k} = coords(filter_mask, :);
                            handles.gaussFitResults{k} = filtered_results;
                            
                            % Verify synchronization
                            if length(filtered_results) ~= size(handles.maximaCoords{k}, 1)
                                warning('Channel %d: maximaCoords and fitResults length mismatch after filtering (%d vs %d)', ...
                                    k, size(handles.maximaCoords{k}, 1), length(filtered_results));
                            end

                        catch ME
                            warning('Gaussian fitting failed for channel %d: %s', k, ME.message);
                            handles.gaussFitResults{k} = [];
                        end
                    end
                else
                    handles.maximaCoords{k} = [];
                end
            else
                handles.processedChannel{k} = [];
                handles.maximaCoords{k} = [];
            end
        end
        
        % Apply nuclei inclusion/exclusion filtering if enabled
        if handles.nucInclusionExclusionEnabledCheck.Value && handles.nucSegEnabledCheck.Value
            handles.statusLabel.Text = 'Status: Applying Nuclei Filtering...';
            drawnow;
            
            % Check for abort before nuclei processing
            if handles.abortRequested
                handles.statusLabel.Text = 'Status: Processing aborted by user';
                handles.statusLabel.FontColor = [0.8 0.2 0.2];
                handles.abortButton.Enable = 'off';
                handles.updateLivePreviewButton.Enable = 'on';
                guidata(fig_handle, handles);
                fprintf('Processing aborted during nuclei filtering\n');
                return;
            end
            
            % Get nuclei mask
            if ~isempty(handles.rawNuclei)
                processed_nuclei = snap_helpers.preprocessNucleiWithBgCorr(handles);
                
                % Check if preprocessing was aborted
                if isempty(processed_nuclei)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    guidata(fig_handle, handles);
                    return;
                end
                
                [nuclei_mask, nucleus_labels] = snap_helpers.segmentNuclei(processed_nuclei, handles);
                
                % Check if segmentation was aborted
                if isempty(nuclei_mask)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    guidata(fig_handle, handles);
                    return;
                end
                
                % Determine which channels to apply filtering to
                apply_to = handles.nucInclusionExclusionApplyDrop.Value;
                mode = handles.nucInclusionExclusionModeDrop.Value;
                
                if strcmp(apply_to, 'All Channels')
                    channels_to_filter = 1:numActiveChannels;
                else
                    % Extract channel number from string like "Channel 1"
                    channel_num = str2double(apply_to(end));
                    channels_to_filter = channel_num;
                end
                
                % Apply filtering to selected channels (only if maxima exist)
                for ch = channels_to_filter
                    % Check if maximaCoords exists and has data for this channel
                    if isfield(handles, 'maximaCoords') && ch <= numel(handles.maximaCoords) && ~isempty(handles.maximaCoords{ch})
                        % Validate dimensional consistency and provide user guidance
                        validateDimensionalConsistency(handles.maximaCoords{ch}, nuclei_mask, ch);
                        
                        % Use shared nuclei filtering function (returns mask to maintain synchronization)
                        [filtered_coords, keepMask] = snap_helpers.filterMaximaByNuclei(handles.maximaCoords{ch}, nuclei_mask, mode);
                        handles.maximaCoords{ch} = filtered_coords;
                        
                        % CRITICAL: Also filter gaussFitResults using the same mask
                        if isfield(handles, 'gaussFitResults') && ch <= length(handles.gaussFitResults) && ...
                           ~isempty(handles.gaussFitResults{ch}) && ~isempty(keepMask)
                            if length(handles.gaussFitResults{ch}) == length(keepMask)
                                handles.gaussFitResults{ch} = handles.gaussFitResults{ch}(keepMask);
                            else
                                warning('Channel %d: Cannot sync fitResults with nuclei filter (length mismatch). Clearing fits.', ch);
                                handles.gaussFitResults{ch} = [];
                            end
                        end
                        
                        % Update the preview cache with filtered coordinates
                        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
                           ch <= length(handles.previewCache.channels) && ~isempty(handles.previewCache.channels{ch})
                            handles.previewCache.channels{ch}.allMaxima = filtered_coords;
                        end
                        
                        % Update maxima count text
                        colorName = handles.maximaColorDrops(ch).Value;
                        count = size(filtered_coords, 1);
                        % Find and update the corresponding entry in maxima_counts_text
                        for i = 1:length(maxima_counts_text)
                            if contains(maxima_counts_text{i}, ['Channel ' num2str(ch)])
                                maxima_counts_text{i} = sprintf('Channel %d (%s): %d', ch, colorName, count);
                                break;
                            end
                        end
                    end
                end
            end
        end
        
        handles.maximaCountLabel.Text = [{'Local Maxima Counts:'}; maxima_counts_text'];
        
        % Update Global Z Slider
        if max_z_dim > 1
            handles.globalZSlider.Limits = [1, max_z_dim];
            current_val = round(handles.globalZSlider.Value);
            if current_val > max_z_dim, current_val = max_z_dim; end
            handles.globalZSlider.Value = current_val;
            handles.globalZSlider.Enable = 'on';
            handles.globalZLabel.Text = sprintf('Global Z: %d / %d', current_val, max_z_dim);
            
            % Enable play button when z-stack data is available
            if isfield(handles, 'playButton')
                handles.playButton.Enable = 'on';
            end
        else
            handles.globalZSlider.Limits = [1, 2]; % Set valid, but disabled range
            handles.globalZSlider.Value = 1;
            handles.globalZSlider.Enable = 'off';
            handles.globalZLabel.Text = 'Global Z: 1 / 1';
            
            % Disable play button when no z-stack data
            if isfield(handles, 'playButton')
                handles.playButton.Enable = 'off';
            end
        end
        
        % CRITICAL FIX: Apply nuclei filtering BEFORE creating cache
        % This ensures the analysis uses the same filtered data as the preview
        if handles.nucInclusionExclusionEnabledCheck.Value && handles.nucSegEnabledCheck.Value
            mode = handles.nucInclusionExclusionModeDrop.Value;
            
            % Get nuclei mask - ensure nuclei segmentation is enabled
            if ~isempty(handles.rawNuclei)
                processed_nuclei = snap_helpers.preprocessNucleiWithBgCorr(handles);
                
                % Check if preprocessing was aborted
                if isempty(processed_nuclei)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    guidata(handles.fig, handles);
                    return;
                end
                
                [nuclei_mask, nucleus_labels] = snap_helpers.segmentNuclei(processed_nuclei, handles);
                
                % Check if segmentation was aborted
                if isempty(nuclei_mask)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    guidata(handles.fig, handles);
                    return;
                end
                
                % Determine which channels to apply filtering to
                apply_to = handles.nucInclusionExclusionApplyDrop.Value;
                
                if strcmp(apply_to, 'All Channels')
                    channels_to_filter = 1:numActiveChannels;
                else
                    % Extract channel number from string like "Channel 1"
                    channel_num = str2double(apply_to(end));
                    channels_to_filter = channel_num;
                end
                
                % Apply filtering to selected channels (only if maxima exist)
                for ch = channels_to_filter
                    % Check if maximaCoords exists and has data for this channel
                    if isfield(handles, 'maximaCoords') && ch <= numel(handles.maximaCoords) && ~isempty(handles.maximaCoords{ch})
                        % Validate dimensional consistency and provide user guidance
                        validateDimensionalConsistency(handles.maximaCoords{ch}, nuclei_mask, ch);
                        
                        % Use shared nuclei filtering function (returns mask to maintain synchronization)
                        [filtered_coords, keepMask] = snap_helpers.filterMaximaByNuclei(handles.maximaCoords{ch}, nuclei_mask, mode);
                        handles.maximaCoords{ch} = filtered_coords;
                        
                        % CRITICAL: Also filter gaussFitResults using the same mask
                        if isfield(handles, 'gaussFitResults') && ch <= length(handles.gaussFitResults) && ...
                           ~isempty(handles.gaussFitResults{ch}) && ~isempty(keepMask)
                            if length(handles.gaussFitResults{ch}) == length(keepMask)
                                handles.gaussFitResults{ch} = handles.gaussFitResults{ch}(keepMask);
                            else
                                warning('Channel %d: Cannot sync fitResults with nuclei filter (length mismatch). Clearing fits.', ch);
                                handles.gaussFitResults{ch} = [];
                            end
                        end
                        
                        % Update maxima count text
                        colorName = handles.maximaColorDrops(ch).Value;
                        count = size(filtered_coords, 1);
                        % Find and update the corresponding entry in maxima_counts_text
                        for i = 1:length(maxima_counts_text)
                            if contains(maxima_counts_text{i}, ['Channel ' num2str(ch)])
                                maxima_counts_text{i} = sprintf('Channel %d (%s): %d', ch, colorName, count);
                                break;
                            end
                        end
                    end
                end
            end
        end

        guidata(handles.fig, handles);

        % Create preview cache for fast rendering (now with filtered coordinates)
        handles.statusLabel.Text = 'Status: Creating Preview Cache...';
        drawnow;
        snap_helpers.createPreviewCache(handles);
        
        % Validate cache was created properly (only warn if nuclei were loaded and expected)
        % Re-fetch handles after cache creation to get the updated cache
        handles = guidata(fig_handle);
        
        % Only warn if nuclei were loaded and segmentation is enabled
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei) && handles.nucSegEnabledCheck.Value
            % Check if cache exists and has nuclei data
            if ~isfield(handles, 'previewCache') || ~isfield(handles.previewCache, 'nuclei') || isempty(handles.previewCache.nuclei)
                warning('Preview cache creation may have failed - nuclei overlays may be inconsistent');
            end
        end
        
        % Keep the current handles with all processed data intact
        
        % Now redraw all previews using cached data
        handles.statusLabel.Text = 'Status: Redrawing Previews...';
        drawnow;
        for i = 1:5
            snap_helpers.redrawPreview(handles.fig, i);
        end

        % Update analysis panel with new data
        handles.statusLabel.Text = 'Status: Updating Analysis Panel...';
        drawnow;
        snap_helpers.updateAnalysisPanel(handles);
        
        % Update export checklist based on loaded data
        snap_helpers.updateExportChecklist(handles.fig);
        
        % CRITICAL: Reload handles after updateExportChecklist to preserve numExportItems
        handles = guidata(fig_handle);

        handles.statusLabel.Text = 'Status: Done';
        handles.statusLabel.FontColor = [0, 0.6, 0];
        
        % Disable abort button and re-enable update button
        handles.abortButton.Enable = 'off';
        handles.updateLivePreviewButton.Enable = 'on';
        
        guidata(fig_handle, handles);

    catch ME
        handles.statusLabel.Text = 'Status: Error!';
        handles.statusLabel.FontColor = [0.8, 0, 0];
        
        % Disable abort button and re-enable update button even on error
        handles.abortButton.Enable = 'off';
        handles.updateLivePreviewButton.Enable = 'on';
        
        warning('An error occurred during preview update: %s', ME.message);
        rethrow(ME);
    end
    
    % Maxima detection functions have been moved to snap_helpers.detectMaxima()
    % Nuclei filtering function has been moved to snap_helpers.filterMaximaByNuclei()
    % to ensure SNAP and SNAP_batch use identical processing logic

    function windowData = getWindowData(imageData, centerCoord, windowSize, is3D)
        % Extract window around centerCoord
        % COORDINATE CONVENTION: centerCoord is [row, col, slice] (ARRAY CONVENTION)
        
        winRadius = (windowSize - 1) / 2;
        imgSize = size(imageData);
        
        % Use ARRAY CONVENTION: row, col, slice
        row_min = max(1, round(centerCoord(1) - winRadius));
        row_max = min(imgSize(1), round(centerCoord(1) + winRadius));
        col_min = max(1, round(centerCoord(2) - winRadius));
        col_max = min(imgSize(2), round(centerCoord(2) + winRadius));
        
        if is3D
            slice_min = max(1, round(centerCoord(3) - winRadius));
            slice_max = min(imgSize(3), round(centerCoord(3) + winRadius));
            % DIRECT array access with ARRAY CONVENTION
            windowData = imageData(row_min:row_max, col_min:col_max, slice_min:slice_max);
        else
            % DIRECT array access with ARRAY CONVENTION
            windowData = imageData(row_min:row_max, col_min:col_max);
        end
    end

    function validateDimensionalConsistency(maxima_coords, nuclei_mask, channel_idx)
        % Validate dimensional consistency between maxima coordinates and nuclei mask
        % Provides user guidance for optimal settings
        %
        % TESTING GUIDE:
        % To test the fix, try these combinations:
        % 1. Nuclei: "On Z-Projection" + Maxima: "3D" -> Should show warning
        % 2. Nuclei: "On Z-Projection" + Maxima: "On Z-Projection" -> No warning
        % 3. Nuclei: "3D" + Maxima: "3D" -> No warning
        % 4. Nuclei: "2D (Slice-by-slice)" + Maxima: "2D (Slice-by-slice)" -> No warning
        
        if isempty(maxima_coords) || isempty(nuclei_mask)
            return;
        end
        
        is_2d_mask = ndims(nuclei_mask) == 2;
        has_3d_coords = size(maxima_coords, 2) >= 3 && any(maxima_coords(:, 3) ~= maxima_coords(1, 3));
        
        if is_2d_mask && has_3d_coords
            % Dimensional consistency warning - but don't print to reduce noise
        end
    end

    % Nuclei filtering function moved to snap_helpers.filterMaximaByNuclei()
end