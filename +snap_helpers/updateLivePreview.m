function handles = updateLivePreview(fig_handle)
    handles = guidata(fig_handle);
    
    if isfield(handles, 'isProcessingPreview') && handles.isProcessingPreview
        return;
    end
    
    % Enable abort button and reset abort flag
    handles.isProcessingPreview = true;
    handles.abortRequested = false;
    handles.abortButton.Enable = 'on';
    handles.updateLivePreviewButton.Enable = 'off'; % Disable update during processing
    guidata(fig_handle, handles);
    
    handles.statusLabel.Text = 'Status: Loading...';
    handles.statusLabel.FontColor = [0.8, 0, 0];
    drawnow;
    try
        numActiveChannels = str2double(handles.numChanDrop.Value);
        handles = normalizeRuntimeCache(handles, numActiveChannels);

        % Clear preview cache before processing new data to ensure fresh results
        snap_helpers.clearPreviewCache(handles);
        
        [handles.rawDIC, handles.runtimeCache.dic] = loadVolumeWithCache( ...
            handles.runtimeCache.dic, handles.dicPathText.Value, 'DIC');
        [handles.rawNuclei, handles.runtimeCache.nuclei] = loadVolumeWithCache( ...
            handles.runtimeCache.nuclei, handles.nucPathText.Value, 'Nuclei');

        handles.rawChannel = cell(1, numActiveChannels);
        handles.processedChannel = resizeCellArray(getFieldOrDefault(handles, 'processedChannel', {}), numActiveChannels);
        handles.maximaCoords = resizeCellArray(getFieldOrDefault(handles, 'maximaCoords', {}), numActiveChannels);
        handles.gaussFitResults = resizeCellArray(getFieldOrDefault(handles, 'gaussFitResults', {}), numActiveChannels);
        for k = 1:numActiveChannels
            chPath = handles.channelPathTexts(k).Value;
            [handles.rawChannel{k}, handles.runtimeCache.channels{k}] = loadVolumeWithCache( ...
                handles.runtimeCache.channels{k}, chPath, sprintf('Ch. %d', k));
        end

        % First, update the controls to populate dropdowns with all possible items
        handles.statusLabel.Text = 'Status: Updating UI...';
        drawnow;
        guidata(fig_handle, handles); % Save handles with raw data before updating controls
        snap_helpers.updateControls(fig_handle);
        handles = guidata(fig_handle);

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

            latest_handles = guidata(fig_handle);
            if latest_handles.abortRequested
                handles = latest_handles;
                handles.statusLabel.Text = 'Status: Processing aborted by user';
                handles.statusLabel.FontColor = [0.8 0.2 0.2];
                handles.abortButton.Enable = 'off';
                handles.updateLivePreviewButton.Enable = 'on';
                handles.isProcessingPreview = false;
                guidata(fig_handle, handles);
                return;
            end

            chCache = handles.runtimeCache.channels{k};
            chPath = handles.channelPathTexts(k).Value;

            if k > numel(handles.rawChannel) || isempty(handles.rawChannel{k})
                handles.processedChannel{k} = [];
                handles.maximaCoords{k} = [];
                handles.gaussFitResults{k} = [];
                chCache.procSignature = '';
                chCache.maximaSignature = '';
                chCache.fitSignature = '';
                chCache.maximaCoordsRaw = [];
                chCache.maximaCoordsFiltered = [];
                chCache.fitResults = [];
                handles.runtimeCache.channels{k} = chCache;
                continue;
            end

            processSig = buildChannelProcessingSignature(handles, k, chPath);
            if isfield(chCache, 'procSignature') && strcmp(chCache.procSignature, processSig) && ...
                    isfield(chCache, 'processed') && ~isempty(chCache.processed)
                handles.processedChannel{k} = chCache.processed;
            else
                handles.statusLabel.Text = ['Status: Processing Ch. ' num2str(k) '...'];
                drawnow;
                handles.processedChannel{k} = snap_helpers.processImage(handles, k, handles.statusLabel);
                if isempty(handles.processedChannel{k})
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    handles.isProcessingPreview = false;
                    guidata(handles.fig, handles);
                    return;
                end
                chCache.processed = handles.processedChannel{k};
                chCache.procSignature = processSig;
                chCache.maximaSignature = '';
                chCache.fitSignature = '';
            end

            handles.gaussFitResults{k} = [];

            if handles.maximaEnabledChecks(k).Value
                maximaSig = buildChannelMaximaSignature(handles, k, processSig);
                if isfield(chCache, 'maximaSignature') && strcmp(chCache.maximaSignature, maximaSig) && ...
                        isfield(chCache, 'maximaCoordsRaw')
                    coords_raw = chCache.maximaCoordsRaw;
                else
                    handles.statusLabel.Text = ['Status: Finding Maxima Ch. ' num2str(k) '...'];
                    drawnow;
                    coords_raw = snap_helpers.detectMaxima(handles.processedChannel{k}, handles, k);
                    chCache.maximaCoordsRaw = coords_raw;
                    chCache.maximaSignature = maximaSig;
                    chCache.fitSignature = '';
                end

                final_coords = coords_raw;
                fit_results = [];

                if handles.gaussFitEnabledChecks(k).Value && ~isempty(coords_raw)
                    fitSig = buildChannelFitSignature(handles, k, maximaSig);
                    if isfield(chCache, 'fitSignature') && strcmp(chCache.fitSignature, fitSig) && ...
                            isfield(chCache, 'maximaCoordsFiltered') && isfield(chCache, 'fitResults')
                        final_coords = chCache.maximaCoordsFiltered;
                        fit_results = chCache.fitResults;
                    else
                        handles.statusLabel.Text = ['Status: Fitting Maxima Ch. ' num2str(k) '...'];
                        drawnow;

                        fitParams.gaussFitVoxelWindowSize = handles.gaussFitVoxelWindowSlider(k).Value;
                        fitParams.gaussFitBgCorrMethod = handles.gaussFitBgCorrMethodDrop(k).Value;
                        fitParams.gaussFitBgCorrWidth = handles.gaussFitBgCorrWidthEdit(k).Value;
                        fitParams.gaussFitPolyDegree = handles.gaussFitPolyDegreeEdit(k).Value;
                        fitParams.gaussFitMethod = handles.gaussFitMethodDrop(k).Value;
                        fitParams.gaussFitMaxIterations = handles.gaussFitMaxIterationsEdit(k).Value;
                        fitParams.gaussFitTolerance = handles.gaussFitToleranceEdit(k).Value;
                        fitParams.gaussFitRadialRadius = handles.gaussFitRadialRadiusEdit(k).Value;
                        fitParams.gaussFitPlotCheck = false;
                        fitParams.xySpacing = handles.xySpacingInputs(k).Value;
                        fitParams.zSpacing = handles.zSpacingInputs(k).Value;
                        fitParams = sanitizeFitParams(fitParams);
                        mode = handles.maximaModeDrops(k).Value;
                        is3D = ~strcmp(mode, 'On Z-Projection');

                        try
                            numCoords = size(coords_raw, 1);
                            chunkSize = 250;
                            numChunks = ceil(numCoords / chunkSize);
                            resultChunks = cell(1, numChunks);
                            chunkCount = 0;
                            abortedDuringFit = false;
                            for startIdx = 1:chunkSize:numCoords
                                latest_handles = guidata(fig_handle);
                                if latest_handles.abortRequested
                                    abortedDuringFit = true;
                                    break;
                                end
                                endIdx = min(startIdx + chunkSize - 1, numCoords);
                                handles.statusLabel.Text = sprintf('Status: Fitting Maxima Ch. %d (%d/%d)...', k, endIdx, numCoords);
                                drawnow;
                                chunkCount = chunkCount + 1;
                                chunkCoords = coords_raw(startIdx:endIdx, :);
                                [resultChunks{chunkCount}, chunkAborted] = snap_helpers.fitGaussians( ...
                                    handles.rawChannel{k}, chunkCoords, fitParams, is3D, @() isAbortRequested(fig_handle));
                                if chunkAborted
                                    abortedDuringFit = true;
                                    break;
                                end
                            end

                            if abortedDuringFit
                                handles = guidata(fig_handle);
                                handles.statusLabel.Text = 'Status: Processing aborted by user';
                                handles.statusLabel.FontColor = [0.8 0.2 0.2];
                                handles.abortButton.Enable = 'off';
                                handles.updateLivePreviewButton.Enable = 'on';
                                handles.isProcessingPreview = false;
                                guidata(fig_handle, handles);
                                return;
                            end

                            results = [resultChunks{1:chunkCount}];
                            [filtered_results, filter_mask] = snap_helpers.applyFitFiltering(results, k, handles);
                            final_coords = coords_raw(filter_mask, :);
                            fit_results = filtered_results;
                            chCache.fitSignature = fitSig;
                            chCache.maximaCoordsFiltered = final_coords;
                            chCache.fitResults = fit_results;
                        catch ME
                            warning('Gaussian fitting failed for channel %d: %s', k, ME.message);
                            final_coords = coords_raw;
                            fit_results = [];
                            chCache.fitSignature = '';
                            chCache.maximaCoordsFiltered = final_coords;
                            chCache.fitResults = fit_results;
                        end
                    end
                else
                    chCache.fitSignature = '';
                    chCache.maximaCoordsFiltered = [];
                    chCache.fitResults = [];
                end

                handles.maximaCoords{k} = final_coords;
                handles.gaussFitResults{k} = fit_results;

                colorName = handles.maximaColorDrops(k).Value;
                count = size(final_coords, 1);
                maxima_counts_text{end+1} = sprintf('Channel %d (%s): %d', k, colorName, count);
            else
                handles.maximaCoords{k} = [];
                handles.gaussFitResults{k} = [];
                chCache.maximaSignature = '';
                chCache.fitSignature = '';
                chCache.maximaCoordsRaw = [];
                chCache.maximaCoordsFiltered = [];
                chCache.fitResults = [];
            end

            handles.runtimeCache.channels{k} = chCache;
        end
        
        % Build nuclei segmentation once and reuse across filtering/cache/export.
        processed_nuclei = [];
        nuclei_mask = [];
        nucleus_labels = [];
        if ~isempty(handles.rawNuclei)
            nucSig = buildNucleiSignature(handles, handles.nucPathText.Value);
            if isfield(handles.runtimeCache, 'nucleiSeg') && ...
                    isfield(handles.runtimeCache.nucleiSeg, 'signature') && ...
                    strcmp(handles.runtimeCache.nucleiSeg.signature, nucSig)
                processed_nuclei = handles.runtimeCache.nucleiSeg.processed;
                nuclei_mask = handles.runtimeCache.nucleiSeg.mask;
                nucleus_labels = handles.runtimeCache.nucleiSeg.labels;
            else
                handles.statusLabel.Text = 'Status: Processing Nuclei...';
                drawnow;
                processed_nuclei = snap_helpers.preprocessNucleiWithBgCorr(handles);
                if isempty(processed_nuclei)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    handles.isProcessingPreview = false;
                    guidata(fig_handle, handles);
                    return;
                end
                [nuclei_mask, nucleus_labels] = snap_helpers.segmentNuclei(processed_nuclei, handles);
                if isempty(nuclei_mask)
                    handles.statusLabel.Text = 'Status: Processing aborted by user';
                    handles.statusLabel.FontColor = [0.8 0.2 0.2];
                    handles.abortButton.Enable = 'off';
                    handles.updateLivePreviewButton.Enable = 'on';
                    handles.isProcessingPreview = false;
                    guidata(fig_handle, handles);
                    return;
                end
                handles.runtimeCache.nucleiSeg.signature = nucSig;
                handles.runtimeCache.nucleiSeg.processed = processed_nuclei;
                handles.runtimeCache.nucleiSeg.mask = nuclei_mask;
                handles.runtimeCache.nucleiSeg.labels = nucleus_labels;
            end
        else
            handles.runtimeCache.nucleiSeg.signature = '';
            handles.runtimeCache.nucleiSeg.processed = [];
            handles.runtimeCache.nucleiSeg.mask = [];
            handles.runtimeCache.nucleiSeg.labels = [];
        end

        handles.processedNuclei = processed_nuclei;
        handles.nucleiMask = nuclei_mask;
        handles.nucleusLabels = nucleus_labels;

        % Apply nuclei inclusion/exclusion filtering if enabled
        if handles.nucInclusionExclusionEnabledCheck.Value && handles.nucSegEnabledCheck.Value && ~isempty(nuclei_mask)
            handles.statusLabel.Text = 'Status: Applying Nuclei Filtering...';
            drawnow;
            apply_to = handles.nucInclusionExclusionApplyDrop.Value;
            mode = handles.nucInclusionExclusionModeDrop.Value;
            if strcmp(apply_to, 'All Channels')
                channels_to_filter = 1:numActiveChannels;
            else
                channel_num = sscanf(apply_to, 'Channel %d');
                if isempty(channel_num) || ~isfinite(channel_num) || channel_num < 1 || channel_num > numActiveChannels
                    warning('Invalid nuclei filter target "%s". Applying to all active channels.', char(string(apply_to)));
                    channels_to_filter = 1:numActiveChannels;
                else
                    channels_to_filter = channel_num;
                end
            end

            for ch = channels_to_filter
                if isfield(handles, 'maximaCoords') && ch <= numel(handles.maximaCoords) && ~isempty(handles.maximaCoords{ch})
                    validateDimensionalConsistency(handles.maximaCoords{ch}, nuclei_mask, ch);
                    [filtered_coords, keepMask] = snap_helpers.filterMaximaByNuclei(handles.maximaCoords{ch}, nuclei_mask, mode);
                    handles.maximaCoords{ch} = filtered_coords;
                    if ch <= length(handles.gaussFitResults) && ~isempty(handles.gaussFitResults{ch}) && ~isempty(keepMask)
                        if length(handles.gaussFitResults{ch}) == length(keepMask)
                            handles.gaussFitResults{ch} = handles.gaussFitResults{ch}(keepMask);
                        else
                            warning('Channel %d: Cannot sync fitResults with nuclei filter (length mismatch). Clearing fits.', ch);
                            handles.gaussFitResults{ch} = [];
                        end
                    end
                    colorName = handles.maximaColorDrops(ch).Value;
                    count = size(filtered_coords, 1);
                    for i = 1:length(maxima_counts_text)
                        if contains(maxima_counts_text{i}, ['Channel ' num2str(ch)])
                            maxima_counts_text{i} = sprintf('Channel %d (%s): %d', ch, colorName, count);
                            break;
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
        handles.isProcessingPreview = false;
        
        guidata(fig_handle, handles);

    catch ME
        handles.statusLabel.Text = 'Status: Error!';
        handles.statusLabel.FontColor = [0.8, 0, 0];
        
        % Disable abort button and re-enable update button even on error
        handles.abortButton.Enable = 'off';
        handles.updateLivePreviewButton.Enable = 'on';
        handles.isProcessingPreview = false;
        guidata(fig_handle, handles);
        
        warning('An error occurred during preview update: %s', ME.message);
        try
            uialert(handles.fig, ME.message, 'Preview Update Error');
        catch
            % Keep non-fatal in environments where uialert is unavailable.
        end
    end
    
    % Maxima detection functions have been moved to snap_helpers.detectMaxima()
    % Nuclei filtering function has been moved to snap_helpers.filterMaximaByNuclei()
    % to ensure SNAP and SNAP_batch use identical processing logic

    function handlesOut = normalizeRuntimeCache(handlesIn, numChannels)
        handlesOut = handlesIn;
        if ~isfield(handlesOut, 'runtimeCache') || ~isstruct(handlesOut.runtimeCache)
            handlesOut.runtimeCache = struct();
        end
        if ~isfield(handlesOut.runtimeCache, 'dic') || ~isstruct(handlesOut.runtimeCache.dic)
            handlesOut.runtimeCache.dic = struct('path', '', 'raw', []);
        end
        if ~isfield(handlesOut.runtimeCache, 'nuclei') || ~isstruct(handlesOut.runtimeCache.nuclei)
            handlesOut.runtimeCache.nuclei = struct('path', '', 'raw', []);
        end
        if ~isfield(handlesOut.runtimeCache, 'nucleiSeg') || ~isstruct(handlesOut.runtimeCache.nucleiSeg)
            handlesOut.runtimeCache.nucleiSeg = struct('signature', '', 'processed', [], 'mask', [], 'labels', []);
        end
        if ~isfield(handlesOut.runtimeCache, 'channels') || ~iscell(handlesOut.runtimeCache.channels)
            handlesOut.runtimeCache.channels = cell(1, numChannels);
        end
        if numel(handlesOut.runtimeCache.channels) < numChannels
            handlesOut.runtimeCache.channels{numChannels} = [];
        end
        for ci = 1:numChannels
            if isempty(handlesOut.runtimeCache.channels{ci}) || ~isstruct(handlesOut.runtimeCache.channels{ci})
                handlesOut.runtimeCache.channels{ci} = struct( ...
                    'path', '', ...
                    'raw', [], ...
                    'processed', [], ...
                    'procSignature', '', ...
                    'maximaCoordsRaw', [], ...
                    'maximaCoordsFiltered', [], ...
                    'fitResults', [], ...
                    'maximaSignature', '', ...
                    'fitSignature', '');
            end
        end
    end

    function [volumeData, cacheEntry] = loadVolumeWithCache(cacheEntry, filePath, label)
        if nargin < 3
            label = 'Image';
        end
        if isempty(filePath)
            volumeData = [];
            cacheEntry.path = '';
            cacheEntry.raw = [];
            cacheEntry.fileDatenum = [];
            cacheEntry.fileBytes = [];
            return;
        end
        fileInfo = dir(filePath);
        fileDatenum = [];
        fileBytes = [];
        if ~isempty(fileInfo)
            fileDatenum = fileInfo(1).datenum;
            fileBytes = fileInfo(1).bytes;
        end

        samePath = isfield(cacheEntry, 'path') && strcmp(cacheEntry.path, filePath);
        sameMeta = isfield(cacheEntry, 'fileDatenum') && isfield(cacheEntry, 'fileBytes') && ...
            isequal(cacheEntry.fileDatenum, fileDatenum) && isequal(cacheEntry.fileBytes, fileBytes);
        if samePath && sameMeta && isfield(cacheEntry, 'raw') && ~isempty(cacheEntry.raw)
            volumeData = cacheEntry.raw;
            return;
        end
        handles.statusLabel.Text = ['Status: Loading ' label '...'];
        drawnow;
        try
            volumeData = tiffreadVolume(filePath);
            cacheEntry.path = filePath;
            cacheEntry.raw = volumeData;
            cacheEntry.fileDatenum = fileDatenum;
            cacheEntry.fileBytes = fileBytes;
        catch ME
            warning('Failed to load %s (%s): %s', label, filePath, ME.message);
            volumeData = [];
            cacheEntry.path = filePath;
            cacheEntry.raw = [];
            cacheEntry.fileDatenum = [];
            cacheEntry.fileBytes = [];
        end
    end

    function outCells = resizeCellArray(inCells, n)
        outCells = cell(1, n);
        if isempty(inCells)
            return;
        end
        copyN = min(numel(inCells), n);
        outCells(1:copyN) = inCells(1:copyN);
    end

    function val = getFieldOrDefault(inStruct, fieldName, defaultVal)
        if isfield(inStruct, fieldName)
            val = inStruct.(fieldName);
        else
            val = defaultVal;
        end
    end

    function sig = buildChannelProcessingSignature(h, ch, chPath)
        s = struct();
        s.path = chPath;
        s.xySpacing = h.xySpacingInputs(ch).Value;
        s.zSpacing = h.zSpacingInputs(ch).Value;
        s.deconvEnabled = h.deconvEnabledChecks(ch).Value;
        s.deconvMethod = h.deconvMethodDrops(ch).Value;
        s.deconvPSFSource = h.deconvPSFSourceDrops(ch).Value;
        s.deconvPSFPath = h.deconvPSFPathTexts(ch).Value;
        s.deconvPSFSigmaXY = h.deconvPSFSigmaXYInputs(ch).Value;
        s.deconvPSFSigmaZ = h.deconvPSFSigmaZInputs(ch).Value;
        s.deconvPSFSizeXY = h.deconvPSFSizeXYInputs(ch).Value;
        s.deconvPSFSizeZ = h.deconvPSFSizeZInputs(ch).Value;
        s.deconvLRIterations = h.deconvLRIterationsInputs(ch).Value;
        s.deconvLRDamping = h.deconvLRDampingInputs(ch).Value;
        s.deconvWienerNSR = h.deconvWienerNSRInputs(ch).Value;
        s.deconvBlindIterations = h.deconvBlindIterationsInputs(ch).Value;
        s.deconvBlindUnderRelax = h.deconvBlindUnderRelaxInputs(ch).Value;
        s.preprocEnabled = h.preprocEnabledChecks(ch).Value;
        s.preprocMode = h.preprocessModeDrops(ch).Value;
        s.preprocProjection = h.preprocessProjectionDrops(ch).Value;
        s.preprocScale = h.preprocessScaleChecks(ch).Value;
        s.preprocMethod = h.preprocMethodDrops(ch).Value;
        s.preprocParam1 = h.preprocParam1Inputs(ch).Value;
        s.preprocParam2 = h.preprocParam2Inputs(ch).Value;
        s.waveletName = h.waveletNameDrops(ch).Value;
        s.preprocParam3 = h.preprocParam3Inputs(ch).Value;
        s.preprocParam4 = h.preprocParam4Inputs(ch).Value;
        s.nlmStrength = h.nlmFilterStrengthInputs(ch).Value;
        s.nlmSearch = h.nlmSearchWindowInputs(ch).Value;
        s.nlmComparison = h.nlmComparisonWindowInputs(ch).Value;
        s.preprocClip = h.preprocClipChecks(ch).Value;
        s.bgEnabled = h.bgCorrEnabledChecks(ch).Value;
        s.bgMode = h.bgCorrModeDrops(ch).Value;
        s.bgProjection = h.bgCorrProjectionDrops(ch).Value;
        s.bgScale = h.bgCorrScaleChecks(ch).Value;
        s.bgMethod = h.bgMethodDrops(ch).Value;
        s.bgParam = h.bgParamInputs(ch).Value;
        s.bgClip = h.bgCorrClipChecks(ch).Value;
        sig = jsonencode(s);
    end

    function sig = buildChannelMaximaSignature(h, ch, processSig)
        s = struct();
        s.processSig = processSig;
        s.maximaEnabled = h.maximaEnabledChecks(ch).Value;
        s.maximaMode = h.maximaModeDrops(ch).Value;
        s.maximaProjection = h.maximaProjectionDrops(ch).Value;
        s.maximaMethod = h.maximaMethodDrops(ch).Value;
        s.maximaNeighborhood = h.maximaNeighborhoodInputs(ch).Value;
        s.maximaScale = h.maximaScaleChecks(ch).Value;
        s.maximaH = h.hMaxInputs(ch).Value;
        s.logSigma = h.logSigmaInputs(ch).Value;
        s.logThreshold = h.logThresholdInputs(ch).Value;
        sig = jsonencode(s);
    end

    function sig = buildChannelFitSignature(h, ch, maximaSig)
        s = struct();
        s.maximaSig = maximaSig;
        s.fitEnabled = h.gaussFitEnabledChecks(ch).Value;
        s.fitMethod = h.gaussFitMethodDrop(ch).Value;
        s.fitWindow = h.gaussFitVoxelWindowSlider(ch).Value;
        s.fitBgMethod = h.gaussFitBgCorrMethodDrop(ch).Value;
        s.fitBgWidth = h.gaussFitBgCorrWidthEdit(ch).Value;
        s.fitPolyDegree = h.gaussFitPolyDegreeEdit(ch).Value;
        s.fitMaxIterations = h.gaussFitMaxIterationsEdit(ch).Value;
        s.fitTolerance = h.gaussFitToleranceEdit(ch).Value;
        s.fitRadialRadius = h.gaussFitRadialRadiusEdit(ch).Value;
        s.xySpacing = h.xySpacingInputs(ch).Value;
        s.zSpacing = h.zSpacingInputs(ch).Value;
        s.filterEnabled = h.fitFilterEnabledChecks(ch).Value;
        s.filterR2Enabled = h.fitFilterRSquaredEnabledChecks(ch).Value;
        s.filterR2Min = h.fitFilterRSquaredMinInputs(ch).Value;
        s.filterR2Max = h.fitFilterRSquaredMaxInputs(ch).Value;
        s.filterSigmaEnabled = h.fitFilterSigmaSumEnabledChecks(ch).Value;
        s.filterSigmaMin = h.fitFilterSigmaSumMinInputs(ch).Value;
        s.filterSigmaMax = h.fitFilterSigmaSumMaxInputs(ch).Value;
        s.filterAmpEnabled = h.fitFilterAmplitudeEnabledChecks(ch).Value;
        s.filterAmpMin = h.fitFilterAmplitudeMinInputs(ch).Value;
        s.filterAmpMax = h.fitFilterAmplitudeMaxInputs(ch).Value;
        s.filterIntensityEnabled = h.fitFilterIntensityEnabledChecks(ch).Value;
        s.filterIntensityMin = h.fitFilterIntensityMinInputs(ch).Value;
        s.filterIntensityMax = h.fitFilterIntensityMaxInputs(ch).Value;
        sig = jsonencode(s);
    end

    function sig = buildNucleiSignature(h, nucPath)
        s = struct();
        s.path = nucPath;
        s.xySpacing = h.nucXYSpacingInput.Value;
        s.zSpacing = h.nucZSpacingInput.Value;
        s.deconvEnabled = h.nucDeconvEnabledCheck.Value;
        s.deconvMethod = h.nucDeconvMethodDrop.Value;
        s.deconvPSFSource = h.nucDeconvPSFSourceDrop.Value;
        s.deconvPSFPath = h.nucDeconvPSFPathText.Value;
        s.deconvPSFSigmaXY = h.nucDeconvPSFSigmaXYInput.Value;
        s.deconvPSFSigmaZ = h.nucDeconvPSFSigmaZInput.Value;
        s.deconvPSFSizeXY = h.nucDeconvPSFSizeXYInput.Value;
        s.deconvPSFSizeZ = h.nucDeconvPSFSizeZInput.Value;
        s.deconvLRIterations = h.nucDeconvLRIterationsInput.Value;
        s.deconvLRDamping = h.nucDeconvLRDampingInput.Value;
        s.deconvWienerNSR = h.nucDeconvWienerNSRInput.Value;
        s.deconvBlindIterations = h.nucDeconvBlindIterationsInput.Value;
        s.deconvBlindUnderRelax = h.nucDeconvBlindUnderRelaxInput.Value;
        s.preprocEnabled = h.nucPreprocEnabledCheck.Value;
        s.preprocMode = h.nucPreprocessModeDrop.Value;
        s.preprocProjection = h.nucPreprocessProjectionDrop.Value;
        s.preprocScale = h.nucPreprocessScaleCheck.Value;
        s.preprocMethod = h.nucPreprocMethodDrop.Value;
        s.preprocParam1 = h.nucPreprocParam1Inputs.Value;
        s.preprocParam2 = h.nucPreprocParam2Inputs.Value;
        s.preprocParam3 = h.nucPreprocParam3Inputs.Value;
        s.preprocParam4 = h.nucPreprocParam4Inputs.Value;
        s.preprocWavelet = h.nucWaveletNameDrop.Value;
        s.preprocNlmStrength = h.nucNlmFilterStrengthInput.Value;
        s.preprocNlmSearch = h.nucNlmSearchWindowInput.Value;
        s.preprocNlmComparison = h.nucNlmComparisonWindowInput.Value;
        s.preprocClip = h.nucPreprocClipChecks.Value;
        s.bgEnabled = h.nucBgCorrEnabledCheck.Value;
        s.bgMode = h.nucBgCorrModeDrop.Value;
        s.bgProjection = h.nucBgCorrProjectionDrop.Value;
        s.bgScale = h.nucBgCorrScaleCheck.Value;
        s.bgMethod = h.nucBgMethodDrop.Value;
        s.bgParam = h.nucBgParamInput.Value;
        s.bgClip = h.nucBgCorrClipChecks.Value;
        s.segEnabled = h.nucSegEnabledCheck.Value;
        s.segMode = h.nucSegModeDrop.Value;
        s.segProjection = h.nucSegProjectionDrop.Value;
        s.segMainMethod = h.nucSegMainMethodDrop.Value;
        s.segSubMethod = h.nucSegSubMethodDrop.Value;
        s.segParam1 = h.nucSegParam1Input.Value;
        s.segAlgorithm = h.nucSegAlgorithmDrop.Value;
        s.segAlgParam1 = h.nucSegAlgParamInput.Value;
        s.segAlgParam2 = h.nucSegAlgParam2Input.Value;
        s.segAlgDefault1 = h.nucSegAlgParamDefaultCheck.Value;
        s.segAlgDefault2 = h.nucSegAlgParam2DefaultCheck.Value;
        s.filterEnabled = h.nucFilterEnabledCheck.Value;
        s.filterSizeEnabled = h.nucFilterSizeEnabledCheck.Value;
        s.filterMinSize = h.nucFilterMinSizeInput.Value;
        s.filterSizeUnit = h.nucFilterSizeUnitDrop.Value;
        s.filterCircularityEnabled = h.nucFilterCircularityEnabledCheck.Value;
        s.filterMinCircularity = h.nucFilterMinCircularityInput.Value;
        s.filterSolidityEnabled = h.nucFilterSolidityEnabledCheck.Value;
        s.filterMinSolidity = h.nucFilterMinSolidityInput.Value;
        s.excludeEdges = h.nucExcludeEdgesCheck.Value;
        sig = jsonencode(s);
    end

    function p = sanitizeFitParams(p)
        p.gaussFitVoxelWindowSize = sanitizeOddInteger(p.gaussFitVoxelWindowSize, 3, 21, 7);
        p.gaussFitBgCorrWidth = sanitizeInteger(p.gaussFitBgCorrWidth, 0, 10, 1);
        p.gaussFitPolyDegree = sanitizeInteger(p.gaussFitPolyDegree, 1, 3, 2);
        p.gaussFitMaxIterations = sanitizeInteger(p.gaussFitMaxIterations, 10, 1000, 200);
        p.gaussFitTolerance = sanitizePositive(p.gaussFitTolerance, 1e-12, 1, 1e-6);
        p.gaussFitRadialRadius = sanitizePositive(p.gaussFitRadialRadius, 0.5, 20, 3);
    end

    function out = sanitizeOddInteger(in, minVal, maxVal, fallback)
        out = sanitizeInteger(in, minVal, maxVal, fallback);
        if mod(out, 2) == 0
            if out >= maxVal
                out = out - 1;
            else
                out = out + 1;
            end
        end
        out = max(minVal, min(maxVal, out));
    end

    function out = sanitizeInteger(in, minVal, maxVal, fallback)
        out = fallback;
        if isnumeric(in) && isscalar(in) && isfinite(in)
            out = round(in);
        end
        out = max(minVal, min(maxVal, out));
    end

    function out = sanitizePositive(in, minVal, maxVal, fallback)
        out = fallback;
        if isnumeric(in) && isscalar(in) && isfinite(in) && in > 0
            out = in;
        end
        out = max(minVal, min(maxVal, out));
    end

    function tf = isAbortRequested(figRef)
        tf = false;
        try
            if isempty(figRef) || ~isvalid(figRef)
                tf = true;
                return;
            end
            h = guidata(figRef);
            tf = isstruct(h) && isfield(h, 'abortRequested') && logical(h.abortRequested);
        catch
            tf = true;
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
