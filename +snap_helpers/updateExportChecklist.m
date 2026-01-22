function updateExportChecklist(fig_handle)
% Update the export checklist based on what data is loaded and processed
    handles = guidata(fig_handle);
    
    % Save current checkbox states before updating (to preserve user selections)
    previousChecked = false(1, length(handles.exportItemChecks));
    for i = 1:length(handles.exportItemChecks)
        % Save current value regardless of visibility
        previousChecked(i) = handles.exportItemChecks(i).Value;
        handles.exportItemChecks(i).Visible = 'off';
    end
    
    itemIdx = 1; % Track which checkbox to use
    numChannels = str2double(handles.numChanDrop.Value);
    
    % --- Parameters (ALWAYS FIRST - Always available regardless of data) ---
    tooltipText = 'Export complete workflow parameters:\n  - All processing settings for deconvolution, preprocessing, background correction\n  - Segmentation parameters and thresholds\n  - Maxima detection settings\n  - Gaussian fitting configuration\n  - Filter thresholds and criteria\n  - Saved as MAT file for reproducibility';
    
    handles.exportItemChecks(itemIdx).Text = 'Parameters';
    handles.exportItemLabels{itemIdx} = 'Parameters';
    handles.exportItemTypes{itemIdx} = 'parameters';
    handles.exportItemChecks(itemIdx).Visible = 'on';
    handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx); % Preserve previous state
    handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
    itemIdx = itemIdx + 1;
    
    % --- Export Items in Requested Order ---
    % Order: Parameters, Nuclei Image, Nuclei Data, Channel Image, Channel Data, Cluster Signal Data
    
    % 1. Nuclei Image
    if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
        % Check what processing is enabled
        processingSteps = {};
        processingDescriptions = {};
        if isfield(handles, 'nucDeconvEnabledCheck') && handles.nucDeconvEnabledCheck.Value
            processingSteps{end+1} = 'deconvolution';
            processingDescriptions{end+1} = handles.nucDeconvMethodDrop.Value;
        end
        if handles.nucPreprocEnabledCheck.Value
            processingSteps{end+1} = 'pre-processing';
            processingDescriptions{end+1} = handles.nucPreprocMethodDrop.Value;
        end
        if isfield(handles, 'nucBgCorrEnabledCheck') && handles.nucBgCorrEnabledCheck.Value
            processingSteps{end+1} = 'background correction';
            processingDescriptions{end+1} = handles.nucBgMethodDrop.Value;
        end
        
        if isempty(processingSteps)
            tooltipText = 'Export raw nuclei image (no processing applied)';
        else
            tooltipText = sprintf('Export nuclei image with: %s', strjoin(processingSteps, ', '));
            for i = 1:length(processingDescriptions)
                tooltipText = sprintf('%s\n  - %s: %s', tooltipText, processingSteps{i}, processingDescriptions{i});
            end
        end
        
        handles.exportItemChecks(itemIdx).Text = 'Nuclei Image';
        handles.exportItemLabels{itemIdx} = 'Nuclei Image';
        handles.exportItemTypes{itemIdx} = 'image_nuclei';
        handles.exportItemChecks(itemIdx).Visible = 'on';
        % Preserve user's previous selection (defaults to false on first appearance)
        handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx);
        handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
        itemIdx = itemIdx + 1;
    end
    
    % 2. Nuclei Data
    if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
        tooltipText = 'Export nuclei segmentation results including:\n  - Binary masks\n  - Nucleus labels and IDs\n  - Morphological properties (area, circularity, centroids)\n  - Connected component information';
        
        handles.exportItemChecks(itemIdx).Text = 'Nuclei Data';
        handles.exportItemLabels{itemIdx} = 'Nuclei Data';
        handles.exportItemTypes{itemIdx} = 'nuclei_data';
        handles.exportItemChecks(itemIdx).Visible = 'on';
        % Preserve user's previous selection (defaults to false on first appearance)
        handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx);
        handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
        itemIdx = itemIdx + 1;
    end
    
    % 3. Channel Images
    hasChannelData = false;
    for ch = 1:numChannels
        if isfield(handles, 'rawChannel') && ch <= length(handles.rawChannel) && ~isempty(handles.rawChannel{ch})
            hasChannelData = true;
            % Check what processing is enabled for this channel
            processingSteps = {};
            processingDescriptions = {};
            if isfield(handles, 'deconvEnabledChecks') && handles.deconvEnabledChecks(ch).Value
                processingSteps{end+1} = 'deconvolution';
                processingDescriptions{end+1} = handles.deconvMethodDrops(ch).Value;
            end
            if handles.preprocEnabledChecks(ch).Value
                processingSteps{end+1} = 'pre-processing';
                processingDescriptions{end+1} = handles.preprocMethodDrops(ch).Value;
            end
            if handles.bgCorrEnabledChecks(ch).Value
                processingSteps{end+1} = 'background correction';
                processingDescriptions{end+1} = handles.bgMethodDrops(ch).Value;
            end
            
            if isempty(processingSteps)
                tooltipText = sprintf('Export raw Channel %d image (no processing applied)', ch);
            else
                tooltipText = sprintf('Export Channel %d image with: %s', ch, strjoin(processingSteps, ', '));
                for i = 1:length(processingDescriptions)
                    tooltipText = sprintf('%s\n  - %s: %s', tooltipText, processingSteps{i}, processingDescriptions{i});
                end
            end
            
            handles.exportItemChecks(itemIdx).Text = sprintf('Channel %d Image', ch);
            handles.exportItemLabels{itemIdx} = sprintf('Channel %d Image', ch);
            handles.exportItemTypes{itemIdx} = sprintf('image_channel_%d', ch);
            handles.exportItemChecks(itemIdx).Visible = 'on';
            % Preserve user's previous selection (defaults to false on first appearance)
            handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx);
            handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
            itemIdx = itemIdx + 1;
        end
    end
    
    % 4. Channel Data
    if hasChannelData
        tooltipText = 'Export channel analysis results including:\n  - Local maxima coordinates (X, Y, Z)\n  - Gaussian fitting parameters (amplitude, sigma, centers)\n  - Fit quality metrics (R-squared, integrated intensity)\n  - Background-corrected values\n  - Per-channel CSV files';
        
        handles.exportItemChecks(itemIdx).Text = 'Channel Data';
        handles.exportItemLabels{itemIdx} = 'Channel Data';
        handles.exportItemTypes{itemIdx} = 'channel_data';
        handles.exportItemChecks(itemIdx).Visible = 'on';
        % Preserve user's previous selection (defaults to false on first appearance)
        handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx);
        handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
        itemIdx = itemIdx + 1;
    end
    
    % 5. Cluster Signal Data (only if nuclei segmentation is enabled AND channels loaded)
    nucleiSegmentationEnabled = isfield(handles, 'nucSegEnabledCheck') && handles.nucSegEnabledCheck.Value;
    if hasChannelData && isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei) && nucleiSegmentationEnabled
        tooltipText = 'Export signal clustering analysis per nucleus:\n  - Number of signals per nucleus (each channel)\n  - Signal density (signals per unit area/volume)\n  - Mean signal properties within each nucleus\n  - Nucleus-to-signal mapping\n  - Cross-channel colocalization data';
        
        handles.exportItemChecks(itemIdx).Text = 'Cluster Signal Data';
        handles.exportItemLabels{itemIdx} = 'Cluster Signal Data';
        handles.exportItemTypes{itemIdx} = 'signal_composition';
        handles.exportItemChecks(itemIdx).Visible = 'on';
        % Preserve user's previous selection (defaults to false on first appearance)
        handles.exportItemChecks(itemIdx).Value = previousChecked(itemIdx);
        handles.exportItemChecks(itemIdx).Tooltip = tooltipText;
        itemIdx = itemIdx + 1;
    end
    
    % Store the number of active export items (always at least 1 for Parameters)
    handles.numExportItems = itemIdx - 1;
    
    guidata(fig_handle, handles);
end
