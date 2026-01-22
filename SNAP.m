function SNAP
    % SNAP - Spot & Nuclei Analysis Pipeline
    % Main function to create the UI and set up callbacks.
    
    % Ensure proper path for package access
    current_dir = fileparts(mfilename('fullpath'));
    if ~contains(path, current_dir)
        addpath(current_dir);
    end
    
    Nmax = 5;
    
    % --- Load Parameters ---
    % Always start with default parameters
    % User can manually load saved parameters using "Load Parameters" button
    lastUsed = snap_helpers.initializeParameters(Nmax);
    
    % --- Force Navigation States to Start at 0 ---
    % Always start all panels in Stage 0 regardless of saved state
    lastUsed.nucNavPanelIndex = 0; % Nuclei navigation
    for k = 1:Nmax
        lastUsed.navPanelIndex{k} = 0; % Channel navigation
    end
    
    % --- Migrate Old Fit Method Names ---
    % Update old fit method names to new names with "Gaussian" suffix
    if isfield(lastUsed, 'gaussFitMethod')
        for k = 1:length(lastUsed.gaussFitMethod)
            if strcmp(lastUsed.gaussFitMethod{k}, '1D (X,Y,Z)')
                lastUsed.gaussFitMethod{k} = '1D (X,Y,Z) Gaussian';
            elseif strcmp(lastUsed.gaussFitMethod{k}, '2D (XY) + 1D (Z)')
                lastUsed.gaussFitMethod{k} = '2D (XY) + 1D (Z) Gaussian';
            end
        end
    end
    
    % --- Create UI ---
    try
        handles = snap_helpers.createUI(Nmax, lastUsed);
    catch ME
        errordlg(['An error occurred during UI creation: ' ME.message], 'UI Creation Error');
        rethrow(ME);
    end
    
    % --- Initial UI State (BEFORE setting callbacks to prevent infinite recursion) ---
    guidata(handles.fig, handles); % Store handles before first call
    snap_helpers.updateControls(handles.fig); % Initialize control visibility and state
    handles = guidata(handles.fig); % Re-fetch handles after update
    
    % Initialize export checklist (makes "Parameters" export available from start)
    snap_helpers.updateExportChecklist(handles.fig);
    handles = guidata(handles.fig);
    
    % --- Initialize abort flag ---
    handles.abortRequested = false;
    
    % --- Set UI Callbacks (AFTER initial state is set) ---
    handles.numChanDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.updateLivePreviewButton.ButtonPushedFcn = @(src,evt) snap_helpers.updateLivePreview(handles.fig);
    handles.abortButton.ButtonPushedFcn = @(src,evt) abortProcessing(handles.fig);
    handles.globalZSlider.ValueChangedFcn = @(src,evt) snap_helpers.synchronizeZSliders(handles.fig);
    
    % Analysis panel callbacks
    handles.analysisCollapseButton.ButtonPushedFcn = @(src,evt) toggleAnalysisPanel(handles.fig);
    
    % Play/Pause button callbacks
    handles.playButton.ButtonPushedFcn = @(src,evt) snap_helpers.startZPlayback(handles.fig);
    handles.pauseButton.ButtonPushedFcn = @(src,evt) snap_helpers.stopZPlayback(handles.fig);
    
    % Export panel callbacks
    handles.exportAllSelectedButton.ButtonPushedFcn = @(src,evt) snap_helpers.exportData(handles.fig, 'selected');
    handles.exportDirButton.ButtonPushedFcn = @(src,evt) snap_helpers.exportData(handles.fig, 'browse');
    handles.exportSelectAllButton.ButtonPushedFcn = @(src,evt) selectAllExportItems(handles.fig);
    
    % Parameters will be saved automatically when the application closes

% Nuclei tab callbacks
handles.nucSegMainMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegSubMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegLocalAlgorithmDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

% Dynamic algorithm parameter callbacks
handles.nucSegAlgParamInput.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegAlgParam2Input.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegAlgParamDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegAlgParam2DefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

% Default checkbox callbacks for algorithm parameters
handles.nucSegBernsenContrastDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegMeanCDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegMedianCDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegMidGreyCDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegNiblackKDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegNiblackCDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegPhansalkarKDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegPhansalkarRDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegSauvolaKDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
handles.nucSegSauvolaRDefaultCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

    % Nuclei pre-processing callbacks
    handles.nucPreprocessModeDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucPreprocMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucPreprocClipChecks.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    
    % Nuclei Non-Local Means parameter callbacks
    handles.nucNlmFilterStrengthInput.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucNlmSearchWindowInput.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucNlmComparisonWindowInput.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    
    % Nuclei Wavelet Denoising parameter callbacks
    handles.nucWaveletNameDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucWaveletLevelInput.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucWaveletThresholdRuleDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucWaveletThresholdMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    
    
    % Nuclei background correction callbacks
    if isfield(handles, 'nucBgCorrModeDrop')
        handles.nucBgCorrModeDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    end
    if isfield(handles, 'nucBgMethodDrop')
        handles.nucBgMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    end
    if isfield(handles, 'nucBgCorrClipChecks')
        handles.nucBgCorrClipChecks.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    end
    
    % Nuclei segmentation mode callbacks
    handles.nucSegModeDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

    % Nuclei panel enable/disable callbacks - instantly update UI state
    if isfield(handles, 'nucDeconvEnabledCheck')
        handles.nucDeconvEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.nucDeconvMethodDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.nucDeconvPSFSourceDrop.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.nucDeconvPSFBrowseButton.ButtonPushedFcn = @(src,evt) browseNucleiPSFFile(handles.fig);
    end
    handles.nucPreprocEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    if isfield(handles, 'nucBgCorrEnabledCheck')
        handles.nucBgCorrEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    end
    handles.nucSegEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucFilterEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucFilterSizeEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucFilterCircularityEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucFilterSolidityEnabledCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    handles.nucExcludeEdgesCheck.ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
    
    % Nuclei navigation callbacks
    handles.nucUpButton.ButtonPushedFcn = @(src,evt) snap_helpers.navigateNucleiPanels(handles.fig, 'up');
    handles.nucDownButton.ButtonPushedFcn = @(src,evt) snap_helpers.navigateNucleiPanels(handles.fig, 'down');

    % File browsing callbacks
    handles.paramLoadButton.ButtonPushedFcn = @(src,evt) loadParameterFile(handles.fig);
    handles.dicBrowseButton.ButtonPushedFcn = @(src,evt) loadFileCallback(handles.fig, 'dic');
    handles.nucBrowseButton.ButtonPushedFcn = @(src,evt) loadFileCallback(handles.fig, 'nuc');

    for i = 1:Nmax
        handles.channelBrowseButtons(i).ButtonPushedFcn = @(src,evt) loadFileCallback(handles.fig, 'channel', i);
        
        % Navigation button callbacks
        handles.upButtons(i).ButtonPushedFcn = @(src,evt) snap_helpers.navigatePanels(handles.fig, i, 'up');
        handles.downButtons(i).ButtonPushedFcn = @(src,evt) snap_helpers.navigatePanels(handles.fig, i, 'down');
        
        % Callbacks that trigger a full UI state update
        if isfield(handles, 'deconvMethodDrops')
            handles.deconvMethodDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
            handles.deconvPSFSourceDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        end
        if isfield(handles, 'deconvPSFBrowseButtons')
            handles.deconvPSFBrowseButtons(i).ButtonPushedFcn = @(src,evt) browsePSFFile(handles.fig, i);
        end
        handles.preprocessModeDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.preprocMethodDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.preprocClipChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        
        % Channel Non-Local Means parameter callbacks
        handles.nlmFilterStrengthInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.nlmSearchWindowInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.nlmComparisonWindowInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        
        % Channel Wavelet Denoising parameter callbacks
        handles.waveletNameDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.waveletLevelInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.waveletThresholdRuleDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.waveletThresholdMethodDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        
        
        handles.bgCorrModeDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.bgMethodDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.bgCorrClipChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.maximaModeDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.maximaMethodDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

        handles.xySpacingInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.zSpacingInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        
        % Callbacks for Local Maxima Fitting controls
        handles.gaussFitBgCorrMethodDrop(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.gaussFitMethodDrop(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

        % Callbacks for Fit Filtering controls
    handles.fitFilterRSquaredEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterRSquaredMinInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterRSquaredMaxInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterSigmaSumEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterSigmaSumMinInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterSigmaSumMaxInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterAmplitudeEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterAmplitudeMinInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterAmplitudeMaxInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterIntensityEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterIntensityMinInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterIntensityMaxInputs(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);

        % Channel panel enable/disable callbacks - instantly update UI state
        if isfield(handles, 'deconvEnabledChecks')
            handles.deconvEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        end
        handles.preprocEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.bgCorrEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.maximaEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.gaussFitEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        handles.fitFilterEnabledChecks(i).ValueChangedFcn = @(src,evt) snap_helpers.updateControls(handles.fig);
        
        % Classification callbacks
        handles.classifyEnabledChecks(i).ValueChangedFcn = @(src,evt) updateClassificationControls(i, handles.fig);
        handles.classifyLoadButtons(i).ButtonPushedFcn = @(src,evt) loadClassifierForChannel(i, handles.fig);
        handles.classifyTrainButtons(i).ButtonPushedFcn = @(src,evt) SNAP_classify();
        handles.classifyApplyButtons(i).ButtonPushedFcn = @(src,evt) applyClassifierToChannel(i, handles.fig);
        handles.classifySelectFeaturesButtons(i).ButtonPushedFcn = @(src,evt) selectClassifierFeatures(i, handles.fig);
    end
    
    % Callbacks for preview controls
    for i = 1:5
        handles.previewContentDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.redrawPreview(handles.fig, i);
        handles.previewModeDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.redrawPreview(handles.fig, i);
        handles.zSliders(i).ValueChangedFcn = @(src,evt) snap_helpers.redrawPreview(handles.fig, i);
        handles.previewProjectionDrops(i).ValueChangedFcn = @(src,evt) snap_helpers.redrawPreview(handles.fig, i);
        handles.brightnessSliders(i).ValueChangedFcn = @(src,evt) snap_helpers.updateBrightness(src, handles.fig, i);
        handles.brightnessResets(i).ButtonPushedFcn = @(src,evt) snap_helpers.resetBrightness(src, handles.fig, i);
    end
    
    % Store final handles after all callbacks are set
    guidata(handles.fig, handles);
    
    % Set the close request function (no auto-save)
    handles.fig.CloseRequestFcn = @(src,evt) closeSNAP(src);
    
    % Add emergency close function (Ctrl+Alt+Q) in case normal close fails
    handles.fig.KeyPressFcn = @(src,evt) emergencyClose(src, evt);
    
end

function closeSNAP(fig_handle)
    % Close SNAP without auto-saving parameters
    % User should manually export parameters if they want to save them
    fprintf('Closing SNAP...\n');
    
    % Set a flag to prevent any processing during close
    try
        handles = guidata(fig_handle);
        handles.abortRequested = true; % Stop any ongoing processing
        handles.fig.Tag = 'closing'; % Mark figure as closing
        guidata(fig_handle, handles);
    catch
        % If we can't set the flag, continue anyway
    end
    
    try
        % Try to get handles for cleanup
        handles = guidata(fig_handle);
    
    % Clean up any active z-playback timer
        if isfield(handles, 'zPlaybackTimer') && isvalid(handles.zPlaybackTimer)
        try
                stop(handles.zPlaybackTimer);
                delete(handles.zPlaybackTimer);
                fprintf('Z-playback timer cleaned up.\n');
            catch
                % Ignore timer cleanup errors
        end
    end
    
        % Clean up any other timers
        timers_list = timerfind;
        if ~isempty(timers_list)
            try
                stop(timers_list);
                delete(timers_list);
        catch
                % Ignore timer cleanup errors
        end
    end
        
        catch
        % If we can't get handles, just continue to close
        fprintf('Could not get handles, proceeding with close.\n');
    end
    
    % Force close the figure immediately
    try
        % Force a UI update before closing
        drawnow;
        delete(fig_handle);
        fprintf('SNAP closed successfully.\n');
    catch ME
        fprintf('Warning: Error closing figure: %s\n', ME.message);
        % Force close using backup methods
        try
            closereq;
        catch
        try
                close all force;
        catch
                % Final attempt - force deletion
        try
                    if isvalid(fig_handle)
                        close(fig_handle, 'force');
        end
        catch
                    fprintf('Unable to close figure. Please close MATLAB manually.\n');
        end
    end
        end
    end
end

function emergencyClose(fig_handle, evt)
    % Emergency close function - bypasses parameter saving
    % Activated with Ctrl+Alt+Q
    if strcmp(evt.Key, 'q') && ~isempty(evt.Modifier) && ...
       length(evt.Modifier) >= 2 && ...
       any(strcmp(evt.Modifier, 'control')) && any(strcmp(evt.Modifier, 'alt'))
        
        fprintf('EMERGENCY CLOSE activated (Ctrl+Alt+Q)!\n');
        fprintf('Bypassing parameter saving and forcing close...\n');
        
        % Clean up timers first
        try
            handles = guidata(fig_handle);
            if isfield(handles, 'zPlaybackTimer') && isvalid(handles.zPlaybackTimer)
                stop(handles.zPlaybackTimer);
                delete(handles.zPlaybackTimer);
            end
        catch
            % Ignore errors during cleanup
        end
        
        % Force close immediately
        try
            delete(fig_handle);
        catch
            close all force;
        end
    end
end


% --- Analysis Panel Helper Functions ---

function toggleAnalysisPanel(fig_handle)
    % Toggles the analysis panel between collapsed and expanded states
    try
        handles = guidata(fig_handle);
        
        % Find the main grid layout
        mainGrid = handles.fig.Children;
        if isa(mainGrid, 'matlab.ui.container.GridLayout')
            grid = mainGrid;
        else
            % If there are multiple children, find the GridLayout
            for i = 1:length(mainGrid)
                if isa(mainGrid(i), 'matlab.ui.container.GridLayout')
                    grid = mainGrid(i);
                    break;
        end
    end
        end
        
        if handles.analysisCollapsed
            % Expand panel
            grid.ColumnWidth{3} = 300;  % Restore full width
            handles.analysisCollapseButton.Text = '◄ Collapse';
            handles.analysisTablePanel.Visible = 'on';
            handles.analysisInfoLabel.Visible = 'on';
            handles.analysisCollapsed = false;
        else
            % Collapse panel - keep button accessible
            grid.ColumnWidth{3} = 50;   % Wider than before to keep button clickable
            handles.analysisCollapseButton.Text = '►';  % Shorter text for narrow space
            handles.analysisTablePanel.Visible = 'off';
            handles.analysisInfoLabel.Visible = 'off';
            handles.analysisCollapsed = true;
        end
        
        guidata(fig_handle, handles);
    catch ME
        warning('Failed to toggle analysis panel: %s', ME.message);
        fprintf('Debug info: Error = %s\n', ME.message);
        if exist('handles', 'var') && isfield(handles, 'analysisCollapsed')
            fprintf('Current state: analysisCollapsed = %d\n', handles.analysisCollapsed);
        end
    end
end

function loadFileCallback(fig, fileType, varargin)
    % Local callback wrapper for loadFile function
    try
        snap_helpers.loadFile(fig, fileType, varargin{:});
    catch ME
        fprintf('Warning: Error loading file with package call: %s\n', ME.message);
        % Try alternative approach if package path fails
        try
            % Get the current file directory and add to path temporarily
            current_dir = fileparts(mfilename('fullpath'));
            helpers_dir = fullfile(current_dir, '+snap_helpers');
            if exist(helpers_dir, 'dir')
                addpath(helpers_dir);
                loadFile(fig, fileType, varargin{:});
                rmpath(helpers_dir);
                fprintf('Successfully loaded file using alternative path method\n');
            else
                error('Cannot find snap_helpers directory at: %s', helpers_dir);
            end
        catch ME2
            fprintf('Error: Failed to load file: %s\n', ME2.message);
            errordlg(sprintf('Failed to load file: %s\n\nPlease ensure you are running MATLAB from the correct directory:\n%s', ME2.message, current_dir), 'File Load Error');
        end
        end
    end

function browsePSFFile(fig_handle, channel_idx)
    % Browse for PSF file for deconvolution
    handles = guidata(fig_handle);
    [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'; '*.*', 'All Files (*.*)'}, 'Select PSF Image');
    if file ~= 0
        handles.deconvPSFPathTexts(channel_idx).Value = fullfile(path, file);
        guidata(fig_handle, handles);
        end
    end
    
function browseNucleiPSFFile(fig_handle)
    % Browse for PSF file for nuclei deconvolution
    handles = guidata(fig_handle);
    [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'; '*.*', 'All Files (*.*)'}, 'Select Nuclei PSF Image');
    if file ~= 0
        handles.nucDeconvPSFPathText.Value = fullfile(path, file);
        guidata(fig_handle, handles);
        end
    end

function selectAllExportItems(fig_handle)
    % Select all visible export items
    handles = guidata(fig_handle);
    for i = 1:handles.numExportItems
        if strcmp(handles.exportItemChecks(i).Visible, 'on')
            handles.exportItemChecks(i).Value = true;
        end
    end
    guidata(fig_handle, handles);
end

function abortProcessing(fig_handle)
    % Request abort of ongoing processing
    handles = guidata(fig_handle);
    handles.abortRequested = true;
    handles.abortButton.Enable = 'off';
    handles.statusLabel.Text = 'Status: Aborting...';
    handles.statusLabel.FontColor = [0.8 0.2 0.2];
    guidata(fig_handle, handles);
    fprintf('ABORT requested by user\n');
end

function loadParameterFile(fig_handle)
    % Load parameters from a saved .mat file and apply them to the UI
    handles = guidata(fig_handle);
    
    % Browse for parameter file
    [file, path] = uigetfile({'*.mat', 'MATLAB Files (*.mat)'}, ...
        'Select Parameter File');
    
    if file == 0
        return; % User cancelled
    end
    
    fullPath = fullfile(path, file);
    handles.paramLoadPathText.Value = fullPath;
    
    try
        % Load the parameter file
        paramData = load(fullPath);
        
        % Extract parameters from the file
        if isfield(paramData, 'lastUsed')
            loadedParams = paramData.lastUsed;
        elseif isfield(paramData, 'batchConfig')
            loadedParams = paramData.batchConfig.parameters;
        elseif isfield(paramData, 'paramData')
            loadedParams = paramData.paramData.parameters;
        else
            uialert(fig_handle, 'Could not find parameters in file. Expected lastUsed, batchConfig, or paramData structure.', 'Invalid Parameter File');
            return;
        end
        
        % Apply parameters to UI (updates handles.lastUsed)
        applyParametersToUI(handles, loadedParams);
        
        % Save handles after updating lastUsed
        guidata(fig_handle, handles);
        
        % Update UI state
        snap_helpers.updateControls(fig_handle);
        
        % Clear any existing preview cache since parameters changed
        snap_helpers.clearPreviewCache(fig_handle);
        
        uialert(fig_handle, sprintf('Parameters loaded successfully from:\n%s', file), ...
            'Parameters Loaded', 'Icon', 'success');
        
    catch ME
        uialert(fig_handle, sprintf('Failed to load parameters:\n%s', ME.message), ...
            'Load Error', 'Icon', 'error');
    end
    
    guidata(fig_handle, handles);
end

function applyParametersToUI(handles, params)
    % Apply loaded parameters to all UI elements
    % This function directly sets UI element values from loaded parameters
    
    try
        % === BASIC SETTINGS ===
        if isfield(params, 'numChannels') && isfield(handles, 'numChanDrop')
            handles.numChanDrop.Value = num2str(params.numChannels);
        end
        
        % === NUCLEI SPACING ===
        if isfield(params, 'nucXYSpacing') && isfield(handles, 'nucXYSpacingInput')
            handles.nucXYSpacingInput.Value = params.nucXYSpacing;
        end
        if isfield(params, 'nucZSpacing') && isfield(handles, 'nucZSpacingInput')
            handles.nucZSpacingInput.Value = params.nucZSpacing;
        end
        
        % === NUCLEI DECONVOLUTION ===
        if isfield(params, 'nucDeconvEnabled') && isfield(handles, 'nucDeconvEnabledCheck')
            handles.nucDeconvEnabledCheck.Value = params.nucDeconvEnabled;
        end
        if isfield(params, 'nucDeconvMethod') && isfield(handles, 'nucDeconvMethodDrop')
            handles.nucDeconvMethodDrop.Value = params.nucDeconvMethod;
        end
        if isfield(params, 'nucDeconvLRIterations') && isfield(handles, 'nucDeconvLRIterationsInput')
            handles.nucDeconvLRIterationsInput.Value = params.nucDeconvLRIterations;
        end
        if isfield(params, 'nucDeconvLRDamping') && isfield(handles, 'nucDeconvLRDampingInput')
            handles.nucDeconvLRDampingInput.Value = params.nucDeconvLRDamping;
        end
        if isfield(params, 'nucDeconvWienerNSR') && isfield(handles, 'nucDeconvWienerNSRInput')
            handles.nucDeconvWienerNSRInput.Value = params.nucDeconvWienerNSR;
        end
        if isfield(params, 'nucDeconvBlindIterations') && isfield(handles, 'nucDeconvBlindIterationsInput')
            handles.nucDeconvBlindIterationsInput.Value = params.nucDeconvBlindIterations;
        end
        if isfield(params, 'nucDeconvBlindUnderRelax') && isfield(handles, 'nucDeconvBlindUnderRelaxInput')
            handles.nucDeconvBlindUnderRelaxInput.Value = params.nucDeconvBlindUnderRelax;
        end
        if isfield(params, 'nucDeconvPSFSource') && isfield(handles, 'nucDeconvPSFSourceDrop')
            handles.nucDeconvPSFSourceDrop.Value = params.nucDeconvPSFSource;
        end
        if isfield(params, 'nucDeconvPSFFilePath') && isfield(handles, 'nucDeconvPSFPathText')
            handles.nucDeconvPSFPathText.Value = params.nucDeconvPSFFilePath;
        end
        if isfield(params, 'nucDeconvPSFSigmaXY') && isfield(handles, 'nucDeconvPSFSigmaXYInput')
            handles.nucDeconvPSFSigmaXYInput.Value = params.nucDeconvPSFSigmaXY;
        end
        if isfield(params, 'nucDeconvPSFSigmaZ') && isfield(handles, 'nucDeconvPSFSigmaZInput')
            handles.nucDeconvPSFSigmaZInput.Value = params.nucDeconvPSFSigmaZ;
        end
        if isfield(params, 'nucDeconvPSFSizeXY') && isfield(handles, 'nucDeconvPSFSizeXYInput')
            handles.nucDeconvPSFSizeXYInput.Value = params.nucDeconvPSFSizeXY;
        end
        if isfield(params, 'nucDeconvPSFSizeZ') && isfield(handles, 'nucDeconvPSFSizeZInput')
            handles.nucDeconvPSFSizeZInput.Value = params.nucDeconvPSFSizeZ;
        end
        
        % === NUCLEI PREPROCESSING ===
        if isfield(params, 'nucPreprocEnabled') && isfield(handles, 'nucPreprocEnabledCheck')
            handles.nucPreprocEnabledCheck.Value = params.nucPreprocEnabled;
        end
        if isfield(params, 'nucPreProcMode') && isfield(handles, 'nucPreprocessModeDrop')
            handles.nucPreprocessModeDrop.Value = params.nucPreProcMode;
        end
        if isfield(params, 'nucPreProcScale') && isfield(handles, 'nucPreprocessScaleCheck')
            handles.nucPreprocessScaleCheck.Value = params.nucPreProcScale;
        end
        if isfield(params, 'nucPreProcProjection') && isfield(handles, 'nucPreprocessProjectionDrop')
            handles.nucPreprocessProjectionDrop.Value = params.nucPreProcProjection;
        end
        if isfield(params, 'nucPreProcMethod') && isfield(handles, 'nucPreprocMethodDrop')
            handles.nucPreprocMethodDrop.Value = params.nucPreProcMethod;
        end
        if isfield(params, 'nucPreprocClipAtZero') && isfield(handles, 'nucPreprocClipChecks')
            handles.nucPreprocClipChecks.Value = params.nucPreprocClipAtZero;
        end
        
        % Nuclei preprocessing method-specific parameters
        if isfield(params, 'nucSmoothGaussianValue') && isfield(handles, 'nucPreprocParam1Inputs')
            handles.nucPreprocParam1Inputs.Value = params.nucSmoothGaussianValue;
        end
        if isfield(params, 'nucWaveletName') && isfield(handles, 'nucWaveletNameDrop')
            handles.nucWaveletNameDrop.Value = params.nucWaveletName;
        end
        if isfield(params, 'nucWaveletLevel') && isfield(handles, 'nucWaveletLevelInput')
            handles.nucWaveletLevelInput.Value = params.nucWaveletLevel;
        end
        if isfield(params, 'nucWaveletThresholdRule') && isfield(handles, 'nucWaveletThresholdRuleDrop')
            handles.nucWaveletThresholdRuleDrop.Value = params.nucWaveletThresholdRule;
        end
        if isfield(params, 'nucWaveletThresholdMethod') && isfield(handles, 'nucWaveletThresholdMethodDrop')
            handles.nucWaveletThresholdMethodDrop.Value = params.nucWaveletThresholdMethod;
        end
        if isfield(params, 'nucNlmFilterStrength') && isfield(handles, 'nucNlmFilterStrengthInput')
            handles.nucNlmFilterStrengthInput.Value = params.nucNlmFilterStrength;
        end
        if isfield(params, 'nucNlmSearchWindow') && isfield(handles, 'nucNlmSearchWindowInput')
            handles.nucNlmSearchWindowInput.Value = params.nucNlmSearchWindow;
        end
        if isfield(params, 'nucNlmComparisonWindow') && isfield(handles, 'nucNlmComparisonWindowInput')
            handles.nucNlmComparisonWindowInput.Value = params.nucNlmComparisonWindow;
        end
        
        % === NUCLEI BACKGROUND CORRECTION ===
        if isfield(params, 'nucBgCorrEnabled') && isfield(handles, 'nucBgCorrEnabledCheck')
            handles.nucBgCorrEnabledCheck.Value = params.nucBgCorrEnabled;
        end
        if isfield(params, 'nucBgCorrMode') && isfield(handles, 'nucBgCorrModeDrop')
            handles.nucBgCorrModeDrop.Value = params.nucBgCorrMode;
        end
        if isfield(params, 'nucBgMethod') && isfield(handles, 'nucBgMethodDrop')
            handles.nucBgMethodDrop.Value = params.nucBgMethod;
        end
        if isfield(params, 'nucBgParam') && isfield(handles, 'nucBgParamInput')
            handles.nucBgParamInput.Value = params.nucBgParam;
        end
        if isfield(params, 'nucBgCorrScale') && isfield(handles, 'nucBgCorrScaleCheck')
            handles.nucBgCorrScaleCheck.Value = params.nucBgCorrScale;
        end
        if isfield(params, 'nucBgCorrProjection') && isfield(handles, 'nucBgCorrProjectionDrop')
            handles.nucBgCorrProjectionDrop.Value = params.nucBgCorrProjection;
        end
        if isfield(params, 'nucBgCorrClipAtZero') && isfield(handles, 'nucBgCorrClipCheck')
            handles.nucBgCorrClipCheck.Value = params.nucBgCorrClipAtZero;
        end
        
        % === NUCLEI SEGMENTATION ===
        if isfield(params, 'nucSegEnabled') && isfield(handles, 'nucSegEnabledCheck')
            handles.nucSegEnabledCheck.Value = params.nucSegEnabled;
        end
        if isfield(params, 'nucSegMode') && isfield(handles, 'nucSegModeDrop')
            handles.nucSegModeDrop.Value = params.nucSegMode;
        end
        if isfield(params, 'nucSegProjection') && isfield(handles, 'nucSegProjectionDrop')
            handles.nucSegProjectionDrop.Value = params.nucSegProjection;
        end
        if isfield(params, 'nucSegMainMethod') && isfield(handles, 'nucSegMainMethodDrop')
            handles.nucSegMainMethodDrop.Value = params.nucSegMainMethod;
        end
        if isfield(params, 'nucSegSubMethod') && isfield(handles, 'nucSegSubMethodDrop')
            handles.nucSegSubMethodDrop.Value = params.nucSegSubMethod;
        end
        if isfield(params, 'nucSegAbsoluteThreshold') && isfield(handles, 'nucSegAbsoluteThresholdInput')
            handles.nucSegAbsoluteThresholdInput.Value = params.nucSegAbsoluteThreshold;
        end
        if isfield(params, 'nucSegStdMultiplier') && isfield(handles, 'nucSegStdMultiplierInput')
            handles.nucSegStdMultiplierInput.Value = params.nucSegStdMultiplier;
        end
        if isfield(params, 'nucSegOffset') && isfield(handles, 'nucSegOffsetInput')
            handles.nucSegOffsetInput.Value = params.nucSegOffset;
        end
        if isfield(params, 'nucSegLocalAlgorithm') && isfield(handles, 'nucSegLocalAlgorithmDrop')
            handles.nucSegLocalAlgorithmDrop.Value = params.nucSegLocalAlgorithm;
        end
        
        % Algorithm-specific parameters - map from named fields back to dynamic inputs
        if isfield(params, 'nucSegLocalAlgorithm') && isfield(handles, 'nucSegLocalAlgorithmDrop')
            algorithm = params.nucSegLocalAlgorithm;
            
            % Map algorithm-specific parameters to dynamic input
            switch algorithm
                case 'Bernsen'
                    if isfield(params, 'nucSegBernsenContrast') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegBernsenContrast;
                    end
                    if isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = 0;
                    end
                case 'Mean'
                    if isfield(params, 'nucSegMeanC') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegMeanC;
                    end
                    if isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = 0;
                    end
                case 'Median'
                    if isfield(params, 'nucSegMedianC') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegMedianC;
                    end
                    if isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = 0;
                    end
                case 'MidGrey'
                    if isfield(params, 'nucSegMidGreyC') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegMidGreyC;
                    end
                    if isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = 0;
                    end
                case 'Niblack'
                    if isfield(params, 'nucSegNiblackK') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegNiblackK;
                    end
                    if isfield(params, 'nucSegNiblackC') && isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = params.nucSegNiblackC;
                    end
                case 'Phansalkar'
                    if isfield(params, 'nucSegPhansalkarK') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegPhansalkarK;
                    end
                    if isfield(params, 'nucSegPhansalkarR') && isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = params.nucSegPhansalkarR;
                    end
                case 'Sauvola'
                    if isfield(params, 'nucSegSauvolaK') && isfield(handles, 'nucSegAlgParamInput')
                        handles.nucSegAlgParamInput.Value = params.nucSegSauvolaK;
                    end
                    if isfield(params, 'nucSegSauvolaR') && isfield(handles, 'nucSegAlgParam2Input')
                        handles.nucSegAlgParam2Input.Value = params.nucSegSauvolaR;
            end
        end
        
            % Apply default checkbox states
            switch algorithm
                case 'Bernsen'
                    if isfield(params, 'nucSegBernsenContrastDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegBernsenContrastDefault;
                    end
                    if isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = true;
                    end
                case 'Mean'
                    if isfield(params, 'nucSegMeanCDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegMeanCDefault;
                    end
                    if isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = true;
                    end
                case 'Median'
                    if isfield(params, 'nucSegMedianCDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegMedianCDefault;
                    end
                    if isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = true;
                    end
                case 'MidGrey'
                    if isfield(params, 'nucSegMidGreyCDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegMidGreyCDefault;
                    end
                    if isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = true;
                    end
                case 'Niblack'
                    if isfield(params, 'nucSegNiblackKDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegNiblackKDefault;
                    end
                    if isfield(params, 'nucSegNiblackCDefault') && isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = params.nucSegNiblackCDefault;
                    end
                case 'Phansalkar'
                    if isfield(params, 'nucSegPhansalkarKDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegPhansalkarKDefault;
                    end
                    if isfield(params, 'nucSegPhansalkarRDefault') && isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = params.nucSegPhansalkarRDefault;
                    end
                case 'Sauvola'
                    if isfield(params, 'nucSegSauvolaKDefault') && isfield(handles, 'nucSegAlgParamDefaultCheck')
                        handles.nucSegAlgParamDefaultCheck.Value = params.nucSegSauvolaKDefault;
                    end
                    if isfield(params, 'nucSegSauvolaRDefault') && isfield(handles, 'nucSegAlgParam2DefaultCheck')
                        handles.nucSegAlgParam2DefaultCheck.Value = params.nucSegSauvolaRDefault;
                    end
            end
        end
        
        if isfield(params, 'nucSegAlgParam1') && isfield(handles, 'nucSegAlgParamInput')
            handles.nucSegAlgParamInput.Value = params.nucSegAlgParam1;
        end
        if isfield(params, 'nucSegAlgParam2') && isfield(handles, 'nucSegAlgParam2Input')
            handles.nucSegAlgParam2Input.Value = params.nucSegAlgParam2;
        end
        if isfield(params, 'nucShowSeg') && isfield(handles, 'nucShowSegCheck')
            handles.nucShowSegCheck.Value = params.nucShowSeg;
        end
        
        % === NUCLEI FILTERING ===
        if isfield(params, 'nucFilterEnabled') && isfield(handles, 'nucFilterEnabledCheck')
            handles.nucFilterEnabledCheck.Value = params.nucFilterEnabled;
        end
        if isfield(params, 'nucFilterSizeEnabled') && isfield(handles, 'nucFilterSizeEnabledCheck')
            handles.nucFilterSizeEnabledCheck.Value = params.nucFilterSizeEnabled;
        end
        if isfield(params, 'nucFilterMinSize') && isfield(handles, 'nucFilterMinSizeInput')
            handles.nucFilterMinSizeInput.Value = params.nucFilterMinSize;
        end
        if isfield(params, 'nucFilterSizeUnit') && isfield(handles, 'nucFilterSizeUnitDrop')
            handles.nucFilterSizeUnitDrop.Value = params.nucFilterSizeUnit;
        end
        if isfield(params, 'nucFilterCircularityEnabled') && isfield(handles, 'nucFilterCircularityEnabledCheck')
            handles.nucFilterCircularityEnabledCheck.Value = params.nucFilterCircularityEnabled;
        end
        if isfield(params, 'nucFilterMinCircularity') && isfield(handles, 'nucFilterMinCircularityInput')
            handles.nucFilterMinCircularityInput.Value = params.nucFilterMinCircularity;
        end
        if isfield(params, 'nucFilterSolidityEnabled') && isfield(handles, 'nucFilterSolidityEnabledCheck')
            handles.nucFilterSolidityEnabledCheck.Value = params.nucFilterSolidityEnabled;
        end
        if isfield(params, 'nucFilterMinSolidity') && isfield(handles, 'nucFilterMinSolidityInput')
            handles.nucFilterMinSolidityInput.Value = params.nucFilterMinSolidity;
        end
        if isfield(params, 'nucExcludeEdges') && isfield(handles, 'nucExcludeEdgesCheck')
            handles.nucExcludeEdgesCheck.Value = params.nucExcludeEdges;
        end
        
        % === NUCLEI INCLUSION/EXCLUSION ===
        if isfield(params, 'nucInclusionExclusionEnabled') && isfield(handles, 'nucInclusionExclusionEnabledCheck')
            handles.nucInclusionExclusionEnabledCheck.Value = params.nucInclusionExclusionEnabled;
        end
        if isfield(params, 'nucInclusionExclusionMode') && isfield(handles, 'nucInclusionExclusionModeDrop')
            handles.nucInclusionExclusionModeDrop.Value = params.nucInclusionExclusionMode;
        end
        if isfield(params, 'nucInclusionExclusionApplyTo') && isfield(handles, 'nucInclusionExclusionApplyDrop')
            handles.nucInclusionExclusionApplyDrop.Value = params.nucInclusionExclusionApplyTo;
        end
        
        % === CHANNEL PARAMETERS ===
        numChannels = params.numChannels;
        for ch = 1:min(numChannels, handles.Nmax)
            try
            % Spacing
            if isfield(params, 'xySpacing') && ch <= length(params.xySpacing) && length(handles.xySpacingInputs) >= ch
                handles.xySpacingInputs(ch).Value = params.xySpacing{ch};
            end
            if isfield(params, 'zSpacing') && ch <= length(params.zSpacing) && length(handles.zSpacingInputs) >= ch
                handles.zSpacingInputs(ch).Value = params.zSpacing{ch};
            end
            
            % Deconvolution
            if isfield(params, 'deconvEnabled') && ch <= length(params.deconvEnabled) && length(handles.deconvEnabledChecks) >= ch
                handles.deconvEnabledChecks(ch).Value = params.deconvEnabled{ch};
            end
            if isfield(params, 'deconvMethod') && ch <= length(params.deconvMethod) && length(handles.deconvMethodDrops) >= ch
                handles.deconvMethodDrops(ch).Value = params.deconvMethod{ch};
            end
            if isfield(params, 'deconvLRIterations') && ch <= length(params.deconvLRIterations) && length(handles.deconvLRIterationsInputs) >= ch
                handles.deconvLRIterationsInputs(ch).Value = params.deconvLRIterations{ch};
            end
            if isfield(params, 'deconvLRDamping') && ch <= length(params.deconvLRDamping) && length(handles.deconvLRDampingInputs) >= ch
                handles.deconvLRDampingInputs(ch).Value = params.deconvLRDamping{ch};
            end
            if isfield(params, 'deconvWienerNSR') && ch <= length(params.deconvWienerNSR) && length(handles.deconvWienerNSRInputs) >= ch
                handles.deconvWienerNSRInputs(ch).Value = params.deconvWienerNSR{ch};
            end
            if isfield(params, 'deconvBlindIterations') && ch <= length(params.deconvBlindIterations) && length(handles.deconvBlindIterationsInputs) >= ch
                handles.deconvBlindIterationsInputs(ch).Value = params.deconvBlindIterations{ch};
            end
            if isfield(params, 'deconvBlindUnderRelax') && ch <= length(params.deconvBlindUnderRelax) && length(handles.deconvBlindUnderRelaxInputs) >= ch
                handles.deconvBlindUnderRelaxInputs(ch).Value = params.deconvBlindUnderRelax{ch};
            end
            if isfield(params, 'deconvPSFSource') && ch <= length(params.deconvPSFSource) && length(handles.deconvPSFSourceDrops) >= ch
                handles.deconvPSFSourceDrops(ch).Value = params.deconvPSFSource{ch};
            end
            if isfield(params, 'deconvPSFFilePath') && ch <= length(params.deconvPSFFilePath) && length(handles.deconvPSFPathTexts) >= ch
                handles.deconvPSFPathTexts(ch).Value = params.deconvPSFFilePath{ch};
            end
            if isfield(params, 'deconvPSFSigmaXY') && ch <= length(params.deconvPSFSigmaXY) && length(handles.deconvPSFSigmaXYInputs) >= ch
                handles.deconvPSFSigmaXYInputs(ch).Value = params.deconvPSFSigmaXY{ch};
            end
            if isfield(params, 'deconvPSFSigmaZ') && ch <= length(params.deconvPSFSigmaZ) && length(handles.deconvPSFSigmaZInputs) >= ch
                handles.deconvPSFSigmaZInputs(ch).Value = params.deconvPSFSigmaZ{ch};
            end
            if isfield(params, 'deconvPSFSizeXY') && ch <= length(params.deconvPSFSizeXY) && length(handles.deconvPSFSizeXYInputs) >= ch
                handles.deconvPSFSizeXYInputs(ch).Value = params.deconvPSFSizeXY{ch};
            end
            if isfield(params, 'deconvPSFSizeZ') && ch <= length(params.deconvPSFSizeZ) && length(handles.deconvPSFSizeZInputs) >= ch
                handles.deconvPSFSizeZInputs(ch).Value = params.deconvPSFSizeZ{ch};
            end
            
            % Preprocessing
            if isfield(params, 'preprocEnabled') && ch <= length(params.preprocEnabled) && length(handles.preprocEnabledChecks) >= ch
                handles.preprocEnabledChecks(ch).Value = params.preprocEnabled{ch};
            end
            if isfield(params, 'preProcMode') && ch <= length(params.preProcMode) && length(handles.preprocessModeDrops) >= ch
                handles.preprocessModeDrops(ch).Value = params.preProcMode{ch};
            end
            if isfield(params, 'preProcProjection') && ch <= length(params.preProcProjection) && length(handles.preprocessProjectionDrops) >= ch
                handles.preprocessProjectionDrops(ch).Value = params.preProcProjection{ch};
            end
            if isfield(params, 'preProcScale') && ch <= length(params.preProcScale) && length(handles.preprocessScaleChecks) >= ch
                handles.preprocessScaleChecks(ch).Value = params.preProcScale{ch};
            end
            if isfield(params, 'preProcMethod') && ch <= length(params.preProcMethod) && length(handles.preprocMethodDrops) >= ch
                handles.preprocMethodDrops(ch).Value = params.preProcMethod{ch};
            end
            if isfield(params, 'preprocClipAtZero') && ch <= length(params.preprocClipAtZero) && length(handles.preprocClipChecks) >= ch
                handles.preprocClipChecks(ch).Value = params.preprocClipAtZero{ch};
            end
            
            % Method-specific preprocessing parameters
            if isfield(params, 'smoothGaussianValues') && ch <= length(params.smoothGaussianValues) && length(handles.gaussInputs) >= ch
                handles.gaussInputs(ch).Value = params.smoothGaussianValues{ch};
            end
            if isfield(params, 'smoothMedianValues') && ch <= length(params.smoothMedianValues) && length(handles.medianInputs) >= ch
                handles.medianInputs(ch).Value = params.smoothMedianValues{ch};
            end
            if isfield(params, 'nlmFilterStrength') && ch <= length(params.nlmFilterStrength) && length(handles.nlmFilterStrengthInputs) >= ch
                handles.nlmFilterStrengthInputs(ch).Value = params.nlmFilterStrength{ch};
            end
            if isfield(params, 'nlmSearchWindow') && ch <= length(params.nlmSearchWindow) && length(handles.nlmSearchWindowInputs) >= ch
                handles.nlmSearchWindowInputs(ch).Value = params.nlmSearchWindow{ch};
            end
            if isfield(params, 'nlmComparisonWindow') && ch <= length(params.nlmComparisonWindow) && length(handles.nlmComparisonWindowInputs) >= ch
                handles.nlmComparisonWindowInputs(ch).Value = params.nlmComparisonWindow{ch};
            end
            if isfield(params, 'waveletName') && ch <= length(params.waveletName) && length(handles.waveletNameDrops) >= ch
                handles.waveletNameDrops(ch).Value = params.waveletName{ch};
            end
            if isfield(params, 'waveletLevel') && ch <= length(params.waveletLevel) && length(handles.waveletLevelInputs) >= ch
                handles.waveletLevelInputs(ch).Value = params.waveletLevel{ch};
            end
            if isfield(params, 'waveletThresholdRule') && ch <= length(params.waveletThresholdRule) && length(handles.waveletThresholdRuleDrops) >= ch
                handles.waveletThresholdRuleDrops(ch).Value = params.waveletThresholdRule{ch};
            end
            if isfield(params, 'waveletThresholdMethod') && ch <= length(params.waveletThresholdMethod) && length(handles.waveletThresholdMethodDrops) >= ch
                handles.waveletThresholdMethodDrops(ch).Value = params.waveletThresholdMethod{ch};
            end
            
            % Background Correction
            if isfield(params, 'bgCorrEnabled') && ch <= length(params.bgCorrEnabled) && length(handles.bgCorrEnabledChecks) >= ch
                handles.bgCorrEnabledChecks(ch).Value = params.bgCorrEnabled{ch};
            end
            if isfield(params, 'bgCorrMode') && ch <= length(params.bgCorrMode) && length(handles.bgCorrModeDrops) >= ch
                handles.bgCorrModeDrops(ch).Value = params.bgCorrMode{ch};
            end
            if isfield(params, 'bgCorrScale') && ch <= length(params.bgCorrScale) && length(handles.bgCorrScaleChecks) >= ch
                handles.bgCorrScaleChecks(ch).Value = params.bgCorrScale{ch};
            end
            if isfield(params, 'bgCorrProjection') && ch <= length(params.bgCorrProjection) && length(handles.bgCorrProjectionDrops) >= ch
                handles.bgCorrProjectionDrops(ch).Value = params.bgCorrProjection{ch};
            end
            if isfield(params, 'bgMethod') && ch <= length(params.bgMethod) && length(handles.bgMethodDrops) >= ch
                handles.bgMethodDrops(ch).Value = params.bgMethod{ch};
            end
            if isfield(params, 'bgParam') && ch <= length(params.bgParam) && length(handles.bgParamInputs) >= ch
                handles.bgParamInputs(ch).Value = params.bgParam{ch};
            end
            if isfield(params, 'bgCorrClipAtZero') && ch <= length(params.bgCorrClipAtZero) && length(handles.bgCorrClipChecks) >= ch
                handles.bgCorrClipChecks(ch).Value = params.bgCorrClipAtZero{ch};
            end
            
            % Maxima Detection
            if isfield(params, 'maximaEnabled') && ch <= length(params.maximaEnabled) && length(handles.maximaEnabledChecks) >= ch
                handles.maximaEnabledChecks(ch).Value = params.maximaEnabled{ch};
            end
            if isfield(params, 'maximaMode') && ch <= length(params.maximaMode) && length(handles.maximaModeDrops) >= ch
                handles.maximaModeDrops(ch).Value = params.maximaMode{ch};
            end
            if isfield(params, 'maximaScale') && ch <= length(params.maximaScale) && length(handles.maximaScaleChecks) >= ch
                handles.maximaScaleChecks(ch).Value = params.maximaScale{ch};
            end
            if isfield(params, 'maximaProjection') && ch <= length(params.maximaProjection) && length(handles.maximaProjectionDrops) >= ch
                handles.maximaProjectionDrops(ch).Value = params.maximaProjection{ch};
            end
            if isfield(params, 'maximaMethod') && ch <= length(params.maximaMethod) && length(handles.maximaMethodDrops) >= ch
                handles.maximaMethodDrops(ch).Value = params.maximaMethod{ch};
            end
            if isfield(params, 'maximaNeighborhoodSize') && ch <= length(params.maximaNeighborhoodSize) && length(handles.maximaNeighborhoodInputs) >= ch
                handles.maximaNeighborhoodInputs(ch).Value = params.maximaNeighborhoodSize{ch};
            end
            if isfield(params, 'hMaxValue') && ch <= length(params.hMaxValue) && length(handles.hMaxInputs) >= ch
                handles.hMaxInputs(ch).Value = params.hMaxValue{ch};
            end
            if isfield(params, 'sigmaValue') && ch <= length(params.sigmaValue) && length(handles.logSigmaInputs) >= ch
                handles.logSigmaInputs(ch).Value = params.sigmaValue{ch};
            end
            if isfield(params, 'peakThresholdValue') && ch <= length(params.peakThresholdValue) && length(handles.logThresholdInputs) >= ch
                handles.logThresholdInputs(ch).Value = params.peakThresholdValue{ch};
            end
            if isfield(params, 'showMaxima') && ch <= length(params.showMaxima) && length(handles.showMaximaChecks) >= ch && isvalid(handles.showMaximaChecks(ch))
                handles.showMaximaChecks(ch).Value = params.showMaxima{ch};
            end
            if isfield(params, 'maximaColor') && ch <= length(params.maximaColor) && length(handles.maximaColorDrops) >= ch && isvalid(handles.maximaColorDrops(ch))
                handles.maximaColorDrops(ch).Value = params.maximaColor{ch};
            end
            if isfield(params, 'displayOnAllPreviews') && ch <= length(params.displayOnAllPreviews) && length(handles.displayOnAllPreviewsChecks) >= ch && isvalid(handles.displayOnAllPreviewsChecks(ch))
                handles.displayOnAllPreviewsChecks(ch).Value = params.displayOnAllPreviews{ch};
            end
            
            % Gaussian Fitting
            if isfield(params, 'gaussFitEnabled') && ch <= length(params.gaussFitEnabled) && length(handles.gaussFitEnabledChecks) >= ch
                handles.gaussFitEnabledChecks(ch).Value = params.gaussFitEnabled{ch};
            end
            if isfield(params, 'gaussFitMethod') && ch <= length(params.gaussFitMethod) && length(handles.gaussFitMethodDrop) >= ch && isvalid(handles.gaussFitMethodDrop(ch))
                handles.gaussFitMethodDrop(ch).Value = params.gaussFitMethod{ch};
            end
            if isfield(params, 'gaussFitVoxelWindowSize') && ch <= length(params.gaussFitVoxelWindowSize) && length(handles.gaussFitVoxelWindowSlider) >= ch && isvalid(handles.gaussFitVoxelWindowSlider(ch))
                handles.gaussFitVoxelWindowSlider(ch).Value = params.gaussFitVoxelWindowSize{ch};
            end
            if isfield(params, 'gaussFitBgCorrMethod') && ch <= length(params.gaussFitBgCorrMethod) && length(handles.gaussFitBgCorrMethodDrop) >= ch && isvalid(handles.gaussFitBgCorrMethodDrop(ch))
                handles.gaussFitBgCorrMethodDrop(ch).Value = params.gaussFitBgCorrMethod{ch};
            end
            if isfield(params, 'gaussFitBgCorrWidth') && ch <= length(params.gaussFitBgCorrWidth) && length(handles.gaussFitBgCorrWidthEdit) >= ch && isvalid(handles.gaussFitBgCorrWidthEdit(ch))
                handles.gaussFitBgCorrWidthEdit(ch).Value = params.gaussFitBgCorrWidth{ch};
            end
            if isfield(params, 'gaussFitPolyDegree') && ch <= length(params.gaussFitPolyDegree) && length(handles.gaussFitPolyDegreeEdit) >= ch && isvalid(handles.gaussFitPolyDegreeEdit(ch))
                handles.gaussFitPolyDegreeEdit(ch).Value = params.gaussFitPolyDegree{ch};
            end
            if isfield(params, 'gaussFitMaxIterations') && ch <= length(params.gaussFitMaxIterations) && length(handles.gaussFitMaxIterationsEdit) >= ch && isvalid(handles.gaussFitMaxIterationsEdit(ch))
                handles.gaussFitMaxIterationsEdit(ch).Value = params.gaussFitMaxIterations{ch};
            end
            if isfield(params, 'gaussFitTolerance') && ch <= length(params.gaussFitTolerance) && length(handles.gaussFitToleranceEdit) >= ch && isvalid(handles.gaussFitToleranceEdit(ch))
                handles.gaussFitToleranceEdit(ch).Value = params.gaussFitTolerance{ch};
            end
            if isfield(params, 'gaussFitRadialRadius') && ch <= length(params.gaussFitRadialRadius) && length(handles.gaussFitRadialRadiusEdit) >= ch && isvalid(handles.gaussFitRadialRadiusEdit(ch))
                handles.gaussFitRadialRadiusEdit(ch).Value = params.gaussFitRadialRadius{ch};
            end
            if isfield(params, 'gaussFitPlotCheck') && ch <= length(params.gaussFitPlotCheck) && length(handles.gaussFitPlotCheck) >= ch && isvalid(handles.gaussFitPlotCheck(ch))
                handles.gaussFitPlotCheck(ch).Value = params.gaussFitPlotCheck{ch};
            end
            
            % Fit Filtering
            if isfield(params, 'fitFilterEnabled') && ch <= length(params.fitFilterEnabled) && length(handles.fitFilterEnabledChecks) >= ch
                handles.fitFilterEnabledChecks(ch).Value = params.fitFilterEnabled{ch};
            end
            if isfield(params, 'fitFilterRSquaredEnabled') && ch <= length(params.fitFilterRSquaredEnabled) && length(handles.fitFilterRSquaredEnabledChecks) >= ch
                handles.fitFilterRSquaredEnabledChecks(ch).Value = params.fitFilterRSquaredEnabled{ch};
            end
            if isfield(params, 'fitFilterRSquaredMin') && ch <= length(params.fitFilterRSquaredMin) && length(handles.fitFilterRSquaredMinInputs) >= ch
                handles.fitFilterRSquaredMinInputs(ch).Value = params.fitFilterRSquaredMin{ch};
            end
            if isfield(params, 'fitFilterRSquaredMax') && ch <= length(params.fitFilterRSquaredMax) && length(handles.fitFilterRSquaredMaxInputs) >= ch
                handles.fitFilterRSquaredMaxInputs(ch).Value = params.fitFilterRSquaredMax{ch};
            end
            if isfield(params, 'fitFilterSigmaSumEnabled') && ch <= length(params.fitFilterSigmaSumEnabled) && length(handles.fitFilterSigmaSumEnabledChecks) >= ch
                handles.fitFilterSigmaSumEnabledChecks(ch).Value = params.fitFilterSigmaSumEnabled{ch};
            end
            if isfield(params, 'fitFilterSigmaSumMin') && ch <= length(params.fitFilterSigmaSumMin) && length(handles.fitFilterSigmaSumMinInputs) >= ch
                handles.fitFilterSigmaSumMinInputs(ch).Value = params.fitFilterSigmaSumMin{ch};
            end
            if isfield(params, 'fitFilterSigmaSumMax') && ch <= length(params.fitFilterSigmaSumMax) && length(handles.fitFilterSigmaSumMaxInputs) >= ch
                handles.fitFilterSigmaSumMaxInputs(ch).Value = params.fitFilterSigmaSumMax{ch};
            end
            if isfield(params, 'fitFilterAmplitudeEnabled') && ch <= length(params.fitFilterAmplitudeEnabled) && length(handles.fitFilterAmplitudeEnabledChecks) >= ch
                handles.fitFilterAmplitudeEnabledChecks(ch).Value = params.fitFilterAmplitudeEnabled{ch};
            end
            if isfield(params, 'fitFilterAmplitudeMin') && ch <= length(params.fitFilterAmplitudeMin) && length(handles.fitFilterAmplitudeMinInputs) >= ch
                handles.fitFilterAmplitudeMinInputs(ch).Value = params.fitFilterAmplitudeMin{ch};
            end
            if isfield(params, 'fitFilterAmplitudeMax') && ch <= length(params.fitFilterAmplitudeMax) && length(handles.fitFilterAmplitudeMaxInputs) >= ch
                handles.fitFilterAmplitudeMaxInputs(ch).Value = params.fitFilterAmplitudeMax{ch};
            end
            if isfield(params, 'fitFilterIntensityEnabled') && ch <= length(params.fitFilterIntensityEnabled) && length(handles.fitFilterIntensityEnabledChecks) >= ch
                handles.fitFilterIntensityEnabledChecks(ch).Value = params.fitFilterIntensityEnabled{ch};
            end
            if isfield(params, 'fitFilterIntensityMin') && ch <= length(params.fitFilterIntensityMin) && length(handles.fitFilterIntensityMinInputs) >= ch
                handles.fitFilterIntensityMinInputs(ch).Value = params.fitFilterIntensityMin{ch};
            end
            if isfield(params, 'fitFilterIntensityMax') && ch <= length(params.fitFilterIntensityMax) && length(handles.fitFilterIntensityMaxInputs) >= ch
                handles.fitFilterIntensityMaxInputs(ch).Value = params.fitFilterIntensityMax{ch};
            end
            catch ME
                % Skip this channel if UI elements are invalid (e.g., GraphicsPlaceholder)
                % This can happen for channels that weren't fully initialized
                % fprintf('Warning: Could not apply parameters for channel %d: %s\n', ch, ME.message);
            end
        end
    
        % Preview settings
        try
            if isfield(params, 'previewContents') && isfield(handles, 'previewContentDrops')
                for i = 1:min(5, length(params.previewContents), length(handles.previewContentDrops))
                    if ~isempty(params.previewContents{i})
                        % Only set if the value is in the Items list
                        if ismember(params.previewContents{i}, handles.previewContentDrops(i).Items)
                            handles.previewContentDrops(i).Value = params.previewContents{i};
                        end
                    end
                end
            end
            if isfield(params, 'previewModes') && isfield(handles, 'previewModeDrops')
                for i = 1:min(5, length(params.previewModes), length(handles.previewModeDrops))
                    if ~isempty(params.previewModes{i})
                        % Only set if the value is in the Items list
                        if ismember(params.previewModes{i}, handles.previewModeDrops(i).Items)
                            handles.previewModeDrops(i).Value = params.previewModes{i};
                        end
                    end
                end
            end
            if isfield(params, 'previewProjections') && isfield(handles, 'previewProjectionDrops')
                for i = 1:min(5, length(params.previewProjections), length(handles.previewProjectionDrops))
                    if ~isempty(params.previewProjections{i})
                        % Only set if the value is in the Items list
                        if ismember(params.previewProjections{i}, handles.previewProjectionDrops(i).Items)
                            handles.previewProjectionDrops(i).Value = params.previewProjections{i};
                        end
                    end
                end
            end
        catch
            % Silently skip preview settings if there's an issue (e.g., value not in Items list)
        end
    
        % Export settings
        if isfield(params, 'exportImageFormat') && isfield(handles, 'exportImageFormatDrop')
            handles.exportImageFormatDrop.Value = params.exportImageFormat;
        end
        
    catch ME
        warning('Error applying parameters to UI: %s', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).file, ME.stack(i).line);
        end
    end
end

% ============================================================================
% CLASSIFICATION HELPER FUNCTIONS
% ============================================================================

function updateClassificationControls(channelIdx, fig)
% Update classification panel controls based on enabled state
    handles = guidata(fig);
    enabled = handles.classifyEnabledChecks(channelIdx).Value;
    
    if enabled
        enableState = 'on';
    else
        enableState = 'off';
    end
    
    handles.classifyLoadButtons(channelIdx).Enable = enableState;
    handles.classifyTrainButtons(channelIdx).Enable = enableState;
    handles.classifyFilterNoiseChecks(channelIdx).Enable = enableState;
    
    % Apply button only if classifier is loaded
    if enabled && ~isempty(handles.classifiers{channelIdx})
        handles.classifyApplyButtons(channelIdx).Enable = 'on';
        handles.classifySelectFeaturesButtons(channelIdx).Enable = 'on';
    else
        handles.classifyApplyButtons(channelIdx).Enable = 'off';
        handles.classifySelectFeaturesButtons(channelIdx).Enable = 'off';
    end
    
    guidata(fig, handles);
end

function loadClassifierForChannel(channelIdx, fig)
% Load a trained classifier for the specified channel
    handles = guidata(fig);
    
    [file, path] = uigetfile('*.mat', 'Load Classifier');
    if file == 0, return; end
    
    % Load classifier with custom expressions and normalization params
    [model, features, featureInfo, trainStats, fittingMethod, ~, success, customExpr, normParams] = ...
        snap_helpers.classification.loadClassifier(fullfile(path, file));
    
    if success
        handles.classifiers{channelIdx} = model;
        handles.classifierFeatures{channelIdx} = features;
        handles.classifierFeatureInfo{channelIdx} = featureInfo;
        handles.classifierTrainStats{channelIdx} = trainStats;
        
        % Store custom expressions and norm params
        if ~isfield(handles, 'classifierCustomExpressions')
            handles.classifierCustomExpressions = cell(1, 10);
        end
        if ~isfield(handles, 'classifierNormParams')
            handles.classifierNormParams = cell(1, 10);
        end
        handles.classifierCustomExpressions{channelIdx} = customExpr;
        handles.classifierNormParams{channelIdx} = normParams;
        
        % Update UI
        nBase = numel(features);
        nCustom = numel(customExpr);
        handles.classifyStatusLabels(channelIdx).Text = sprintf('Loaded (%s)', fittingMethod);
        handles.classifyStatusLabels(channelIdx).FontColor = [0.2 0.6 0.2];
        
        if nCustom > 0
            handles.classifyFeatureCountLabels(channelIdx).Text = sprintf('%d base + %d expr', nBase, nCustom);
        else
            handles.classifyFeatureCountLabels(channelIdx).Text = sprintf('%d features', nBase);
        end
        handles.classifyFeatureCountLabels(channelIdx).FontColor = [0.2 0.6 0.2];
        handles.classifyApplyButtons(channelIdx).Enable = 'on';
        handles.classifySelectFeaturesButtons(channelIdx).Enable = 'on';
        
        % Check fitting method compatibility
        currentMethod = handles.gaussFitMethodDrop(channelIdx).Value;
        if ~contains(currentMethod, fittingMethod, 'IgnoreCase', true) && ...
           ~contains(fittingMethod, currentMethod, 'IgnoreCase', true)
            warndlg(sprintf('Warning: Classifier was trained on "%s" but current method is "%s".\nResults may be unreliable.', ...
                fittingMethod, currentMethod), 'Method Mismatch');
        end
        
        guidata(fig, handles);
    else
        uialert(fig, 'Failed to load classifier', 'Load Error');
    end
end

function applyClassifierToChannel(channelIdx, fig)
% Apply loaded classifier to current fit results
    handles = guidata(fig);
    
    if isempty(handles.classifiers{channelIdx})
        uialert(fig, 'No classifier loaded for this channel', 'No Classifier');
        return;
    end
    
    % Get current fit results from preview cache
    if ~isfield(handles, 'previewCache') || isempty(handles.previewCache)
        uialert(fig, 'No fit results available. Run processing first.', 'No Data');
        return;
    end
    
    if ~isfield(handles.previewCache, 'channels') || channelIdx > numel(handles.previewCache.channels)
        uialert(fig, 'No data for this channel', 'No Data');
        return;
    end
    
    channelData = handles.previewCache.channels(channelIdx);
    if ~isfield(channelData, 'fitResults') || isempty(channelData.fitResults)
        uialert(fig, 'No fit results for this channel. Enable fitting first.', 'No Fits');
        return;
    end
    
    fitResults = channelData.fitResults;
    
    % Build feature matrix with custom expressions
    features = handles.classifierFeatures{channelIdx};
    featureInfo = handles.classifierFeatureInfo{channelIdx};
    
    % Get custom expressions if available
    customExpr = struct('name', {}, 'expression', {});
    if isfield(handles, 'classifierCustomExpressions') && ...
       numel(handles.classifierCustomExpressions) >= channelIdx && ...
       ~isempty(handles.classifierCustomExpressions{channelIdx})
        customExpr = handles.classifierCustomExpressions{channelIdx};
    end
    
    [X, featureNames, validMask] = snap_helpers.classification.buildFeatureMatrix(...
        fitResults, features, featureInfo, customExpr);
    
    % Get normalization parameters
    normParams = struct('mu', [], 'sigma', [], 'standardized', false);
    if isfield(handles, 'classifierNormParams') && ...
       numel(handles.classifierNormParams) >= channelIdx && ...
       ~isempty(handles.classifierNormParams{channelIdx})
        normParams = handles.classifierNormParams{channelIdx};
    end
    
    % Apply classifier with normalization
    [predictions, ~, confidence] = snap_helpers.classification.applyClassifier(...
        handles.classifiers{channelIdx}, X, featureNames, featureNames, normParams);
    
    % Store results
    handles.previewCache.channels(channelIdx).classifications = struct();
    handles.previewCache.channels(channelIdx).classifications.predictions = predictions;
    handles.previewCache.channels(channelIdx).classifications.confidence = confidence;
    handles.previewCache.channels(channelIdx).classifications.validMask = validMask;
    
    % Count results
    nReal = sum(predictions == 1 & validMask);
    nNoise = sum(predictions == 0 & validMask);
    nInvalid = sum(~validMask);
    
    % Update filter mask if auto-filter is enabled
    if handles.classifyFilterNoiseChecks(channelIdx).Value
        if isfield(channelData, 'filterMask')
            oldMask = channelData.filterMask;
            newMask = oldMask;
            % Only keep spots that are: passing old filter AND (predicted real OR invalid)
            validIdx = find(oldMask);
            for i = 1:numel(validIdx)
                if validMask(validIdx(i)) && predictions(validIdx(i)) == 0
                    newMask(validIdx(i)) = false;
                end
            end
            handles.previewCache.channels(channelIdx).filterMask = newMask;
            
            nFiltered = sum(oldMask) - sum(newMask);
            msgStr = sprintf('Classification complete:\n- Real: %d\n- Noise: %d (filtered)\n- Invalid: %d', ...
                nReal, nNoise, nInvalid);
        else
            msgStr = sprintf('Classification complete:\n- Real: %d\n- Noise: %d\n- Invalid: %d\n\nEnable fit filtering to auto-remove noise.', ...
                nReal, nNoise, nInvalid);
        end
    else
        msgStr = sprintf('Classification complete:\n- Real: %d\n- Noise: %d\n- Invalid features: %d', ...
            nReal, nNoise, nInvalid);
    end
    
    guidata(fig, handles);
    
    % Refresh preview
    snap_helpers.redrawPreview(fig, channelIdx);
    
    uialert(fig, msgStr, 'Classification Complete', 'Icon', 'success');
end

function selectClassifierFeatures(channelIdx, fig)
% Open feature selection UI for this channel
    handles = guidata(fig);
    
    fittingMethod = handles.gaussFitMethodDrop(channelIdx).Value;
    
    % Determine if 3D
    has3D = false;
    if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'channels') && ...
       channelIdx <= numel(handles.previewCache.channels) && ...
       isfield(handles.previewCache.channels(channelIdx), 'processed')
        has3D = size(handles.previewCache.channels(channelIdx).processed, 3) > 1;
    end
    
    previousSelection = handles.classifierFeatures{channelIdx};
    if isempty(previousSelection)
        previousSelection = {};
    end
    
    [selected, cancelled] = snap_helpers.classification.featureSelectionUI(...
        fittingMethod, has3D, false, previousSelection);
    
    if ~cancelled
        handles.classifierFeatures{channelIdx} = selected;
        handles.classifyFeatureCountLabels(channelIdx).Text = sprintf('%d features', numel(selected));
        handles.classifyFeatureCountLabels(channelIdx).FontColor = [0.2 0.6 0.2];
        guidata(fig, handles);
    end
end