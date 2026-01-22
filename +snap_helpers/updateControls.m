function updateControls(fig_handle)
% Central function to update the entire UI state based on current selections.
    handles = guidata(fig_handle);
    
    % Prevent infinite recursion by checking if we're already updating
    if isfield(handles, 'isUpdatingControls') && handles.isUpdatingControls
        return;
    end
    
    % Set flag to prevent recursion
    handles.isUpdatingControls = true;
    guidata(fig_handle, handles);
    
    numChannels = str2double(handles.numChanDrop.Value);
    
    % --- Part 1: Update Nuclei Tab Controls ---
    % Update nuclei segmentation method controls
    updateNucleiSegmentationControls(handles);
    
    % Pre-processing controls visibility
    nuc_mode = handles.nucPreprocessModeDrop.Value;
    is_nuc_projection = contains(nuc_mode, 'Projection');
    is_nuc_3d = strcmp(nuc_mode, '3D');
    
    handles.nucPreprocessProjectionDrop.Visible = is_nuc_projection;
    handles.nucPreprocessScaleCheck.Visible = is_nuc_3d;
    
    % Wavelet controls are only for 2D modes
    nuc_wavelet_controls = findall(handles.fig, 'Tag', 'nuc_preproc_param_ctrl');
    set(nuc_wavelet_controls, 'Enable', matlab.lang.OnOffSwitchState(~is_nuc_3d));
    % Note: The actual method selection and parameter visibility is now handled by the method dropdown
    
    % Enable/disable parameter inputs based on method selection
    % Note: The nuclei preprocessing now uses the same dynamic system as fluorescent channels
    % The actual enabling/disabling is handled by the method dropdown and updatePreprocControls function
    
    % Update nuclei deconvolution controls based on method selection
    if isfield(handles, 'nucDeconvMethodDrop')
        updateNucleiDeconvControls(handles);
    end
    
    % Update nuclei preprocessing controls based on method selection
    updateNucleiPreprocControls(handles);

    % Background correction controls visibility
    if isfield(handles, 'nucBgCorrModeDrop')
        nuc_bg_mode = handles.nucBgCorrModeDrop.Value;
        is_nuc_bg_projection = contains(nuc_bg_mode, 'Projection');
        is_nuc_bg_3d = strcmp(nuc_bg_mode, '3D');

        if isfield(handles, 'nucBgCorrProjectionDrop')
            handles.nucBgCorrProjectionDrop.Visible = is_nuc_bg_projection;
        end
        if isfield(handles, 'nucBgCorrScaleCheck')
            handles.nucBgCorrScaleCheck.Visible = is_nuc_bg_3d;
        end

        % Update nuclei background correction controls based on method selection
        updateNucleiBgCorrControls(handles);
    end

    % Segmentation projection dropdown visibility
    handles.nucSegProjectionDrop.Visible = strcmp(handles.nucSegModeDrop.Value, 'On Z-Projection');
    
        % Update Nuclei Scale Checkbox State
    updateNucleiScaleCheckboxState(handles);
    
    % --- Part 1a: Update Nuclei Panel Enabled States ---
    if isfield(handles, 'nucDeconvPanel')
        updatePanelEnabledState(handles.nucDeconvPanel, handles.nucDeconvEnabledCheck);
    end
    updatePanelEnabledState(handles.nucPreprocPanel, handles.nucPreprocEnabledCheck);
    if isfield(handles, 'nucBgPanel') && isfield(handles, 'nucBgCorrEnabledCheck')
        updatePanelEnabledState(handles.nucBgPanel, handles.nucBgCorrEnabledCheck);
    end
    updatePanelEnabledState(handles.segPanel, handles.nucSegEnabledCheck);
    updatePanelEnabledState(handles.nucFilterPanel, handles.nucFilterEnabledCheck);
    
    % Update segmentation controls visibility
    updateNucleiSegmentationControls(handles);
    
    % Update filter controls
    updateNucleiFilterControls(handles);

    % --- Part 2: Enable/disable channel tabs ---
    for k = 1:handles.Nmax
        is_active = k <= numChannels;
        set(findall(handles.channelTabs(k), '-property', 'Enable'), 'Enable', matlab.lang.OnOffSwitchState(is_active));
    end
    
    % --- Part 3: Update preview content dropdowns ---
    new_items = {'None', 'DIC', 'Nuclei'};
    for k = 1:numChannels
        new_items{end+1} = ['Channel ' num2str(k)];
    end
    
    for i = 1:5
        current_value = handles.previewContentDrops(i).Value;
        handles.previewContentDrops(i).Items = unique(new_items, 'stable');
        if ~ismember(current_value, handles.previewContentDrops(i).Items)
            handles.previewContentDrops(i).Value = 'None';
        else
            handles.previewContentDrops(i).Value = current_value;
        end
    end

    % --- Part 4: Update controls within each active channel tab ---
    for k = 1:handles.Nmax
        % Update panel enabled states FIRST to set the baseline
        if isfield(handles, 'deconvPanels')
            updatePanelEnabledState(handles.deconvPanels(k), handles.deconvEnabledChecks(k));
        end
        updatePanelEnabledState(handles.preprocPanels(k), handles.preprocEnabledChecks(k));
        updatePanelEnabledState(handles.bgPanels(k), handles.bgCorrEnabledChecks(k));
        updatePanelEnabledState(handles.maximaPanels(k), handles.maximaEnabledChecks(k));
        updatePanelEnabledState(handles.gaussFitPanels(k), handles.gaussFitEnabledChecks(k));
        updatePanelEnabledState(handles.fitFilterPanels(k), handles.fitFilterEnabledChecks(k));

        % Update visibility based on processing modes (3D vs. Projection)
        updateProcessingModeControls(k, 'preprocess', handles);
        updateProcessingModeControls(k, 'bg', handles);
        updateProcessingModeControls(k, 'maxima', handles);
        
        % Update method-specific controls
        if isfield(handles, 'deconvMethodDrops')
            updateDeconvControls(k, handles);
        end
        updatePreprocControls(k, handles);
        updateMaximaControls(k, handles);
        updateDisplayOnAllPreviewsControls(k, handles);
        
        % Update controls for the Local Maxima Fitting panel
        updateGaussianFitControls(k, handles);
        
        % Update fit filtering controls
        updateFitFilterControls(k, handles);
        
        % Enable/disable background parameter based on method
        is_bg_none = strcmp(handles.bgMethodDrops(k).Value, 'None');
        handles.bgParamInputs(k).Enable = ~is_bg_none;
        
        % Update scale checkbox usability
        updateScaleCheckboxState(k, handles);
    end
    
    % Clear recursion flag before saving
    handles.isUpdatingControls = false;
    guidata(handles.fig, handles);
    enforcePreviewModes(handles); % Enforce preview modes after all updates
end

% --- Nested Helper Functions ---

function updatePreprocControls(channel_idx, handles)
    % Manages the dynamic parameter display for pre-processing methods.
    if ~handles.preprocEnabledChecks(channel_idx).Value
        return;
    end
    
    method = handles.preprocMethodDrops(channel_idx).Value;
    is_3d_mode = strcmp(handles.preprocessModeDrops(channel_idx).Value, '3D');
    
    % Handle 3D mode restrictions first
    if is_3d_mode
        if strcmp(method, 'Wavelet Denoising')
            handles.preprocMethodDrops(channel_idx).Value = 'None';
            method = 'None';
        end
        % Disable the wavelet option in the dropdown
        handles.preprocMethodDrops(channel_idx).Items = {'None', 'Gaussian', 'Median', 'Non-Local Means'};
    else
        handles.preprocMethodDrops(channel_idx).Items = {'None', 'Gaussian', 'Median', 'Non-Local Means', 'Wavelet Denoising'};
    end
    
    % Get handles to dynamic parameter controls
    param1Label = handles.preprocParam1Labels(channel_idx);
    param1Input = handles.preprocParam1Inputs(channel_idx);
    param2Label = handles.preprocParam2Labels(channel_idx);
    param2Input = handles.preprocParam2Inputs(channel_idx);
    param3Label = handles.preprocParam3Labels(channel_idx);
    param3Input = handles.preprocParam3Inputs(channel_idx);
    param4Label = handles.preprocParam4Labels(channel_idx);
    param4Input = handles.preprocParam4Inputs(channel_idx);
    waveletDropdown = handles.waveletNameDrops(channel_idx);
    
    % Configure controls based on selected method
    switch method
        case 'None'
            % Hide all parameter controls
            param1Label.Visible = 'off';
            param1Input.Visible = 'off';
            param2Label.Visible = 'off';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            
        case 'Gaussian'
            % Show: "Sigma: [input] (units)"
            param1Label.Text = 'Sigma:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            param2Label.Text = '(units)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            
        case 'Median'
            % Show: "Size: [input] (units)"
            param1Label.Text = 'Size:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            param2Label.Text = '(units)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            
        case 'Non-Local Means'
            % Show all three NLM parameters:
            % Row 4: "Filter Strength: [input] (0.0-1.0)"
            % Row 5: "Search Window: [input]"
            % Row 6: "Comparison Window: [input]"
            param1Label.Text = 'Filter Strength:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'off'; % Use specific NLM controls instead
            param2Label.Text = '(0.0-1.0)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Text = 'Search Window:';
            param3Label.Visible = 'on';
            param3Input.Visible = 'off'; % Use specific NLM controls instead
            param4Label.Text = 'Comparison Window:';
            param4Label.Visible = 'on';
            param4Input.Visible = 'off'; % Use specific NLM controls instead
            waveletDropdown.Visible = 'off';
            
            % Show NLM-specific controls
            handles.nlmFilterStrengthInputs(channel_idx).Visible = 'on';
            handles.nlmSearchWindowInputs(channel_idx).Visible = 'on';
            handles.nlmComparisonWindowInputs(channel_idx).Visible = 'on';
            
        case 'Wavelet Denoising'
            % Show all wavelet parameters:
            % Row 4: "Wavelet: [dropdown] Level: [input]"
            % Row 5: "Rule: [dropdown]"
            % Row 6: "Method: [dropdown]"
            param1Label.Text = 'Wavelet:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'off';
            param2Label.Text = 'Level:';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
            param3Label.Text = 'Threshold Rule:';
            param3Label.Visible = 'on';
            param3Input.Visible = 'on';
            param4Label.Text = 'Threshold Method:';
            param4Label.Visible = 'on';
            param4Input.Visible = 'on';
            waveletDropdown.Visible = 'on';
            % Hide NLM-specific controls
            handles.nlmFilterStrengthInputs(channel_idx).Visible = 'off';
            handles.nlmSearchWindowInputs(channel_idx).Visible = 'off';
            handles.nlmComparisonWindowInputs(channel_idx).Visible = 'off';
            
    end
    
    % Ensure NLM controls are hidden for all methods except Non-Local Means
    if ~strcmp(method, 'Non-Local Means')
        handles.nlmFilterStrengthInputs(channel_idx).Visible = 'off';
        handles.nlmSearchWindowInputs(channel_idx).Visible = 'off';
        handles.nlmComparisonWindowInputs(channel_idx).Visible = 'off';
    end
    
end

function updateGaussianFitControls(channel_idx, handles)
    % This function robustly manages the state of the Local Maxima Fitting panel's controls.

    % First, check if the whole panel is enabled. If not, do nothing.
    if ~handles.gaussFitEnabledChecks(channel_idx).Value
        return;
    end

    % --- Polynomial Degree Controls ---
    bg_method = handles.gaussFitBgCorrMethodDrop(channel_idx).Value;
    is_poly = strcmp(bg_method, 'Local Polynomial Fitting');
    
    % Find controls ONLY within the current channel's tab
    poly_label = findall(handles.channelTabs(channel_idx), 'Tag', ['poly_degree_label_' num2str(channel_idx)]);
    poly_edit = findall(handles.channelTabs(channel_idx), 'Tag', ['poly_degree_edit_' num2str(channel_idx)]);
    
    set([poly_label; poly_edit], 'Enable', matlab.lang.OnOffSwitchState(is_poly));


    % --- Dropdown Item Logic for 2D/3D ---
    fit_drop = handles.gaussFitMethodDrop(channel_idx);
    is_3d = strcmp(handles.zSpacingInputs(channel_idx).Enable, 'on');
    if ~is_3d
        if any(strcmp(fit_drop.Value, {'2D (XY) + 1D (Z) Gaussian', '2D (XY) + 1D (Z)', '3D Gaussian', 'Distorted 3D Gaussian'}))
            fit_drop.Value = '1D (X,Y,Z) Gaussian';
        end
        fit_drop.Items = {'1D (X,Y,Z) Gaussian', 'Radial Symmetry'};
    else
        fit_drop.Items = {'1D (X,Y,Z) Gaussian', '2D (XY) + 1D (Z) Gaussian', '3D Gaussian', 'Distorted 3D Gaussian', 'Radial Symmetry'};
    end
    
    % --- Radial Radius Visibility ---
    % Show radial radius input only when Radial Symmetry is selected
    is_radial = strcmp(fit_drop.Value, 'Radial Symmetry');
    if isfield(handles, 'gaussFitRadialRadiusEdit') && length(handles.gaussFitRadialRadiusEdit) >= channel_idx
        % Find the label for radial radius
        radial_radius_label = findall(handles.channelTabs(channel_idx), 'Tag', 'radial_radius_label');
        if ~isempty(radial_radius_label)
            radial_radius_label.Visible = matlab.lang.OnOffSwitchState(is_radial);
        end
        handles.gaussFitRadialRadiusEdit(channel_idx).Visible = matlab.lang.OnOffSwitchState(is_radial);
    end
end

function updateProcessingModeControls(channel_idx, category, handles)
    % First, check if the whole panel is enabled. If not, do nothing.
    if strcmp(category, 'preprocess') && ~handles.preprocEnabledChecks(channel_idx).Value
        return;
    elseif strcmp(category, 'bg') && ~handles.bgCorrEnabledChecks(channel_idx).Value
        return;
    elseif strcmp(category, 'maxima') && ~handles.maximaEnabledChecks(channel_idx).Value
        return;
    end
    
    if strcmp(category, 'bg')
        category = 'bgCorr'; % Correct the category name
    end
    mode_drop = handles.([category 'ModeDrops'])(channel_idx);
    proj_drop = handles.([category 'ProjectionDrops'])(channel_idx);
    
    is_projection = contains(mode_drop.Value, 'Projection');
    is_3d_mode = strcmp(mode_drop.Value, '3D');
    
    proj_drop.Visible = is_projection;
    
    if ~strcmp(category, 'maxima') % Scale check was removed from maxima
        scale_check = handles.([category 'ScaleChecks'])(channel_idx);
        scale_check.Visible = is_3d_mode;
    end
    
    if strcmp(category, 'preprocess')
        wavelet_controls = findall(handles.channelTabs(channel_idx), 'Tag', 'wavelet_control');
        set(wavelet_controls, 'Enable', matlab.lang.OnOffSwitchState(~is_3d_mode));
        if is_3d_mode, handles.waveletDenoiseChecks(channel_idx).Value = false; end
    end
end

function updateMaximaControls(channel_idx, handles)
    % First, check if the whole panel is enabled. If not, do nothing.
    if ~handles.maximaEnabledChecks(channel_idx).Value
        return;
    end
    
    selected_method = handles.maximaMethodDrops(channel_idx).Value;
    
    is_hmax = strcmp(selected_method, 'Extended Maxima');
    is_log = strcmp(selected_method, 'Laplacian of Gaussian');
    
    % Get handles to the specific controls that can be enabled/disabled
    hmax_panel_controls = findall(handles.hMaxPanel(channel_idx), '-property', 'Enable');
    log_panel_controls = findall(handles.logPanel(channel_idx), '-property', 'Enable');

    % Enable/disable the controls within those panels, but keep them visible.
    set(hmax_panel_controls, 'Enable', matlab.lang.OnOffSwitchState(is_hmax));
    set(log_panel_controls, 'Enable', matlab.lang.OnOffSwitchState(is_log));
end

function updateWaveletInputState(channel_idx, handles)
    % This function is now obsolete as its logic is handled by updatePreprocControls
end

function updateScaleCheckboxState(channel_idx, handles)
    xy_val = handles.xySpacingInputs(channel_idx).Value;
    z_val = handles.zSpacingInputs(channel_idx).Value;
    is_3d_image = strcmp(handles.zSpacingInputs(channel_idx).Enable, 'on');
    
    can_scale_base = is_3d_image && (xy_val > 0) && (z_val > 0);
    
    % Check state of Pre-processing panel's scale checkbox
    % Scale is only available for 3D mode since anisotropic spacing only matters for 3D operations
    preproc_mode = handles.preprocessModeDrops(channel_idx).Value;
    is_preproc_3d = strcmp(preproc_mode, '3D');
    preproc_can_scale = can_scale_base && handles.preprocEnabledChecks(channel_idx).Value && is_preproc_3d;
    set(handles.preprocessScaleChecks(channel_idx), 'Enable', matlab.lang.OnOffSwitchState(preproc_can_scale));
    if ~preproc_can_scale, handles.preprocessScaleChecks(channel_idx).Value = false; end

    % Check state of Background Correction panel's scale checkbox
    % Scale is only available for 3D mode since anisotropic spacing only matters for 3D operations
    bg_mode = handles.bgCorrModeDrops(channel_idx).Value;
    is_bg_3d = strcmp(bg_mode, '3D');
    bg_can_scale = can_scale_base && handles.bgCorrEnabledChecks(channel_idx).Value && is_bg_3d;
    set(handles.bgCorrScaleChecks(channel_idx), 'Enable', matlab.lang.OnOffSwitchState(bg_can_scale));
    if ~bg_can_scale, handles.bgCorrScaleChecks(channel_idx).Value = false; end
    
    % Check state of Maxima Detection panel's scale checkbox
    % Scale is only available for 3D mode since anisotropic spacing only matters for 3D operations
    maxima_mode = handles.maximaModeDrops(channel_idx).Value;
    is_maxima_3d = strcmp(maxima_mode, '3D');
    maxima_can_scale = can_scale_base && handles.maximaEnabledChecks(channel_idx).Value && is_maxima_3d;
    set(handles.maximaScaleChecks(channel_idx), 'Enable', matlab.lang.OnOffSwitchState(maxima_can_scale));
    if ~maxima_can_scale, handles.maximaScaleChecks(channel_idx).Value = false; end
end

function updateNucleiScaleCheckboxState(handles)
    is_3d_image = isfield(handles, 'nucZSpacingInput') && strcmp(handles.nucZSpacingInput.Enable, 'on');
    xy_val = handles.nucXYSpacingInput.Value;
    z_val = handles.nucZSpacingInput.Value;
    
    % Scale is only available for 3D mode since anisotropic spacing only matters for 3D operations
    preproc_mode = handles.nucPreprocessModeDrop.Value;
    is_preproc_3d = strcmp(preproc_mode, '3D');
    
    can_scale_preproc = is_3d_image && (xy_val > 0) && (z_val > 0) && is_preproc_3d;
    
    handles.nucPreprocessScaleCheck.Enable = can_scale_preproc;
    if ~can_scale_preproc
        handles.nucPreprocessScaleCheck.Value = false;
    end
    
    % Handle background correction scaling if controls exist
    if isfield(handles, 'nucBgCorrModeDrop') && isfield(handles, 'nucBgCorrScaleCheck')
        bg_mode = handles.nucBgCorrModeDrop.Value;
        is_bg_3d = strcmp(bg_mode, '3D');
        
        can_scale_bg = is_3d_image && (xy_val > 0) && (z_val > 0) && is_bg_3d;
        
        handles.nucBgCorrScaleCheck.Enable = can_scale_bg;
        if ~can_scale_bg
            handles.nucBgCorrScaleCheck.Value = false;
        end
    end
end

function enforcePreviewModes(handles)
% Allow all preview modes for nuclei - preview mode should be independent of segmentation mode
    % Always enable all nuclei preview modes regardless of segmentation settings
    % Users should be free to view nuclei in any preview mode they prefer
    for i = 1:5
        if contains(handles.previewContentDrops(i).Value, 'Nuclei')
            handles.previewModeDrops(i).Enable = 'on';
        end
    end
end

function updatePanelEnabledState(panel_handle, checkbox_handle)
    % Enables or disables all uicontrols within a panel based on a checkbox.
    is_enabled = checkbox_handle.Value;
    
    % Find all UI components within the panel, excluding the checkbox and title.
    outer_grid = checkbox_handle.Parent.Parent; % The checkbox is in a title_grid, which is in the outer_grid
    children = outer_grid.Children;
    
    % The inner_grid is the other child of the outer_grid (that isn't the title_grid)
    inner_grid = children(children ~= checkbox_handle.Parent);

    if ~isempty(inner_grid)
        controls_to_toggle = findall(inner_grid, '-property', 'Enable');
        set(controls_to_toggle, 'Enable', matlab.lang.OnOffSwitchState(is_enabled));
    end

    % Also update the visual state of the panel (e.g., title color)
    title_label = findobj(checkbox_handle.Parent, 'Type', 'uilabel');
    if is_enabled
        title_label.FontColor = 'k'; % Black for enabled
        panel_handle.BorderColor = [0.15 0.15 0.15];
    else
        title_label.FontColor = [0.5 0.5 0.5]; % Gray for disabled
        panel_handle.BorderColor = [0.8 0.8 0.8];
    end
end


function updateInclusionExclusionControls(handles)
    % Updates the visibility and state of inclusion/exclusion controls
    % Apply to Channels and Exclude Edges are dependent on segmentation being enabled
    is_segmentation_enabled = handles.nucSegEnabledCheck.Value;
    
    % Enable/disable Apply to Channels controls based on segmentation enabled state
    handles.nucInclusionExclusionApplyLabel.Enable = matlab.lang.OnOffSwitchState(is_segmentation_enabled);
    handles.nucInclusionExclusionApplyDrop.Enable = matlab.lang.OnOffSwitchState(is_segmentation_enabled);
    
    % Update the apply to dropdown based on number of active channels
    numChannels = str2double(handles.numChanDrop.Value);
    current_items = handles.nucInclusionExclusionApplyDrop.Items;
    
    % Create new items list
    new_items = {'All Channels'};
    for k = 1:numChannels
        new_items{end+1} = ['Channel ' num2str(k)];
    end
    
    % Update items if they've changed
    if ~isequal(current_items, new_items)
        current_value = handles.nucInclusionExclusionApplyDrop.Value;
        handles.nucInclusionExclusionApplyDrop.Items = new_items;
        
        % Try to keep the same value, or default to 'All Channels'
        if ismember(current_value, new_items)
            handles.nucInclusionExclusionApplyDrop.Value = current_value;
        else
            handles.nucInclusionExclusionApplyDrop.Value = 'All Channels';
        end
    end
end

function updateNucleiBgCorrControls(handles)
    % Updates the nuclei background correction controls based on method selection
    % Similar to updatePreprocControls but for nuclei background correction
    
    if ~isfield(handles, 'nucBgCorrEnabledCheck') || ~isfield(handles, 'nucBgParamInput')
        return;
    end
    
    if ~handles.nucBgCorrEnabledCheck.Value
        % If background correction is disabled, disable all parameter inputs
        handles.nucBgParamInput.Enable = 'off';
        return;
    end
    
    if ~isfield(handles, 'nucBgMethodDrop') || ~isfield(handles, 'nucBgCorrModeDrop')
        return;
    end
    
    method = handles.nucBgMethodDrop.Value;
    mode = handles.nucBgCorrModeDrop.Value;
    is_3d_mode = strcmp(mode, '3D');
    
    % Update parameter visibility based on method
    switch method
        case 'None'
            % Hide parameter input
            handles.nucBgParamInput.Visible = 'off';
            
        case {'Gaussian', 'Rolling-ball', 'Top-hat'}
            % Show parameter input
            handles.nucBgParamInput.Visible = 'on';
            handles.nucBgParamInput.Enable = 'on';
            
            % Update tooltip based on method
            if strcmp(method, 'Gaussian')
                handles.nucBgParamInput.Tooltip = 'Gaussian: Sigma. Units are µm if "Scale" is checked (3D), otherwise pixels/voxels.';
            else
                handles.nucBgParamInput.Tooltip = 'Rolling-ball/Top-hat: Radius. Units are µm if "Scale" is checked (3D), otherwise pixels/voxels.';
            end
    end
end

function updateNucleiPreprocControls(handles)
    % Manages the dynamic parameter display for nuclei pre-processing methods.
    % Exactly matches the fluorescent channel preprocessing logic.
    if ~handles.nucPreprocEnabledCheck.Value
        return;
    end
    
    method = handles.nucPreprocMethodDrop.Value;
    is_3d_mode = strcmp(handles.nucPreprocessModeDrop.Value, '3D');
    
    % Handle 3D mode restrictions first
    if is_3d_mode
        if strcmp(method, 'Wavelet Denoising')
            handles.nucPreprocMethodDrop.Value = 'None';
            method = 'None';
        end
        % Disable the wavelet option in the dropdown
        handles.nucPreprocMethodDrop.Items = {'None', 'Gaussian', 'Median', 'Non-Local Means'};
    else
        handles.nucPreprocMethodDrop.Items = {'None', 'Gaussian', 'Median', 'Non-Local Means', 'Wavelet Denoising'};
    end
    
    % Get handles to dynamic parameter controls
    param1Label = handles.nucPreprocParam1Labels;
    param1Input = handles.nucPreprocParam1Inputs;
    param2Label = handles.nucPreprocParam2Labels;
    param2Input = handles.nucPreprocParam2Inputs;
    param3Label = handles.nucPreprocParam3Labels;
    param3Input = handles.nucPreprocParam3Inputs;
    param4Label = handles.nucPreprocParam4Labels;
    param4Input = handles.nucPreprocParam4Inputs;
    waveletDropdown = handles.nucWaveletNameDrop;
    
    % Configure controls based on selected method
    switch method
        case 'None'
            % Hide all parameter controls
            param1Label.Visible = 'off';
            param1Input.Visible = 'off';
            param2Label.Visible = 'off';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            % Hide NLM-specific controls
            handles.nucNlmFilterStrengthInput.Visible = 'off';
            handles.nucNlmSearchWindowInput.Visible = 'off';
            handles.nucNlmComparisonWindowInput.Visible = 'off';
            
        case 'Gaussian'
            % Show: "Sigma: [input] (units)"
            param1Label.Text = 'Sigma:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            param2Label.Text = '(units)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            
        case 'Median'
            % Show: "Size: [input] (units)"
            param1Label.Text = 'Size:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            param2Label.Text = '(units)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Visible = 'off';
            param3Input.Visible = 'off';
            param4Label.Visible = 'off';
            param4Input.Visible = 'off';
            waveletDropdown.Visible = 'off';
            
        case 'Non-Local Means'
            % Show all three NLM parameters:
            % Row 4: "Filter Strength: [input] (0.0-1.0)"
            % Row 5: "Search Window: [input]"
            % Row 6: "Comparison Window: [input]"
            param1Label.Text = 'Filter Strength:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'off'; % Use specific NLM controls instead
            param2Label.Text = '(0.0-1.0)';
            param2Label.Visible = 'on';
            param2Input.Visible = 'off';
            param3Label.Text = 'Search Window:';
            param3Label.Visible = 'on';
            param3Input.Visible = 'off'; % Use specific NLM controls instead
            param4Label.Text = 'Comparison Window:';
            param4Label.Visible = 'on';
            param4Input.Visible = 'off'; % Use specific NLM controls instead
            waveletDropdown.Visible = 'off';
            
            % Show NLM-specific controls
            handles.nucNlmFilterStrengthInput.Visible = 'on';
            handles.nucNlmSearchWindowInput.Visible = 'on';
            handles.nucNlmComparisonWindowInput.Visible = 'on';
            
        case 'Wavelet Denoising'
            % Show: "Wavelet: [dropdown] Level: [input]"
            %       "Rule: [dropdown] Method: [dropdown]"
            param1Label.Text = 'Wavelet:';
            param1Label.Visible = 'on';
            param1Input.Visible = 'off'; % Hide numeric input
            waveletDropdown.Visible = 'on'; % Show wavelet dropdown
            
            param2Label.Text = 'Level:';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
            
            param3Label.Text = 'Threshold Rule:';
            param3Label.Visible = 'on';
            param3Input.Visible = 'on';
            
            param4Label.Text = 'Threshold Method:';
            param4Label.Visible = 'on';
            param4Input.Visible = 'on';
            % Hide NLM-specific controls
            handles.nucNlmFilterStrengthInput.Visible = 'off';
            handles.nucNlmSearchWindowInput.Visible = 'off';
            handles.nucNlmComparisonWindowInput.Visible = 'off';
            
    end
    
    % Ensure NLM controls are hidden for all methods except Non-Local Means
    if ~strcmp(method, 'Non-Local Means')
        handles.nucNlmFilterStrengthInput.Visible = 'off';
        handles.nucNlmSearchWindowInput.Visible = 'off';
        handles.nucNlmComparisonWindowInput.Visible = 'off';
    end
    
end

function updateDisplayOnAllPreviewsControls(channel_idx, handles)
    % Updates the visibility and state of the "Display on All Previews" checkbox
    % This checkbox should only be enabled if maxima detection is enabled for this channel
    
    is_maxima_enabled = handles.maximaEnabledChecks(channel_idx).Value;
    handles.displayOnAllPreviewsChecks(channel_idx).Enable = matlab.lang.OnOffSwitchState(is_maxima_enabled);
end

function updateNucleiFilterControls(handles)
    % Updates the nuclei filter controls based on nuclei data dimensionality
    
    if ~handles.nucFilterEnabledCheck.Value
        return;
    end
    
    % Enable/disable size filter controls based on size filter checkbox
    is_size_filter_enabled = handles.nucFilterSizeEnabledCheck.Value;
    handles.nucFilterMinSizeInput.Enable = matlab.lang.OnOffSwitchState(is_size_filter_enabled);
    handles.nucFilterSizeUnitDrop.Enable = matlab.lang.OnOffSwitchState(is_size_filter_enabled);
    
    % Check if nuclei data is 3D by looking at Z spacing input enable state
    is_3d_image = isfield(handles, 'nucZSpacingInput') && strcmp(handles.nucZSpacingInput.Enable, 'on');
    
    % Get current value to preserve it if valid
    current_value = handles.nucFilterSizeUnitDrop.Value;
    
    if is_3d_image
        % 3D data: show voxels and microns^3
        available_units = {'voxels', 'microns^3'};
    else
        % 2D data: show pixels and microns^2
        available_units = {'pixels', 'microns^2'};
    end
    
    % Update the dropdown items
    handles.nucFilterSizeUnitDrop.Items = available_units;
    
    % Try to keep the current value if it's still valid, otherwise default to first option
    if ismember(current_value, available_units)
        handles.nucFilterSizeUnitDrop.Value = current_value;
    else
        handles.nucFilterSizeUnitDrop.Value = available_units{1};
    end
    
    % Enable/disable circularity input based on circularity checkbox
    is_circularity_enabled = handles.nucFilterCircularityEnabledCheck.Value;
    handles.nucFilterMinCircularityInput.Enable = matlab.lang.OnOffSwitchState(is_circularity_enabled);
    
    % Enable/disable solidity input based on solidity checkbox
    if isfield(handles, 'nucFilterSolidityEnabledCheck') && isvalid(handles.nucFilterSolidityEnabledCheck)
        is_solidity_enabled = handles.nucFilterSolidityEnabledCheck.Value;
        handles.nucFilterMinSolidityInput.Enable = matlab.lang.OnOffSwitchState(is_solidity_enabled);
    end
end

function updateNucleiSegmentationControls(handles)
    % Update visibility of nuclei segmentation method controls using dynamic parameter system
    if ~isfield(handles, 'nucSegMainMethodDrop')
        return;
    end
    
    main_method = handles.nucSegMainMethodDrop.Value;
    
    % Determine what should be visible
    is_absolute = strcmp(main_method, 'Absolute');
    is_mean_or_median = strcmp(main_method, 'Mean') || strcmp(main_method, 'Median');
    is_auto_local = strcmp(main_method, 'Auto Local Threshold');
    
    % Show/hide sub-method dropdown (for Mean/Median only)
    if isfield(handles, 'nucSegSubMethodLabel')
        handles.nucSegSubMethodLabel.Visible = matlab.lang.OnOffSwitchState(is_mean_or_median);
    end
    if isfield(handles, 'nucSegSubMethodDrop')
        handles.nucSegSubMethodDrop.Visible = matlab.lang.OnOffSwitchState(is_mean_or_median);
    end
    
    % Update dynamic parameter controls based on method selection
    updateSegmentationParameterControls(handles, main_method);
    
    % Show/hide Auto Local Threshold controls
    if isfield(handles, 'nucSegAlgorithmLabel')
        handles.nucSegAlgorithmLabel.Visible = matlab.lang.OnOffSwitchState(is_auto_local);
    end
    if isfield(handles, 'nucSegAlgorithmDrop')
        handles.nucSegAlgorithmDrop.Visible = matlab.lang.OnOffSwitchState(is_auto_local);
    end
    
    % Note: nucSegParam2Label and nucSegParam2Input removed - redundant with algorithm-specific parameters
    
    % Show/hide algorithm-specific parameters
    if is_auto_local && isfield(handles, 'nucSegAlgorithmDrop')
        algorithm = handles.nucSegAlgorithmDrop.Value;
        updateAlgorithmParameterControls(handles, algorithm);
    else
        % Hide algorithm-specific parameters when not using Auto Local Threshold
        if isfield(handles, 'nucSegAlgParamLabel')
            handles.nucSegAlgParamLabel.Visible = 'off';
        end
        if isfield(handles, 'nucSegAlgParamInput')
            handles.nucSegAlgParamInput.Visible = 'off';
        end
        if isfield(handles, 'nucSegAlgParamDefaultCheck')
            handles.nucSegAlgParamDefaultCheck.Visible = 'off';
        end
        % Hide second parameter controls
        if isfield(handles, 'nucSegAlgParam2Label')
            handles.nucSegAlgParam2Label.Visible = 'off';
        end
        if isfield(handles, 'nucSegAlgParam2Input')
            handles.nucSegAlgParam2Input.Visible = 'off';
        end
        if isfield(handles, 'nucSegAlgParam2DefaultCheck')
            handles.nucSegAlgParam2DefaultCheck.Visible = 'off';
        end
    end
    
    % Update inclusion/exclusion controls
    updateInclusionExclusionControls(handles);
end

function updateSegmentationParameterControls(handles, main_method)
    % Update the dynamic parameter controls based on selected method
    if ~isfield(handles, 'nucSegParam1Label') || ~isfield(handles, 'nucSegParam1Input')
        return;
    end
    
    switch main_method
        case 'Absolute'
            handles.nucSegParam1Label.Text = 'Threshold:';
            handles.nucSegParam1Input.Tooltip = 'Direct intensity threshold. Pixels above this value = nuclei. Higher values = smaller nuclei regions.';
            handles.nucSegParam1Input.Value = handles.nucSegAbsoluteInput.Value;
            handles.nucSegParam1Input.Visible = 'on';
            handles.nucSegParam2Label.Visible = 'off';
            handles.nucSegParam2Input.Visible = 'off';
            
        case {'Mean', 'Median'}
            if isfield(handles, 'nucSegSubMethodDrop')
                sub_method = handles.nucSegSubMethodDrop.Value;
                if strcmp(sub_method, 'Std Multiplier')
                    handles.nucSegParam1Label.Text = 'Std×:';
                    handles.nucSegParam1Input.Tooltip = 'Standard deviation multiplier. Threshold = mean/median + k×std. Higher k = less sensitive. Range: 1.5-3.0';
                    handles.nucSegParam1Input.Value = handles.nucSegStdMultiplierInput.Value;
                else % Absolute Offset
                    handles.nucSegParam1Label.Text = 'Offset:';
                    handles.nucSegParam1Input.Tooltip = 'Fixed intensity offset. Threshold = mean/median + offset. Use when nuclei are consistently X units brighter than background.';
                    handles.nucSegParam1Input.Value = handles.nucSegOffsetInput.Value;
                end
            end
            handles.nucSegParam1Input.Visible = 'on';
            handles.nucSegParam2Label.Visible = 'off';
            handles.nucSegParam2Input.Visible = 'off';
            
        case 'Auto Local Threshold'
            handles.nucSegParam1Label.Text = 'Radius:';
            handles.nucSegParam1Input.Tooltip = 'Local neighborhood radius (pixels). Larger = smoother thresholds. Range: 7-50';
            % Note: nucSegLocalRadiusInput removed - radius now handled by algorithm-specific parameters
            handles.nucSegParam1Input.Visible = 'on';
            handles.nucSegParam2Label.Text = '(pixels)';
            handles.nucSegParam2Label.Visible = 'on';
            handles.nucSegParam2Input.Visible = 'off';
    end
end

function updateAlgorithmParameterControls(handles, algorithm)
    % Update algorithm-specific parameter controls
    if ~isfield(handles, 'nucSegAlgParamLabel') || ~isfield(handles, 'nucSegAlgParamInput')
        return;
    end
    
    % Default values for each algorithm (first parameter)
    algorithm_defaults = containers.Map({
        'Bernsen', 'Mean', 'Median', 'MidGrey', 'Niblack', 'Otsu', 'Phansalkar', 'Sauvola'
    }, {
        15, 0, 0, 0, 0.2, 0, 0.25, 0.5
    });
    
    % Default values for dual-parameter algorithms (second parameter)
    algorithm_defaults2 = containers.Map({
        'Niblack', 'Phansalkar', 'Sauvola'
    }, {
        0, 0.5, 128
    });
    
    % Parameter names for each algorithm
    algorithm_names = containers.Map({
        'Bernsen', 'Mean', 'Median', 'MidGrey', 'Niblack', 'Otsu', 'Phansalkar', 'Sauvola'
    }, {
        'Contrast', 'C Value', 'C Value', 'C Value', 'K Value', 'N/A', 'K Value', 'K Value'
    });
    
    % Second parameter names for dual-parameter algorithms
    algorithm_names2 = containers.Map({
        'Niblack', 'Phansalkar', 'Sauvola'
    }, {
        'C Value', 'R Value', 'R Value'
    });
    
    % Check if this is a dual-parameter algorithm
    is_dual_param = isKey(algorithm_defaults2, algorithm);
    
    if isKey(algorithm_names, algorithm)
        handles.nucSegAlgParamLabel.Text = algorithm_names(algorithm);
        handles.nucSegAlgParamLabel.Visible = 'on';
        
        if strcmp(algorithm, 'Otsu')
            handles.nucSegAlgParamInput.Visible = 'off';
            handles.nucSegAlgParamDefaultCheck.Visible = 'off';
            % Hide second parameter for Otsu
            if isfield(handles, 'nucSegAlgParam2Label')
                handles.nucSegAlgParam2Label.Visible = 'off';
                handles.nucSegAlgParam2Input.Visible = 'off';
                handles.nucSegAlgParam2DefaultCheck.Visible = 'off';
            end
        else
            % Only set default value if default checkbox is checked, otherwise preserve user value
            if handles.nucSegAlgParamDefaultCheck.Value
                handles.nucSegAlgParamInput.Value = algorithm_defaults(algorithm);
            end
            handles.nucSegAlgParamInput.Visible = 'on';
            handles.nucSegAlgParamDefaultCheck.Visible = 'on';
            
            % Handle second parameter for dual-parameter algorithms
            if is_dual_param && isfield(handles, 'nucSegAlgParam2Label')
                handles.nucSegAlgParam2Label.Text = algorithm_names2(algorithm);
                handles.nucSegAlgParam2Label.Visible = 'on';
                % Only set default value if default checkbox is checked, otherwise preserve user value
                if handles.nucSegAlgParam2DefaultCheck.Value
                    handles.nucSegAlgParam2Input.Value = algorithm_defaults2(algorithm);
                end
                handles.nucSegAlgParam2Input.Visible = 'on';
                handles.nucSegAlgParam2DefaultCheck.Visible = 'on';
            elseif isfield(handles, 'nucSegAlgParam2Label')
                % Hide second parameter for single-parameter algorithms
                handles.nucSegAlgParam2Label.Visible = 'off';
                handles.nucSegAlgParam2Input.Visible = 'off';
                handles.nucSegAlgParam2DefaultCheck.Visible = 'off';
            end
        end
    else
        handles.nucSegAlgParamLabel.Visible = 'off';
        handles.nucSegAlgParamInput.Visible = 'off';
        handles.nucSegAlgParamDefaultCheck.Visible = 'off';
        % Hide second parameter
        if isfield(handles, 'nucSegAlgParam2Label')
            handles.nucSegAlgParam2Label.Visible = 'off';
            handles.nucSegAlgParam2Input.Visible = 'off';
            handles.nucSegAlgParam2DefaultCheck.Visible = 'off';
        end
    end
end

function updateFitFilterControls(channel_idx, handles)
    % Updates the fit filtering controls based on enabled state and fitting method
    % Throws error if radial symmetry fitting is used
    
    if ~handles.fitFilterEnabledChecks(channel_idx).Value
        % If fit filtering is disabled, disable all parameter inputs
        handles.fitFilterRSquaredMinInputs(channel_idx).Enable = 'off';
        handles.fitFilterRSquaredMaxInputs(channel_idx).Enable = 'off';
        handles.fitFilterSigmaSumMinInputs(channel_idx).Enable = 'off';
        handles.fitFilterSigmaSumMaxInputs(channel_idx).Enable = 'off';
        handles.fitFilterAmplitudeMinInputs(channel_idx).Enable = 'off';
        handles.fitFilterAmplitudeMaxInputs(channel_idx).Enable = 'off';
        handles.fitFilterIntensityMinInputs(channel_idx).Enable = 'off';
        handles.fitFilterIntensityMaxInputs(channel_idx).Enable = 'off';
        return;
    end
    
    % Check if radial symmetry fitting is used - throw error if so
    fitting_method = handles.gaussFitMethodDrop(channel_idx).Value;
    if strcmp(fitting_method, 'Radial Symmetry')
        error('Fit filtering is not supported with Radial Symmetry fitting. Please use Gaussian fitting instead.');
    end
    
    % Enable/disable individual filter controls based on their enabled state
    r_squared_enabled = handles.fitFilterRSquaredEnabledChecks(channel_idx).Value;
    handles.fitFilterRSquaredMinInputs(channel_idx).Enable = r_squared_enabled;
    handles.fitFilterRSquaredMaxInputs(channel_idx).Enable = r_squared_enabled;
    
    sigma_sum_enabled = handles.fitFilterSigmaSumEnabledChecks(channel_idx).Value;
    handles.fitFilterSigmaSumMinInputs(channel_idx).Enable = sigma_sum_enabled;
    handles.fitFilterSigmaSumMaxInputs(channel_idx).Enable = sigma_sum_enabled;
    
    amplitude_enabled = handles.fitFilterAmplitudeEnabledChecks(channel_idx).Value;
    handles.fitFilterAmplitudeMinInputs(channel_idx).Enable = amplitude_enabled;
    handles.fitFilterAmplitudeMaxInputs(channel_idx).Enable = amplitude_enabled;
    
    intensity_enabled = handles.fitFilterIntensityEnabledChecks(channel_idx).Value;
    handles.fitFilterIntensityMinInputs(channel_idx).Enable = intensity_enabled;
    handles.fitFilterIntensityMaxInputs(channel_idx).Enable = intensity_enabled;
end

function updateDeconvControls(channel_idx, handles)
    % Manages the dynamic parameter display for deconvolution methods.
    if ~handles.deconvEnabledChecks(channel_idx).Value
        return;
    end
    
    method = handles.deconvMethodDrops(channel_idx).Value;
    psf_source = handles.deconvPSFSourceDrops(channel_idx).Value;
    
    % Get handles to dynamic labels
    param1Label = handles.deconvParam1Labels(channel_idx);
    param2Label = handles.deconvParam2Labels(channel_idx);
    param1Input = handles.deconvLRIterationsInputs(channel_idx);
    param2Input = handles.deconvLRDampingInputs(channel_idx);
    
    % Update labels and visibility based on selected method
    switch method
        case 'Lucy-Richardson'
            param1Label.Text = 'Iterations:';
            param1Input.Tooltip = 'Number of iterations for Lucy-Richardson deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Text = 'Damping:';
            param2Input.Tooltip = 'Damping parameter (0 = no damping).';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
            
        case 'Wiener'
            param1Label.Text = 'NSR:';
            param1Input.Tooltip = 'Noise-to-Signal Ratio for Wiener deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Visible = 'off';
            param2Input.Visible = 'off';
            
        case 'Blind'
            param1Label.Text = 'Iterations:';
            param1Input.Tooltip = 'Number of iterations for blind deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Text = 'Under-Relax:';
            param2Input.Tooltip = 'Under-relaxation parameter for blind deconvolution.';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
    end
    
    % Show/hide PSF controls based on PSF source
    psf_file_controls = findall(handles.fig, 'Tag', ['deconv_psf_file_' num2str(channel_idx)]);
    psf_gen_controls = findall(handles.fig, 'Tag', ['deconv_psf_gen_' num2str(channel_idx)]);
    
    if strcmp(psf_source, 'Load File')
        set(psf_file_controls, 'Visible', 'on');
        set(psf_gen_controls, 'Visible', 'off');
    else % 'Generate'
        set(psf_file_controls, 'Visible', 'off');
        set(psf_gen_controls, 'Visible', 'on');
    end
end

function updateNucleiDeconvControls(handles)
    % Manages the dynamic parameter display for nuclei deconvolution methods.
    if ~handles.nucDeconvEnabledCheck.Value
        return;
    end
    
    method = handles.nucDeconvMethodDrop.Value;
    psf_source = handles.nucDeconvPSFSourceDrop.Value;
    
    % Get handles to dynamic labels
    param1Label = handles.nucDeconvParam1Label;
    param2Label = handles.nucDeconvParam2Label;
    param1Input = handles.nucDeconvLRIterationsInput;
    param2Input = handles.nucDeconvLRDampingInput;
    
    % Update labels and visibility based on selected method
    switch method
        case 'Lucy-Richardson'
            param1Label.Text = 'Iterations:';
            param1Input.Tooltip = 'Number of iterations for Lucy-Richardson deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Text = 'Damping:';
            param2Input.Tooltip = 'Damping parameter (0 = no damping).';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
            
        case 'Wiener'
            param1Label.Text = 'NSR:';
            param1Input.Tooltip = 'Noise-to-Signal Ratio for Wiener deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Visible = 'off';
            param2Input.Visible = 'off';
            
        case 'Blind'
            param1Label.Text = 'Iterations:';
            param1Input.Tooltip = 'Number of iterations for blind deconvolution.';
            param1Label.Visible = 'on';
            param1Input.Visible = 'on';
            
            param2Label.Text = 'Under-Relax:';
            param2Input.Tooltip = 'Under-relaxation parameter for blind deconvolution.';
            param2Label.Visible = 'on';
            param2Input.Visible = 'on';
    end
    
    % Show/hide PSF controls based on PSF source
    psf_file_controls = findall(handles.fig, 'Tag', 'nuc_deconv_psf_file');
    psf_gen_controls = findall(handles.fig, 'Tag', 'nuc_deconv_psf_gen');
    
    if strcmp(psf_source, 'Load File')
        set(psf_file_controls, 'Visible', 'on');
        set(psf_gen_controls, 'Visible', 'off');
    else % 'Generate'
        set(psf_file_controls, 'Visible', 'off');
        set(psf_gen_controls, 'Visible', 'on');
    end
end