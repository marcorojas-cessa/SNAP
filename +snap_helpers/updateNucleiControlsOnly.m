function updateNucleiControlsOnly(fig_handle)
% Updates only nuclei-related controls for faster response when nuclei checkboxes change
    handles = guidata(fig_handle);
    
    % Prevent infinite recursion
    if isfield(handles, 'isUpdatingNucleiControls') && handles.isUpdatingNucleiControls
        return;
    end
    
    % Set flag to prevent recursion
    handles.isUpdatingNucleiControls = true;
    guidata(fig_handle, handles);
    
    try
        % Update only nuclei segmentation controls
        updateNucleiSegmentationControls(handles);
        
        % Update nuclei preprocessing controls
        updateNucleiPreprocControls(handles);
        
        % Update nuclei deconvolution controls
        if isfield(handles, 'nucDeconvMethodDrop')
            updateNucleiDeconvControls(handles);
        end
        
        % Update nuclei background correction controls
        if isfield(handles, 'nucBgCorrModeDrop')
            updateNucleiBgControls(handles);
        end
        
        % Update nuclei filtering controls
        updateNucleiFilteringControls(handles);
        
    catch ME
        warning('Error updating nuclei controls: %s', ME.message);
    end
    
    % Clear the flag
    handles.isUpdatingNucleiControls = false;
    guidata(fig_handle, handles);
end

function updateNucleiSegmentationControls(handles)
% Update nuclei segmentation method controls
    if isfield(handles, 'nucSegMethodDrop')
        method = handles.nucSegMethodDrop.Value;
        is_3d = strcmp(method, '3D Local Threshold');
        
        % Enable/disable 3D-specific controls
        if isfield(handles, 'nucSeg3DParam1Input')
            handles.nucSeg3DParam1Input.Enable = matlab.lang.OnOffSwitchState(is_3d);
        end
        if isfield(handles, 'nucSeg3DParam2Input')
            handles.nucSeg3DParam2Input.Enable = matlab.lang.OnOffSwitchState(is_3d);
        end
    end
end

function updateNucleiPreprocControls(handles)
% Update nuclei preprocessing controls
    if isfield(handles, 'nucPreprocessModeDrop')
        mode = handles.nucPreprocessModeDrop.Value;
        is_projection = contains(mode, 'Projection');
        is_3d = strcmp(mode, '3D');
        
        % Update projection dropdown visibility
        if isfield(handles, 'nucPreprocessProjectionDrop')
            handles.nucPreprocessProjectionDrop.Visible = is_projection;
        end
        
        % Update 3D-specific controls
        if isfield(handles, 'nucPreprocessScaleCheck')
            handles.nucPreprocessScaleCheck.Visible = is_3d;
        end
    end
end

function updateNucleiDeconvControls(handles)
% Update nuclei deconvolution controls
    if isfield(handles, 'nucDeconvMethodDrop')
        method = handles.nucDeconvMethodDrop.Value;
        
        % Update parameter visibility based on method
        if isfield(handles, 'nucDeconvParam1Input')
            handles.nucDeconvParam1Input.Visible = true;
        end
        if isfield(handles, 'nucDeconvParam2Input')
            handles.nucDeconvParam2Input.Visible = strcmp(method, 'Blind');
        end
    end
end

function updateNucleiBgControls(handles)
% Update nuclei background correction controls
    if isfield(handles, 'nucBgCorrModeDrop')
        mode = handles.nucBgCorrModeDrop.Value;
        is_projection = contains(mode, 'Projection');
        is_3d = strcmp(mode, '3D');
        
        % Update projection dropdown visibility
        if isfield(handles, 'nucBgCorrProjectionDrop')
            handles.nucBgCorrProjectionDrop.Visible = is_projection;
        end
        
        % Update 3D-specific controls
        if isfield(handles, 'nucBgCorrScaleCheck')
            handles.nucBgCorrScaleCheck.Visible = is_3d;
        end
    end
end

function updateNucleiFilteringControls(handles)
% Update nuclei filtering controls
    % Update filter parameter visibility based on enabled filters
    if isfield(handles, 'nucFilterSizeEnabledCheck')
        size_enabled = handles.nucFilterSizeEnabledCheck.Value;
        if isfield(handles, 'nucFilterMinSizeInput')
            handles.nucFilterMinSizeInput.Enable = matlab.lang.OnOffSwitchState(size_enabled);
        end
        if isfield(handles, 'nucFilterMaxSizeInput')
            handles.nucFilterMaxSizeInput.Enable = matlab.lang.OnOffSwitchState(size_enabled);
        end
    end
    
    if isfield(handles, 'nucFilterCircularityEnabledCheck')
        circ_enabled = handles.nucFilterCircularityEnabledCheck.Value;
        if isfield(handles, 'nucFilterMinCircularityInput')
            handles.nucFilterMinCircularityInput.Enable = matlab.lang.OnOffSwitchState(circ_enabled);
        end
    end
    
    if isfield(handles, 'nucFilterSolidityEnabledCheck')
        sol_enabled = handles.nucFilterSolidityEnabledCheck.Value;
        if isfield(handles, 'nucFilterMinSolidityInput')
            handles.nucFilterMinSolidityInput.Enable = matlab.lang.OnOffSwitchState(sol_enabled);
        end
    end
end
