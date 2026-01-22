function advanceZFrame(fig_handle)
% Advances to the next z-frame in the playback sequence
    
    try
        handles = guidata(fig_handle);
        if isempty(handles)
            return;
        end
        
        % Get current z position
        current_z = round(handles.globalZSlider.Value);
        max_z = handles.globalZSlider.Limits(2);
        min_z = handles.globalZSlider.Limits(1);
        
        % Advance to next frame (loop back to min_z if at end)
        if current_z >= max_z
            next_z = min_z;  % Loop back to start
        else
            next_z = current_z + 1;
        end
        
        % Ensure next_z is within valid range
        next_z = max(min_z, min(max_z, next_z));
        
        % Update global z slider
        handles.globalZSlider.Value = next_z;
        
        % Update global z label
        handles.globalZLabel.Text = sprintf('Global Z: %d / %d', next_z, max_z);
        
        % Trigger z-slider synchronization to update all previews
        snap_helpers.synchronizeZSliders(fig_handle);
        
        % Force a drawnow to ensure UI updates are visible
        drawnow;
        
    catch ME
        warning('Error advancing z-frame: %s', ME.message);
        % Stop playback on error
        if isfield(handles, 'zPlaybackTimer') && isvalid(handles.zPlaybackTimer)
            stop(handles.zPlaybackTimer);
            delete(handles.zPlaybackTimer);
            handles = rmfield(handles, 'zPlaybackTimer');
            
            % Reset button states
            handles.playButton.Enable = 'on';
            handles.pauseButton.Enable = 'off';
            handles.playButton.Text = '▶ Play';
            handles.pauseButton.Text = '⏸ Pause';
            
            % Update status
            handles.statusLabel.Text = 'Status: Playback stopped due to error';
            
            guidata(fig_handle, handles);
        end
    end
end
