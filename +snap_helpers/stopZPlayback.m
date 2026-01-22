function stopZPlayback(fig_handle)
% Stops automatic playback through z-frames and cleans up the timer
    
    handles = guidata(fig_handle);
    if isempty(handles)
        return;
    end
    
    % Stop and delete the timer if it exists
    if isfield(handles, 'zPlaybackTimer') && isvalid(handles.zPlaybackTimer)
        stop(handles.zPlaybackTimer);
        delete(handles.zPlaybackTimer);
        handles = rmfield(handles, 'zPlaybackTimer');
    end
    
    % Update button states
    handles.pauseButton.Enable = 'off';
    handles.pauseButton.Text = '⏸ Pause';
    
    % Only re-enable play button if there's valid z-stack data
    if handles.globalZSlider.Enable && handles.globalZSlider.Limits(2) > 1
        handles.playButton.Enable = 'on';
        handles.playButton.Text = '▶ Play';
    else
        handles.playButton.Enable = 'off';
        handles.playButton.Text = '▶ Play';
    end
    
    % Update status
    handles.statusLabel.Text = 'Status: Idle';
    
    % Store handles
    guidata(fig_handle, handles);
end
