function startZPlayback(fig_handle)
% Starts automatic playback through z-frames with a timer
    
    handles = guidata(fig_handle);
    if isempty(handles)
        return;
    end
    
    % Check if the global z slider is enabled and has valid limits
    if ~handles.globalZSlider.Enable || handles.globalZSlider.Limits(2) <= 1
        % No z-stack data available - use the same logic as updateLivePreview
        max_z = 1;
        if isfield(handles, 'rawDIC') && ~isempty(handles.rawDIC)
            max_z = max(max_z, size(handles.rawDIC, 3));
        end
        
        if isfield(handles, 'rawNuclei') && ~isempty(handles.rawNuclei)
            max_z = max(max_z, size(handles.rawNuclei, 3));
        end
        
        for k = 1:handles.Nmax
            if isfield(handles, 'rawChannels') && length(handles.rawChannels) >= k && ~isempty(handles.rawChannels{k})
                max_z = max(max_z, size(handles.rawChannels{k}, 3));
            end
        end
        
        if max_z <= 1
            warndlg('No z-stack data available for playback. Please load some data and update previews first.', 'Playback Error');
            return;
        end
        
        % Update the global z slider with the detected z-stack size
        handles.globalZSlider.Limits = [1, max_z];
        handles.globalZSlider.Enable = 'on';
        handles.globalZLabel.Text = sprintf('Global Z: 1 / %d', max_z);
    end
    
    % Create a timer for z-frame playback
    if isfield(handles, 'zPlaybackTimer') && isvalid(handles.zPlaybackTimer)
        stop(handles.zPlaybackTimer);
        delete(handles.zPlaybackTimer);
    end
    
    % Create new timer
    handles.zPlaybackTimer = timer('ExecutionMode', 'fixedRate', ...
                                   'Period', 0.3, ... % ~3.3 FPS - slower for better UI responsiveness
                                   'TimerFcn', @(src, event) snap_helpers.advanceZFrame(fig_handle));
    
    % Start the timer
    start(handles.zPlaybackTimer);
    
    % Update button states
    handles.playButton.Enable = 'off';
    handles.pauseButton.Enable = 'on';
    handles.playButton.Text = '▶ Playing';
    handles.pauseButton.Text = '⏸ Pause';
    
    % Update status
    handles.statusLabel.Text = 'Status: Playing Z-frames';
    
    % Store handles
    guidata(fig_handle, handles);
end
