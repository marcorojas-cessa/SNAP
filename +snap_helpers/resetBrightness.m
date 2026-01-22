function resetBrightness(button, fig_handle, preview_idx)
% Resets the associated brightness slider to 1 and updates the axes.
    handles = guidata(fig_handle);

    % Use the preview_idx to find the correct slider
    slider = handles.brightnessSliders(preview_idx);

    % Set slider value to default
    slider.Value = 1;
    
    % Call the update function to apply the change
    snap_helpers.updateBrightness(slider, handles.fig, preview_idx);
end