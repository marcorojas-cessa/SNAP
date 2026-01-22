function synchronizeZSliders(fig)
% Synchronizes the Z-position of all preview windows set to 'Z-Stack' mode.

    handles = guidata(fig);
    
    % Get the value from the global slider
    globalZValue = round(handles.globalZSlider.Value);
    
    % Update the global slider's label
    maxZ = handles.globalZSlider.Limits(2);
    handles.globalZLabel.Text = sprintf('Global Z: %d / %d', globalZValue, maxZ);
    
    % Iterate through all preview panels
    for i = 1:5
        % Check if the preview is in Z-Stack mode
        if strcmp(handles.previewModeDrops(i).Value, 'Z-Stack')
            
            % Update the individual slider's value
            handles.zSliders(i).Value = globalZValue;
            
            % Update the z-label for this preview
            handles.zLabels(i).Text = sprintf('Z: %d / %d', globalZValue, maxZ);
            
            % Redraw the preview to show the new Z-slice
            snap_helpers.redrawPreview(fig, i);
        end
    end
    
    % Store updated handles
    guidata(fig, handles);
end
