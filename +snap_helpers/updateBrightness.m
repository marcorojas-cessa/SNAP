function updateBrightness(slider, fig_handle, preview_idx)
% Adjusts the CLim of the associated axes based on the slider value.
    handles = guidata(fig_handle);

    % Use the preview_idx to get the correct axes
    ax = handles.previewAxes(preview_idx);

    % Check if the axis has an image and original CLim data stored
    if ~isempty(ax.Children) && isfield(ax.UserData, 'OriginalCLim') && ~isempty(ax.UserData.OriginalCLim)
        originalCLim = ax.UserData.OriginalCLim;
        minValue = originalCLim(1);
        originalRange = originalCLim(2) - minValue;
        
        % Adjust the upper limit based on the slider value.
        % A value of 1 is normal. < 1 makes it brighter, > 1 makes it darker.
        newMaxValue = minValue + (originalRange / slider.Value);
        
        % Ensure the new max is not less than the min
        if newMaxValue <= minValue
            newMaxValue = minValue + eps; % Add a small value
        end

        ax.CLim = [minValue, newMaxValue];
    end
end