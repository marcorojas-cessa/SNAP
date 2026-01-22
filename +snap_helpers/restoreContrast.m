function restoreContrast(ax, axRef)
% Restores contrast for the given axes to the full min/max of its data.
if isempty(axRef.Children) || ~isprop(axRef.Children, 'CData'), return; end
img = axRef.Children.CData;
if isempty(img) || ~isnumeric(img) || ~isreal(img) || ~all(isfinite(img(:))), return; end

minVal = min(img(:));
maxVal = max(img(:));
if minVal < maxVal
    ax.CLim = [minVal, maxVal];
else
    ax.CLim = [minVal, minVal + 1]; % Handle constant image
end
drawnow;
end