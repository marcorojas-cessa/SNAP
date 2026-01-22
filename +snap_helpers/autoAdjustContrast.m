function autoAdjustContrast(ax)
% Auto-adjusts contrast for the given axes using stretchlim.
if isempty(ax.Children) || ~isprop(ax.Children, 'CData'), return; end
img = ax.Children.CData;
if isempty(img) || ~isnumeric(img) || ~isreal(img) || ~all(isfinite(img(:))), return; end

img_double = double(img);
if all(img_double(:) == img_double(1)) % Handle constant image
    val = img_double(1);
    ax.CLim = [val, val + 1];
    return;
end

try
    limits = stretchlim(img_double);
    if limits(1) < limits(2)
        ax.CLim = [limits(1), limits(2)];
    end
catch
    % Fallback for safety
    minVal = min(img_double(:));
    maxVal = max(img_double(:));
    if minVal < maxVal
        ax.CLim = [minVal, maxVal];
    end
end
drawnow;
end