function displayGaussianFitPlot(handles, channel_idx, maxima_idx)
% Creates a new figure and displays detailed plots of a Gaussian fit.
% channel_idx: The channel number for color selection

fprintf('\n=== DISPLAY GAUSSIAN FIT CALLED ===\n');
fprintf('Channel: %d, Maxima: %d\n', channel_idx, maxima_idx);

% --- 1. Retrieve Data ---
if isempty(handles.gaussFitResults{channel_idx})
    disp('No Gaussian fit results available for this channel.');
    return;
end

if maxima_idx > length(handles.gaussFitResults{channel_idx})
    disp('Maximum index out of range.');
    return;
end

fit_result = handles.gaussFitResults{channel_idx}(maxima_idx);
if isempty(fit_result) || ~isfield(fit_result, 'amplitude')
    disp('No valid fit result found for this maximum.');
    disp(['Channel: ' num2str(channel_idx) ', Maxima: ' num2str(maxima_idx)]);
    disp(['Results available: ' num2str(length(handles.gaussFitResults{channel_idx}))]);
    if ~isempty(fit_result)
        disp('Available fields:');
        disp(fieldnames(fit_result));
    end
    return;
end

% Check if this is a valid fit result based on the fitting method
is_valid_fit = false;
if isfield(fit_result, 'fitMethod') && strcmp(fit_result.fitMethod, 'Radial Symmetry')
    % For radial symmetry, check if we have center coordinates (even if NaN, still allow visualization)
    if isfield(fit_result, 'center_x') && isfield(fit_result, 'center_y')
        is_valid_fit = true;
    end
else
    % For Gaussian fits, check if amplitude is not NaN
    if ~isnan(fit_result.amplitude)
        is_valid_fit = true;
    end
end

if ~is_valid_fit
    disp('No valid fit result found for this maximum.');
    disp(['Channel: ' num2str(channel_idx) ', Maxima: ' num2str(maxima_idx)]);
    disp(['Results available: ' num2str(length(handles.gaussFitResults{channel_idx}))]);
    if ~isempty(fit_result)
        disp('Available fields:');
        disp(fieldnames(fit_result));
    end
    return;
end

% Get channel color for plotting
colorName = handles.maximaColorDrops(channel_idx).Value;
colorMap = containers.Map({'Red','Green','Blue','Yellow','Magenta','Cyan'}, {'r','g','b','y','m','c'});
channelColor = colorMap(colorName);

% Create a new figure for the plot
fig = uifigure('Name', ['Gaussian Fit: Ch ' num2str(channel_idx) ', Maxima ' num2str(maxima_idx)], 'Position', [200 200 1000 800]);
mainGrid = uigridlayout(fig, [1 1]);

% --- 2. Plot Based on Fit Method ---
fprintf('Fit Method: %s\n', fit_result.fitMethod);
fprintf('===================================\n\n');

switch fit_result.fitMethod
    case {'1D (X,Y,Z)', '1D (X,Y,Z) Gaussian'}
        plot_1D_fits(mainGrid, fit_result);
    case {'2D (XY) + 1D (Z)', '2D (XY) + 1D (Z) Gaussian'}
        plot_2D_plus_1D_fits(mainGrid, fit_result);
    case {'3D Gaussian', 'Distorted 3D Gaussian'}
        plot_3D_fits(mainGrid, fit_result);
    case 'Radial Symmetry'
        plot_radial_symmetry_fits(mainGrid, fit_result, channelColor);
    otherwise
        uilabel(mainGrid, 'Text', 'Plotting for this fit method is not implemented.');
end

end

% --- Plotting Helper Functions ---

function plot_1D_fits(parent, result)
    
    fprintf('\n\n*** 1D FIT VISUALIZATION STARTED ***\n\n');
    
    data = result.rawDataWindow;
    
    % Create a grid layout with 1 row and 4 columns (3 plots + 1 info panel)
    grid = uigridlayout(parent, [1, 4]);
    grid.ColumnWidth = {'1x', '1x', '1x', 300};
    grid.RowHeight = {'1x'};
    
    % Plot X Profile
    ax_x = uiaxes(grid);
    ax_x.Layout.Column = 1;
    
    % Use the exact same X profile that was used during fitting
    if isfield(result, 'x_profile') && ~isempty(result.x_profile)
        x_profile = result.x_profile;
    else
        % Fallback to extracting profile if not stored (shouldn't happen for 1D fits)
        if size(data, 3) > 1  % 3D data
            center_y = round(size(data, 1)/2);
            center_z = round(size(data, 3)/2);
            x_profile = squeeze(data(center_y, :, center_z));
        else  % 2D data
            center_y = round(size(data, 1)/2);
            x_profile = squeeze(data(center_y, :));
        end
    end
    
    x_coords = 1:length(x_profile);
    plot(ax_x, x_coords, x_profile, 'b.', 'MarkerSize', 10);
    hold(ax_x, 'on');
    x_fine = linspace(1, length(x_profile), 200);
    
    % Use the exact fitted amplitude from the 1D Gaussian fit
    if isfield(result, 'amplitude_x') && ~isnan(result.amplitude_x)
        x_fitted_amplitude = result.amplitude_x;
    else
        % Fallback to profile peak if fitted amplitude not available
        x_fitted_amplitude = max(x_profile);
    end
    
    fit_x = x_fitted_amplitude * exp(-(x_fine - result.center_x).^2 / (2 * result.sigma_x^2));
    plot(ax_x, x_fine, fit_x, 'r-', 'LineWidth', 2);
    
    % Add center line for X at the FITTED Gaussian center
    % The fitted center represents where the Gaussian fit determined the peak should be
    % Use exact fitted center value for sub-pixel precision (no rounding)
    center_x_exact = result.center_x;
    if center_x_exact >= 1 && center_x_exact <= length(x_profile)
        % Use full range to ensure line goes through the center point
        plot(ax_x, [center_x_exact, center_x_exact], [0, max(x_profile)*1.1], 'g--', 'LineWidth', 2);
    end
    
    hold(ax_x, 'off');
    title(ax_x, 'X Profile Fit', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_x, 'X (pixels)', 'FontSize', 12);
    ylabel(ax_x, 'Intensity', 'FontSize', 12);
    % Legend created after all plot elements (Data, Fit, Center) are added
    legend(ax_x, 'Data', 'Fit', 'Center', 'Location', 'best');
    ax_x.GridAlpha = 0.3;
    ax_x.GridLineStyle = '-';
    
    % Plot Y Profile
    ax_y = uiaxes(grid);
    ax_y.Layout.Column = 2;
    
    % Use the exact same Y profile that was used during fitting
    if isfield(result, 'y_profile') && ~isempty(result.y_profile)
        y_profile = result.y_profile;
    else
        % Fallback to extracting profile if not stored (shouldn't happen for 1D fits)
        if size(data, 3) > 1  % 3D data
            center_x = round(size(data, 2)/2);
            center_z = round(size(data, 3)/2);
            y_profile = squeeze(data(:, center_x, center_z));
        else  % 2D data
            center_x = round(size(data, 2)/2);
            y_profile = squeeze(data(:, center_x));
        end
    end
    
    y_coords = 1:length(y_profile);
    plot(ax_y, y_coords, y_profile, 'b.', 'MarkerSize', 10);
    hold(ax_y, 'on');
    y_fine = linspace(1, length(y_profile), 200);
    
    % Use the exact fitted amplitude from the 1D Gaussian fit
    if isfield(result, 'amplitude_y') && ~isnan(result.amplitude_y)
        y_fitted_amplitude = result.amplitude_y;
    else
        % Fallback to profile peak if fitted amplitude not available
        y_fitted_amplitude = max(y_profile);
    end
    
    fit_y = y_fitted_amplitude * exp(-(y_fine - result.center_y).^2 / (2 * result.sigma_y^2));
    plot(ax_y, y_fine, fit_y, 'r-', 'LineWidth', 2);
    
    % Add center line for Y at the FITTED Gaussian center
    % The fitted center represents where the Gaussian fit determined the peak should be
    % Use exact fitted center value for sub-pixel precision (no rounding)
    center_y_exact = result.center_y;
    if center_y_exact >= 1 && center_y_exact <= length(y_profile)
        % Use full range to ensure line goes through the center point
        plot(ax_y, [center_y_exact, center_y_exact], [0, max(y_profile)*1.1], 'g--', 'LineWidth', 2);
    end
    
    hold(ax_y, 'off');
    title(ax_y, 'Y Profile Fit', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_y, 'Y (pixels)', 'FontSize', 12);
    ylabel(ax_y, 'Intensity', 'FontSize', 12);
    % Legend created after all plot elements (Data, Fit, Center) are added
    legend(ax_y, 'Data', 'Fit', 'Center', 'Location', 'best');
    ax_y.GridAlpha = 0.3;
    ax_y.GridLineStyle = '-';
    
    % Plot Z Profile
    ax_z = uiaxes(grid);
    ax_z.Layout.Column = 3;
    if ndims(data) == 3
        % Use the exact same Z profile that was used during fitting
        if isfield(result, 'z_profile') && ~isempty(result.z_profile)
            z_profile = result.z_profile;
        else
            % Fallback to extracting profile if not stored (shouldn't happen for 1D fits)
            center_x = round(size(data, 2)/2);
            center_y = round(size(data, 1)/2);
            z_profile = squeeze(data(center_y, center_x, :));
        end
        
        z_coords = 1:length(z_profile);
        plot(ax_z, z_coords, z_profile, 'b.', 'MarkerSize', 10);
        hold(ax_z, 'on');
        z_fine = linspace(1, length(z_profile), 200);
        
        % Use the exact fitted amplitude from the 1D Gaussian fit
        if isfield(result, 'amplitude_z') && ~isnan(result.amplitude_z)
            z_fitted_amplitude = result.amplitude_z;
        else
            % Fallback to profile peak if fitted amplitude not available
            z_fitted_amplitude = max(z_profile);
        end
        
        fit_z = z_fitted_amplitude * exp(-(z_fine - result.center_z).^2 / (2 * result.sigma_z^2));
        plot(ax_z, z_fine, fit_z, 'r-', 'LineWidth', 2);
        
            % Add center line for Z at the FITTED Gaussian center
            % The fitted center represents where the Gaussian fit determined the peak should be
            % Use exact fitted center value for sub-pixel precision (no rounding)
            center_z_exact = result.center_z;
            if center_z_exact >= 1 && center_z_exact <= length(z_profile)
                % Use full range to ensure line goes through the center point
                plot(ax_z, [center_z_exact, center_z_exact], [0, max(z_profile)*1.1], 'g--', 'LineWidth', 2);
            end
        
        hold(ax_z, 'off');
        title(ax_z, 'Z Profile Fit', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel(ax_z, 'Z (pixels)', 'FontSize', 12);
        ylabel(ax_z, 'Intensity', 'FontSize', 12);
        % Legend created after all plot elements (Data, Fit, Center) are added
        legend(ax_z, 'Data', 'Fit', 'Center', 'Location', 'best');
    else
        text(ax_z, 0.5, 0.5, 'No Z dimension', 'HorizontalAlignment', 'center', 'FontSize', 14);
        title(ax_z, 'Z Profile Fit', 'FontSize', 14, 'FontWeight', 'bold');
    end
    ax_z.GridAlpha = 0.3;
    ax_z.GridLineStyle = '-';
    
    % Create info panel
    info_panel = uipanel(grid);
    info_panel.Layout.Column = 4;
    info_panel.Title = 'Fit Parameters';
    info_panel.FontSize = 12;
    info_panel.FontWeight = 'bold';
    
    % Set proper layout properties for the panel
    info_panel.Layout.Row = 1;
    info_panel.Layout.Column = 4;
    
    % Add fit details with safe field access
    info_parts = {};
    
    % Fit Method
    if isfield(result, 'fitMethod')
        info_parts{end+1} = sprintf('Fit Method: %s', result.fitMethod);
    else
        info_parts{end+1} = 'Fit Method: Unknown';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Amplitude
    if isfield(result, 'amplitude') && ~isnan(result.amplitude)
        info_parts{end+1} = sprintf('Amplitude: %.2f', result.amplitude);
    else
        info_parts{end+1} = 'Amplitude: N/A';
    end
    
    % Centers
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        info_parts{end+1} = sprintf('Center X: %.2f', result.center_x);
    else
        info_parts{end+1} = 'Center X: N/A';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        info_parts{end+1} = sprintf('Center Y: %.2f', result.center_y);
    else
        info_parts{end+1} = 'Center Y: N/A';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        info_parts{end+1} = sprintf('Center Z: %.2f', result.center_z);
    else
        info_parts{end+1} = 'Center Z: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Sigmas
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        info_parts{end+1} = sprintf('Sigma X: %.2f', result.sigma_x);
    else
        info_parts{end+1} = 'Sigma X: N/A';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        info_parts{end+1} = sprintf('Sigma Y: %.2f', result.sigma_y);
    else
        info_parts{end+1} = 'Sigma Y: N/A';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        info_parts{end+1} = sprintf('Sigma Z: %.2f', result.sigma_z);
    else
        info_parts{end+1} = 'Sigma Z: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % R-squared
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared);
        elseif length(result.r_squared) >= 3
            info_parts{end+1} = sprintf('R² X: %.3f', result.r_squared(1));
            info_parts{end+1} = sprintf('R² Y: %.3f', result.r_squared(2));
            info_parts{end+1} = sprintf('R² Z: %.3f', result.r_squared(3));
        else
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared(1));
        end
    else
        info_parts{end+1} = 'R²: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Background
    if isfield(result, 'background') && ~isnan(result.background)
        info_parts{end+1} = sprintf('Background: %.2f', result.background);
    else
        info_parts{end+1} = 'Background: N/A';
    end
    
    % Integrated Intensity
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        info_parts{end+1} = sprintf('Integrated Intensity: %.2f', result.integratedIntensity);
    else
        info_parts{end+1} = 'Integrated Intensity: N/A';
    end
    
    % Create a table to display fit parameters
    % Prepare data for the table
    param_names = {};
    param_values = {};
    param_units = {};
    
    % Add fit parameters to the table
    if isfield(result, 'fitMethod')
        param_names{end+1} = 'Fit Method';
        param_values{end+1} = result.fitMethod;
        param_units{end+1} = '';
    end
    
    % Add equation information
    param_names{end+1} = '--- EQUATIONS ---';
    param_values{end+1} = '--- USED ---';
    param_units{end+1} = '---';
    
    param_names{end+1} = '1D Gaussian X';
    param_values{end+1} = 'f(x) = A_x*exp(-((x-μ_x)²)/(2σ_x²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Y';
    param_values{end+1} = 'f(y) = A_y*exp(-((y-μ_y)²)/(2σ_y²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Z';
    param_values{end+1} = 'f(z) = A_z*exp(-((z-μ_z)²)/(2σ_z²))';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Background';
    param_values{end+1} = 'Pre-subtracted (no offset)';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Constraints';
    param_values{end+1} = 'A>0, σ<window_size/2';
    param_units{end+1} = '';
    
    % Check if individual X, Y, Z amplitudes are available (for 1D fits)
    if isfield(result, 'amplitude_x') && isfield(result, 'amplitude_y') && isfield(result, 'amplitude_z')
        param_names{end+1} = 'Amplitude X';
        param_values{end+1} = sprintf('%.2f', result.amplitude_x);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Y';
        param_values{end+1} = sprintf('%.2f', result.amplitude_y);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Z';
        param_values{end+1} = sprintf('%.2f', result.amplitude_z);
        param_units{end+1} = 'counts';
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        param_names{end+1} = 'Amplitude';
        param_values{end+1} = sprintf('%.2f', result.amplitude);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        param_names{end+1} = 'Center X';
        param_values{end+1} = sprintf('%.2f', result.center_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        param_names{end+1} = 'Center Y';
        param_values{end+1} = sprintf('%.2f', result.center_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        param_names{end+1} = 'Center Z';
        param_values{end+1} = sprintf('%.2f', result.center_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        param_names{end+1} = 'Sigma X';
        param_values{end+1} = sprintf('%.2f', result.sigma_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        param_names{end+1} = 'Sigma Y';
        param_values{end+1} = sprintf('%.2f', result.sigma_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        param_names{end+1} = 'Sigma Z';
        param_values{end+1} = sprintf('%.2f', result.sigma_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared);
            param_units{end+1} = '';
        elseif length(result.r_squared) >= 3
            param_names{end+1} = 'R² X';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
            param_names{end+1} = 'R² Y';
            param_values{end+1} = sprintf('%.3f', result.r_squared(2));
            param_units{end+1} = '';
            param_names{end+1} = 'R² Z';
            param_values{end+1} = sprintf('%.3f', result.r_squared(3));
            param_units{end+1} = '';
        else
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
        end
    end
    
    if isfield(result, 'background') && ~isnan(result.background)
        param_names{end+1} = 'Background';
        param_values{end+1} = sprintf('%.2f', result.background);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        param_names{end+1} = 'Integrated Intensity';
        param_values{end+1} = sprintf('%.2f', result.integratedIntensity);
        param_units{end+1} = 'counts';
    end
    
    % Create the table first
    param_table = uitable(info_panel, ...
        'Data', [param_names', param_values', param_units'], ...
        'ColumnName', {'Parameter', 'Value', 'Unit'}, ...
        'ColumnWidth', {120, 100, 80}, ...
        'FontSize', 10, ...
        'RowName', []);
    
    % Try to set layout properties with error handling
    try
        param_table.Layout.Row = 1;
        param_table.Layout.Column = 1;
    catch ME
    end
    
    % Final visualization summary
    fprintf('\n=== VISUALIZATION SUMMARY ===\n');
    fprintf('Displayed amplitudes in graphs:\n');
    if exist('x_fitted_amplitude', 'var')
        fprintf('  X: %.2e\n', x_fitted_amplitude);
    end
    if exist('y_fitted_amplitude', 'var')
        fprintf('  Y: %.2e\n', y_fitted_amplitude);
    end
    if exist('z_fitted_amplitude', 'var')
        fprintf('  Z: %.2e\n', z_fitted_amplitude);
    end
    fprintf('Reported amplitudes in parameters table:\n');
    if isfield(result, 'amplitude_x')
        fprintf('  X: %.2e\n', result.amplitude_x);
    end
    if isfield(result, 'amplitude_y')
        fprintf('  Y: %.2e\n', result.amplitude_y);
    end
    if isfield(result, 'amplitude_z')
        fprintf('  Z: %.2e\n', result.amplitude_z);
    end
    fprintf('=== END VISUALIZATION DEBUG ===\n\n');
end

function plot_2D_plus_1D_fits(parent, result)
    
    data = result.rawDataWindow;
    
    % Create a grid layout with 1 row and 3 columns (2 plots + 1 info panel)
    grid = uigridlayout(parent, [1, 3]);
    grid.ColumnWidth = {'1x', '1x', 300};
    grid.RowHeight = {'1x'};
    
    % Plot 2D XY Fit
    ax_xy = uiaxes(grid);
    ax_xy.Layout.Column = 1;
    
    % Use background-corrected data for plotting
    % The background field contains the local background that was subtracted during fitting
    if isfield(result, 'background') && ~isnan(result.background)
        xy_proj = sum(data, 3) - result.background;
        fprintf('XY projection: Applied local background correction of %.2f\n', result.background);
    else
    xy_proj = sum(data, 3);
        fprintf('XY projection: No background correction available\n');
    end
    
    imagesc(ax_xy, xy_proj);
    axis(ax_xy, 'image');
    colormap(ax_xy, 'gray');
    hold(ax_xy, 'on');
    [Y, X] = meshgrid(1:size(xy_proj, 2), 1:size(xy_proj, 1));
    coords = [X(:), Y(:)];
    p = [result.amplitude, result.center_x, result.center_y, result.sigma_x, result.sigma_y, result.rho_xy];
    fit_vals = model_2d_distorted_for_plot(p, coords);
    fit_grid = reshape(fit_vals, size(xy_proj));
    
    % Note: No background correction needed here - fit parameters are already based on background-corrected data
    fprintf('2D contour: Using fitted values without additional background correction\n');
    
    % Validate fit_grid before plotting contour
    if any(~isfinite(fit_grid(:)))
        warning('Non-finite values detected in fit grid. Attempting to clean data.');
        fit_grid(~isfinite(fit_grid)) = 0;
    end
    
    % Debug: Print fit grid statistics
    fprintf('Fit grid stats: min=%.2e, max=%.2e, mean=%.2e, any_nan=%d, any_inf=%d\n', ...
            min(fit_grid(:)), max(fit_grid(:)), mean(fit_grid(:)), ...
            any(isnan(fit_grid(:))), any(isinf(fit_grid(:))));
    
    % Only plot contour if we have valid data
    if all(isfinite(fit_grid(:))) && max(fit_grid(:)) > 0
        try
            % Calculate contour levels based on amplitude
            if isfield(result, 'amplitude') && ~isnan(result.amplitude) && result.amplitude > 0
                contour_levels = [result.amplitude/2, result.amplitude*0.8];
                fprintf('2D contour: Using levels [%.2e, %.2e]\n', contour_levels(1), contour_levels(2));
            else
                % Fallback to percentage of max
                max_val = max(fit_grid(:));
                contour_levels = [max_val/2, max_val*0.8];
            end
            
            fprintf('Plotting contour with levels: [%.2e, %.2e]\n', contour_levels(1), contour_levels(2));
            contour(ax_xy, Y, X, fit_grid, contour_levels, 'r-', 'LineWidth', 2);
            fprintf('Contour plotted successfully\n');
        catch ME
            warning('Contour plotting failed: %s. Skipping contour.', ME.message);
            fprintf('Contour error details: %s\n', ME.message);
        end
    else
        warning('Invalid fit data detected. Skipping contour plot.');
        fprintf('Fit grid validation failed: isfinite=%d, max>0=%d\n', ...
                all(isfinite(fit_grid(:))), max(fit_grid(:)) > 0);
    end
    
    % Add center marker (large and clearly visible) at the ORIGINAL local maxima center
    if isfield(result, 'originalMaximaCoords')
        % Convert global coordinates to local window coordinates
        local_x = result.originalMaximaCoords(1) - (size(data, 2) - 1) / 2;
        local_y = result.originalMaximaCoords(2) - (size(data, 1) - 1) / 2;
        local_x_round = round(local_x);
        local_y_round = round(local_y);
        
        if local_x_round >= 1 && local_x_round <= size(xy_proj, 2) && ...
           local_y_round >= 1 && local_y_round <= size(xy_proj, 1)
            % Plot center marker: plot(x, y) where x=Y (left-right), y=X (up-down)
            % CORRECTED: X is up-down (row), Y is left-right (column)
            plot(ax_xy, local_y_round, local_x_round, 'r+', 'MarkerSize', 20, 'LineWidth', 3);
            fprintf('2D center marker: Placed at original local maxima (%.1f, %.1f) -> local (%.1f, %.1f)\n', ...
                    result.originalMaximaCoords(1), result.originalMaximaCoords(2), local_x_round, local_y_round);
            
            % Add center lines at ORIGINAL local maxima coordinates (passing through centroid)
            % X line: vertical line at x=local_y_round (left-right position)
            plot(ax_xy, [local_y_round, local_y_round], [1, size(xy_proj, 1)], 'g--', 'LineWidth', 2);
            % Y line: horizontal line at y=local_x_round (up-down position)
            plot(ax_xy, [1, size(xy_proj, 2)], [local_x_round, local_x_round], 'g--', 'LineWidth', 2);
        else
            fprintf('2D center marker: Original coordinates (%.1f, %.1f) map to local (%.1f, %.1f) - out of bounds\n', ...
                    result.originalMaximaCoords(1), result.originalMaximaCoords(2), local_x_round, local_y_round);
            % Fallback to fitted center
            % CORRECTED: X is up-down (row), Y is left-right (column)
            plot(ax_xy, result.center_y, result.center_x, 'r+', 'MarkerSize', 20, 'LineWidth', 3);
            % X line: vertical line at x=result.center_y (left-right position)
            plot(ax_xy, [result.center_y, result.center_y], [1, size(xy_proj, 1)], 'g--', 'LineWidth', 2);
            % Y line: horizontal line at y=result.center_x (up-down position)
            plot(ax_xy, [1, size(xy_proj, 2)], [result.center_x, result.center_x], 'g--', 'LineWidth', 2);
        end
    else
        % Fallback to fitted center if original coordinates not available
        % CORRECTED: X is up-down (row), Y is left-right (column)
        plot(ax_xy, result.center_y, result.center_x, 'r+', 'MarkerSize', 20, 'LineWidth', 3);
        % X line: vertical line at x=result.center_y (left-right position)
        plot(ax_xy, [result.center_y, result.center_y], [1, size(xy_proj, 1)], 'g--', 'LineWidth', 2);
        % Y line: horizontal line at y=result.center_x (up-down position)
        plot(ax_xy, [1, size(xy_proj, 2)], [result.center_x, result.center_x], 'g--', 'LineWidth', 2);
        fprintf('2D center marker: Using fitted center (%.1f, %.1f) as fallback\n', result.center_x, result.center_y);
    end
    
    hold(ax_xy, 'off');
    title(ax_xy, '2D Fit on Z-Projection', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_xy, 'X (pixels)', 'FontSize', 12);
    ylabel(ax_xy, 'Y (pixels)', 'FontSize', 12);
    colorbar(ax_xy);

    % Plot Z Profile
    ax_z = uiaxes(grid);
    ax_z.Layout.Column = 2;
    
    % Use background-corrected data for plotting
    % The background field contains the local background that was subtracted during fitting
    if isfield(result, 'background') && ~isnan(result.background)
        % Apply background correction to data BEFORE extracting Z profile
        data_bg_corrected = data - result.background;
        % For 2D+1D fitting, Z profile should be extracted at the ORIGINAL local maxima X,Y coordinates
        % NOT at the fitted X,Y center, because we want the precise Z position of the original maximum
        if isfield(result, 'originalMaximaCoords')
            % Convert global coordinates to local window coordinates
            local_x = result.originalMaximaCoords(1) - (size(data, 2) - 1) / 2;
            local_y = result.originalMaximaCoords(2) - (size(data, 1) - 1) / 2;
            local_x_round = round(local_x);
            local_y_round = round(local_y);
            
            if local_x_round >= 1 && local_x_round <= size(data_bg_corrected, 2) && ...
               local_y_round >= 1 && local_y_round <= size(data_bg_corrected, 1)
                z_profile = squeeze(data_bg_corrected(local_y_round, local_x_round, :));
                fprintf('Z profile (2D+1D): Extracted at ORIGINAL local maxima X,Y (%.1f, %.1f), applied per-pixel background correction of %.2f, profile range=[%.2f, %.2f]\n', ...
                        result.originalMaximaCoords(1), result.originalMaximaCoords(2), result.background, min(z_profile), max(z_profile));
            else
                % Fallback to summing if coordinates are out of bounds
                z_profile = squeeze(sum(sum(data_bg_corrected, 1), 2));
                fprintf('Z profile (2D+1D): Original coordinates out of bounds, using summed profile, applied per-pixel background correction of %.2f, profile range=[%.2f, %.2f]\n', ...
                        result.background, min(z_profile), max(z_profile));
            end
        else
            % Fallback to fitted center if original coordinates not available
            % CORRECTED: X is up-down (row), Y is left-right (column)
            center_x_round = round(result.center_x);
            center_y_round = round(result.center_y);
            if center_x_round >= 1 && center_x_round <= size(data_bg_corrected, 1) && ...
               center_y_round >= 1 && center_y_round <= size(data_bg_corrected, 2)
                % Use correct coordinate order: data(row, column, :) where row=X, column=Y
                z_profile = squeeze(data_bg_corrected(center_x_round, center_y_round, :));
                fprintf('Z profile (2D+1D): Original coordinates not available, extracted at fitted center (%.1f, %.1f), applied per-pixel background correction of %.2f, profile range=[%.2f, %.2f]\n', ...
                        result.center_x, result.center_y, result.background, min(z_profile), max(z_profile));
            else
                z_profile = squeeze(sum(sum(data_bg_corrected, 1), 2));
                fprintf('Z profile (2D+1D): Both original and fitted coordinates out of bounds, using summed profile, applied per-pixel background correction of %.2f, profile range=[%.2f, %.2f]\n', ...
                        result.background, min(z_profile), max(z_profile));
            end
        end
    else
        % Extract Z profile without background correction
        if isfield(result, 'originalMaximaCoords')
            local_x = result.originalMaximaCoords(1) - (size(data, 2) - 1) / 2;
            local_y = result.originalMaximaCoords(2) - (size(data, 1) - 1) / 2;
            local_x_round = round(local_x);
            local_y_round = round(local_y);
            
            if local_x_round >= 1 && local_x_round <= size(data, 2) && ...
               local_y_round >= 1 && local_y_round <= size(data, 1)
                z_profile = squeeze(data(local_y_round, local_x_round, :));
                fprintf('Z profile (2D+1D): Extracted at ORIGINAL local maxima X,Y (%.1f, %.1f), no background correction, profile range=[%.2f, %.2f]\n', ...
                        result.originalMaximaCoords(1), result.originalMaximaCoords(2), min(z_profile), max(z_profile));
            else
    z_profile = squeeze(sum(sum(data, 1), 2));
                fprintf('Z profile (2D+1D): Original coordinates out of bounds, using summed profile, no background correction, profile range=[%.2f, %.2f]\n', ...
                        min(z_profile), max(z_profile));
            end
        else
            % Fallback to fitted center if original coordinates not available
            % CORRECTED: X is up-down (row), Y is left-right (column)
            center_x_round = round(result.center_x);
            center_y_round = round(result.center_y);
            if center_x_round >= 1 && center_x_round <= size(data, 1) && ...
               center_y_round >= 1 && center_y_round <= size(data, 2)
                % Use correct coordinate order: data(row, column, :) where row=X, column=Y
                z_profile = squeeze(data(center_x_round, center_y_round, :));
                fprintf('Z profile (2D+1D): Original coordinates not available, extracted at fitted center (%.1f, %.1f), no background correction, profile range=[%.2f, %.2f]\n', ...
                        result.center_x, result.center_y, min(z_profile), max(z_profile));
            else
                z_profile = squeeze(sum(sum(data, 1), 2));
                fprintf('Z profile (2D+1D): Both original and fitted coordinates out of bounds, using summed profile, no background correction, profile range=[%.2f, %.2f]\n', ...
                        min(z_profile), max(z_profile));
            end
        end
    end
    
    z_coords = 1:length(z_profile);
    plot(ax_z, z_coords, z_profile, 'b.', 'MarkerSize', 10);
    hold(ax_z, 'on');
    z_fine = linspace(1, length(z_profile), 200);
    
    % For 2D+1D fitting, the Z profile amplitude should match the actual Z profile peak
    % The result.amplitude is from the 2D XY fit, but we need the Z profile amplitude
    if isfield(result, 'background') && ~isnan(result.background)
        % Use the actual peak of the background-corrected Z profile
        z_profile_amplitude = max(z_profile);
        fprintf('Z profile (2D+1D): Using actual Z profile peak amplitude: %.2e\n', z_profile_amplitude);
    else
        % Use the actual peak of the raw Z profile
        z_profile_amplitude = max(z_profile);
        fprintf('Z profile (2D+1D): Using actual Z profile peak amplitude: %.2e\n', z_profile_amplitude);
    end
    
    fit_z = z_profile_amplitude * exp(-(z_fine - result.center_z).^2 / (2 * result.sigma_z^2));
    fprintf('Z fit (2D+1D): Z profile amplitude=%.2e, center_z=%.2f, sigma_z=%.2f, max_fit=%.2e\n', ...
            z_profile_amplitude, result.center_z, result.sigma_z, max(fit_z));
    plot(ax_z, z_fine, fit_z, 'r-', 'LineWidth', 2);
    
    % Add center line for Z at the FITTED Gaussian center
    % The fitted center represents where the Gaussian fit determined the peak should be
    % Use exact fitted center value for sub-pixel precision (no rounding)
    center_z_exact = result.center_z;
    if center_z_exact >= 1 && center_z_exact <= length(z_profile)
        % Use full range to ensure line goes through the center point
        plot(ax_z, [center_z_exact, center_z_exact], [0, max(z_profile)*1.1], 'g--', 'LineWidth', 2);
        fprintf('Z center line (2D+1D): Placed at fitted Gaussian center Z coordinate %.3f (sub-pixel precision)\n', ...
                center_z_exact);
    else
        fprintf('Z center line (2D+1D): Fitted center Z coordinate %.3f is out of bounds [1, %d]\n', ...
                center_z_exact, length(z_profile));
    end
    
    hold(ax_z, 'off');
    title(ax_z, 'Z Profile Fit', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_z, 'Z (pixels)', 'FontSize', 12);
    ylabel(ax_z, 'Intensity', 'FontSize', 12);
    % Legend created after all plot elements (Data, Fit, Center) are added
    legend(ax_z, 'Data', 'Fit', 'Center', 'Location', 'best');
    ax_z.GridAlpha = 0.3;
    ax_z.GridLineStyle = '-';
    
    % Create info panel
    info_panel = uipanel(grid);
    info_panel.Layout.Column = 3;
    info_panel.Title = 'Fit Parameters';
    info_panel.FontSize = 12;
    info_panel.FontWeight = 'bold';
    
    % Set proper layout properties for the panel
    info_panel.Layout.Row = 1;
    info_panel.Layout.Column = 3;
    
    % Debug: Verify the info panel was created
    fprintf('Info panel created: %s\n', class(info_panel));
    fprintf('Info panel title: %s\n', info_panel.Title);
    fprintf('Info panel column: %d\n', info_panel.Layout.Column);
    
    % Add fit details with safe field access
    % First, let's see what fields are actually available
    available_fields = fieldnames(result);
    
    % Build info text safely
    info_parts = {};
    
    % Fit Method
    if isfield(result, 'fitMethod')
        info_parts{end+1} = sprintf('Fit Method: %s', result.fitMethod);
    else
        info_parts{end+1} = 'Fit Method: Unknown';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Amplitude
    if isfield(result, 'amplitude') && ~isnan(result.amplitude)
        info_parts{end+1} = sprintf('Amplitude: %.2f', result.amplitude);
    else
        info_parts{end+1} = 'Amplitude: N/A';
    end
    
    % Centers
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        info_parts{end+1} = sprintf('Center X: %.2f', result.center_x);
    else
        info_parts{end+1} = 'Center X: N/A';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        info_parts{end+1} = sprintf('Center Y: %.2f', result.center_y);
    else
        info_parts{end+1} = 'Center Y: N/A';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        info_parts{end+1} = sprintf('Center Z: %.2f', result.center_z);
    else
        info_parts{end+1} = 'Center Z: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Sigmas
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        info_parts{end+1} = sprintf('Sigma X: %.2f', result.sigma_x);
    else
        info_parts{end+1} = 'Sigma X: N/A';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        info_parts{end+1} = sprintf('Sigma Y: %.2f', result.sigma_y);
    else
        info_parts{end+1} = 'Sigma Y: N/A';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        info_parts{end+1} = sprintf('Sigma Z: %.2f', result.sigma_z);
    else
        info_parts{end+1} = 'Sigma Z: N/A';
    end
    
    % Rho values (for debugging)
    if isfield(result, 'rho_xy') && ~isnan(result.rho_xy)
        info_parts{end+1} = sprintf('Rho XY: %.3f', result.rho_xy);
    else
        info_parts{end+1} = 'Rho XY: N/A';
    end
    
    if isfield(result, 'rho_xz') && ~isnan(result.rho_xz)
        info_parts{end+1} = sprintf('Rho XZ: %.3f', result.rho_xz);
    else
        info_parts{end+1} = 'Rho XZ: N/A';
    end
    
    if isfield(result, 'rho_yz') && ~isnan(result.rho_yz)
        info_parts{end+1} = sprintf('Rho YZ: %.3f', result.rho_yz);
    else
        info_parts{end+1} = 'Rho YZ: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % R-squared
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared);
        elseif length(result.r_squared) >= 2
            info_parts{end+1} = sprintf('R² XY: %.3f', result.r_squared(1));
            info_parts{end+1} = sprintf('R² Z: %.3f', result.r_squared(2));
        else
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared(1));
        end
    else
        info_parts{end+1} = 'R²: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Background
    if isfield(result, 'background') && ~isnan(result.background)
        info_parts{end+1} = sprintf('Background: %.2f', result.background);
    else
        info_parts{end+1} = 'Background: N/A';
    end
    
    % Integrated Intensity
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        info_parts{end+1} = sprintf('Integrated Intensity: %.2f', result.integratedIntensity);
    else
        info_parts{end+1} = 'Integrated Intensity: N/A';
    end
    
    % Available fields (for debugging)
    info_parts{end+1} = '';
    info_parts{end+1} = 'Available fields:';
    for i = 1:min(10, length(available_fields))
        info_parts{end+1} = sprintf('  %s', available_fields{i});
    end
    if length(available_fields) > 10
        info_parts{end+1} = sprintf('  ... and %d more', length(available_fields) - 10);
    end
    
    % Create a table to display fit parameters
    % Prepare data for the table
    param_names = {};
    param_values = {};
    param_units = {};
    
    % Add fit parameters to the table
    if isfield(result, 'fitMethod')
        param_names{end+1} = 'Fit Method';
        param_values{end+1} = result.fitMethod;
        param_units{end+1} = '';
    end
    
    % Add equation information
    param_names{end+1} = '--- EQUATIONS ---';
    param_values{end+1} = '--- USED ---';
    param_units{end+1} = '---';
    
    param_names{end+1} = '1D Gaussian X';
    param_values{end+1} = 'f(x) = A_x*exp(-((x-μ_x)²)/(2σ_x²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Y';
    param_values{end+1} = 'f(y) = A_y*exp(-((y-μ_y)²)/(2σ_y²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Z';
    param_values{end+1} = 'f(z) = A_z*exp(-((z-μ_z)²)/(2σ_z²))';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Background';
    param_values{end+1} = 'Pre-subtracted (no offset)';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Constraints';
    param_values{end+1} = 'A>0, σ<window_size/2';
    param_units{end+1} = '';
    
    % Check if individual X, Y, Z amplitudes are available (for 1D fits)
    if isfield(result, 'amplitude_x') && isfield(result, 'amplitude_y') && isfield(result, 'amplitude_z')
        param_names{end+1} = 'Amplitude X';
        param_values{end+1} = sprintf('%.2f', result.amplitude_x);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Y';
        param_values{end+1} = sprintf('%.2f', result.amplitude_y);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Z';
        param_values{end+1} = sprintf('%.2f', result.amplitude_z);
        param_units{end+1} = 'counts';
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        param_names{end+1} = 'Amplitude';
        param_values{end+1} = sprintf('%.2f', result.amplitude);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        param_names{end+1} = 'Center X';
        param_values{end+1} = sprintf('%.2f', result.center_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        param_names{end+1} = 'Center Y';
        param_values{end+1} = sprintf('%.2f', result.center_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        param_names{end+1} = 'Center Z';
        param_values{end+1} = sprintf('%.2f', result.center_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        param_names{end+1} = 'Sigma X';
        param_values{end+1} = sprintf('%.2f', result.sigma_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        param_names{end+1} = 'Sigma Y';
        param_values{end+1} = sprintf('%.2f', result.sigma_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        param_names{end+1} = 'Sigma Z';
        param_values{end+1} = sprintf('%.2f', result.sigma_z);
        param_units{end+1} = 'pixels';
    end
    
    % Rho values (for debugging)
    if isfield(result, 'rho_xy') && ~isnan(result.rho_xy)
        param_names{end+1} = 'Rho XY';
        param_values{end+1} = sprintf('%.3f', result.rho_xy);
        param_units{end+1} = '';
    end
    
    if isfield(result, 'rho_xz') && ~isnan(result.rho_xz)
        param_names{end+1} = 'Rho XZ';
        param_values{end+1} = sprintf('%.3f', result.rho_xz);
        param_units{end+1} = '';
    end
    
    if isfield(result, 'rho_yz') && ~isnan(result.rho_yz)
        param_names{end+1} = 'Rho YZ';
        param_values{end+1} = sprintf('%.3f', result.rho_yz);
        param_units{end+1} = '';
    end
    
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared);
            param_units{end+1} = '';
        elseif length(result.r_squared) >= 2
            param_names{end+1} = 'R² XY';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
            param_names{end+1} = 'R² Z';
            param_values{end+1} = sprintf('%.3f', result.r_squared(2));
            param_units{end+1} = '';
        else
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
        end
    end
    
    if isfield(result, 'background') && ~isnan(result.background)
        param_names{end+1} = 'Background';
        param_values{end+1} = sprintf('%.2f', result.background);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        param_names{end+1} = 'Integrated Intensity';
        param_values{end+1} = sprintf('%.2f', result.integratedIntensity);
        param_units{end+1} = 'counts';
    end
    
    % Create the table first
    param_table = uitable(info_panel, ...
        'Data', [param_names', param_values', param_units'], ...
        'ColumnName', {'Parameter', 'Value', 'Unit'}, ...
        'ColumnWidth', {120, 100, 80}, ...
        'FontSize', 10, ...
        'RowName', []);
    
    % Debug: Check if table was created successfully
    fprintf('Table created: %s, empty: %d\n', class(param_table), isempty(param_table));
    
    % Try to set layout properties with error handling
    try
        param_table.Layout.Row = 1;
        param_table.Layout.Column = 1;
        fprintf('Layout properties set successfully\n');
    catch ME
    end
    
    fprintf('Fit parameters table created successfully\n');
end


function plot_3D_fits(parent, result)
    % Debug: Print result structure
    fprintf('=== 3D Plot Debug ===\n');
    fprintf('Result fields: %s\n', strjoin(fieldnames(result), ', '));
    if isfield(result, 'amplitude')
        fprintf('Amplitude: %.2e\n', result.amplitude);
    end
    if isfield(result, 'center_x')
        fprintf('Center X: %.2f\n', result.center_x);
    end
    if isfield(result, 'center_y')
        fprintf('Center Y: %.2f\n', result.center_y);
    end
    if isfield(result, 'center_z')
        fprintf('Center Z: %.2f\n', result.center_z);
    end
    if isfield(result, 'sigma_x')
        fprintf('Sigma X: %.2f\n', result.sigma_x);
    end
    if isfield(result, 'sigma_y')
        fprintf('Sigma Y: %.2f\n', result.sigma_y);
    end
    if isfield(result, 'sigma_z')
        fprintf('Sigma Z: %.2f\n', result.sigma_z);
    end
    fprintf('=======================\n');
    
    data = result.rawDataWindow;
    sz = size(data);
    
    % Ensure data is double type for slice function
    if ~isa(data, 'double')
        data = double(data);
        fprintf('3D plot: Converted data from %s to double\n', class(result.rawDataWindow));
    end
    
    % Create a grid layout with 1 row and 2 columns (3D plot + info panel)
    grid = uigridlayout(parent, [1, 2]);
    grid.ColumnWidth = {'1x', 300};
    grid.RowHeight = {'1x'};
    
    % Create the 3D plot
    ax = uiaxes(grid);
    ax.Layout.Column = 1;
    
    % --- Create Grid ---
    [Y, X, Z] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3)); % Note: X/Y swap for plotting
    
    % --- Plot Orthogonal Slices of Background-Corrected Data ---
    center_x = round(result.center_x);
    center_y = round(result.center_y);
    center_z = round(result.center_z);
    
    % Clamp centers to be within the window bounds
    center_x = max(1, min(sz(2), center_x));
    center_y = max(1, min(sz(1), center_y));
    center_z = max(1, min(sz(3), center_z));
    
    % Use background-corrected data for plotting
    % The background field contains the local background that was subtracted during fitting
    if isfield(result, 'background') && ~isnan(result.background)
        plot_data = data - result.background;
        fprintf('3D plot: Applied local background correction of %.2f\n', result.background);
    else
        plot_data = data;
        fprintf('3D plot: No background correction available\n');
    end
    
    % Ensure plot_data is double type for slice function
    if ~isa(plot_data, 'double')
        plot_data = double(plot_data);
        fprintf('3D plot: Converted plot_data to double\n');
    end
    
    % Validate plot_data before using slice function
    if any(~isfinite(plot_data(:)))
        warning('Non-finite values detected in plot_data. Attempting to clean data.');
        plot_data(~isfinite(plot_data)) = 0;
    end
    
    try
        slice(ax, Y, X, Z, plot_data, center_x, center_y, center_z);
        fprintf('3D slice plotted successfully\n');
    catch ME
        warning('3D slice plotting failed: %s. Attempting alternative visualization.', ME.message);
        % Fallback: just show the center point and lines without the slice
        fprintf('3D slice error details: %s\n', ME.message);
    end
    
    hold(ax, 'on');
    shading(ax, 'flat');
    colormap(ax, 'gray');
    
    % --- Generate and Plot FWHM Isosurface of the Fit ---
    p = [result.amplitude, result.center_x, result.center_y, result.center_z, ...
         result.sigma_x, result.sigma_y, result.sigma_z, ...
         result.rho_xy, result.rho_xz, result.rho_yz, ...
         result.alpha_x, result.alpha_y, result.alpha_z];
         
    fit_model = @(coords) snap_helpers.model_3d_skewed(p, coords);
    if strcmp(result.fitMethod, 'Distorted 3D Gaussian')
        fit_model = @(coords) snap_helpers.model_3d_distorted(p(1:10), coords);
    elseif strcmp(result.fitMethod, '3D Gaussian')
        % A simple model is needed if fitGaussians doesn't have it
        fit_model = @(coords) p(1) * exp(-((coords(:,1)-p(2)).^2/(2*p(5)^2) + (coords(:,2)-p(3)).^2/(2*p(6)^2) + (coords(:,3)-p(4)).^2/(2*p(7)^2)));
    end
    
    coords_vec = [X(:), Y(:), Z(:)];
    fit_values = fit_model(coords_vec);
    fit_grid = reshape(fit_values, sz);
    
    % Note: No background correction needed here - fit parameters are already based on background-corrected data
    fprintf('3D isosurface: Using fitted values without additional background correction\n');
    
    % Calculate FWHM level
    fwhm_level = result.amplitude / 2;
    fprintf('3D FWHM: Using level of %.2e\n', fwhm_level);
    
    % Validate fit_grid before creating isosurface
    if any(~isfinite(fit_grid(:)))
        warning('Non-finite values detected in fit_grid. Attempting to clean data.');
        fit_grid(~isfinite(fit_grid)) = 0;
    end
    
    try
    fv = isosurface(Y, X, Z, fit_grid, fwhm_level);
        if ~isempty(fv.vertices)
            patch_obj = patch(ax, fv);
            patch_obj.FaceColor = 'red';
            patch_obj.EdgeColor = 'none';
            patch_obj.FaceAlpha = 0.4;
            fprintf('3D isosurface plotted successfully\n');
        else
            fprintf('3D isosurface: No vertices found at FWHM level %.2e\n', fwhm_level);
        end
    catch ME
        warning('3D isosurface plotting failed: %s. Skipping isosurface.', ME.message);
        fprintf('3D isosurface error details: %s\n', ME.message);
    end
    
    % --- Plot Centerline/Point (large and clearly visible) ---
    % Ensure center lines pass through the actual Gaussian centroid
    % CORRECTED: X is up-down (row), Y is left-right (column)
    % X center line: vertical line at X position (up-down)
    plot3(ax, [result.center_y, result.center_y], [result.center_x, result.center_x], [1, sz(3)], 'g--', 'LineWidth', 3);
    % Y center line: horizontal line at Y position (left-right)  
    plot3(ax, [result.center_y, result.center_y], [1, sz(1)], [result.center_z, result.center_z], 'g--', 'LineWidth', 3);
    % Z center line: depth line at Z position
    plot3(ax, [1, sz(2)], [result.center_x, result.center_x], [result.center_z, result.center_z], 'g--', 'LineWidth', 3);
    % Center point: CORRECTED coordinate order
    plot3(ax, result.center_y, result.center_x, result.center_z, 'go', 'MarkerSize', 25, 'LineWidth', 4, 'MarkerFaceColor', 'g');
    
    hold(ax, 'off');
    
    % --- Finalize Plot Appearance ---
    title(ax, ['3D Fit (FWHM shell) - ' result.fitMethod], 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax, 'X (pixels)', 'FontSize', 12);
    ylabel(ax, 'Y (pixels)', 'FontSize', 12);
    zlabel(ax, 'Z (pixels)', 'FontSize', 12);
    axis(ax, 'equal', 'tight');
    view(ax, 3);
    ax.GridAlpha = 0.3;
    ax.GridLineStyle = '-';
    rotate3d(ax, 'on');
    
    % Correct axis direction to match image convention (Y increases downwards)
    set(ax, 'YDir', 'reverse');
    
    % Create info panel
    info_panel = uipanel(grid);
    info_panel.Layout.Column = 2;
    info_panel.Title = 'Fit Parameters';
    info_panel.FontSize = 12;
    info_panel.FontWeight = 'bold';
    
    % Set proper layout properties for the panel
    info_panel.Layout.Row = 1;
    info_panel.Layout.Column = 2;
    
    % Add fit details with safe field access
    info_parts = {};
    
    % Fit Method
    if isfield(result, 'fitMethod')
        info_parts{end+1} = sprintf('Fit Method: %s', result.fitMethod);
    else
        info_parts{end+1} = 'Fit Method: Unknown';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Amplitude
    if isfield(result, 'amplitude') && ~isnan(result.amplitude)
        info_parts{end+1} = sprintf('Amplitude: %.2f', result.amplitude);
    else
        info_parts{end+1} = 'Amplitude: N/A';
    end
    
    % Centers
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        info_parts{end+1} = sprintf('Center X: %.2f', result.center_x);
    else
        info_parts{end+1} = 'Center X: N/A';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        info_parts{end+1} = sprintf('Center Y: %.2f', result.center_y);
    else
        info_parts{end+1} = 'Center Y: N/A';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        info_parts{end+1} = sprintf('Center Z: %.2f', result.center_z);
    else
        info_parts{end+1} = 'Center Z: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Sigmas
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        info_parts{end+1} = sprintf('Sigma X: %.2f', result.sigma_x);
    else
        info_parts{end+1} = 'Sigma X: N/A';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        info_parts{end+1} = sprintf('Sigma Y: %.2f', result.sigma_y);
    else
        info_parts{end+1} = 'Sigma Y: N/A';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        info_parts{end+1} = sprintf('Sigma Z: %.2f', result.sigma_z);
    else
        info_parts{end+1} = 'Sigma Z: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % R-squared
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared);
        else
            info_parts{end+1} = sprintf('R²: %.3f', result.r_squared(1));
        end
    else
        info_parts{end+1} = 'R²: N/A';
    end
    
    info_parts{end+1} = ''; % Empty line
    
    % Background
    if isfield(result, 'background') && ~isnan(result.background)
        info_parts{end+1} = sprintf('Background: %.2f', result.background);
    else
        info_parts{end+1} = 'Background: N/A';
    end
    
    % Integrated Intensity
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        info_parts{end+1} = sprintf('Integrated Intensity: %.2f', result.integratedIntensity);
    else
        info_parts{end+1} = 'Integrated Intensity: N/A';
    end
    
    % Create a table to display fit parameters
    % Prepare data for the table
    param_names = {};
    param_values = {};
    param_units = {};
    
    % Add fit parameters to the table
    if isfield(result, 'fitMethod')
        param_names{end+1} = 'Fit Method';
        param_values{end+1} = result.fitMethod;
        param_units{end+1} = '';
    end
    
    % Add equation information
    param_names{end+1} = '--- EQUATIONS ---';
    param_values{end+1} = '--- USED ---';
    param_units{end+1} = '---';
    
    param_names{end+1} = '1D Gaussian X';
    param_values{end+1} = 'f(x) = A_x*exp(-((x-μ_x)²)/(2σ_x²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Y';
    param_values{end+1} = 'f(y) = A_y*exp(-((y-μ_y)²)/(2σ_y²))';
    param_units{end+1} = '';
    
    param_names{end+1} = '1D Gaussian Z';
    param_values{end+1} = 'f(z) = A_z*exp(-((z-μ_z)²)/(2σ_z²))';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Background';
    param_values{end+1} = 'Pre-subtracted (no offset)';
    param_units{end+1} = '';
    
    param_names{end+1} = 'Constraints';
    param_values{end+1} = 'A>0, σ<window_size/2';
    param_units{end+1} = '';
    
    % Check if individual X, Y, Z amplitudes are available (for 1D fits)
    if isfield(result, 'amplitude_x') && isfield(result, 'amplitude_y') && isfield(result, 'amplitude_z')
        param_names{end+1} = 'Amplitude X';
        param_values{end+1} = sprintf('%.2f', result.amplitude_x);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Y';
        param_values{end+1} = sprintf('%.2f', result.amplitude_y);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Z';
        param_values{end+1} = sprintf('%.2f', result.amplitude_z);
        param_units{end+1} = 'counts';
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        param_names{end+1} = 'Amplitude';
        param_values{end+1} = sprintf('%.2f', result.amplitude);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        param_names{end+1} = 'Center X';
        param_values{end+1} = sprintf('%.2f', result.center_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        param_names{end+1} = 'Center Y';
        param_values{end+1} = sprintf('%.2f', result.center_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        param_names{end+1} = 'Center Z';
        param_values{end+1} = sprintf('%.2f', result.center_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        param_names{end+1} = 'Sigma X';
        param_values{end+1} = sprintf('%.2f', result.sigma_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        param_names{end+1} = 'Sigma Y';
        param_values{end+1} = sprintf('%.2f', result.sigma_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        param_names{end+1} = 'Sigma Z';
        param_values{end+1} = sprintf('%.2f', result.sigma_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared);
            param_units{end+1} = '';
        else
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
        end
    end
    
    if isfield(result, 'background') && ~isnan(result.background)
        param_names{end+1} = 'Background';
        param_values{end+1} = sprintf('%.2f', result.background);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        param_names{end+1} = 'Integrated Intensity';
        param_values{end+1} = sprintf('%.2f', result.integratedIntensity);
        param_units{end+1} = 'counts';
    end
    
    % Create the table first
    param_table = uitable(info_panel, ...
        'Data', [param_names', param_values', param_units'], ...
        'ColumnName', {'Parameter', 'Value', 'Unit'}, ...
        'ColumnWidth', {120, 100, 80}, ...
        'FontSize', 10, ...
        'RowName', []);
    
    % Debug: Check if table was created successfully
    fprintf('Table created: %s, empty: %d\n', class(param_table), isempty(param_table));
    
    % Try to set layout properties with error handling
    try
        param_table.Layout.Row = 1;
        param_table.Layout.Column = 1;
        fprintf('Layout properties set successfully\n');
    catch ME
    end
    
    fprintf('Fit parameters table created successfully\n');
end


function intensity = model_2d_distorted_for_plot(p, xy)
    % This is a local copy of the model from fitGaussians for plotting
    % It's slightly different from the fitting version (rho vs theta)
    % For now, we assume rho_xy is close to theta for visualization.
    
    % Debug: Print input parameters
    fprintf('Model input: A=%.2e, mu=[%.2f,%.2f], sig=[%.2f,%.2f], rho=%.2f\n', ...
            p(1), p(2), p(3), p(4), p(5), p(6));
    
    A = p(1);
    mu = p(2:3);
    sig = p(4:5);
    
    % Handle rho_xy safely - clamp to valid range for asin
    rho_xy = p(6);
    if isnan(rho_xy) || isinf(rho_xy)
        fprintf('Warning: rho_xy=%.2f is invalid, using 0 (no rotation)\n', rho_xy);
        theta = 0; % No rotation for invalid rho_xy
    elseif abs(rho_xy) > 1
        fprintf('Warning: rho_xy=%.2f clamped to ±1\n', rho_xy);
        rho_xy = sign(rho_xy); % Clamp to ±1
        theta = asin(rho_xy);
    else
        theta = asin(rho_xy); % Now safe to call asin
    end
    fprintf('Calculated theta: %.2f radians\n', theta);
    
    x_centered = xy - mu;
    
    x_rot = x_centered(:,1)*cos(theta) + x_centered(:,2)*sin(theta);
    y_rot = -x_centered(:,1)*sin(theta) + x_centered(:,2)*cos(theta);
    
    % Ensure we don't get non-finite values
    intensity = A * exp(-((x_rot.^2 / (2*sig(1)^2)) + (y_rot.^2 / (2*sig(2)^2))));
    
    % Replace any non-finite values with 0
    intensity(~isfinite(intensity)) = 0;
    
    % If the rotated model produced all zeros, fall back to simple 2D Gaussian
    if max(intensity) == 0
        fprintf('Warning: Rotated model produced all zeros, using simple 2D Gaussian\n');
        x_centered = xy - mu;
        intensity = A * exp(-((x_centered(:,1).^2 / (2*sig(1)^2)) + (x_centered(:,2).^2 / (2*sig(2)^2))));
        intensity(~isfinite(intensity)) = 0;
    end
    
    % Debug: Print output statistics
    fprintf('Model output: min=%.2e, max=%.2e, mean=%.2e, any_nan=%d, any_inf=%d\n', ...
            min(intensity), max(intensity), mean(intensity), ...
            any(isnan(intensity)), any(isinf(intensity)));
end

function plot_radial_symmetry_fits(parent, result, channelColor)
    % Debug: Print result structure
    fprintf('=== Radial Symmetry Plot Debug ===\n');
    fprintf('Result fields: %s\n', strjoin(fieldnames(result), ', '));
    if isfield(result, 'amplitude')
        fprintf('Amplitude: %.2e\n', result.amplitude);
    end
    if isfield(result, 'center_x')
        fprintf('Center X: %.2f\n', result.center_x);
    end
    if isfield(result, 'center_y')
        fprintf('Center Y: %.2f\n', result.center_y);
    end
    if isfield(result, 'center_z')
        fprintf('Center Z: %.2f\n', result.center_z);
    end
    if isfield(result, 'sigma_x')
        fprintf('Sigma X: %.2f\n', result.sigma_x);
    end
    if isfield(result, 'sigma_y')
        fprintf('Sigma Y: %.2f\n', result.sigma_y);
    end
    if isfield(result, 'sigma_z')
        fprintf('Sigma Z: %.2f\n', result.sigma_z);
    end
    fprintf('=====================================\n');
    
    data = result.rawDataWindow;
    % Determine if this is 3D based on multiple criteria
    is_3d = false;
    
    % Method 1: Check if data has multiple Z slices
    if ndims(data) == 3 && size(data, 3) > 1
        is_3d = true;
    end
    
    % Method 2: Check if center_z exists and is not NaN (indicating 3D fit)
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        is_3d = true;
    end
    
    % Method 3: Check if fitMethod indicates 3D processing
    if isfield(result, 'fitMethod')
        method = result.fitMethod;
        if contains(method, '3D') || strcmp(method, 'Radial Symmetry') 
            % For radial symmetry, we need to check if it was processed in 3D mode
            % This can be inferred from having a valid center_z
            if strcmp(method, 'Radial Symmetry') && isfield(result, 'center_z') && ~isnan(result.center_z)
                is_3d = true;
            elseif contains(method, '3D')
                is_3d = true;
            end
        end
    end
    
    % 3D detection successful - removed debug output
    
    if is_3d
        % 3D Radial Symmetry - show 3D vector map visualization
        show_3d_vector_map(parent, result, data, channelColor);
    else
        % 2D Radial Symmetry - show 2D vector map visualization
        show_2d_vector_map(parent, result, data, channelColor);
    end
end

function plot_2d_radial_symmetry(parent, result, data, channelColor)
    % Create a grid layout with 1 row and 3 columns (2D image with vector overlay + radial profile + info panel)
    grid = uigridlayout(parent, [1, 3]);
    grid.ColumnWidth = {'1x', '1x', 300};
    grid.RowHeight = {'1x'};
    
    % Plot 2D Image with Center Marker
    ax_2d = uiaxes(grid);
    ax_2d.Layout.Column = 1;
    
    % Use background-corrected data for plotting (ensure double precision)
    if isfield(result, 'background') && ~isnan(result.background)
        plot_data = double(data - result.background);
        fprintf('2D radial: Applied local background correction of %.2f\n', result.background);
    else
        plot_data = double(data);
        fprintf('2D radial: No background correction available\n');
    end
    
    imagesc(ax_2d, plot_data);
    axis(ax_2d, 'image');
    colormap(ax_2d, 'gray');
    hold(ax_2d, 'on');
    
    % Add gradient vector overlay
    [rows, cols] = size(plot_data);
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Calculate gradients (ensure data is double precision)
    [Gx, Gy] = gradient(double(plot_data));
    
    % Subsample for cleaner visualization
    step = max(1, floor(min(rows, cols) / 20));
    X_sub = X(1:step:end, 1:step:end);
    Y_sub = Y(1:step:end, 1:step:end);
    Gx_sub = Gx(1:step:end, 1:step:end);
    Gy_sub = Gy(1:step:end, 1:step:end);
    
    % Normalize gradients for consistent arrow lengths
    magnitude = sqrt(Gx_sub.^2 + Gy_sub.^2);
    magnitude(magnitude == 0) = 1; % Avoid division by zero
    Gx_norm = Gx_sub ./ magnitude;
    Gy_norm = Gy_sub ./ magnitude;
    
    % Plot the vector field overlay
    quiver(ax_2d, X_sub, Y_sub, Gx_norm, Gy_norm, 0.5, 'Color', channelColor, 'LineWidth', 1);
    
    % Add center marker at the detected radial symmetry center using channel color
    if ~isnan(result.center_x) && ~isnan(result.center_y)
        plot(ax_2d, result.center_x, result.center_y, [channelColor 'o'], 'MarkerSize', 20, 'LineWidth', 4, 'MarkerFaceColor', channelColor);
        
        % Add center label with better positioning
        text(ax_2d, result.center_x + 3, result.center_y + 3, 'Center', 'Color', channelColor, 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', channelColor);
    end
    
    % Add center lines using channel color
    if ~isnan(result.center_x) && ~isnan(result.center_y)
        plot(ax_2d, [result.center_x, result.center_x], [1, size(plot_data, 1)], [channelColor '--'], 'LineWidth', 2);
        plot(ax_2d, [1, size(plot_data, 2)], [result.center_y, result.center_y], [channelColor '--'], 'LineWidth', 2);
    end
    
    hold(ax_2d, 'off');
    title(ax_2d, '2D Radial Symmetry with Vector Map', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_2d, 'X (pixels)', 'FontSize', 12);
    ylabel(ax_2d, 'Y (pixels)', 'FontSize', 12);
    colorbar(ax_2d);
    
    % Plot Radial Profile
    ax_radial = uiaxes(grid);
    ax_radial.Layout.Column = 2;
    
    % Calculate radial profile around the center
    [rows, cols] = size(plot_data);
    center_x = result.center_x;
    center_y = result.center_y;
    
    % Create coordinate grids
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Calculate distances from center
    distances = sqrt((X - center_x).^2 + (Y - center_y).^2);
    
    % Create radial profile
    max_radius = min([center_x-1, cols-center_x, center_y-1, rows-center_y]);
    if max_radius < 1
        max_radius = 1;
    end
    
    radial_profile = zeros(1, max_radius);
    for r = 1:max_radius
        % Find pixels at distance r
        mask = (distances >= r-0.5) & (distances < r+0.5);
        if any(mask(:))
            radial_profile(r) = mean(plot_data(mask));
        end
    end
    
    % Plot radial profile
    plot(ax_radial, 1:max_radius, radial_profile, 'b.-', 'MarkerSize', 8, 'LineWidth', 2);
    hold(ax_radial, 'on');
    
    % Overlay estimated Gaussian fit
    if ~isnan(result.sigma_x) && result.sigma_x > 0
        r_fine = linspace(1, max_radius, 200);
        center_val = plot_data(round(center_y), round(center_x));
        gaussian_fit = center_val * exp(-(r_fine.^2) / (2 * result.sigma_x^2));
        plot(ax_radial, r_fine, gaussian_fit, 'r-', 'LineWidth', 2);
        legend(ax_radial, 'Radial Profile', 'Gaussian Fit', 'Location', 'best');
    else
        legend(ax_radial, 'Radial Profile', 'Location', 'best');
    end
    
    hold(ax_radial, 'off');
    title(ax_radial, 'Radial Intensity Profile', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_radial, 'Distance from Center (pixels)', 'FontSize', 12);
    ylabel(ax_radial, 'Intensity', 'FontSize', 12);
    ax_radial.GridAlpha = 0.3;
    ax_radial.GridLineStyle = '-';
    
    % Create info panel
    create_radial_symmetry_info_panel(grid, result, 3);
end

function plot_3d_radial_symmetry(parent, result, data, channelColor)
    % Create a grid layout with 1 row and 3 columns (XY projection with vector overlay + Z profile + info panel)
    grid = uigridlayout(parent, [1, 3]);
    grid.ColumnWidth = {'1x', '1x', 300};
    grid.RowHeight = {'1x'};
    
    % Plot XY Projection with Center Marker
    ax_xy = uiaxes(grid);
    ax_xy.Layout.Column = 1;
    
    % Use background-corrected data for plotting (ensure double precision)
    if isfield(result, 'background') && ~isnan(result.background)
        plot_data = double(data - result.background);
        fprintf('3D radial: Applied local background correction of %.2f\n', result.background);
    else
        plot_data = double(data);
        fprintf('3D radial: No background correction available\n');
    end
    
    % Create XY projection (ensure data is double precision)
    xy_proj = double(max(plot_data, [], 3));
    
    imagesc(ax_xy, xy_proj);
    axis(ax_xy, 'image');
    colormap(ax_xy, 'gray');
    hold(ax_xy, 'on');
    
    % Add gradient vector overlay
    [rows, cols] = size(xy_proj);
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Calculate gradients on the XY projection (ensure data is double precision)
    [Gx, Gy] = gradient(double(xy_proj));
    
    % Subsample for cleaner visualization
    step = max(1, floor(min(rows, cols) / 20));
    X_sub = X(1:step:end, 1:step:end);
    Y_sub = Y(1:step:end, 1:step:end);
    Gx_sub = Gx(1:step:end, 1:step:end);
    Gy_sub = Gy(1:step:end, 1:step:end);
    
    % Normalize gradients for consistent arrow lengths
    magnitude = sqrt(Gx_sub.^2 + Gy_sub.^2);
    magnitude(magnitude == 0) = 1; % Avoid division by zero
    Gx_norm = Gx_sub ./ magnitude;
    Gy_norm = Gy_sub ./ magnitude;
    
    % Plot the vector field overlay
    quiver(ax_xy, X_sub, Y_sub, Gx_norm, Gy_norm, 0.5, 'Color', channelColor, 'LineWidth', 1);
    
    % Add center marker at the detected radial symmetry center using channel color
    if ~isnan(result.center_x) && ~isnan(result.center_y)
        plot(ax_xy, result.center_x, result.center_y, [channelColor 'o'], 'MarkerSize', 20, 'LineWidth', 4, 'MarkerFaceColor', channelColor);
        
        % Add center label with better positioning
        text(ax_xy, result.center_x + 3, result.center_y + 3, 'Center', 'Color', channelColor, 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', channelColor);
    end
    
    % Add center lines using channel color
    if ~isnan(result.center_x) && ~isnan(result.center_y)
        plot(ax_xy, [result.center_x, result.center_x], [1, size(xy_proj, 1)], [channelColor '--'], 'LineWidth', 2);
        plot(ax_xy, [1, size(xy_proj, 2)], [result.center_y, result.center_y], [channelColor '--'], 'LineWidth', 2);
    end
    
    hold(ax_xy, 'off');
    title(ax_xy, '3D Radial Symmetry with Vector Map (XY)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_xy, 'X (pixels)', 'FontSize', 12);
    ylabel(ax_xy, 'Y (pixels)', 'FontSize', 12);
    colorbar(ax_xy);
    
    % Plot Z Profile
    ax_z = uiaxes(grid);
    ax_z.Layout.Column = 2;
    
    % Extract Z profile at the center
    center_x = round(result.center_x);
    center_y = round(result.center_y);
    center_z = round(result.center_z);
    
    if center_x >= 1 && center_x <= size(plot_data, 2) && ...
       center_y >= 1 && center_y <= size(plot_data, 1) && ...
       center_z >= 1 && center_z <= size(plot_data, 3)
        z_profile = squeeze(plot_data(center_y, center_x, :));
    else
        % Fallback to center slice
        z_profile = squeeze(plot_data(round(size(plot_data,1)/2), round(size(plot_data,2)/2), :));
    end
    
    % Plot Z profile
    plot(ax_z, 1:length(z_profile), z_profile, 'b.-', 'MarkerSize', 8, 'LineWidth', 2);
    hold(ax_z, 'on');
    
    % Overlay estimated Gaussian fit
    if ~isnan(result.sigma_z) && result.sigma_z > 0
        z_fine = linspace(1, length(z_profile), 200);
        center_val = z_profile(center_z);
        gaussian_fit = center_val * exp(-((z_fine - center_z).^2) / (2 * result.sigma_z^2));
        plot(ax_z, z_fine, gaussian_fit, 'r-', 'LineWidth', 2);
        legend(ax_z, 'Z Profile', 'Gaussian Fit', 'Location', 'best');
    else
        legend(ax_z, 'Z Profile', 'Location', 'best');
    end
    
    % Add center line using channel color
    plot(ax_z, [center_z, center_z], [0, max(z_profile)*1.1], [channelColor '--'], 'LineWidth', 2);
    
    hold(ax_z, 'off');
    title(ax_z, 'Z Profile at Center', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax_z, 'Z (pixels)', 'FontSize', 12);
    ylabel(ax_z, 'Intensity', 'FontSize', 12);
    ax_z.GridAlpha = 0.3;
    ax_z.GridLineStyle = '-';
    
    % Create info panel
    create_radial_symmetry_info_panel(grid, result, 3);
end

function create_radial_symmetry_info_panel(parent, result, column)
    % Create info panel for radial symmetry results
    info_panel = uipanel(parent);
    info_panel.Layout.Column = column;
    info_panel.Title = 'Radial Symmetry Parameters';
    info_panel.FontSize = 12;
    info_panel.FontWeight = 'bold';
    
    % Set proper layout properties for the panel
    info_panel.Layout.Row = 1;
    info_panel.Layout.Column = column;
    
    % Prepare data for the table
    param_names = {};
    param_values = {};
    param_units = {};
    
    % Add fit parameters to the table
    param_names{end+1} = 'Fit Method';
    param_values{end+1} = 'Radial Symmetry';
    param_units{end+1} = '';
    
    % Check if individual X, Y, Z amplitudes are available (for 1D fits)
    if isfield(result, 'amplitude_x') && isfield(result, 'amplitude_y') && isfield(result, 'amplitude_z')
        param_names{end+1} = 'Amplitude X';
        param_values{end+1} = sprintf('%.2f', result.amplitude_x);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Y';
        param_values{end+1} = sprintf('%.2f', result.amplitude_y);
        param_units{end+1} = 'counts';
        
        param_names{end+1} = 'Amplitude Z';
        param_values{end+1} = sprintf('%.2f', result.amplitude_z);
        param_units{end+1} = 'counts';
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        param_names{end+1} = 'Amplitude';
        param_values{end+1} = sprintf('%.2f', result.amplitude);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        param_names{end+1} = 'Center X';
        param_values{end+1} = sprintf('%.2f', result.center_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        param_names{end+1} = 'Center Y';
        param_values{end+1} = sprintf('%.2f', result.center_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        param_names{end+1} = 'Center Z';
        param_values{end+1} = sprintf('%.2f', result.center_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        param_names{end+1} = 'Sigma X';
        param_values{end+1} = sprintf('%.2f', result.sigma_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        param_names{end+1} = 'Sigma Y';
        param_values{end+1} = sprintf('%.2f', result.sigma_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        param_names{end+1} = 'Sigma Z';
        param_values{end+1} = sprintf('%.2f', result.sigma_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'r_squared')
        if isscalar(result.r_squared)
            param_names{end+1} = 'R² (Symmetry Score)';
            param_values{end+1} = sprintf('%.3f', result.r_squared);
            param_units{end+1} = '';
        else
            param_names{end+1} = 'R² (Symmetry Score)';
            param_values{end+1} = sprintf('%.3f', result.r_squared(1));
            param_units{end+1} = '';
        end
    end
    
    if isfield(result, 'background') && ~isnan(result.background)
        param_names{end+1} = 'Background';
        param_values{end+1} = sprintf('%.2f', result.background);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'integratedIntensity') && ~isnan(result.integratedIntensity)
        param_names{end+1} = 'Integrated Intensity';
        param_values{end+1} = sprintf('%.2f', result.integratedIntensity);
        param_units{end+1} = 'counts';
    end
    
    % Create the table
    param_table = uitable(info_panel, ...
        'Data', [param_names', param_values', param_units'], ...
        'ColumnName', {'Parameter', 'Value', 'Unit'}, ...
        'ColumnWidth', {140, 100, 80}, ...
        'FontSize', 10, ...
        'RowName', []);
    
    % Set layout properties
    try
        param_table.Layout.Row = 1;
        param_table.Layout.Column = 1;
    catch ME
        fprintf('Warning: Could not set layout properties: %s\n', ME.message);
    end
    
    fprintf('Radial symmetry parameters table created successfully\n');
end

function show_3d_vector_map(parent, result, data, channelColor)
    % Show 3D vector map visualization in separate figure
    
    % Create a simple fitParams structure with the necessary fields
    fitParams = struct();
    fitParams.gaussFitRadialRadius = 5; % Default radius for visualization
    
    % Get center coordinates in the format expected by the visualization
    center = [result.center_y, result.center_x, result.center_z]; % [y, x, z] format
    
    % Get spacing information if available, otherwise use defaults
    xy_spacing = 1; % Default
    z_spacing = 1;  % Default
    
    try
        % Call the existing 3D vector visualization function (creates its own figure)
        snap_helpers.visualizeRadialSymmetry3D(data, center, fitParams, ...
            'Title', sprintf('3D Radial Symmetry Vector Map - Score: %.3f', result.r_squared), ...
            'XYSpacing', xy_spacing, 'ZSpacing', z_spacing, ...
            'VectorDensity', 2, 'ShowScoreMap', true);
        
        fprintf('3D vector map visualization displayed successfully\n');
    catch ME
        warning('3D vector map visualization failed: %s', ME.message);
        
        % Fallback to the original 3D radial symmetry plot
        plot_3d_radial_symmetry(parent, result, data, channelColor);
    end
end

function show_2d_vector_map(parent, result, data, channelColor)
    % Show 2D vector map visualization in separate figure
    
    % Create a simple fitParams structure with the necessary fields
    fitParams = struct();
    fitParams.gaussFitRadialRadius = 5; % Default radius for visualization
    
    % Get center coordinates in the format expected by the visualization
    center = [result.center_y, result.center_x]; % [y, x] format
    
    % Get spacing information if available, otherwise use defaults
    xy_spacing = 1; % Default
    
    try
        % Call the existing 2D vector visualization function (creates its own figure)
        snap_helpers.visualizeRadialSymmetry2D(data, center, fitParams, ...
            'Title', sprintf('2D Radial Symmetry Vector Map - Score: %.3f', result.r_squared), ...
            'XYSpacing', xy_spacing, ...
            'VectorDensity', 3, 'ShowScoreMap', true);
        
        fprintf('2D vector map visualization displayed successfully\n');
    catch ME
        warning('2D vector map visualization failed: %s', ME.message);
        
        % Fallback to the original 2D radial symmetry plot
        plot_2d_radial_symmetry(parent, result, data, channelColor);
    end
end




