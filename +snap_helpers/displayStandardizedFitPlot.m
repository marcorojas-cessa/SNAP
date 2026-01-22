function displayStandardizedFitPlot(handles, channel_idx, maxima_idx)
% Creates a standardized three-part fitting visualization window
% Part 1: Fit graphs (left section)
% Part 2: Background information (middle section) 
% Part 3: Numerical information (right section)
%
% COORDINATE CONVENTION: 
%   Data stored in ARRAY CONVENTION [row, col, slice]
%   Display uses CARTESIAN CONVENTION [x, y, z] = [col, row, slice]
%   
%   For visualization:
%   - X-axis (horizontal) = column dimension
%   - Y-axis (vertical) = row dimension  
%   - Z-axis (depth) = slice dimension
%
%   NOTE: Fit results (center_x, center_y, center_z) are LOCAL coordinates
%   within the fitting window, stored in ARRAY CONVENTION

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
    return;
end

% Check if this is a valid fit result
is_valid_fit = false;
if isfield(fit_result, 'fitMethod') && strcmp(fit_result.fitMethod, 'Radial Symmetry')
    if isfield(fit_result, 'center_x') && isfield(fit_result, 'center_y')
        is_valid_fit = true;
    end
else
    if ~isnan(fit_result.amplitude)
        is_valid_fit = true;
    end
end

if ~is_valid_fit
    disp('No valid fit result found for this maximum.');
    return;
end

% Get channel color for plotting
colorName = handles.maximaColorDrops(channel_idx).Value;
colorMap = containers.Map({'Red','Green','Blue','Yellow','Magenta','Cyan'}, {'r','g','b','y','m','c'});
channelColor = colorMap(colorName);

% Create standardized figure with three sections
fig = uifigure('Name', ['Fit Analysis: Ch ' num2str(channel_idx) ', Maxima ' num2str(maxima_idx)], ...
               'Position', [100 100 1400 800]);

% Main grid: 3 columns for the three sections
mainGrid = uigridlayout(fig, [1 3]);
mainGrid.ColumnWidth = {'1x', '1x', 350}; % Fit graphs, Background, Numerical info
mainGrid.Padding = [10 10 10 10];

% --- SECTION 1: FIT GRAPHS (Left) ---
fitSection = uipanel(mainGrid, 'Title', 'Fit Analysis', 'FontWeight', 'bold', 'FontSize', 14);
fitSection.Layout.Column = 1;
createFitGraphs(fitSection, fit_result, channelColor);

% --- SECTION 2: BACKGROUND INFORMATION (Middle) ---
bgSection = uipanel(mainGrid, 'Title', 'Background Analysis', 'FontWeight', 'bold', 'FontSize', 14);
bgSection.Layout.Column = 2;
createBackgroundGraphs(bgSection, fit_result, channelColor);

% --- SECTION 3: NUMERICAL INFORMATION (Right) ---
infoSection = uipanel(mainGrid, 'Title', 'Fit Parameters & Statistics', 'FontWeight', 'bold', 'FontSize', 14);
infoSection.Layout.Column = 3;
createNumericalInfo(infoSection, fit_result, channel_idx, maxima_idx);

end

function createFitGraphs(parent, result, channelColor)
% Creates the fit visualization graphs based on the fitting method
    
    switch result.fitMethod
        case {'1D (X,Y,Z)', '1D (X,Y,Z) Gaussian'}
            create1DFitGraphs(parent, result, channelColor);
        case {'2D (XY) + 1D (Z)', '2D (XY) + 1D (Z) Gaussian'}
            create2DPlus1DFitGraphs(parent, result, channelColor);
        case {'3D Gaussian', 'Distorted 3D Gaussian'}
            create3DFitGraphs(parent, result, channelColor);
        case 'Radial Symmetry'
            createRadialSymmetryGraphs(parent, result, channelColor);
        otherwise
            % Fallback for unknown methods
            grid = uigridlayout(parent, [1 1]);
            uilabel(grid, 'Text', ['Fit visualization for ' result.fitMethod ' is not implemented.'], ...
                   'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end

function create1DFitGraphs(parent, result, channelColor)
% Creates 1D fit graphs (X, Y, Z profiles)
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col, result.center_z = slice (ARRAY CONVENTION)
%   Profiles show intensity along each dimension
%   NO coordinate conversion needed for 1D profiles (just plotting intensity vs position)
    
    data = result.rawDataWindow;
    
    % Create grid for three 1D plots
    grid = uigridlayout(parent, [2 2]);
    grid.RowHeight = {'1x', '1x'};
    grid.ColumnWidth = {'1x', '1x'};
    
    % Apply background correction if available and ensure double precision
    if isfield(result, 'background') && ~isnan(result.background)
        data_corrected = double(data) - double(result.background);
    else
        data_corrected = double(data);
    end
    
    % X Profile
    ax_x = uiaxes(grid);
    ax_x.Layout.Row = 1;
    ax_x.Layout.Column = 1;
    plotProfile(ax_x, data_corrected, result, 'X', channelColor);
    
    % Y Profile  
    ax_y = uiaxes(grid);
    ax_y.Layout.Row = 1;
    ax_y.Layout.Column = 2;
    plotProfile(ax_y, data_corrected, result, 'Y', channelColor);
    
    % Z Profile (if 3D data)
    if ndims(data) == 3 && size(data, 3) > 1
        ax_z = uiaxes(grid);
        ax_z.Layout.Row = 2;
        ax_z.Layout.Column = [1 2];
        plotProfile(ax_z, data_corrected, result, 'Z', channelColor);
    else
        % Show 2D summary
        ax_summary = uiaxes(grid);
        ax_summary.Layout.Row = 2;
        ax_summary.Layout.Column = [1 2];
        imagesc(ax_summary, data_corrected);
        axis(ax_summary, 'image');
        colormap(ax_summary, 'gray');
        title(ax_summary, '2D Data Window', 'FontWeight', 'bold');
        colorbar(ax_summary);
    end
    
    % Print visualization summary to Command Window
    fprintf('\n=== 1D FIT VISUALIZATION SUMMARY ===\n');
    fprintf('Available fields in result:\n');
    disp(fieldnames(result));
    fprintf('\nFitted amplitudes:\n');
    if isfield(result, 'amplitude_x') && ~isnan(result.amplitude_x)
        fprintf('  X: %.2e\n', result.amplitude_x);
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        fprintf('  X: %.2e\n', result.amplitude);
    end
    if isfield(result, 'amplitude_y') && ~isnan(result.amplitude_y)
        fprintf('  Y: %.2e\n', result.amplitude_y);
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        fprintf('  Y: %.2e\n', result.amplitude);
    end
    if isfield(result, 'amplitude_z') && ~isnan(result.amplitude_z)
        fprintf('  Z: %.2e\n', result.amplitude_z);
    elseif isfield(result, 'amplitude') && ~isnan(result.amplitude)
        fprintf('  Z: %.2e\n', result.amplitude);
    end
    fprintf('\nR² values:\n');
    if isfield(result, 'r_squared_x')
        fprintf('  X: %.4f\n', result.r_squared_x);
    elseif isfield(result, 'r_squared') && length(result.r_squared) >= 1
        fprintf('  X: %.4f\n', result.r_squared(1));
    end
    if isfield(result, 'r_squared_y')
        fprintf('  Y: %.4f\n', result.r_squared_y);
    elseif isfield(result, 'r_squared') && length(result.r_squared) >= 2
        fprintf('  Y: %.4f\n', result.r_squared(2));
    end
    if isfield(result, 'r_squared_z')
        fprintf('  Z: %.4f\n', result.r_squared_z);
    elseif isfield(result, 'r_squared') && length(result.r_squared) >= 3
        fprintf('  Z: %.4f\n', result.r_squared(3));
    end
    fprintf('====================================\n\n');
end

function create2DPlus1DFitGraphs(parent, result, channelColor)
% Creates 2D+1D fit graphs (XY contour + Z profile)
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col, result.center_z = slice (ARRAY CONVENTION)
%   For 2D plot: Convert to CARTESIAN (X=col, Y=row)
%   For Z profile: No conversion needed
    
    data = result.rawDataWindow;
    
    % Create grid with 2D plot on top (larger) and Z profile below (smaller)
    grid = uigridlayout(parent, [2 1]);
    grid.RowHeight = {'2x', '1x'}; % 2D plot gets 2/3, Z profile gets 1/3
    
    % Apply background correction if available and ensure double precision
    if isfield(result, 'background') && ~isnan(result.background)
        data_corrected = double(data) - double(result.background);
    else
        data_corrected = double(data);
    end
    
    % 2D XY Plot with contours (larger, on top)
    ax_xy = uiaxes(grid);
    ax_xy.Layout.Row = 1;
    
    % Use the stored XY slice that was actually fitted, if available
    if isfield(result, 'xy_slice') && ~isempty(result.xy_slice)
        xy_proj = result.xy_slice;
    else
        % Fallback to center slice
        center_z_idx = round(result.center_z);
        center_z_idx = max(1, min(size(data_corrected, 3), center_z_idx));
        xy_proj = data_corrected(:, :, center_z_idx);
    end
    imagesc(ax_xy, xy_proj);
    axis(ax_xy, 'image');
    colormap(ax_xy, 'gray');
    hold(ax_xy, 'on');
    
    % Add fit contours with validation
    try
        [Y, X] = meshgrid(1:size(xy_proj, 2), 1:size(xy_proj, 1));
        coords = [X(:), Y(:)];
        
        % Validate fit parameters before creating model
        if isfield(result, 'rho_xy') && ~isnan(result.rho_xy)
            p = [result.amplitude, result.center_x, result.center_y, result.sigma_x, result.sigma_y, result.rho_xy];
        else
            p = [result.amplitude, result.center_x, result.center_y, result.sigma_x, result.sigma_y, 0];
        end
        
        % Check if parameters are valid
        if any(isnan(p)) || any(~isfinite(p)) || result.amplitude <= 0 || result.sigma_x <= 0 || result.sigma_y <= 0
            fprintf('Skipping contour plot: Invalid fit parameters\n');
        else
            fit_vals = model_2d_gaussian(p, coords);
            fit_grid = reshape(fit_vals, size(xy_proj));
            
            % Check if fit_grid has variation (not constant)
            if max(fit_grid(:)) > min(fit_grid(:)) && max(fit_grid(:)) > 0
                contour_levels = [result.amplitude/2, result.amplitude*0.8];
                % Filter out levels that are too small or invalid
                contour_levels = contour_levels(contour_levels > max(fit_grid(:))*0.01);
                
                if ~isempty(contour_levels)
                    contour(ax_xy, Y, X, fit_grid, contour_levels, 'r-', 'LineWidth', 2);
                else
                    fprintf('Skipping contour plot: No valid contour levels\n');
                end
            else
                fprintf('Skipping contour plot: Fit surface has no variation (constant values)\n');
            end
        end
    catch ME
        fprintf('Contour plotting failed: %s\n', ME.message);
    end
    
    % Add center marker (CARTESIAN CONVERSION for display)
    % result.center_x = row, result.center_y = col (ARRAY CONVENTION)
    % For plot: X=col, Y=row (CARTESIAN)
    plot(ax_xy, result.center_y, result.center_x, 'r+', 'MarkerSize', 15, 'LineWidth', 3);
    
    hold(ax_xy, 'off');
    title(ax_xy, '2D Fit (XY Projection)', 'FontWeight', 'bold');
    xlabel(ax_xy, 'X (pixels)');
    ylabel(ax_xy, 'Y (pixels)');
    colorbar(ax_xy);
    
    % Z Profile (smaller, below)
    ax_z = uiaxes(grid);
    ax_z.Layout.Row = 2;
    plotProfile(ax_z, data_corrected, result, 'Z', channelColor);
    
    % Print visualization summary to Command Window
    fprintf('\n=== 2D+1D FIT VISUALIZATION SUMMARY ===\n');
    fprintf('Amplitudes:\n');
    if isfield(result, 'amplitude_x') && ~isnan(result.amplitude_x)
        fprintf('  X: %.2e\n', result.amplitude_x);
    end
    if isfield(result, 'amplitude_y') && ~isnan(result.amplitude_y)
        fprintf('  Y: %.2e\n', result.amplitude_y);
    end
    if isfield(result, 'amplitude_z') && ~isnan(result.amplitude_z)
        fprintf('  Z: %.2e\n', result.amplitude_z);
    end
    fprintf('\nR² Values:\n');
    if isfield(result, 'r_squared')
        dimension_labels = {'X', 'Y', 'Z'};
        for i = 1:length(result.r_squared)
            if i <= length(dimension_labels)
                fprintf('  %s: %.4f\n', dimension_labels{i}, result.r_squared(i));
            end
        end
    end
    fprintf('=======================================\n\n');
end

function create3DFitGraphs(parent, result, channelColor)
% Creates 3D fit visualization with isosurface shell
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col, result.center_z = slice (ARRAY CONVENTION)
%   For plot3: Convert to CARTESIAN (X=col, Y=row, Z=slice)
%   Slice displays data with proper pixel-centered coordinates
    
    % Use the actual fitted data if available, otherwise fall back to raw data with background correction
    if isfield(result, 'fittedDataWindow')
        plot_data = double(result.fittedDataWindow);
    elseif isfield(result, 'background') && ~isnan(result.background)
        plot_data = double(result.rawDataWindow) - double(result.background);
    else
        plot_data = double(result.rawDataWindow);
    end
    
    sz = size(plot_data);  % sz = [rows, cols, slices] = [Y, X, Z]
    
    % Create single 3D plot
    grid = uigridlayout(parent, [1 1]);
    ax = uiaxes(grid);
    
    try
        % The slice planes should be at the CENTER of the fitting window
        % The window is extracted symmetrically around the maxima, so the center
        % is simply the middle of the window dimensions
        center_x = ceil(sz(2) / 2);  % Middle column
        center_y = ceil(sz(1) / 2);  % Middle row
        center_z = ceil(sz(3) / 2);  % Middle slice
        
        % Create coordinate grids with pixel centers at integers
        % For 7×7×7 data to display with centers at [1,2,3,4,5,6,7],
        % we need grid points at pixel edges [0.5, 1.5, 2.5, ..., 7.5] (8 points)
        x_edges = 0.5:1:(sz(2)+0.5);
        y_edges = 0.5:1:(sz(1)+0.5);
        z_edges = 0.5:1:(sz(3)+0.5);
        
        [X_grid, Y_grid, Z_grid] = meshgrid(x_edges, y_edges, z_edges);
        
        % Pad data by 1 on the 'post' (end) side only to get 8×8×8 from 7×7×7
        % Each data value represents the pixel value at that grid edge
        plot_data_extended = padarray(plot_data, [1 1 1], 'replicate', 'post');
        
        % Display orthogonal slice planes through the data at the fit center
        h = slice(ax, X_grid, Y_grid, Z_grid, plot_data_extended, center_x, center_y, center_z);
        
        % Style the slices
        for i = 1:length(h)
            set(h(i), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
        
        % Set colormap to grayscale (black and white)
        colormap(ax, 'gray');
        
        % Generate and display the fitted Gaussian isosurface
        % Use swapped X<->Y coordinates (no Z flip)
        fitted_gaussian = generate3DGaussian_version(result, sz, 'swapped');
        
        % For isosurface, use pixel-centered coordinates (not extended)
        [X_iso, Y_iso, Z_iso] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
        
        % Create isosurface at half-maximum level
        iso_level = max(fitted_gaussian(:)) * 0.5;
        
        % Plot the fitted Gaussian as a semi-transparent isosurface
        if iso_level > 0 && iso_level < max(fitted_gaussian(:))
            fv = isosurface(X_iso, Y_iso, Z_iso, fitted_gaussian, iso_level);
            
            if ~isempty(fv.vertices) && ~isempty(fv.faces)
                hold(ax, 'on');
                patch(ax, 'Vertices', fv.vertices, 'Faces', fv.faces, ...
                      'FaceColor', 'cyan', 'EdgeColor', 'none', ...
                      'FaceAlpha', 0.5);
                hold(ax, 'off');
            end
        end
        
        % Add center point marker (CARTESIAN CONVERSION for display)
        % result: center_x=row, center_y=col, center_z=slice (ARRAY CONVENTION)
        % For plot3: X=col, Y=row, Z=slice (CARTESIAN)
        hold(ax, 'on');
        plot3(ax, result.center_y, result.center_x, result.center_z, 'go', ...
              'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'g');
        hold(ax, 'off');
        
        % Set axis labels
        xlabel(ax, 'X (columns →)');
        ylabel(ax, 'Y (rows →)');
        zlabel(ax, 'Z (slices →)');
        
        % Set axis limits so pixel edges are at half-integers (centers at integers)
        xlim(ax, [0.5, sz(2) + 0.5]);
        ylim(ax, [0.5, sz(1) + 0.5]);
        zlim(ax, [0.5, sz(3) + 0.5]);
        
        % Ticks are already at pixel centers (integers 1, 2, 3, ...)
        ax.XTick = 1:sz(2);
        ax.YTick = 1:sz(1);
        ax.ZTick = 1:sz(3);
        
        % Set title and view
        title(ax, ['3D Fit - ' result.fitMethod], 'FontWeight', 'bold');
        view(ax, 3);
        
        % Enable grid using property syntax (compatible with uiaxes)
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.ZGrid = 'on';
        
        % Add colorbar
        colorbar(ax);
        
    catch ME
        % Fallback to simple text if 3D plotting fails
        fprintf('ERROR: 3D visualization failed: %s\n', ME.message);
        fprintf('ERROR: Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
        cla(ax);
        text(ax, 0.5, 0.5, sprintf('3D visualization failed: %s', ME.message), ...
             'HorizontalAlignment', 'center');
    end
    
    % Print visualization summary to Command Window
    fprintf('\n=== 3D FIT VISUALIZATION SUMMARY ===\n');
    fprintf('Fit Method: %s\n', result.fitMethod);
    fprintf('Amplitude: %.2e\n', result.amplitude);
    fprintf('Center (X,Y,Z): (%.2f, %.2f, %.2f)\n', result.center_x, result.center_y, result.center_z);
    fprintf('Sigma (X,Y,Z): (%.2f, %.2f, %.2f)\n', result.sigma_x, result.sigma_y, result.sigma_z);
    fprintf('R² value: ');
    if isfield(result, 'r_squared')
        if isnumeric(result.r_squared) && ~isnan(result.r_squared)
            fprintf('%.4f\n', result.r_squared);
        end
    end
    fprintf('====================================\n\n');
end

function fitted_gaussian = generate3DGaussian_version(result, sz, version)
% Generate a 3D Gaussian based on fit parameters
% sz = [rows, cols, slices] = [Y_size, X_size, Z_size]
% 
% COORDINATE TRANSFORMATION FOR VISUALIZATION:
% The fitting coordinate system doesn't match the visual display, so we need to:
% 1. Swap X and Y coordinates
% 2. Flip Z coordinates (Z_visual = sz(3) + 1 - Z_fit)
    
    % Create coordinate grids matching image indexing
    [X, Y, Z] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
    
    % Extract fit parameters and apply coordinate transformations based on version
    A = result.amplitude;
    
    switch version
        case 'original'
            % No transformation
            mu_x = result.center_x;
            mu_y = result.center_y;
            mu_z = result.center_z;
            sigma_x = result.sigma_x;
            sigma_y = result.sigma_y;
            sigma_z = result.sigma_z;
            
        case 'swapped'
            % Swap X and Y only
            mu_x = result.center_y;
            mu_y = result.center_x;
            mu_z = result.center_z;
            sigma_x = result.sigma_y;
            sigma_y = result.sigma_x;
            sigma_z = result.sigma_z;
            
        case 'swapped_flipped'
            % Swap X and Y, and flip Z
            mu_x = result.center_y;
            mu_y = result.center_x;
            mu_z = sz(3) + 1 - result.center_z;
            sigma_x = result.sigma_y;
            sigma_y = result.sigma_x;
            sigma_z = result.sigma_z;
    end
    
    % Check for distortion/rotation parameters
    if isfield(result, 'rho_xy') && ~isnan(result.rho_xy)
        % Distorted Gaussian with correlation
        rho_xy = result.rho_xy;
        rho_xz = result.rho_xz;
        rho_yz = result.rho_yz;
        
        % Construct covariance matrix
        Sigma = [sigma_x^2,              rho_xy*sigma_x*sigma_y,   rho_xz*sigma_x*sigma_z;
                 rho_xy*sigma_x*sigma_y,   sigma_y^2,               rho_yz*sigma_y*sigma_z;
                 rho_xz*sigma_x*sigma_z,   rho_yz*sigma_y*sigma_z,   sigma_z^2];
        
        % Check if covariance matrix is valid
        [~, flag] = chol(Sigma);
        if flag == 0
            % Valid covariance matrix - use multivariate Gaussian
            Sigma_inv = inv(Sigma);
            
            % Vectorize computation
            X_centered = X(:) - mu_x;
            Y_centered = Y(:) - mu_y;
            Z_centered = Z(:) - mu_z;
            coords = [X_centered, Y_centered, Z_centered];
            
            % Compute Gaussian
            exponent = -0.5 * sum((coords * Sigma_inv) .* coords, 2);
            fitted_gaussian = A * exp(exponent);
            fitted_gaussian = reshape(fitted_gaussian, sz);
        else
            % Invalid covariance - fall back to simple Gaussian
            fitted_gaussian = A * exp(-((X - mu_x).^2 / (2*sigma_x^2) + ...
                                        (Y - mu_y).^2 / (2*sigma_y^2) + ...
                                        (Z - mu_z).^2 / (2*sigma_z^2)));
        end
    else
        % Simple 3D Gaussian (no rotation/distortion)
        fitted_gaussian = A * exp(-((X - mu_x).^2 / (2*sigma_x^2) + ...
                                    (Y - mu_y).^2 / (2*sigma_y^2) + ...
                                    (Z - mu_z).^2 / (2*sigma_z^2)));
    end
end

function createRadialSymmetryGraphs(parent, result, channelColor)
% Creates radial symmetry visualization
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col, result.center_z = slice (ARRAY CONVENTION)
%   For display: Convert to CARTESIAN (X=col, Y=row, Z=slice)
    
    % Use the actual fitted data if available, otherwise fall back to raw data
    if isfield(result, 'fittedDataWindow')
        data = result.fittedDataWindow;
    else
        data = result.rawDataWindow;
    end
    
    % Determine if 3D or 2D
    is_3d = (ndims(data) == 3 && size(data, 3) > 1) || ...
            (isfield(result, 'center_z') && ~isnan(result.center_z));
    
    if is_3d
        createRadialSymmetry3D(parent, result, data, channelColor);
    else
        createRadialSymmetry2D(parent, result, data, channelColor);
    end
end

function createRadialSymmetry2D(parent, result, data, channelColor)
% Creates 2D radial symmetry visualization
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col (ARRAY CONVENTION)
%   For plot: Convert to CARTESIAN (X=col, Y=row)
    
    grid = uigridlayout(parent, [1 2]);
    grid.ColumnWidth = {'1x', '1x'};
    
    % Data is already background-corrected if fittedDataWindow was used
    plot_data = double(data);
    
    % Image with center marker
    ax_img = uiaxes(grid);
    ax_img.Layout.Column = 1;
    
    imagesc(ax_img, plot_data);
    axis(ax_img, 'image');
    colormap(ax_img, 'gray');
    hold(ax_img, 'on');
    
    % Add center marker (CARTESIAN CONVERSION for display)
    % result.center_x = row, result.center_y = col (ARRAY CONVENTION)
    % For plot: X=col, Y=row (CARTESIAN)
    plot(ax_img, result.center_y, result.center_x, [channelColor 'o'], ...
         'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', channelColor);
    
    hold(ax_img, 'off');
    title(ax_img, '2D Radial Symmetry', 'FontWeight', 'bold');
    xlabel(ax_img, 'X (pixels)');
    ylabel(ax_img, 'Y (pixels)');
    colorbar(ax_img);
    
    % Radial profile
    ax_radial = uiaxes(grid);
    ax_radial.Layout.Column = 2;
    
    % Calculate and plot radial profile
    [rows, cols] = size(plot_data);
    [X, Y] = meshgrid(1:cols, 1:rows);
    distances = sqrt((X - result.center_x).^2 + (Y - result.center_y).^2);
    
    max_radius = min([result.center_x-1, cols-result.center_x, result.center_y-1, rows-result.center_y]);
    max_radius = max(1, max_radius);
    
    radial_profile = zeros(1, max_radius);
    for r = 1:max_radius
        mask = (distances >= r-0.5) & (distances < r+0.5);
        if any(mask(:))
            radial_profile(r) = mean(plot_data(mask));
        end
    end
    
    plot(ax_radial, 1:max_radius, radial_profile, 'b.-', 'MarkerSize', 8, 'LineWidth', 2);
    title(ax_radial, 'Radial Profile', 'FontWeight', 'bold');
    xlabel(ax_radial, 'Distance (pixels)');
    ylabel(ax_radial, 'Intensity');
    ax_radial.XGrid = 'on';
    ax_radial.YGrid = 'on';
end

function createRadialSymmetry3D(parent, result, data, channelColor)
% Creates 3D radial symmetry visualization with slices and gradient vectors
%
% COORDINATE CONVENTION:
%   result.center_x = row, result.center_y = col, result.center_z = slice (ARRAY CONVENTION)
%   For plot3: Convert to CARTESIAN (X=col, Y=row, Z=slice)
    
    % Data is already background-corrected if fittedDataWindow was used
    plot_data = double(data);
    
    sz = size(plot_data);  % sz = [rows, cols, slices] = [Y, X, Z]
    
    % Create single 3D plot
    grid = uigridlayout(parent, [1 1]);
    ax = uiaxes(grid);
    
    try
        % The center is at the detected radial symmetry position
        % result: center_x=row, center_y=col, center_z=slice (ARRAY CONVENTION)
        % Convert to CARTESIAN for display: X=col, Y=row, Z=slice
        center_x = result.center_y;  % col → X (horizontal)
        center_y = result.center_x;  % row → Y (vertical)
        center_z = result.center_z;  % slice → Z (depth)
        
        % Round to nearest pixel for slicing
        slice_x = round(center_x);
        slice_y = round(center_y);
        slice_z = round(center_z);
        
        % Ensure slices are within bounds
        slice_x = max(1, min(sz(2), slice_x));
        slice_y = max(1, min(sz(1), slice_y));
        slice_z = max(1, min(sz(3), slice_z));
        
        % Create coordinate grids with pixel centers at integers
        x_edges = 0.5:1:(sz(2)+0.5);
        y_edges = 0.5:1:(sz(1)+0.5);
        z_edges = 0.5:1:(sz(3)+0.5);
        
        [X_grid, Y_grid, Z_grid] = meshgrid(x_edges, y_edges, z_edges);
        
        % Pad data by 1 on the 'post' side only to match grid (8×8×8 from 7×7×7)
        plot_data_extended = padarray(plot_data, [1 1 1], 'replicate', 'post');
        
        % Display orthogonal slice planes through the data at the detected center
        h = slice(ax, X_grid, Y_grid, Z_grid, plot_data_extended, slice_x, slice_y, slice_z);
        
        % Style the slices
        for i = 1:length(h)
            set(h(i), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
        
        % Set colormap to grayscale
        colormap(ax, 'gray');
        
        % Calculate 3D gradient for vector overlay
        [Gx, Gy, Gz] = gradient(plot_data);
        
        % Subsample the gradient vectors for cleaner visualization
        step = max(2, floor(min(sz) / 6));
        
        % Create 3D grid of gradient vectors throughout the volume
        hold(ax, 'on');
        
        [X_3d, Y_3d, Z_3d] = meshgrid(1:step:sz(2), 1:step:sz(1), 1:step:sz(3));
        Gx_3d = Gx(1:step:sz(1), 1:step:sz(2), 1:step:sz(3));
        Gy_3d = Gy(1:step:sz(1), 1:step:sz(2), 1:step:sz(3));
        Gz_3d = Gz(1:step:sz(1), 1:step:sz(2), 1:step:sz(3));
        
        % Normalize for consistent arrow lengths
        mag_3d = sqrt(Gx_3d.^2 + Gy_3d.^2 + Gz_3d.^2);
        mag_3d(mag_3d == 0) = 1;
        
        % Plot 3D gradient vector field with thicker arrows
        quiver3(ax, X_3d, Y_3d, Z_3d, ...
               Gx_3d./mag_3d, Gy_3d./mag_3d, Gz_3d./mag_3d, 0.6, ...
               'Color', 'cyan', 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
        
        % Add center point marker (already in CARTESIAN after conversion above)
        plot3(ax, center_x, center_y, center_z, 'o', ...
              'Color', channelColor, 'MarkerSize', 15, 'LineWidth', 3, ...
              'MarkerFaceColor', channelColor);
        
        hold(ax, 'off');
        
        % Set axis labels
        xlabel(ax, 'X (columns →)');
        ylabel(ax, 'Y (rows →)');
        zlabel(ax, 'Z (slices →)');
        
        % Set axis limits so pixel edges are at half-integers (centers at integers)
        xlim(ax, [0.5, sz(2) + 0.5]);
        ylim(ax, [0.5, sz(1) + 0.5]);
        zlim(ax, [0.5, sz(3) + 0.5]);
        
        % Ticks are already at pixel centers (integers 1, 2, 3, ...)
        ax.XTick = 1:sz(2);
        ax.YTick = 1:sz(1);
        ax.ZTick = 1:sz(3);
        
        % Set title and view
        title(ax, '3D Radial Symmetry with Gradient Vectors', 'FontWeight', 'bold');
        view(ax, 3);
        
        % Enable grid
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.ZGrid = 'on';
        
        % Add colorbar
        colorbar(ax);
        
    catch ME
        % Fallback to simple text if 3D plotting fails
        fprintf('ERROR: 3D radial symmetry visualization failed: %s\n', ME.message);
        cla(ax);
        text(ax, 0.5, 0.5, sprintf('3D visualization failed: %s', ME.message), ...
             'HorizontalAlignment', 'center');
    end
end

function plotProfile(ax, data, result, dimension, channelColor)
% Helper function to plot 1D profiles with fits
    
    switch dimension
        case 'X'
            % Use the stored profile that was actually fitted, if available
            if isfield(result, 'x_profile') && ~isempty(result.x_profile)
                profile = result.x_profile;
            elseif isfield(result, 'originalMaximaCoords')
                local_y = round(result.originalMaximaCoords(2) - (size(data, 1) - 1) / 2);
                local_z = round(result.originalMaximaCoords(3) - (size(data, 3) - 1) / 2);
                local_y = max(1, min(size(data, 1), local_y));
                local_z = max(1, min(size(data, 3), local_z));
                profile = squeeze(data(local_y, :, local_z));
            else
                profile = squeeze(sum(sum(data, 1), 3));
            end
            center = result.center_x;
            sigma = result.sigma_x;
            % Use fitted amplitude if available, otherwise fall back to profile max
            if isfield(result, 'amplitude_x') && ~isnan(result.amplitude_x)
                amplitude = result.amplitude_x;
            else
                amplitude = max(profile);
            end
            
        case 'Y'
            % Use the stored profile that was actually fitted, if available
            if isfield(result, 'y_profile') && ~isempty(result.y_profile)
                profile = result.y_profile;
            elseif isfield(result, 'originalMaximaCoords')
                local_x = round(result.originalMaximaCoords(1) - (size(data, 2) - 1) / 2);
                local_z = round(result.originalMaximaCoords(3) - (size(data, 3) - 1) / 2);
                local_x = max(1, min(size(data, 2), local_x));
                local_z = max(1, min(size(data, 3), local_z));
                profile = squeeze(data(:, local_x, local_z));
            else
                profile = squeeze(sum(sum(data, 2), 3));
            end
            center = result.center_y;
            sigma = result.sigma_y;
            % Use fitted amplitude if available, otherwise fall back to profile max
            if isfield(result, 'amplitude_y') && ~isnan(result.amplitude_y)
                amplitude = result.amplitude_y;
            else
                amplitude = max(profile);
            end
            
        case 'Z'
            % Use the stored profile that was actually fitted, if available
            if isfield(result, 'z_profile') && ~isempty(result.z_profile)
                profile = result.z_profile;
            elseif isfield(result, 'originalMaximaCoords')
                local_x = round(result.originalMaximaCoords(1) - (size(data, 2) - 1) / 2);
                local_y = round(result.originalMaximaCoords(2) - (size(data, 1) - 1) / 2);
                local_x = max(1, min(size(data, 2), local_x));
                local_y = max(1, min(size(data, 1), local_y));
                profile = squeeze(data(local_y, local_x, :));
            else
                profile = squeeze(sum(sum(data, 1), 2));
            end
            center = result.center_z;
            sigma = result.sigma_z;
            % Use fitted amplitude if available, otherwise fall back to profile max
            if isfield(result, 'amplitude_z') && ~isnan(result.amplitude_z)
                amplitude = result.amplitude_z;
            else
                amplitude = max(profile);
            end
    end
    
    coords = 1:length(profile);
    plot(ax, coords, double(profile), 'b.', 'MarkerSize', 10);
    hold(ax, 'on');
    
    % Plot fit curve
    if ~isnan(sigma) && sigma > 0
        fine_coords = linspace(1, length(profile), 200);
        amplitude = double(amplitude);
        center = double(center);
        sigma = double(sigma);
        fit_curve = amplitude * exp(-(fine_coords - center).^2 / (2 * sigma^2));
        plot(ax, fine_coords, fit_curve, 'r-', 'LineWidth', 2);
        
        % Add center line
        if center >= 1 && center <= length(profile)
            plot(ax, [center, center], [0, double(max(profile))*1.1], 'g--', 'LineWidth', 2);
        end
        
        legend(ax, 'Data', 'Fit', 'Center', 'Location', 'best');
    else
        legend(ax, 'Data', 'Location', 'best');
    end
    
    hold(ax, 'off');
    title(ax, [dimension ' Profile Fit'], 'FontWeight', 'bold');
    xlabel(ax, [dimension ' (pixels)']);
    ylabel(ax, 'Intensity');
    ax.XGrid = 'on';
    ax.YGrid = 'on';
end

function createBackgroundGraphs(parent, result, channelColor)
% Creates background analysis visualization
    
    grid = uigridlayout(parent, [2 1]);
    grid.RowHeight = {'1x', '1x'};
    
    % Background correction visualization
    ax_bg = uiaxes(grid);
    ax_bg.Layout.Row = 1;
    
    data = result.rawDataWindow;
    
    if isfield(result, 'background') && ~isnan(result.background)
        % Show before/after background correction
        % Use fittedDataWindow for "after" if available to ensure consistency
        if isfield(result, 'fittedDataWindow')
            fitted_data = result.fittedDataWindow;
        else
            fitted_data = data - result.background;
        end
        
        if ndims(data) == 3
            before = max(data, [], 3);
            after = max(fitted_data, [], 3);
        else
            before = data;
            after = fitted_data;
        end
        
        % Create side-by-side comparison
        combined = [before, after];
        imagesc(ax_bg, combined);
        axis(ax_bg, 'image');
        colormap(ax_bg, 'gray');
        
        % Add dividing line
        hold(ax_bg, 'on');
        plot(ax_bg, [size(before, 2) + 0.5, size(before, 2) + 0.5], ...
             [1, size(before, 1)], 'r--', 'LineWidth', 2);
        hold(ax_bg, 'off');
        
        title(ax_bg, sprintf('Background Correction (%.2f)', result.background), 'FontWeight', 'bold');
        xlabel(ax_bg, 'Before | After');
        colorbar(ax_bg);
    else
        % No background correction applied
        if ndims(data) == 3
            plot_data = max(data, [], 3);
        else
            plot_data = data;
        end
        
        imagesc(ax_bg, plot_data);
        axis(ax_bg, 'image');
        colormap(ax_bg, 'gray');
        title(ax_bg, 'No Background Correction', 'FontWeight', 'bold');
        colorbar(ax_bg);
    end
    
    % Background method information
    ax_method = uiaxes(grid);
    ax_method.Layout.Row = 2;
    
    % Create background method summary
    if isfield(result, 'backgroundMethod')
        method_text = result.backgroundMethod;
    else
        method_text = 'Unknown';
    end
    
    if isfield(result, 'background') && ~isnan(result.background)
        bg_value = result.background;
        bg_text = sprintf('Background Value: %.2f', bg_value);
    else
        bg_value = 0;
        bg_text = 'No background correction';
    end
    
    % Simple bar chart showing background level
    bar(ax_method, 1, bg_value, 'FaceColor', channelColor);
    title(ax_method, 'Background Level', 'FontWeight', 'bold');
    ylabel(ax_method, 'Background Intensity');
    ax_method.XTickLabel = {bg_text};
    ax_method.XGrid = 'on';
    ax_method.YGrid = 'on';
end

function createNumericalInfo(parent, result, channel_idx, maxima_idx)
% Creates comprehensive numerical information display
    
    grid = uigridlayout(parent, [2 1]);
    grid.RowHeight = {'1x', 'fit'};
    
    % Main parameters table
    param_table = uitable(grid);
    param_table.Layout.Row = 1;
    
    % Collect all parameters
    param_names = {};
    param_values = {};
    param_units = {};
    
    % Basic identification
    param_names{end+1} = 'Channel';
    param_values{end+1} = num2str(channel_idx);
    param_units{end+1} = '';
    
    param_names{end+1} = 'Maxima Index';
    param_values{end+1} = num2str(maxima_idx);
    param_units{end+1} = '';
    
    param_names{end+1} = 'Fit Method';
    param_values{end+1} = result.fitMethod;
    param_units{end+1} = '';
    
    % Fit parameters
    % Check if individual X, Y, Z amplitudes are available (for Gaussian fits)
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
        % Only show amplitude if it's not NaN (skip for radial symmetry)
        param_names{end+1} = 'Amplitude';
        param_values{end+1} = sprintf('%.2f', result.amplitude);
        param_units{end+1} = 'counts';
    end
    
    if isfield(result, 'center_x') && ~isnan(result.center_x)
        param_names{end+1} = 'Center X';
        param_values{end+1} = sprintf('%.3f', result.center_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_y') && ~isnan(result.center_y)
        param_names{end+1} = 'Center Y';
        param_values{end+1} = sprintf('%.3f', result.center_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'center_z') && ~isnan(result.center_z)
        param_names{end+1} = 'Center Z';
        param_values{end+1} = sprintf('%.3f', result.center_z);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_x') && ~isnan(result.sigma_x)
        param_names{end+1} = 'Sigma X';
        param_values{end+1} = sprintf('%.3f', result.sigma_x);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_y') && ~isnan(result.sigma_y)
        param_names{end+1} = 'Sigma Y';
        param_values{end+1} = sprintf('%.3f', result.sigma_y);
        param_units{end+1} = 'pixels';
    end
    
    if isfield(result, 'sigma_z') && ~isnan(result.sigma_z)
        param_names{end+1} = 'Sigma Z';
        param_values{end+1} = sprintf('%.3f', result.sigma_z);
        param_units{end+1} = 'pixels';
    end
    
    % For Distorted 3D Gaussian, show rho (correlation) values
    if strcmp(result.fitMethod, 'Distorted 3D Gaussian')
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
    end
    
    param_names{end+1} = '---';
    param_values{end+1} = '---';
    param_units{end+1} = '---';
    
    % Fit window and coordinate information
    if isfield(result, 'fitWindowDimensions') && ~isempty(result.fitWindowDimensions)
        if length(result.fitWindowDimensions) == 3
            param_names{end+1} = 'Fit Window Size';
            param_values{end+1} = sprintf('%d×%d×%d', result.fitWindowDimensions(1), result.fitWindowDimensions(2), result.fitWindowDimensions(3));
            param_units{end+1} = 'pixels';
        else
            param_names{end+1} = 'Fit Window Size';
            param_values{end+1} = sprintf('%d×%d', result.fitWindowDimensions(1), result.fitWindowDimensions(2));
            param_units{end+1} = 'pixels';
        end
    end
    
    if isfield(result, 'originalMaximaCoords') && ~isempty(result.originalMaximaCoords)
        % originalMaximaCoords is stored in ARRAY CONVENTION [row, col, slice]
        % Convert to Cartesian for display [x, y, z] = [col, row, slice]
        if length(result.originalMaximaCoords) >= 3
            row = result.originalMaximaCoords(1);
            col = result.originalMaximaCoords(2);
            slice = result.originalMaximaCoords(3);
            param_names{end+1} = 'Maxima (array)';
            param_values{end+1} = sprintf('[%.1f, %.1f, %.1f]', row, col, slice);
            param_units{end+1} = '[row,col,slice]';
            param_names{end+1} = 'Maxima (display)';
            param_values{end+1} = sprintf('[%.1f, %.1f, %.1f]', col, row, slice);
            param_units{end+1} = '[x,y,z]';
        elseif length(result.originalMaximaCoords) >= 2
            row = result.originalMaximaCoords(1);
            col = result.originalMaximaCoords(2);
            param_names{end+1} = 'Maxima (array)';
            param_values{end+1} = sprintf('[%.1f, %.1f]', row, col);
            param_units{end+1} = '[row,col]';
            param_names{end+1} = 'Maxima (display)';
            param_values{end+1} = sprintf('[%.1f, %.1f]', col, row);
            param_units{end+1} = '[x,y]';
        end
    end
    
    if isfield(result, 'localMaximaInWindow') && ~isempty(result.localMaximaInWindow)
        % localMaximaInWindow is in ARRAY CONVENTION [row, col, slice]
        if length(result.localMaximaInWindow) >= 3
            row = result.localMaximaInWindow(1);
            col = result.localMaximaInWindow(2);
            slice = result.localMaximaInWindow(3);
            param_names{end+1} = 'Maxima in Window';
            param_values{end+1} = sprintf('[%.1f, %.1f, %.1f]', row, col, slice);
            param_units{end+1} = '[row,col,slice]';
        elseif length(result.localMaximaInWindow) >= 2
            row = result.localMaximaInWindow(1);
            col = result.localMaximaInWindow(2);
            param_names{end+1} = 'Maxima in Window';
            param_values{end+1} = sprintf('[%.1f, %.1f]', row, col);
            param_units{end+1} = '[row,col]';
        end
    end
    
    if isfield(result, 'globalFitCenter') && ~isempty(result.globalFitCenter) && ~all(isnan(result.globalFitCenter))
        % globalFitCenter is in ARRAY CONVENTION [row, col, slice]
        % Convert to Cartesian for display
        if length(result.globalFitCenter) >= 3 && ~isnan(result.globalFitCenter(3))
            row = result.globalFitCenter(1);
            col = result.globalFitCenter(2);
            slice = result.globalFitCenter(3);
            param_names{end+1} = 'Fit Center (array)';
            param_values{end+1} = sprintf('[%.3f, %.3f, %.3f]', row, col, slice);
            param_units{end+1} = '[row,col,slice]';
            param_names{end+1} = 'Fit Center (display)';
            param_values{end+1} = sprintf('[%.3f, %.3f, %.3f]', col, row, slice);
            param_units{end+1} = '[x,y,z]';
        elseif length(result.globalFitCenter) >= 2 && ~isnan(result.globalFitCenter(1)) && ~isnan(result.globalFitCenter(2))
            row = result.globalFitCenter(1);
            col = result.globalFitCenter(2);
            param_names{end+1} = 'Fit Center (array)';
            param_values{end+1} = sprintf('[%.3f, %.3f]', row, col);
            param_units{end+1} = '[row,col]';
            param_names{end+1} = 'Fit Center (display)';
            param_values{end+1} = sprintf('[%.3f, %.3f]', col, row);
            param_units{end+1} = '[x,y]';
        end
    end
    
    param_names{end+1} = '---';
    param_values{end+1} = '---';
    param_units{end+1} = '---';
    
    % Quality metrics
    % For Radial Symmetry, show quality score instead of R²
    if strcmp(result.fitMethod, 'Radial Symmetry')
        if isfield(result, 'radialSymmetryScore') && ~isnan(result.radialSymmetryScore)
            param_names{end+1} = 'Symmetry Score';
            param_values{end+1} = sprintf('%.4f', result.radialSymmetryScore);
            param_units{end+1} = '';
        end
        if isfield(result, 'radialSymmetryQuality') && ~isnan(result.radialSymmetryQuality)
            param_names{end+1} = 'Quality Score';
            param_values{end+1} = sprintf('%.4f', result.radialSymmetryQuality);
            param_units{end+1} = '';
        end
    elseif isfield(result, 'r_squared') && ~all(isnan(result.r_squared))
        % For Gaussian fits, show R²
        if isscalar(result.r_squared)
            param_names{end+1} = 'R²';
            param_values{end+1} = sprintf('%.4f', result.r_squared);
            param_units{end+1} = '';
        else
            % Use proper dimension labels for R²
            if any(strcmp(result.fitMethod, {'1D (X,Y,Z)', '1D (X,Y,Z) Gaussian', '2D (XY) + 1D (Z)', '2D (XY) + 1D (Z) Gaussian'}))
                % For 1D X,Y,Z fits and 2D+1D fits: label as X, Y, Z
                dimension_labels = {'X', 'Y', 'Z'};
                for i = 1:length(result.r_squared)
                    if i <= length(dimension_labels)
                        param_names{end+1} = sprintf('R² %s', dimension_labels{i});
                    else
                        param_names{end+1} = sprintf('R² (%d)', i);
                    end
                    param_values{end+1} = sprintf('%.4f', result.r_squared(i));
                    param_units{end+1} = '';
                end
            else
                % Generic fallback
                for i = 1:length(result.r_squared)
                    param_names{end+1} = sprintf('R² (%d)', i);
                    param_values{end+1} = sprintf('%.4f', result.r_squared(i));
                    param_units{end+1} = '';
                end
            end
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
    
    % Set table data
    param_table.Data = [param_names', param_values', param_units'];
    param_table.ColumnName = {'Parameter', 'Value', 'Unit'};
    param_table.ColumnWidth = {150, 120, 80};
    param_table.RowName = {};
    param_table.FontSize = 10;
    
    % Summary statistics
    summary_panel = uipanel(grid);
    summary_panel.Layout.Row = 2;
    summary_panel.Title = 'Fit Quality Summary';
    summary_panel.FontWeight = 'bold';
    
    summary_grid = uigridlayout(summary_panel, [3 1]);
    
    % Quality assessment
    if strcmp(result.fitMethod, 'Radial Symmetry')
        % For Radial Symmetry, use quality score
        if isfield(result, 'radialSymmetryQuality')
            quality_val = result.radialSymmetryQuality;
            if quality_val > 0.25
                quality_text = sprintf('High Quality (Score = %.3f)', quality_val);
                quality_color = [0, 0.7, 0]; % Green
            elseif quality_val > 0.15
                quality_text = sprintf('Medium Quality (Score = %.3f)', quality_val);
                quality_color = [0.8, 0.6, 0]; % Orange
            else
                quality_text = sprintf('Low Quality (Score = %.3f)', quality_val);
                quality_color = [0.8, 0, 0]; % Red
            end
        else
            quality_text = 'Quality: Unknown';
            quality_color = [0.5, 0.5, 0.5]; % Gray
        end
    elseif isfield(result, 'r_squared')
        % For Gaussian fits, use R²
        r2_val = result.r_squared(1);
        if r2_val > 0.9
            quality_text = sprintf('Excellent Fit (R² = %.3f)', r2_val);
            quality_color = [0, 0.7, 0]; % Green
        elseif r2_val > 0.7
            quality_text = sprintf('Good Fit (R² = %.3f)', r2_val);
            quality_color = [0.8, 0.6, 0]; % Orange
        else
            quality_text = sprintf('Poor Fit (R² = %.3f)', r2_val);
            quality_color = [0.8, 0, 0]; % Red
        end
    else
        quality_text = 'Fit Quality: Unknown';
        quality_color = [0.5, 0.5, 0.5]; % Gray
    end
    
    quality_label = uilabel(summary_grid, 'Text', quality_text, 'FontWeight', 'bold', ...
                           'FontColor', quality_color, 'HorizontalAlignment', 'center');
    quality_label.Layout.Row = 1;
    
    % Method-specific info
    method_info = sprintf('Method: %s', result.fitMethod);
    method_label = uilabel(summary_grid, 'Text', method_info, 'HorizontalAlignment', 'center');
    method_label.Layout.Row = 2;
    
    % Data window info
    data_size = size(result.rawDataWindow);
    if length(data_size) == 3
        window_info = sprintf('Data Window: %dx%dx%d pixels', data_size(1), data_size(2), data_size(3));
    else
        window_info = sprintf('Data Window: %dx%d pixels', data_size(1), data_size(2));
    end
    window_label = uilabel(summary_grid, 'Text', window_info, 'HorizontalAlignment', 'center');
    window_label.Layout.Row = 3;
end

function intensity = model_2d_gaussian(p, xy)
% Simple 2D Gaussian model for contour plotting
    % Convert all parameters to double for safety
    A = double(p(1));
    mu_x = double(p(2));
    mu_y = double(p(3));
    sig_x = double(p(4));
    sig_y = double(p(5));
    
    % Validate parameters
    if A <= 0 || sig_x <= 0 || sig_y <= 0
        intensity = zeros(size(xy, 1), 1);
        return;
    end
    
    if length(p) > 5
        rho = double(p(6));
        % Clamp rho to valid range for asin
        if abs(rho) > 1
            rho = sign(rho) * 0.99; % Slightly less than 1 to avoid numerical issues
        end
        if abs(rho) < 1e-10
            theta = 0; % Avoid numerical issues for very small rho
        else
            theta = asin(rho);
        end
    else
        theta = 0;
    end
    
    % Ensure xy is double
    xy = double(xy);
    
    x_centered = xy(:,1) - mu_x;
    y_centered = xy(:,2) - mu_y;
    
    if abs(theta) < 1e-10
        % No rotation case - more numerically stable
        intensity = A * exp(-((x_centered.^2 / (2*sig_x^2)) + (y_centered.^2 / (2*sig_y^2))));
    else
        % Rotation case
        cos_theta = cos(theta);
        sin_theta = sin(theta);
        
        x_rot = x_centered * cos_theta + y_centered * sin_theta;
        y_rot = -x_centered * sin_theta + y_centered * cos_theta;
        
        intensity = A * exp(-((x_rot.^2 / (2*sig_x^2)) + (y_rot.^2 / (2*sig_y^2))));
    end
    
    % Clean up any non-finite values
    intensity(~isfinite(intensity)) = 0;
    
    % Ensure we have some variation in the output
    if max(intensity) == min(intensity)
        fprintf('Warning: 2D Gaussian model produced constant output\n');
    end
end
