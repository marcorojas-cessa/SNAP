function visualizeRadialSymmetry2D(data, center, fitParams, varargin)
% VISUALIZERADIALSYMMETRY2D Create 2D vector map visualization for radial symmetry fits
%
% INPUTS:
%   data - 2D image data around the detected spot
%   center - [y, x] coordinates of the detected center (in image coordinates)
%   fitParams - fitting parameters structure
%   varargin - optional parameters:
%       'Title' - figure title (default: 'Radial Symmetry 2D Vector Map')
%       'VectorScale' - scaling factor for vector display (default: auto)
%       'ShowScoreMap' - show score map as overlay (default: true)
%       'VectorDensity' - subsample factor for vectors (default: 3)
%       'ColorMap' - colormap for intensity (default: 'hot')
%       'XYSpacing' - XY pixel spacing in microns (default: 1)

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'Title', 'Radial Symmetry 2D Vector Map', @ischar);
    addParameter(p, 'VectorScale', [], @isnumeric);
    addParameter(p, 'ShowScoreMap', true, @islogical);
    addParameter(p, 'VectorDensity', 3, @isnumeric);
    addParameter(p, 'ColorMap', 'hot', @ischar);
    addParameter(p, 'XYSpacing', 1, @isnumeric);
    addParameter(p, 'Parent', [], @(x) isempty(x) || isa(x, 'matlab.ui.container.Panel') || isa(x, 'matlab.ui.Figure'));
    parse(p, varargin{:});
    
    title_str = p.Results.Title;
    vector_scale = p.Results.VectorScale;
    show_score_map = p.Results.ShowScoreMap;
    vector_density = p.Results.VectorDensity;
    color_map = p.Results.ColorMap;
    xy_spacing = p.Results.XYSpacing;
    parent = p.Results.Parent;
    
    [rows, cols] = size(data);
    
    % Calculate 2D gradients
    [gx, gy] = gradient(double(data));
    
    % Scale gradients by pixel spacing
    gx = gx / xy_spacing;
    gy = gy / xy_spacing;
    
    % Calculate gradient magnitudes
    grad_mag = sqrt(gx.^2 + gy.^2);
    
    % Create coordinate grids (scaled by spacing)
    [X, Y] = meshgrid((1:cols) * xy_spacing, (1:rows) * xy_spacing);
    
    % Convert center coordinates to physical units
    center_phys = [center(1) * xy_spacing, center(2) * xy_spacing];
    
    % Calculate radial symmetry score map for visualization
    score_map = calculateScoreMap2D(data, gx, gy, fitParams.gaussFitRadialRadius);
    
    % Create figure with subplots
    % Use provided parent or create figure with subplots
    if ~isempty(parent)
        % Use the provided parent panel
        fig = ancestor(parent, 'figure');
        % Clear any existing content in the parent
        delete(parent.Children);
        % Disable AutoResizeChildren to allow subplots
        parent.AutoResizeChildren = 'off';
    else
        % Create new figure
        fig = figure('Name', title_str, 'Position', [100, 100, 1200, 800]);
        parent = fig;
    end
    
    % Subplot 1: Main vector field visualization
    subplot(2, 3, [1, 2, 4, 5], 'Parent', parent);
    
    % Show intensity as background
    imagesc((1:cols) * xy_spacing, (1:rows) * xy_spacing, data);
    colormap(gca, 'gray');
    hold on;
    
    % Show score map overlay if requested
    if show_score_map && ~isempty(score_map)
        % Normalize score map for overlay
        score_normalized = score_map / max(score_map(:));
        
        % Create alpha mask for score map
        alpha_mask = score_normalized;
        alpha_mask(alpha_mask < 0.1) = 0; % Hide very low scores
        
        % Overlay score map with transparency
        h_overlay = imagesc((1:cols) * xy_spacing, (1:rows) * xy_spacing, score_map);
        colormap(h_overlay.Parent, color_map);
        alpha(h_overlay, 0.6);
    end
    
    % Subsample vectors for cleaner visualization
    step = vector_density;
    [x_sub, y_sub] = meshgrid(1:step:cols, 1:step:rows);
    
    % Convert to physical coordinates
    x_sub_phys = x_sub * xy_spacing;
    y_sub_phys = y_sub * xy_spacing;
    
    % Get subsampled gradient components
    gx_sub = gx(1:step:rows, 1:step:cols);
    gy_sub = gy(1:step:rows, 1:step:cols);
    grad_mag_sub = grad_mag(1:step:rows, 1:step:cols);
    
    % Filter out weak gradients for cleaner visualization
    threshold = 0.15 * max(grad_mag_sub(:));
    strong_gradients = grad_mag_sub > threshold;
    
    % Auto-scale vectors if not specified
    if isempty(vector_scale)
        max_grad = max(grad_mag_sub(strong_gradients));
        max_distance = max([range(x_sub_phys(:)), range(y_sub_phys(:))]);
        vector_scale = max_distance / (15 * max_grad);
    end
    
    % Create 2D quiver plot
    quiver(x_sub_phys(strong_gradients), y_sub_phys(strong_gradients), ...
           gx_sub(strong_gradients) * vector_scale, gy_sub(strong_gradients) * vector_scale, ...
           0, 'Color', [0.2, 0.8, 1.0], 'LineWidth', 1.5, 'MaxHeadSize', 0.8);
    
    % Add convergence vectors from sample points to center
    addConvergenceVectors2D(x_sub_phys, y_sub_phys, gx_sub, gy_sub, grad_mag_sub, ...
                           center_phys, strong_gradients, fitParams.gaussFitRadialRadius, vector_scale);
    
    % Mark the detected center
    scatter(center_phys(2), center_phys(1), 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Formatting
    xlabel('X (μm)', 'FontSize', 12);
    ylabel('Y (μm)', 'FontSize', 12);
    title('2D Gradient Vector Field & Radial Symmetry', 'FontSize', 14);
    axis equal tight;
    grid on;
    
    hold off;
    
    % Subplot 2: Score Map
    subplot(2, 3, 3, 'Parent', parent);
    if ~isempty(score_map)
        imagesc((1:cols) * xy_spacing, (1:rows) * xy_spacing, score_map);
        colormap(gca, color_map);
        colorbar;
        hold on;
        scatter(center_phys(2), center_phys(1), 100, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        hold off;
        xlabel('X (μm)');
        ylabel('Y (μm)');
        title('Radial Symmetry Score Map');
        axis equal tight;
    end
    
    % Subplot 3: Radial Profile
    subplot(2, 3, 6, 'Parent', parent);
    plotRadialProfile2D(data, center, xy_spacing, grad_mag);
    
    % Add comprehensive information panel
    info_ax = axes('Position', [0.02, 0.02, 0.30, 0.35], 'Visible', 'off', 'Parent', parent);
    
    % Title
    text(0, 0.95, 'Radial Symmetry Analysis (2D)', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k');
    
    % Method explanation
    text(0, 0.88, 'Method:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.84, 'Radial symmetry detection analyzes how', 'FontSize', 9, 'Color', 'k');
    text(0, 0.80, 'image gradients converge to find spot centers.', 'FontSize', 9, 'Color', 'k');
    
    % Legend
    text(0, 0.72, 'Visualization Elements:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.68, '→ Gradient vectors (cyan arrows)', 'Color', [0.2, 0.8, 1.0], 'FontSize', 9);
    text(0, 0.64, '● Detected center (red circle)', 'Color', 'r', 'FontSize', 9);
    text(0, 0.60, '--- Convergence paths (orange)', 'Color', [0.8, 0.4, 0.2], 'FontSize', 9);
    if show_score_map
        text(0, 0.56, '█ Score map overlay (hot colors)', 'Color', 'r', 'FontSize', 9);
    end
    
    % Algorithm explanation
    text(0, 0.48, 'How it works:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.44, '1. Calculate 2D gradients at each pixel', 'FontSize', 9, 'Color', 'k');
    text(0, 0.40, '2. Project gradients along their direction', 'FontSize', 9, 'Color', 'k');
    text(0, 0.36, '3. Find where most gradients converge', 'FontSize', 9, 'Color', 'k');
    text(0, 0.32, '4. Refine center to sub-pixel precision', 'FontSize', 9, 'Color', 'k');
    
    % Quality interpretation
    text(0, 0.24, 'Quality Score Interpretation:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.20, '0.8-1.0: Excellent radial symmetry', 'FontSize', 9, 'Color', 'g');
    text(0, 0.16, '0.6-0.8: Good radial symmetry', 'FontSize', 9, 'Color', [0.8, 0.6, 0]);
    text(0, 0.12, '0.4-0.6: Moderate radial symmetry', 'FontSize', 9, 'Color', [0.8, 0.4, 0]);
    text(0, 0.08, '0.0-0.4: Poor radial symmetry', 'FontSize', 9, 'Color', 'r');
    
    text(0, 0.02, 'Higher scores = more gradients converge to center', 'FontSize', 8, 'Color', [0.5, 0.5, 0.5]);
end

function score_map = calculateScoreMap2D(data, gx, gy, radius)
    % Calculate radial symmetry score map for visualization
    [rows, cols] = size(data);
    score_map = zeros(rows, cols);
    
    for r = 1:rows
        for c = 1:cols
            if gx(r,c) == 0 && gy(r,c) == 0
                continue;
            end
            
            grad_mag = sqrt(gx(r,c)^2 + gy(r,c)^2);
            if grad_mag == 0, continue; end
            
            grad_unit_x = gx(r,c) / grad_mag;
            grad_unit_y = gy(r,c) / grad_mag;
            
            for d = 1:radius
                center_x = c + d * grad_unit_x;
                center_y = r + d * grad_unit_y;
                
                if center_x >= 1 && center_x <= cols && center_y >= 1 && center_y <= rows
                    center_idx_x = round(center_x);
                    center_idx_y = round(center_y);
                    
                    score_map(center_idx_y, center_idx_x) = ...
                        score_map(center_idx_y, center_idx_x) + grad_mag;
                end
            end
        end
    end
end

function addConvergenceVectors2D(x_sub, y_sub, gx_sub, gy_sub, grad_mag_sub, center_phys, strong_gradients, radius, vector_scale)
    % Add vectors showing convergence to center from selected points
    
    % Select a few representative points
    indices = find(strong_gradients);
    if length(indices) > 15
        % Randomly sample 15 points for clarity
        sample_idx = indices(randperm(length(indices), 15));
    else
        sample_idx = indices;
    end
    
    for i = 1:length(sample_idx)
        idx = sample_idx(i);
        
        % Current position
        pos = [x_sub(idx), y_sub(idx)];
        
        % Gradient direction (unit vector)
        grad_mag = sqrt(gx_sub(idx)^2 + gy_sub(idx)^2);
        if grad_mag > 0
            grad_unit = [gx_sub(idx), gy_sub(idx)] / grad_mag;
            
            % Project along gradient for visualization
            end_pos = pos + radius * grad_unit * vector_scale;
            
            % Draw convergence line
            plot([pos(1), end_pos(1)], [pos(2), end_pos(2)], ...
                 '--', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 1.5, 'Alpha', 0.8);
        end
    end
end

function plotRadialProfile2D(data, center, xy_spacing, grad_mag)
    % Plot radial intensity and gradient profiles
    
    [rows, cols] = size(data);
    center_y = round(center(1));
    center_x = round(center(2));
    
    % Ensure center is within bounds
    center_y = max(1, min(rows, center_y));
    center_x = max(1, min(cols, center_x));
    
    % Create radial distance map
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Physical distances
    dx = (X - center_x) * xy_spacing;
    dy = (Y - center_y) * xy_spacing;
    
    radial_dist = sqrt(dx.^2 + dy.^2);
    
    % Create radial bins
    max_dist = min([center_x, center_y, cols-center_x, rows-center_y]) * xy_spacing;
    dist_bins = 0:xy_spacing:max_dist;
    
    intensity_profile = zeros(size(dist_bins));
    gradient_profile = zeros(size(dist_bins));
    
    for i = 1:length(dist_bins)-1
        mask = radial_dist >= dist_bins(i) & radial_dist < dist_bins(i+1);
        if any(mask(:))
            intensity_profile(i) = mean(data(mask));
            gradient_profile(i) = mean(grad_mag(mask));
        end
    end
    
    % Plot profiles
    yyaxis left;
    plot(dist_bins(1:end-1), intensity_profile(1:end-1), 'b-', 'LineWidth', 2);
    ylabel('Intensity', 'Color', 'b');
    
    yyaxis right;
    plot(dist_bins(1:end-1), gradient_profile(1:end-1), 'r--', 'LineWidth', 2);
    ylabel('Gradient Magnitude', 'Color', 'r');
    
    xlabel('Radial Distance (μm)');
    title('Radial Intensity & Gradient Profiles');
    grid on;
    legend('Intensity', 'Gradient', 'Location', 'best');
end
