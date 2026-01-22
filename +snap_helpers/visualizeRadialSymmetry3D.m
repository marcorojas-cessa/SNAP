function visualizeRadialSymmetry3D(data, center, fitParams, varargin)
% VISUALIZERADIALSYMMETRY3D Create 3D vector map visualization for radial symmetry fits
%
% INPUTS:
%   data - 3D image data around the detected spot
%   center - [y, x, z] coordinates of the detected center (in image coordinates)
%   fitParams - fitting parameters structure
%   varargin - optional parameters:
%       'Title' - figure title (default: 'Radial Symmetry 3D Vector Map')
%       'VectorScale' - scaling factor for vector display (default: auto)
%       'ShowScoreMap' - show score map as volume rendering (default: true)
%       'VectorDensity' - subsample factor for vectors (default: 2)
%       'ColorMap' - colormap for intensity (default: 'hot')
%       'XYSpacing' - XY pixel spacing in microns (default: 1)
%       'ZSpacing' - Z pixel spacing in microns (default: 1)

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'Title', 'Radial Symmetry 3D Vector Map', @ischar);
    addParameter(p, 'VectorScale', [], @isnumeric);
    addParameter(p, 'ShowScoreMap', true, @islogical);
    addParameter(p, 'VectorDensity', 2, @isnumeric);
    addParameter(p, 'ColorMap', 'hot', @ischar);
    addParameter(p, 'XYSpacing', 1, @isnumeric);
    addParameter(p, 'ZSpacing', 1, @isnumeric);
    addParameter(p, 'Parent', [], @(x) isempty(x) || isa(x, 'matlab.ui.container.Panel') || isa(x, 'matlab.ui.Figure'));
    parse(p, varargin{:});
    
    title_str = p.Results.Title;
    vector_scale = p.Results.VectorScale;
    show_score_map = p.Results.ShowScoreMap;
    vector_density = p.Results.VectorDensity;
    color_map = p.Results.ColorMap;
    xy_spacing = p.Results.XYSpacing;
    z_spacing = p.Results.ZSpacing;
    parent = p.Results.Parent;
    
    [rows, cols, slices] = size(data);
    
    % Calculate 3D gradients
    [gx, gy, gz] = gradient(double(data));
    
    % Scale gradients by pixel spacing for proper 3D geometry
    gx = gx / xy_spacing;
    gy = gy / xy_spacing;
    gz = gz / z_spacing;
    
    % Calculate gradient magnitudes
    grad_mag = sqrt(gx.^2 + gy.^2 + gz.^2);
    
    % Create coordinate grids (scaled by spacing)
    [X, Y, Z] = meshgrid((1:cols) * xy_spacing, (1:rows) * xy_spacing, (1:slices) * z_spacing);
    
    % Convert center coordinates to physical units
    center_phys = [center(1) * xy_spacing, center(2) * xy_spacing, center(3) * z_spacing];
    
    % Calculate radial symmetry score map for visualization
    score_map = calculateScoreMap3D(data, gx, gy, gz, fitParams.gaussFitRadialRadius);
    
    % Create figure with multiple subplots
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
        fig = figure('Name', title_str, 'Position', [100, 100, 1400, 800]);
        parent = fig;
    end
    
    % Subplot 1: 3D Vector Field
    subplot(2, 3, [1, 2, 4, 5], 'Parent', parent);
    
    % Subsample vectors for cleaner visualization
    step = vector_density;
    [x_sub, y_sub, z_sub] = meshgrid(1:step:cols, 1:step:rows, 1:step:slices);
    
    % Convert to physical coordinates
    x_sub_phys = x_sub * xy_spacing;
    y_sub_phys = y_sub * xy_spacing;
    z_sub_phys = z_sub * z_spacing;
    
    % Get subsampled gradient components
    gx_sub = gx(1:step:rows, 1:step:cols, 1:step:slices);
    gy_sub = gy(1:step:rows, 1:step:cols, 1:step:slices);
    gz_sub = gz(1:step:rows, 1:step:cols, 1:step:slices);
    grad_mag_sub = grad_mag(1:step:rows, 1:step:cols, 1:step:slices);
    
    % Filter out weak gradients for cleaner visualization
    threshold = 0.1 * max(grad_mag_sub(:));
    strong_gradients = grad_mag_sub > threshold;
    
    % Create 3D quiver plot
    quiver3(x_sub_phys(strong_gradients), y_sub_phys(strong_gradients), z_sub_phys(strong_gradients), ...
           gx_sub(strong_gradients), gy_sub(strong_gradients), gz_sub(strong_gradients), ...
           'Color', [0.2, 0.6, 1.0], 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    
    hold on;
    
    % Show volume rendering of score map if requested
    if show_score_map && ~isempty(score_map)
        % Create isosurface at 50% of maximum score
        iso_value = 0.5 * max(score_map(:));
        if iso_value > 0
            try
                [faces, vertices] = isosurface(X, Y, Z, score_map, iso_value);
                if ~isempty(faces)
                    patch('Vertices', vertices, 'Faces', faces, ...
                          'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                end
            catch
                % If isosurface fails, show slice planes instead
                showSlicePlanes(X, Y, Z, score_map, center_phys, color_map);
            end
        end
    end
    
    % Mark the detected center
    scatter3(center_phys(2), center_phys(1), center_phys(3), 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Add convergence vectors from sample points to center
    addConvergenceVectors(x_sub_phys, y_sub_phys, z_sub_phys, gx_sub, gy_sub, gz_sub, ...
                         grad_mag_sub, center_phys, strong_gradients, fitParams.gaussFitRadialRadius);
    
    % Formatting
    xlabel('X (μm)', 'FontSize', 12);
    ylabel('Y (μm)', 'FontSize', 12);
    zlabel('Z (μm)', 'FontSize', 12);
    title('3D Gradient Vector Field & Radial Symmetry', 'FontSize', 14);
    grid on;
    axis equal;
    view(45, 30);
    
    % Add lighting for better 3D perception
    camlight('headlight');
    lighting gouraud;
    
    hold off;
    
    % Subplot 2: Score Map Cross-sections
    subplot(2, 3, 3, 'Parent', parent);
    if ~isempty(score_map)
        % XY slice at center Z
        center_z_idx = round(center(3));
        center_z_idx = max(1, min(slices, center_z_idx));
        
        imagesc((1:cols) * xy_spacing, (1:rows) * xy_spacing, score_map(:, :, center_z_idx));
        colormap(gca, color_map);
        colorbar;
        hold on;
        scatter(center_phys(2), center_phys(1), 100, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        hold off;
        xlabel('X (μm)');
        ylabel('Y (μm)');
        title(sprintf('Score Map XY (Z=%.2f μm)', center_phys(3)));
        axis equal tight;
    end
    
    % Subplot 3: Radial Profile
    subplot(2, 3, 6, 'Parent', parent);
    plotRadialProfile(data, center, xy_spacing, z_spacing, grad_mag);
    
    % Add comprehensive information panel
    info_ax = axes('Position', [0.02, 0.02, 0.30, 0.35], 'Visible', 'off', 'Parent', parent);
    
    % Title
    text(0, 0.95, 'Radial Symmetry Analysis', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k');
    
    % Method explanation
    text(0, 0.88, 'Method:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.84, 'Radial symmetry detection analyzes how', 'FontSize', 9, 'Color', 'k');
    text(0, 0.80, 'image gradients converge to find spot centers.', 'FontSize', 9, 'Color', 'k');
    
    % Legend
    text(0, 0.72, 'Visualization Elements:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.68, '→ Gradient vectors (blue arrows)', 'Color', [0.2, 0.6, 1.0], 'FontSize', 9);
    text(0, 0.64, '● Detected center (red sphere)', 'Color', 'r', 'FontSize', 9);
    text(0, 0.60, '--- Convergence paths (orange)', 'Color', [0.8, 0.4, 0.2], 'FontSize', 9);
    if show_score_map
        text(0, 0.56, '█ Score isosurface (red surface)', 'Color', 'r', 'FontSize', 9);
    end
    
    % Algorithm explanation
    text(0, 0.48, 'How it works:', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'k');
    text(0, 0.44, '1. Calculate 3D gradients at each voxel', 'FontSize', 9, 'Color', 'k');
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
    
    % Auto-scale vectors if not specified
    if isempty(vector_scale)
        % Calculate appropriate vector scaling
        max_grad = max(grad_mag_sub(strong_gradients));
        max_distance = max([range(x_sub_phys(:)), range(y_sub_phys(:)), range(z_sub_phys(:))]);
        vector_scale = max_distance / (10 * max_grad);
    end
    
    % Apply vector scaling
    children = get(gca, 'Children');
    for i = 1:length(children)
        if strcmp(get(children(i), 'Type'), 'line') && ...
           strcmp(get(children(i), 'Marker'), 'none')
            % This is likely a quiver line - scale it
            try
                set(children(i), 'XData', get(children(i), 'XData') * vector_scale, ...
                                 'YData', get(children(i), 'YData') * vector_scale, ...
                                 'ZData', get(children(i), 'ZData') * vector_scale);
            catch
                % Ignore if scaling fails
            end
        end
    end
end

function score_map = calculateScoreMap3D(data, gx, gy, gz, radius)
    % Calculate radial symmetry score map for visualization
    [rows, cols, slices] = size(data);
    score_map = zeros(rows, cols, slices);
    
    for r = 1:rows
        for c = 1:cols
            for z = 1:slices
                if gx(r,c,z) == 0 && gy(r,c,z) == 0 && gz(r,c,z) == 0
                    continue;
                end
                
                grad_mag = sqrt(gx(r,c,z)^2 + gy(r,c,z)^2 + gz(r,c,z)^2);
                if grad_mag == 0, continue; end
                
                grad_unit_x = gx(r,c,z) / grad_mag;
                grad_unit_y = gy(r,c,z) / grad_mag;
                grad_unit_z = gz(r,c,z) / grad_mag;
                
                for d = 1:radius
                    center_x = c + d * grad_unit_x;
                    center_y = r + d * grad_unit_y;
                    center_z = z + d * grad_unit_z;
                    
                    if center_x >= 1 && center_x <= cols && center_y >= 1 && center_y <= rows && center_z >= 1 && center_z <= slices
                        center_idx_x = round(center_x);
                        center_idx_y = round(center_y);
                        center_idx_z = round(center_z);
                        
                        score_map(center_idx_y, center_idx_x, center_idx_z) = ...
                            score_map(center_idx_y, center_idx_x, center_idx_z) + grad_mag;
                    end
                end
            end
        end
    end
end

function showSlicePlanes(X, Y, Z, score_map, center_phys, color_map)
    % Show orthogonal slice planes through the center
    [rows, cols, slices] = size(score_map);
    
    % Find center indices
    center_y_idx = round(center_phys(1) / (Y(2,1,1) - Y(1,1,1)));
    center_x_idx = round(center_phys(2) / (X(1,2,1) - X(1,1,1)));
    center_z_idx = round(center_phys(3) / (Z(1,1,2) - Z(1,1,1)));
    
    center_y_idx = max(1, min(rows, center_y_idx));
    center_x_idx = max(1, min(cols, center_x_idx));
    center_z_idx = max(1, min(slices, center_z_idx));
    
    % XY plane
    if slices > 1
        slice(X, Y, Z, score_map, [], [], center_phys(3));
    end
    
    % XZ plane
    if rows > 1
        slice(X, Y, Z, score_map, [], center_phys(1), []);
    end
    
    % YZ plane
    if cols > 1
        slice(X, Y, Z, score_map, center_phys(2), [], []);
    end
    
    colormap(color_map);
    shading interp;
    alpha(0.6);
end

function addConvergenceVectors(x_sub, y_sub, z_sub, gx_sub, gy_sub, gz_sub, grad_mag_sub, center_phys, strong_gradients, radius)
    % Add vectors showing convergence to center from selected points
    
    % Select a few representative points
    indices = find(strong_gradients);
    if length(indices) > 20
        % Randomly sample 20 points for clarity
        sample_idx = indices(randperm(length(indices), 20));
    else
        sample_idx = indices;
    end
    
    for i = 1:length(sample_idx)
        idx = sample_idx(i);
        
        % Current position
        pos = [x_sub(idx), y_sub(idx), z_sub(idx)];
        
        % Gradient direction (unit vector)
        grad_mag = sqrt(gx_sub(idx)^2 + gy_sub(idx)^2 + gz_sub(idx)^2);
        if grad_mag > 0
            grad_unit = [gx_sub(idx), gy_sub(idx), gz_sub(idx)] / grad_mag;
            
            % Project along gradient for visualization
            end_pos = pos + radius * grad_unit;
            
            % Draw convergence line
            plot3([pos(1), end_pos(1)], [pos(2), end_pos(2)], [pos(3), end_pos(3)], ...
                  '--', 'Color', [0.8, 0.4, 0.2, 0.7], 'LineWidth', 1);
        end
    end
end

function plotRadialProfile(data, center, xy_spacing, z_spacing, grad_mag)
    % Plot radial intensity and gradient profiles
    
    [rows, cols, slices] = size(data);
    center_y = round(center(1));
    center_x = round(center(2));
    center_z = round(center(3));
    
    % Ensure center is within bounds
    center_y = max(1, min(rows, center_y));
    center_x = max(1, min(cols, center_x));
    center_z = max(1, min(slices, center_z));
    
    % Create radial distance map
    [X, Y, Z] = meshgrid(1:cols, 1:rows, 1:slices);
    
    % Physical distances
    dx = (X - center_x) * xy_spacing;
    dy = (Y - center_y) * xy_spacing;
    dz = (Z - center_z) * z_spacing;
    
    radial_dist = sqrt(dx.^2 + dy.^2 + dz.^2);
    
    % Create radial bins
    max_dist = min([center_x, center_y, center_z, cols-center_x, rows-center_y, slices-center_z]) * xy_spacing;
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
