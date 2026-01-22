function fitResults = fitGaussians(imageData, maximaCoords, fitParams, is3D)
%fitGaussians Performs Gaussian fitting on local maxima.
%   This function takes image data and a list of local maxima and applies
%   a Gaussian fitting workflow based on the provided parameters.
%
% COORDINATE CONVENTION: This function uses ARRAY CONVENTION
%   Input coordinates are [row, col, slice]
%   - row: 1st dimension (vertical index)
%   - col: 2nd dimension (horizontal index)
%   - slice: 3rd dimension (depth index)
%
% Array access: imageData(row, col, slice) - DIRECT indexing
% Fit results store LOCAL coordinates within the fitting window (also [row, col, slice])
%
% Inputs:
%   imageData    - Image array (2D or 3D)
%   maximaCoords - Nx3 array of [row, col, slice] (ARRAY CONVENTION)
%   fitParams    - Structure with fitting parameters
%   is3D         - Boolean indicating 3D data
%
% Outputs:
%   fitResults   - Structure array with fit parameters and coordinate info



    % Pre-allocate results structure
    numMaxima = size(maximaCoords, 1);
    fitResults = struct('amplitude', cell(1, numMaxima), ...
                        'center_x', cell(1, numMaxima), 'center_y', cell(1, numMaxima), 'center_z', cell(1, numMaxima), ...
                        'sigma_x', cell(1, numMaxima), 'sigma_y', cell(1, numMaxima), 'sigma_z', cell(1, numMaxima), ...
                        'rho_xy', cell(1, numMaxima), 'rho_xz', cell(1, numMaxima), 'rho_yz', cell(1, numMaxima), ...
                        'alpha_x', cell(1, numMaxima), 'alpha_y', cell(1, numMaxima), 'alpha_z', cell(1, numMaxima), ...
                        'background', cell(1, numMaxima), 'integratedIntensity', cell(1, numMaxima), 'r_squared', cell(1, numMaxima), ...
                        'originalMaximaCoords', cell(1, numMaxima), 'rawDataWindow', cell(1, numMaxima), 'fitMethod', cell(1, numMaxima), ...
                        'fitWindowDimensions', cell(1, numMaxima), 'fitWindowOrigin', cell(1, numMaxima), ...
                        'localMaximaInWindow', cell(1, numMaxima), 'globalFitCenter', cell(1, numMaxima));

imgSize = size(imageData);

for i = 1:numMaxima
    % Initialize all results to NaN to handle cases where a field is not calculated
    fitResults(i).amplitude = NaN;
    fitResults(i).amplitude_x = NaN; fitResults(i).amplitude_y = NaN; fitResults(i).amplitude_z = NaN;
    fitResults(i).amplitude_xy = NaN;  % For 2D+1D fits
    fitResults(i).center_x = NaN; fitResults(i).center_y = NaN; fitResults(i).center_z = NaN;
    fitResults(i).sigma_x = NaN; fitResults(i).sigma_y = NaN; fitResults(i).sigma_z = NaN;
    fitResults(i).rho_xy = NaN; fitResults(i).rho_xz = NaN; fitResults(i).rho_yz = NaN;
    fitResults(i).alpha_x = NaN; fitResults(i).alpha_y = NaN; fitResults(i).alpha_z = NaN;
    fitResults(i).background = NaN; fitResults(i).integratedIntensity = NaN; fitResults(i).r_squared = NaN;
    fitResults(i).radialSymmetryScore = NaN;  % For radial symmetry methods
    fitResults(i).originalMaximaCoords = maximaCoords(i, :);
    fitResults(i).fitMethod = fitParams.gaussFitMethod;
    % Initialize background fields
    fitResults(i).background_method = '';
    fitResults(i).background_coeffs = [];
    fitResults(i).background_degree = NaN;
    fitResults(i).background_roi_width = fitParams.gaussFitBgCorrWidth;
    fitResults(i).background_roi_bounds = [];  % Will store window bounds in global coords

    % --- 1. Extract Voxel/Pixel Window ---
    % centerCoord is [row, col, slice] (ARRAY CONVENTION)
    centerCoord = round(maximaCoords(i, :));
    winRadius = (fitParams.gaussFitVoxelWindowSize - 1) / 2;

    % Define window boundaries, clamping to image dimensions
    % Use ARRAY CONVENTION: row, col, slice
    row_min = max(1, round(centerCoord(1) - winRadius));
    row_max = min(imgSize(1), round(centerCoord(1) + winRadius));
    col_min = max(1, round(centerCoord(2) - winRadius));
    col_max = min(imgSize(2), round(centerCoord(2) + winRadius));
    
    if is3D
        slice_min = max(1, round(centerCoord(3) - winRadius));
        slice_max = min(imgSize(3), round(centerCoord(3) + winRadius));
        % DIRECT array access with ARRAY CONVENTION
        windowData = imageData(row_min:row_max, col_min:col_max, slice_min:slice_max);
        
        % Store window information (ARRAY CONVENTION)
        fitResults(i).fitWindowDimensions = size(windowData);  % [height, width, depth]
        fitResults(i).fitWindowOrigin = [row_min, col_min, slice_min];  % [row, col, slice] of origin
        % Local coordinates of maxima within the window [row, col, slice]
        fitResults(i).localMaximaInWindow = [centerCoord(1) - row_min + 1, centerCoord(2) - col_min + 1, centerCoord(3) - slice_min + 1];
        % Store window bounds in global image coordinates (for background ROI reference)
        fitResults(i).background_roi_bounds = [row_min, row_max; col_min, col_max; slice_min, slice_max];
    else
        % DIRECT array access with ARRAY CONVENTION
        windowData = imageData(row_min:row_max, col_min:col_max);
        
        % Store window information (ARRAY CONVENTION)
        fitResults(i).fitWindowDimensions = size(windowData);  % [height, width]
        fitResults(i).fitWindowOrigin = [row_min, col_min];  % [row, col] of origin
        % Local coordinates of maxima within the window [row, col]
        fitResults(i).localMaximaInWindow = [centerCoord(1) - row_min + 1, centerCoord(2) - col_min + 1];
        % Store window bounds in global image coordinates (for background ROI reference)
        fitResults(i).background_roi_bounds = [row_min, row_max; col_min, col_max];
    end

    % --- 2. Perform Local Background Correction ---
    bgWidth = fitParams.gaussFitBgCorrWidth;
    
    % Define larger surrounding window for background calculation
    % Use ARRAY CONVENTION: row, col, slice
    bg_row_min = max(1, round(centerCoord(1) - winRadius - bgWidth));
    bg_row_max = min(imgSize(1), round(centerCoord(1) + winRadius + bgWidth));
    bg_col_min = max(1, round(centerCoord(2) - winRadius - bgWidth));
    bg_col_max = min(imgSize(2), round(centerCoord(2) + winRadius + bgWidth));

    if is3D
        bg_slice_min = max(1, round(centerCoord(3) - winRadius - bgWidth));
        bg_slice_max = min(imgSize(3), round(centerCoord(3) + winRadius + bgWidth));
        % DIRECT array access with ARRAY CONVENTION
        surroundingData = imageData(bg_row_min:bg_row_max, bg_col_min:bg_col_max, bg_slice_min:bg_slice_max);
    else
        % DIRECT array access with ARRAY CONVENTION
        surroundingData = imageData(bg_row_min:bg_row_max, bg_col_min:bg_col_max);
    end
    
    % Create a mask to exclude the central window from background calculation
    mask = true(size(surroundingData));
    mask_col_start = col_min - bg_col_min + 1;
    mask_col_end = mask_col_start + (col_max - col_min);
    mask_row_start = row_min - bg_row_min + 1;
    mask_row_end = mask_row_start + (row_max - row_min);

    if is3D
        mask_slice_start = slice_min - bg_slice_min + 1;
        mask_slice_end = mask_slice_start + (slice_max - slice_min);
        mask(round(mask_row_start):round(mask_row_end), round(mask_col_start):round(mask_col_end), round(mask_slice_start):round(mask_slice_end)) = false;
    else
        mask(round(mask_row_start):round(mask_row_end), round(mask_col_start):round(mask_col_end)) = false;
    end
    
    backgroundCorrectedWindow = windowData; % Initialize
    
    switch fitParams.gaussFitBgCorrMethod
        case 'Mean Surrounding Subtraction'
            backgroundVoxels = surroundingData(mask);
            meanBackground = mean(backgroundVoxels, 'all');
            backgroundCorrectedWindow = double(windowData) - meanBackground;
            % Store mean background
            fitResults(i).background = meanBackground;
            fitResults(i).background_method = 'Mean Surrounding Subtraction';
            
        case {'Local Plane Fitting', 'Local Polynomial Fitting'}
            % For plane or polynomial fitting, we solve a linear system.
            % First, get coordinates and values of the background voxels.
            [bg_y, bg_x, bg_z] = ind2sub(size(surroundingData), find(mask));
            bg_vals = double(surroundingData(mask));
            
            % Construct the design matrix A based on the model
            if is3D
                if strcmp(fitParams.gaussFitBgCorrMethod, 'Local Plane Fitting')
                    A = [bg_x, bg_y, bg_z, ones(size(bg_x))];
                else % Polynomial
                    deg = fitParams.gaussFitPolyDegree;
                    A = [bg_x, bg_y, bg_z, bg_x.^2, bg_y.^2, bg_z.^2, bg_x.*bg_y, bg_x.*bg_z, bg_y.*bg_z];
                    if deg == 3
                        A = [A, bg_x.^3, bg_y.^3, bg_z.^3, bg_x.^2.*bg_y, bg_x.*bg_y.^2, bg_x.^2.*bg_z, bg_x.*bg_z.^2, bg_y.^2.*bg_z, bg_y.*bg_z.^2, bg_x.*bg_y.*bg_z];
                    end
                    A = [A, ones(size(bg_x))];
                end
            else % 2D
                if strcmp(fitParams.gaussFitBgCorrMethod, 'Local Plane Fitting')
                    A = [bg_x, bg_y, ones(size(bg_x))];
                else % Polynomial
                    deg = fitParams.gaussFitPolyDegree;
                    A = [bg_x, bg_y, bg_x.^2, bg_y.^2, bg_x.*bg_y];
                    if deg == 3
                        A = [A, bg_x.^3, bg_y.^3, bg_x.^2.*bg_y, bg_x.*bg_y.^2];
                    end
                    A = [A, ones(size(bg_x))];
                end
            end
            
            % Solve for model coefficients
            coeffs = A \ bg_vals;
            
            % Now, calculate the background within the actual window
            [win_y, win_x, win_z] = ind2sub(size(windowData), 1:numel(windowData));
            win_y = win_y'; win_x = win_x'; win_z = win_z'; % Ensure column vectors
            
            % Build design matrix for window
            if is3D
                if strcmp(fitParams.gaussFitBgCorrMethod, 'Local Plane Fitting')
                    A_win = [win_x, win_y, win_z, ones(size(win_x))];
                else
                    deg = fitParams.gaussFitPolyDegree;
                    A_win = [win_x, win_y, win_z, win_x.^2, win_y.^2, win_z.^2, win_x.*win_y, win_x.*win_z, win_y.*win_z];
                    if deg == 3
                        A_win = [A_win, win_x.^3, win_y.^3, win_z.^3, win_x.^2.*win_y, win_x.*win_y.^2, win_x.^2.*win_z, win_x.*win_z.^2, win_y.^2.*win_z, win_y.*win_z.^2, win_x.*win_y.*win_z];
                    end
                    A_win = [A_win, ones(size(win_x))];
                end
            else % 2D
                if strcmp(fitParams.gaussFitBgCorrMethod, 'Local Plane Fitting')
                    A_win = [win_x, win_y, ones(size(win_x))];
                else
                    deg = fitParams.gaussFitPolyDegree;
                    A_win = [win_x, win_y, win_x.^2, win_y.^2, win_x.*win_y];
                    if deg == 3
                        A_win = [A_win, win_x.^3, win_y.^3, win_x.^2.*win_y, win_x.*win_y.^2];
                    end
                    A_win = [A_win, ones(size(win_x))];
                end
            end
            
            % Predict background at all window locations
            fitted_bg = A_win * coeffs;
            fitted_bg_reshaped = reshape(fitted_bg, size(windowData));
            
            % Check for NaN or Inf
            if any(isnan(fitted_bg(:))) || any(isinf(fitted_bg(:)))
                warning('Polynomial/plane fitting produced invalid values (NaN/Inf). Falling back to mean background subtraction.');
                meanBackground = mean(double(surroundingData(mask)), 'all');
                backgroundCorrectedWindow = double(windowData) - meanBackground;
                % Store mean background with fallback indicator
                fitResults(i).background = meanBackground;
                fitResults(i).background_method = 'Mean (Fallback)';
            else
                backgroundCorrectedWindow = double(windowData) - fitted_bg_reshaped;
                % Store polynomial/plane background information
                fitResults(i).background = mean(fitted_bg, 'all');  % Mean for reference
                fitResults(i).background_method = fitParams.gaussFitBgCorrMethod;
                fitResults(i).background_coeffs = coeffs;
                if strcmp(fitParams.gaussFitBgCorrMethod, 'Local Polynomial Fitting')
                    fitResults(i).background_degree = fitParams.gaussFitPolyDegree;
                else
                    fitResults(i).background_degree = 1;  % Plane is degree 1
                end
            end
    end
    
    % Calculate integrated intensity after background correction
    fitResults(i).integratedIntensity = sum(backgroundCorrectedWindow, 'all');

    % --- 3. Fit the Gaussian Model ---
    % Calculate the actual position of the maxima within the extracted window
    % This accounts for boundary clipping and ensures correct initial guesses
    % ARRAY CONVENTION: center_x=row, center_y=col, center_z=slice within window
    if is3D
        actual_center_row_in_window = centerCoord(1) - row_min + 1;
        actual_center_col_in_window = centerCoord(2) - col_min + 1;
        actual_center_slice_in_window = centerCoord(3) - slice_min + 1;
        % NOTE: window_info uses legacy naming (center_x, center_y, center_z) for compatibility with fitting functions
        % but these actually represent [row, col, slice] in ARRAY CONVENTION
        window_info = struct('center_x', actual_center_row_in_window, 'center_y', actual_center_col_in_window, 'center_z', actual_center_slice_in_window);
    else
        actual_center_row_in_window = centerCoord(1) - row_min + 1;
        actual_center_col_in_window = centerCoord(2) - col_min + 1;
        % NOTE: window_info uses legacy naming (center_x, center_y) for compatibility with fitting functions
        % but these actually represent [row, col] in ARRAY CONVENTION
        window_info = struct('center_x', actual_center_row_in_window, 'center_y', actual_center_col_in_window);
    end
    
    switch fitParams.gaussFitMethod
        case {'1D (X,Y,Z)', '1D (X,Y,Z) Gaussian'}
            [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z, r_squared_vals, amplitude_x, amplitude_y, amplitude_z, x_profile, y_profile, z_profile] = fit_1D_gaussians(backgroundCorrectedWindow, fitParams, is3D, window_info);
            fitResults(i).amplitude = amplitude;
            fitResults(i).center_x = center_x;
            fitResults(i).center_y = center_y;
            fitResults(i).center_z = center_z;
            fitResults(i).sigma_x = sigma_x;
            fitResults(i).sigma_y = sigma_y;
            fitResults(i).sigma_z = sigma_z;
            fitResults(i).r_squared = r_squared_vals;
            % Store individual amplitudes for 1D visualization
            fitResults(i).amplitude_x = amplitude_x;
            fitResults(i).amplitude_y = amplitude_y;
            fitResults(i).amplitude_z = amplitude_z;
            % Store the actual fitted profiles for visualization consistency
            fitResults(i).x_profile = x_profile;
            fitResults(i).y_profile = y_profile;
            fitResults(i).z_profile = z_profile;
            
        case {'2D (XY) + 1D (Z)', '2D (XY) + 1D (Z) Gaussian'}
            if is3D
                [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z, r_squared_vals, amplitude_x, amplitude_y, amplitude_z, z_profile, xy_slice] = fit_2D_plus_1D_gaussians(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = amplitude;  % 2D XY amplitude
                fitResults(i).amplitude_xy = amplitude_x;  % XY amplitude from 2D fit
                fitResults(i).amplitude_z = amplitude_z;  % 1D Z amplitude
                fitResults(i).center_x = center_x;
                fitResults(i).center_y = center_y;
                fitResults(i).center_z = center_z;
                fitResults(i).sigma_x = sigma_x;
                fitResults(i).sigma_y = sigma_y;
                fitResults(i).sigma_z = sigma_z;
                fitResults(i).r_squared = r_squared_vals;
                % Store the actual fitted data for visualization consistency
                fitResults(i).z_profile = z_profile;
                fitResults(i).xy_slice = xy_slice;
            else
                % For 2D data, perform a single 2D Gaussian fit
                [fit_xy, r_squared_xy] = fit_2D_distorted(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = fit_xy(1);
                fitResults(i).center_x = fit_xy(2);
                fitResults(i).center_y = fit_xy(3);
                fitResults(i).sigma_x = fit_xy(4);
                fitResults(i).sigma_y = fit_xy(5);
                fitResults(i).center_z = NaN;
                fitResults(i).sigma_z = NaN;
                fitResults(i).r_squared = r_squared_xy;
            end
            
        case 'Distorted 3D Gaussian'
            if is3D
                [fit_params, r_squared] = fit_3D_distorted(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = fit_params(1);
                fitResults(i).center_x = fit_params(2);
                fitResults(i).center_y = fit_params(3);
                fitResults(i).center_z = fit_params(4);
                fitResults(i).sigma_x = fit_params(5);
                fitResults(i).sigma_y = fit_params(6);
                fitResults(i).sigma_z = fit_params(7);
                fitResults(i).rho_xy = fit_params(8);
                fitResults(i).rho_xz = fit_params(9);
                fitResults(i).rho_yz = fit_params(10);
                fitResults(i).r_squared = r_squared;
                % Store the actual fitted data for visualization consistency
                fitResults(i).fittedDataWindow = backgroundCorrectedWindow;
            else
                warning('Distorted 3D Gaussian not applicable for 2D data. Using 2D distorted fit.');
                [fit_2d, r_squared_2d] = fit_2D_distorted(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = fit_2d(1);
                fitResults(i).center_x = fit_2d(2);
                fitResults(i).center_y = fit_2d(3);
                fitResults(i).sigma_x = fit_2d(4);
                fitResults(i).sigma_y = fit_2d(5);
                fitResults(i).rho_xy = fit_2d(6);
                fitResults(i).r_squared = r_squared_2d;
            end

        case '3D Gaussian'
            if is3D
                [fit_params, r_squared] = fit_3D(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = fit_params(1);
                fitResults(i).center_x = fit_params(2);
                fitResults(i).center_y = fit_params(3);
                fitResults(i).center_z = fit_params(4);
                fitResults(i).sigma_x = fit_params(5);
                fitResults(i).sigma_y = fit_params(6);
                fitResults(i).sigma_z = fit_params(7);
                fitResults(i).r_squared = r_squared;
                % Store the actual fitted data for visualization consistency
                fitResults(i).fittedDataWindow = backgroundCorrectedWindow;
            else
                warning('3D Gaussian not applicable for 2D data. Using 2D fit.');
                [fit_2d, r_squared_2d] = fit_2D(backgroundCorrectedWindow, fitParams, window_info);
                fitResults(i).amplitude = fit_2d(1);
                fitResults(i).center_x = fit_2d(2);
                fitResults(i).center_y = fit_2d(3);
                fitResults(i).sigma_x = fit_2d(4);
                fitResults(i).sigma_y = fit_2d(5);
                fitResults(i).r_squared = r_squared_2d;
            end
            
        case 'Radial Symmetry'
            [fit_params, r_squared, quality_metrics] = fit_radial_symmetry(backgroundCorrectedWindow, fitParams, is3D, window_info);
            fitResults(i).amplitude = NaN;  % No amplitude for radial symmetry
            fitResults(i).center_x = fit_params(1);
            fitResults(i).center_y = fit_params(2);
            if is3D
                fitResults(i).center_z = fit_params(3);
                fitResults(i).sigma_x = fit_params(4);
                fitResults(i).sigma_y = fit_params(5);
                fitResults(i).sigma_z = fit_params(6);
            else
                fitResults(i).center_z = NaN;
                fitResults(i).sigma_x = fit_params(3);
                fitResults(i).sigma_y = fit_params(4);
                fitResults(i).sigma_z = NaN;
            end
            fitResults(i).r_squared = NaN;  % No R² for radial symmetry
            fitResults(i).radialSymmetryScore = r_squared;  % Store raw score
            fitResults(i).radialSymmetryQuality = quality_metrics.normalized_quality;  % Store quality
            % Store the actual fitted data for visualization consistency
            fitResults(i).fittedDataWindow = backgroundCorrectedWindow;
            
        otherwise
            warning('Unknown fitting method: %s. Using 1D fitting as fallback.', fitParams.gaussFitMethod);
            [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z, r_squared_vals, amplitude_x, amplitude_y, amplitude_z, x_profile, y_profile, z_profile] = fit_1D_gaussians(backgroundCorrectedWindow, fitParams, is3D, window_info);
            fitResults(i).amplitude = amplitude;
            fitResults(i).center_x = center_x;
            fitResults(i).center_y = center_y;
            fitResults(i).center_z = center_z;
            fitResults(i).sigma_x = sigma_x;
            fitResults(i).sigma_y = sigma_y;
            fitResults(i).sigma_z = sigma_z;
            fitResults(i).r_squared = r_squared_vals;
            % Store individual amplitudes for 1D visualization
            fitResults(i).amplitude_x = amplitude_x;
            fitResults(i).amplitude_y = amplitude_y;
            fitResults(i).amplitude_z = amplitude_z;
            % Store the actual fitted profiles for visualization consistency
            fitResults(i).x_profile = x_profile;
            fitResults(i).y_profile = y_profile;
            fitResults(i).z_profile = z_profile;
    end
    
    % --- 4. Calculate Global Sub-Pixel Fit Center Coordinates ---
    % IMPORTANT: Fit results store center_x, center_y, center_z in LOCAL WINDOW coordinates
    % These are in ARRAY CONVENTION within the window: center_x=row, center_y=col, center_z=slice
    % Convert local fit center (within window) to global image coordinates
    if ~isnan(fitResults(i).center_x) && ~isnan(fitResults(i).center_y)
        if is3D && ~isnan(fitResults(i).center_z)
            % 3D: Global = window origin + [local_row, local_col, local_slice] - 1
            % Result is [global_row, global_col, global_slice] in ARRAY CONVENTION
            fitResults(i).globalFitCenter = fitResults(i).fitWindowOrigin + ...
                                            [fitResults(i).center_x, fitResults(i).center_y, fitResults(i).center_z] - 1;
        else
            % 2D: Global = window origin + [local_row, local_col] - 1
            % Result is [global_row, global_col, NaN] in ARRAY CONVENTION
            fitResults(i).globalFitCenter = [fitResults(i).fitWindowOrigin(1:2) + ...
                                             [fitResults(i).center_x, fitResults(i).center_y] - 1, NaN];
        end
    else
        % If fit failed, use original maxima coordinates (already in ARRAY CONVENTION)
        fitResults(i).globalFitCenter = maximaCoords(i, :);
    end
    
end  % End of maxima loop

% --- 5. Optional Plotting ---
if fitParams.gaussFitPlotCheck
    % Placeholder for plotting logic
    warning('Plotting is not yet implemented.');
end

end


% --- Helper function for 1D Gaussian Fitting ---
function [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z, r_squared_vals, amplitude_x, amplitude_y, amplitude_z, x_profile, y_profile, z_profile] = fit_1D_gaussians(windowData, fitParams, is3D, window_info)

    winSize = size(windowData);
    
    
    % --- Fit X dimension ---
    % Extract 1D slice along X axis through the center Y,Z coordinates
    center_y_idx = round(window_info.center_y);
    if is3D
        center_z_idx = round(window_info.center_z);
        x_profile = squeeze(windowData(center_y_idx, :, center_z_idx)); % 1D slice along X
    else
        x_profile = squeeze(windowData(center_y_idx, :)); % 1D slice along X
    end
    x_coords = (1:length(x_profile))';
    
    try
        [fit_x, r_squared_x] = fit_1D(x_profile(:), x_coords, fitParams, window_info.center_x);
    catch ME
        warning('X dimension fitting failed: %s', ME.message);
        fit_x = [NaN, NaN, NaN];
        r_squared_x = NaN;
    end
    
    % --- Fit Y dimension ---
    % Extract 1D slice along Y axis through the center X,Z coordinates
    center_x_idx = round(window_info.center_x);
    if is3D
        y_profile = squeeze(windowData(:, center_x_idx, center_z_idx)); % 1D slice along Y
    else
        y_profile = squeeze(windowData(:, center_x_idx)); % 1D slice along Y
    end
    y_coords = (1:length(y_profile))';
    
    try
        [fit_y, r_squared_y] = fit_1D(y_profile(:), y_coords, fitParams, window_info.center_y);
    catch ME
        warning('Y dimension fitting failed: %s', ME.message);
        fit_y = [NaN, NaN, NaN];
        r_squared_y = NaN;
    end
    
    % --- Fit Z dimension (if 3D) ---
    if is3D
        % Extract 1D slice along Z axis through the center X,Y coordinates
        z_profile = squeeze(windowData(center_y_idx, center_x_idx, :)); % 1D slice along Z
        z_coords = (1:length(z_profile))';
        
        try
            [fit_z, r_squared_z] = fit_1D(z_profile(:), z_coords, fitParams, window_info.center_z);
            center_z = fit_z(2);
            sigma_z = fit_z(3);
            r_squared_z_val = r_squared_z;
        catch ME
            warning('Z dimension fitting failed: %s', ME.message);
            center_z = NaN;
            sigma_z = NaN;
            r_squared_z_val = NaN;
        end
    else
        center_z = NaN;
        sigma_z = NaN;
        r_squared_z_val = NaN;
        z_profile = [];  % No Z profile for 2D data
    end
    
    % Use the brightest amplitude from the dimensional fits
    amplitude = max([fit_x(1), fit_y(1)]);
    center_x = fit_x(2);
    sigma_x = fit_x(3);
    center_y = fit_y(2);
    sigma_y = fit_y(3);
    
    % Store individual amplitudes for 1D visualization
    % Now using fitted amplitudes from 1D slices (not projections)
    amplitude_x = fit_x(1);
    amplitude_y = fit_y(1);
    if is3D
        amplitude_z = fit_z(1);
    else
        amplitude_z = NaN;
    end
    
    r_squared_vals = [r_squared_x, r_squared_y, r_squared_z_val];

end

% --- Helper function for 2D + 1D Gaussian Fitting ---
function [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z, r_squared_vals, amplitude_x, amplitude_y, amplitude_z, z_profile, xy_slice] = fit_2D_plus_1D_gaussians(windowData, fitParams, window_info)
    
    % --- Fit Z dimension (1D) - use a 1D SLICE through the center, not a sum ---
    center_x_idx = round(window_info.center_x);
    center_y_idx = round(window_info.center_y);
    center_x_idx = max(1, min(size(windowData, 2), center_x_idx));
    center_y_idx = max(1, min(size(windowData, 1), center_y_idx));
    
    % Extract 1D slice through the center coordinates
    z_profile = squeeze(windowData(center_y_idx, center_x_idx, :));
    z_coords = (1:length(z_profile))';
    [fit_z, r_squared_z] = fit_1D(z_profile(:), z_coords, fitParams, window_info.center_z);
    amplitude_z = fit_z(1);  % Store the Z amplitude separately
    center_z = fit_z(2);
    sigma_z = fit_z(3);

    % --- Fit XY dimension (2D) - use a 2D SLICE at the center Z, not a sum ---
    center_z_idx = round(window_info.center_z);
    center_z_idx = max(1, min(size(windowData, 3), center_z_idx));
    
    % Extract 2D slice at the center Z coordinate
    xy_slice = windowData(:, :, center_z_idx);
    
    % Create 2D window_info for XY slice
    xy_window_info = struct('center_x', window_info.center_x, 'center_y', window_info.center_y);
    [fit_xy, r_squared_xy] = fit_2D_distorted(xy_slice, fitParams, xy_window_info);
    
    amplitude = fit_xy(1);  % This is the 2D XY amplitude
    amplitude_x = amplitude;  % Store separately for consistency with 1D fits
    amplitude_y = amplitude;  % Store separately for consistency with 1D fits
    center_x = fit_xy(2);
    center_y = fit_xy(3);
    sigma_x = fit_xy(4);
    sigma_y = fit_xy(5);
    % theta = fit_xy(6); % Rotation angle is available if needed

    r_squared_vals = [r_squared_xy, r_squared_xy, r_squared_z];  % [R²_X, R²_Y, R²_Z]
end

function [fit_params, r_squared] = fit_3D(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    % data is indexed as (row, col, slice)
    % We need coords that match this indexing
    [rows, cols, slices] = size(data);
    [COL_grid, ROW_grid, SLICE_grid] = meshgrid(1:cols, 1:rows, 1:slices);
    % coords(:,1) = row indices, coords(:,2) = col indices, coords(:,3) = slice indices
    coords = [ROW_grid(:), COL_grid(:), SLICE_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    % Objective function for 3D Gaussian (no rotation)
    gaussian_model_3d = @(p, xyz) p(1) * exp(-((xyz(:,1)-p(2)).^2 / (2*p(5)^2) + ...
                                              (xyz(:,2)-p(3)).^2 / (2*p(6)^2) + ...
                                              (xyz(:,3)-p(4)).^2 / (2*p(7)^2)));
    objective = @(p, xyz, v) gaussian_model_3d(p, xyz) - v;

    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % Actual position of maxima in window
    center_y_guess = window_info.center_y;  % Actual position of maxima in window  
    center_z_guess = window_info.center_z;  % Actual position of maxima in window
    % ARRAY CONVENTION: center_x=row, center_y=col, center_z=slice
    amp_guess = data(round(center_x_guess), round(center_y_guess), round(center_z_guess));  % Direct access [row,col,slice]
    sigma_x_guess = size(data, 1) / 4;  % Quarter of window height (row dimension)
    sigma_y_guess = size(data, 2) / 4;  % Quarter of window width (col dimension)
    sigma_z_guess = size(data, 3) / 4;  % Quarter of window depth (slice dimension)
    
    % p = [amplitude, center_x, center_y, center_z, sigma_x, sigma_y, sigma_z]
    % where center_x=row, center_y=col, center_z=slice in ARRAY CONVENTION
    p0 = [amp_guess, center_x_guess, center_y_guess, center_z_guess, sigma_x_guess, sigma_y_guess, sigma_z_guess];
    
    % Bounds - amplitude from 0 to max, centers within window (1 to dimension size), sigma up to half window
    max_amplitude = max(data(:));
    if max_amplitude <= 0
        max_amplitude = 1;  % Fallback if all data is negative/zero
    end
    lb = [0.001, 1, 1, 1, 0.1, 0.1, 0.1];  % Amplitude > 0, centers at window edges, small sigma
    ub = [max_amplitude, size(data,1), size(data,2), size(data,3), size(data,1)/2, size(data,2)/2, size(data,3)/2];  % Max sigma = half window size
    
    % Solver options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    % Run solver
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    % Calculate R-squared
    v_fit = gaussian_model_3d(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - v_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function [fit_params, r_squared] = fit_2D(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    % data is indexed as (row, col)
    [rows, cols] = size(data);
    [COL_grid, ROW_grid] = meshgrid(1:cols, 1:rows);
    % coords(:,1) = row indices, coords(:,2) = col indices
    coords = [ROW_grid(:), COL_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    % Objective function for 2D Gaussian (with rotation)
    gaussian_model_2d = @(p, xy) p(1) * exp(-(((xy(:,1)-p(2))*cos(p(6)) + (xy(:,2)-p(3))*sin(p(6))).^2 / (2*p(4)^2) + ...
                                           ((xy(:,2)-p(3))*cos(p(6)) - (xy(:,1)-p(2))*sin(p(6))).^2 / (2*p(5)^2)));
    objective = @(p, xy, z) gaussian_model_2d(p, xy) - z;

    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % Actual position of maxima in window
    center_y_guess = window_info.center_y;  % Actual position of maxima in window
    amp_guess = data(round(center_y_guess), round(center_x_guess));  % Value at actual center
    sigma_x_guess = size(data, 2) / 4;  % Quarter of window width
    sigma_y_guess = size(data, 1) / 4;  % Quarter of window height
    theta_guess = 0;
    
    p0 = [amp_guess, center_x_guess, center_y_guess, sigma_x_guess, sigma_y_guess, theta_guess];
    
    % Bounds
    % Bounds - amplitude from 0 to max of background-corrected data
    max_amplitude = max(data(:));
    % Ensure max_amplitude is positive and reasonable
    if max_amplitude <= 0
        max_amplitude = 1;  % Fallback if all data is negative/zero
    end
    lb = [0.001, 1, 1, 0.1, 0.1, -pi];  % Amplitude > 0 (small positive), centers at window edges
    ub = [max_amplitude, size(data, 2), size(data, 1), size(data, 2)/2, size(data, 1)/2, pi];  % Max sigma = half window size
    
    % Solver options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    % Run solver
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    % Calculate R-squared
    z_fit = gaussian_model_2d(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - z_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function [fit_params, r_squared] = fit_1D(data, coords, fitParams, actual_center)
    % Objective function for 1D Gaussian (no vertical shift)
    gaussian_model = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2 * p(3)^2));
    objective = @(p, x, y) gaussian_model(p, x) - y;
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    % Initial guesses - use actual maxima position and window-based sigma
    center_guess = actual_center;  % Actual position of maxima in window
    
    % Fix amplitude guess: find the maximum value in the 1D profile
    amp_guess = max(data);  % Use peak of 1D profile, not coordinate indexing
    % Ensure amplitude guess is positive
    if amp_guess <= 0
        amp_guess = 0.1;  % Fallback positive value
    end
    
    sigma_guess = length(coords) / 4;  % Quarter of window length
    if sigma_guess < 0.1, sigma_guess = 0.1; end  % Minimum sigma
    
    p0 = [amp_guess, center_guess, sigma_guess];
    
    % Bounds
    % Bounds - amplitude from 0 to max of background-corrected data
    max_amplitude = max(data(:));
    % Ensure max_amplitude is positive and reasonable
    if max_amplitude <= 0
        max_amplitude = 1;  % Fallback if all data is negative/zero
    end
    lb = [0.001, min(coords), 0.1];  % Amplitude > 0 (small positive), centers at window edges
    ub = [max_amplitude, max(coords), length(coords)/2];  % Max sigma = half window size
    
    % Solver options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
    
    % Run solver
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    % Calculate R-squared
    y_fit = gaussian_model(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - y_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end


% --- Helper function for Distorted 3D Gaussian Fitting ---
function [fit_params, r_squared] = fit_3D_distorted(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    [rows, cols, slices] = size(data);
    [COL_grid, ROW_grid, SLICE_grid] = meshgrid(1:cols, 1:rows, 1:slices);
    % coords(:,1) = row indices, coords(:,2) = col indices, coords(:,3) = slice indices
    coords = [ROW_grid(:), COL_grid(:), SLICE_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    % Objective function for distorted 3D Gaussian
    objective = @(p, xyz, v) model_3d_distorted(p, xyz) - v;

    % p = [A, mu_x, mu_y, mu_z, sig_x, sig_y, sig_z, rho_xy, rho_xz, rho_yz]
    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % row in window (ARRAY CONVENTION)
    center_y_guess = window_info.center_y;  % col in window (ARRAY CONVENTION)
    center_z_guess = window_info.center_z;  % slice in window (ARRAY CONVENTION)
    % ARRAY CONVENTION: Direct access [row, col, slice]
    amp_guess = data(round(center_x_guess), round(center_y_guess), round(center_z_guess));
    sigma_x_guess = size(data, 1) / 4;  % Quarter of row dimension
    sigma_y_guess = size(data, 2) / 4;  % Quarter of col dimension
    sigma_z_guess = size(data, 3) / 4;  % Quarter of slice dimension
    p0 = [amp_guess, center_x_guess, center_y_guess, center_z_guess, sigma_x_guess, sigma_y_guess, sigma_z_guess, 0, 0, 0];
    
    % Bounds (ARRAY CONVENTION)
    max_amplitude = max(data(:));
    if max_amplitude <= 0
        max_amplitude = 1;
    end
    lb = [0.001, 1, 1, 1, 0.1, 0.1, 0.1, -1, -1, -1];
    ub = [max_amplitude, size(data,1), size(data,2), size(data,3), size(data,1)/2, size(data,2)/2, size(data,3)/2, 1, 1, 1];
    
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    % Calculate R-squared
    v_fit = model_3d_distorted(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - v_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function intensity = model_3d_distorted(p, xyz)
    % p = [A, mu_x, mu_y, mu_z, sig_x, sig_y, sig_z, rho_xy, rho_xz, rho_yz]
    A = p(1);
    mu = p(2:4);
    sig = p(5:7);
    rho = p(8:10);
    
    % Construct covariance matrix
    Sigma = [sig(1)^2,              rho(1)*sig(1)*sig(2),   rho(2)*sig(1)*sig(3);
             rho(1)*sig(1)*sig(2),   sig(2)^2,               rho(3)*sig(2)*sig(3);
             rho(2)*sig(1)*sig(3),   rho(3)*sig(2)*sig(3),   sig(3)^2];
    
    % Ensure covariance matrix is positive semi-definite
    [~, flag] = chol(Sigma);
    if flag ~= 0
        intensity = zeros(size(xyz, 1), 1); % Return zero if matrix is invalid
        return;
    end
    
    x_centered = xyz - mu;
    
    intensity = A * exp(-0.5 * sum((x_centered / Sigma) .* x_centered, 2));
end


% --- Helper function for Skewed 3D Gaussian Fitting ---
function [fit_params, r_squared] = fit_3D_skewed(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    [rows, cols, slices] = size(data);
    [COL_grid, ROW_grid, SLICE_grid] = meshgrid(1:cols, 1:rows, 1:slices);
    % coords(:,1) = row indices, coords(:,2) = col indices, coords(:,3) = slice indices
    coords = [ROW_grid(:), COL_grid(:), SLICE_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    objective = @(p, xyz, v) model_3d_skewed(p, xyz) - v;

    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % row in window (ARRAY CONVENTION)
    center_y_guess = window_info.center_y;  % col in window (ARRAY CONVENTION)
    center_z_guess = window_info.center_z;  % slice in window (ARRAY CONVENTION)
    % ARRAY CONVENTION: Direct access [row, col, slice]
    amp_guess = data(round(center_x_guess), round(center_y_guess), round(center_z_guess));
    sigma_x_guess = size(data, 1) / 4;  % Quarter of row dimension
    sigma_y_guess = size(data, 2) / 4;  % Quarter of col dimension
    sigma_z_guess = size(data, 3) / 4;  % Quarter of slice dimension
    
    % p = [A, mu_x..z, sig_x..z, rho_xy..yz, alpha_x..z]
    % where mu_x=row, mu_y=col, mu_z=slice in ARRAY CONVENTION
    p0 = [amp_guess, center_x_guess, center_y_guess, center_z_guess, sigma_x_guess, sigma_y_guess, sigma_z_guess, 0, 0, 0, 0, 0, 0];
    
    % Bounds (ARRAY CONVENTION)
    max_amplitude = max(data(:));
    if max_amplitude <= 0
        max_amplitude = 1;
    end
    lb = [0.001, 1, 1, 1, 0.1, 0.1, 0.1, -1, -1, -1, -inf, -inf, -inf];
    ub = [max_amplitude, size(data,1), size(data,2), size(data,3), size(data,1)/2, size(data,2)/2, size(data,3)/2, 1, 1, 1, inf, inf, inf];
    
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    v_fit = model_3d_skewed(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - v_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function intensity = model_3d_skewed(p, xyz)
    % p_distorted = [A, mu_x..z, sig_x..z, rho_xy..yz]
    p_distorted = p(1:10);
    alpha = p(11:13);
    mu = p(2:4);
    
    % Get the value of the un-skewed distorted Gaussian
    g_val = model_3d_distorted(p_distorted, xyz);
    
    % Calculate the skew term (CDF)
    x_centered = xyz - mu;
    skew_term = 0.5 * (1 + erf((x_centered * alpha') / sqrt(2)));
    
    intensity = 2 * g_val .* skew_term;
end

% --- Helper function for Skewed 2D Gaussian Fitting ---
function [fit_params, r_squared] = fit_2D_skewed(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    [rows, cols] = size(data);
    [COL_grid, ROW_grid] = meshgrid(1:cols, 1:rows);
    % coords(:,1) = row indices, coords(:,2) = col indices
    coords = [ROW_grid(:), COL_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    objective = @(p, xy, v) model_2d_skewed(p, xy) - v;

    % p = [A, mu_x, mu_y, sig_x, sig_y, rho_xy, alpha_x, alpha_y]
    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % row in window (ARRAY CONVENTION)
    center_y_guess = window_info.center_y;  % col in window (ARRAY CONVENTION)
    % ARRAY CONVENTION: Direct access [row, col]
    amp_guess = data(round(center_x_guess), round(center_y_guess));
    sigma_x_guess = size(data, 1) / 4;  % Quarter of row dimension
    sigma_y_guess = size(data, 2) / 4;  % Quarter of col dimension
    
    p0 = [amp_guess, center_x_guess, center_y_guess, sigma_x_guess, sigma_y_guess, 0, 0, 0];
    
    % Bounds (ARRAY CONVENTION)
    max_amplitude = max(data(:));
    if max_amplitude <= 0
        max_amplitude = 1;
    end
    lb = [0.001, 1, 1, 0.1, 0.1, -1, -inf, -inf];
    ub = [max_amplitude, size(data,1), size(data,2), size(data,1)/2, size(data,2)/2, 1, inf, inf];
    
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    v_fit = model_2d_skewed(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - v_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function intensity = model_2d_skewed(p, xy)
    p_distorted = p(1:6);
    alpha = p(7:8);
    mu = p(2:3);
    
    g_val = model_2d_distorted(p_distorted, xy);
    
    x_centered = xy - mu;
    skew_term = 0.5 * (1 + erf((x_centered * alpha') / sqrt(2)));
    
    intensity = 2 * g_val .* skew_term;
end

function [fit_params, r_squared] = fit_2D_distorted(data, fitParams, window_info)
    % ARRAY CONVENTION: Create coordinate grids for fitting
    [rows, cols] = size(data);
    [COL_grid, ROW_grid] = meshgrid(1:cols, 1:rows);
    % coords(:,1) = row indices, coords(:,2) = col indices
    coords = [ROW_grid(:), COL_grid(:)];
    
    % Ensure data types are double for LSQNONLIN
    data = double(data);
    coords = double(coords);
    
    % Objective function for 2D Gaussian (with rotation)
    objective = @(p, xy, z) model_2d_distorted(p, xy) - z;

    % Initial guesses - use actual maxima position and window-based sigma
    center_x_guess = window_info.center_x;  % row in window (ARRAY CONVENTION)
    center_y_guess = window_info.center_y;  % col in window (ARRAY CONVENTION)
    % ARRAY CONVENTION: Direct access [row, col]
    amp_guess = data(round(center_x_guess), round(center_y_guess));
    sigma_x_guess = size(data, 1) / 4;  % Quarter of row dimension
    sigma_y_guess = size(data, 2) / 4;  % Quarter of col dimension
    theta_guess = 0;
    
    % p = [A, mu_x, mu_y, sig_x, sig_y, theta]
    % where mu_x=row, mu_y=col, sig_x=row_sigma, sig_y=col_sigma (ARRAY CONVENTION)
    p0 = [amp_guess, center_x_guess, center_y_guess, sigma_x_guess, sigma_y_guess, theta_guess];
    
    % Bounds - amplitude from 0 to max of background-corrected data
    max_amplitude = max(data(:));
    if max_amplitude <= 0
        max_amplitude = 1;  % Fallback
    end
    lb = [0.001, 1, 1, 0.1, 0.1, -pi];
    ub = [max_amplitude, size(data,1), size(data,2), size(data,1)/2, size(data,2)/2, pi];
    
    % Solver options
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxIterations', fitParams.gaussFitMaxIterations, ...
                           'FunctionTolerance', fitParams.gaussFitTolerance);
                           
    % Run solver
    fit_params = lsqnonlin(@(p) objective(p, coords, data(:)), p0, lb, ub, options);
    
    % Calculate R-squared
    z_fit = model_2d_distorted(fit_params, coords);
    ss_total = sum((data(:) - mean(data(:))).^2);
    ss_residual = sum((data(:) - z_fit).^2);
    r_squared = 1 - (ss_residual / ss_total);
end

function intensity = model_2d_distorted(p, xy)
    % p = [A, mu_x, mu_y, sig_x, sig_y, theta]
    A = p(1);
    mu = p(2:3);
    sig = p(4:5);
    theta = p(6);
    
    x_centered = xy - mu;
    
    x_rot = x_centered(:,1)*cos(theta) + x_centered(:,2)*sin(theta);
    y_rot = -x_centered(:,1)*sin(theta) + x_centered(:,2)*cos(theta);
    
    intensity = A * exp(-((x_rot.^2 / (2*sig(1)^2)) + (y_rot.^2 / (2*sig(2)^2))));
end

% --- Helper function for Radial Symmetry Fitting ---
function [fit_params, r_squared, quality_metrics] = fit_radial_symmetry(data, fitParams, is3D, window_info)
    % Radial Symmetry method for spot detection
    % This method finds the center of radial symmetry by analyzing gradient vectors
    % Always performs symmetry detection regardless of score thresholds
    
    % Get parameters
    radius = fitParams.gaussFitRadialRadius;
    
    if is3D
        [center, score, quality_metrics] = find_radial_symmetry_3d(data, radius);
        
        % ARRAY CONVENTION: center is already [row, col, slice]
        % Return as [center_x, center_y, center_z] where center_x=row, center_y=col, center_z=slice
        % Sigma values are NaN for radial symmetry (not applicable)
        % No amplitude for radial symmetry
        fit_params = [center(1), center(2), center(3), NaN, NaN, NaN];
        
        fprintf('3D Radial symmetry center at [%.1f, %.1f, %.1f] (row,col,slice) with score %.2f (quality: %.3f)\n', ...
                center(1), center(2), center(3), score, quality_metrics.normalized_quality);
        

    else
        [center, score, quality_metrics] = find_radial_symmetry_2d(data, radius);
        
        % ARRAY CONVENTION: center is already [row, col]
        % Return as [center_x, center_y] where center_x=row, center_y=col
        % Sigma values are NaN for radial symmetry (not applicable)
        % No amplitude for radial symmetry
        fit_params = [center(1), center(2), NaN, NaN];
        
        fprintf('2D Radial symmetry center at [%.1f, %.1f] (row,col) with score %.2f (quality: %.3f)\n', ...
                center(1), center(2), score, quality_metrics.normalized_quality);
        

    end
    
    % Return the raw score (not normalized quality) as r_squared for compatibility
    r_squared = score;
end

function [center, score, quality_metrics] = find_radial_symmetry_2d(data, radius)
    % 2D Radial Symmetry Transform
    % Based on the method described in Parthasarathy, R. (2012)
    % Always finds a center regardless of score - no thresholding
    
    [rows, cols] = size(data);
    
    % Calculate gradients
    [gx, gy] = gradient(double(data));
    
    % Initialize score map
    score_map = zeros(rows, cols);
    gradient_count_map = zeros(rows, cols); % Track how many gradients contribute to each point
    
    % For each pixel, calculate radial symmetry score
    total_gradients = 0;
    contributing_gradients = 0;
    
    for r = 1:rows
        for c = 1:cols
            if gx(r,c) == 0 && gy(r,c) == 0
                continue; % Skip zero gradient points
            end
            
            total_gradients = total_gradients + 1;
            
            % Calculate gradient magnitude and direction
            grad_mag = sqrt(gx(r,c)^2 + gy(r,c)^2);
            if grad_mag == 0
                continue;
            end
            
            % Unit gradient vector (pointing towards increasing intensity)
            grad_unit_x = gx(r,c) / grad_mag;
            grad_unit_y = gy(r,c) / grad_mag;
            
            % Calculate potential center positions along gradient direction
            for d = 1:radius
                center_x = c + d * grad_unit_x;
                center_y = r + d * grad_unit_y;
                
                % Check if center is within image bounds
                if center_x >= 1 && center_x <= cols && center_y >= 1 && center_y <= rows
                    center_idx_x = round(center_x);
                    center_idx_y = round(center_y);
                    
                    % Add to score map (weighted by gradient magnitude)
                    score_map(center_idx_y, center_idx_x) = score_map(center_idx_y, center_idx_x) + grad_mag;
                    gradient_count_map(center_idx_y, center_idx_x) = gradient_count_map(center_idx_y, center_idx_x) + 1;
                    contributing_gradients = contributing_gradients + 1;
                end
            end
        end
    end
    
    % Always find the maximum score - no thresholding
    [max_score_val, max_idx] = max(score_map(:));
    [max_y, max_x] = ind2sub(size(score_map), max_idx);
    
    % Normalize the score: divide by the sum of all gradient magnitudes that contributed
    % This makes the score represent the fraction of gradient "votes" that converged
    total_gradient_magnitude = sum(sqrt(gx(:).^2 + gy(:).^2));
    if total_gradient_magnitude > 0
        normalized_score = max_score_val / total_gradient_magnitude;
    else
        normalized_score = 0;
    end
    
    % Refine to sub-pixel precision using parabolic interpolation
    [refined_y, refined_x] = refine_center_2d(score_map, max_y, max_x);
    
    center = [refined_y, refined_x];
    score = normalized_score;
    
    % Create simple quality metrics - the normalized score IS the quality
    quality_metrics = struct();
    quality_metrics.normalized_quality = normalized_score;
    quality_metrics.raw_score = max_score_val;
    quality_metrics.total_gradient_magnitude = total_gradient_magnitude;
    quality_metrics.convergence_ratio = contributing_gradients / max(1, total_gradients);
end

function [center, score, quality_metrics] = find_radial_symmetry_3d(data, radius)
    % 3D Radial Symmetry Transform
    % Always finds a center regardless of score - no thresholding
    
    [rows, cols, slices] = size(data);
    
    % Calculate gradients
    [gx, gy, gz] = gradient(double(data));
    
    % Initialize score map
    score_map = zeros(rows, cols, slices);
    gradient_count_map = zeros(rows, cols, slices); % Track how many gradients contribute to each point
    
    % For each voxel, calculate radial symmetry score
    total_gradients = 0;
    contributing_gradients = 0;
    
    for r = 1:rows
        for c = 1:cols
            for z = 1:slices
                if gx(r,c,z) == 0 && gy(r,c,z) == 0 && gz(r,c,z) == 0
                    continue; % Skip zero gradient points
                end
                
                total_gradients = total_gradients + 1;
                
                % Calculate gradient magnitude and direction
                grad_mag = sqrt(gx(r,c,z)^2 + gy(r,c,z)^2 + gz(r,c,z)^2);
                if grad_mag == 0
                    continue;
                end
                
                % Unit gradient vector (pointing towards increasing intensity)
                grad_unit_x = gx(r,c,z) / grad_mag;
                grad_unit_y = gy(r,c,z) / grad_mag;
                grad_unit_z = gz(r,c,z) / grad_mag;
                
                % Calculate potential center positions along gradient direction
                for d = 1:radius
                    center_x = c + d * grad_unit_x;
                    center_y = r + d * grad_unit_y;
                    center_z = z + d * grad_unit_z;
                    
                    % Check if center is within image bounds
                    if center_x >= 1 && center_x <= cols && center_y >= 1 && center_y <= rows && center_z >= 1 && center_z <= slices
                        center_idx_x = round(center_x);
                        center_idx_y = round(center_y);
                        center_idx_z = round(center_z);
                        
                        % Add to score map (weighted by gradient magnitude)
                        score_map(center_idx_y, center_idx_x, center_idx_z) = score_map(center_idx_y, center_idx_x, center_idx_z) + grad_mag;
                        gradient_count_map(center_idx_y, center_idx_x, center_idx_z) = gradient_count_map(center_idx_y, center_idx_x, center_idx_z) + 1;
                        contributing_gradients = contributing_gradients + 1;
                    end
                end
            end
        end
    end
    
    % Always find the maximum score - no thresholding
    [max_score_val, max_idx] = max(score_map(:));
    [max_y, max_x, max_z] = ind2sub(size(score_map), max_idx);
    
    % Normalize the score: divide by the sum of all gradient magnitudes that contributed
    % This makes the score represent the fraction of gradient "votes" that converged
    total_gradient_magnitude = sum(sqrt(gx(:).^2 + gy(:).^2 + gz(:).^2));
    if total_gradient_magnitude > 0
        normalized_score = max_score_val / total_gradient_magnitude;
    else
        normalized_score = 0;
    end
    
    % Refine to sub-pixel precision using parabolic interpolation
    [refined_y, refined_x, refined_z] = refine_center_3d(score_map, max_y, max_x, max_z);
    
    center = [refined_y, refined_x, refined_z];
    score = normalized_score;
    
    % Create simple quality metrics - the normalized score IS the quality
    quality_metrics = struct();
    quality_metrics.normalized_quality = normalized_score;
    quality_metrics.raw_score = max_score_val;
    quality_metrics.total_gradient_magnitude = total_gradient_magnitude;
    quality_metrics.convergence_ratio = contributing_gradients / max(1, total_gradients);
end

% --- Quality Metrics Calculation Functions ---

function quality_metrics = calculate_quality_metrics_2d(score_map, gradient_count_map, center, total_gradients, contributing_gradients, max_score)
    % Calculate comprehensive quality metrics for 2D radial symmetry detection
    
    [rows, cols] = size(score_map);
    center_y = center(1);
    center_x = center(2);
    
    % 1. Peak sharpness - how well-defined is the maximum?
    % Calculate ratio of max score to mean score in surrounding region
    window_size = 3;
    y_start = max(1, center_y - window_size);
    y_end = min(rows, center_y + window_size);
    x_start = max(1, center_x - window_size);
    x_end = min(cols, center_x + window_size);
    
    surrounding_region = score_map(y_start:y_end, x_start:x_end);
    surrounding_region(surrounding_region == max_score) = []; % Remove the peak itself
    
    if ~isempty(surrounding_region) && mean(surrounding_region) > 0
        peak_sharpness = max_score / mean(surrounding_region);
    else
        peak_sharpness = 1;
    end
    
    % 2. Gradient convergence - what fraction of gradients contribute to the detection?
    convergence_ratio = contributing_gradients / max(1, total_gradients);
    
    % 3. Spatial concentration - how concentrated are the contributing gradients around the center?
    gradient_concentration = gradient_count_map(center_y, center_x) / max(1, max(gradient_count_map(:)));
    
    % 4. Score relative to image intensity
    max_intensity = max(score_map(:));
    relative_score = max_score / max(1, max_intensity);
    
    % 5. Background uniformity - how uniform is the score map background?
    background_scores = score_map(score_map > 0 & score_map < max_score * 0.5);
    if length(background_scores) > 1
        background_uniformity = 1 / (1 + std(background_scores) / mean(background_scores));
    else
        background_uniformity = 0.5;
    end
    
    % Combine metrics into a single quality score (0 to 1)
    % Weight the different components
    weights = [0.3, 0.25, 0.2, 0.15, 0.1]; % peak_sharpness, convergence, concentration, relative_score, uniformity
    
    % Normalize individual metrics to 0-1 range
    normalized_sharpness = min(1, max(0, (peak_sharpness - 1) / (peak_sharpness + 9))); % Sharpness > 1 is good, asymptotic to 1
    normalized_convergence = convergence_ratio; % Already 0-1
    normalized_concentration = gradient_concentration; % Already 0-1
    normalized_relative = relative_score; % Already 0-1, don't scale up
    normalized_uniformity = background_uniformity; % Already 0-1
    
    normalized_quality = weights(1) * normalized_sharpness + ...
                        weights(2) * normalized_convergence + ...
                        weights(3) * normalized_concentration + ...
                        weights(4) * normalized_relative + ...
                        weights(5) * normalized_uniformity;
    
    % Store all metrics
    quality_metrics = struct();
    quality_metrics.normalized_quality = normalized_quality;
    quality_metrics.peak_sharpness = peak_sharpness;
    quality_metrics.convergence_ratio = convergence_ratio;
    quality_metrics.gradient_concentration = gradient_concentration;
    quality_metrics.relative_score = relative_score;
    quality_metrics.background_uniformity = background_uniformity;
end

% --- Sub-pixel Refinement Functions ---

function [refined_y, refined_x] = refine_center_2d(score_map, peak_y, peak_x)
    % Refine center location to sub-pixel precision using parabolic interpolation
    
    [rows, cols] = size(score_map);
    
    % Initialize with integer coordinates
    refined_y = peak_y;
    refined_x = peak_x;
    
    % X-direction refinement (columns)
    if peak_x > 1 && peak_x < cols
        % Get three points around the peak
        left = score_map(peak_y, peak_x - 1);
        center = score_map(peak_y, peak_x);
        right = score_map(peak_y, peak_x + 1);
        
        % Parabolic interpolation to find sub-pixel maximum
        denom = 2 * (left - 2*center + right);
        if abs(denom) > eps % Avoid division by zero
            dx = (left - right) / denom;
            % Clamp the offset to reasonable bounds
            dx = max(-0.5, min(0.5, dx));
            refined_x = peak_x + dx;
        end
    end
    
    % Y-direction refinement (rows)
    if peak_y > 1 && peak_y < rows
        % Get three points around the peak
        top = score_map(peak_y - 1, peak_x);
        center = score_map(peak_y, peak_x);
        bottom = score_map(peak_y + 1, peak_x);
        
        % Parabolic interpolation to find sub-pixel maximum
        denom = 2 * (top - 2*center + bottom);
        if abs(denom) > eps % Avoid division by zero
            dy = (top - bottom) / denom;
            % Clamp the offset to reasonable bounds
            dy = max(-0.5, min(0.5, dy));
            refined_y = peak_y + dy;
        end
    end
end

function [refined_y, refined_x, refined_z] = refine_center_3d(score_map, peak_y, peak_x, peak_z)
    % Refine center location to sub-pixel precision using parabolic interpolation
    
    [rows, cols, slices] = size(score_map);
    
    % Initialize with integer coordinates
    refined_y = peak_y;
    refined_x = peak_x;
    refined_z = peak_z;
    
    % X-direction refinement (columns)
    if peak_x > 1 && peak_x < cols
        % Get three points around the peak
        left = score_map(peak_y, peak_x - 1, peak_z);
        center = score_map(peak_y, peak_x, peak_z);
        right = score_map(peak_y, peak_x + 1, peak_z);
        
        % Parabolic interpolation to find sub-pixel maximum
        denom = 2 * (left - 2*center + right);
        if abs(denom) > eps % Avoid division by zero
            dx = (left - right) / denom;
            % Clamp the offset to reasonable bounds
            dx = max(-0.5, min(0.5, dx));
            refined_x = peak_x + dx;
        end
    end
    
    % Y-direction refinement (rows)
    if peak_y > 1 && peak_y < rows
        % Get three points around the peak
        top = score_map(peak_y - 1, peak_x, peak_z);
        center = score_map(peak_y, peak_x, peak_z);
        bottom = score_map(peak_y + 1, peak_x, peak_z);
        
        % Parabolic interpolation to find sub-pixel maximum
        denom = 2 * (top - 2*center + bottom);
        if abs(denom) > eps % Avoid division by zero
            dy = (top - bottom) / denom;
            % Clamp the offset to reasonable bounds
            dy = max(-0.5, min(0.5, dy));
            refined_y = peak_y + dy;
        end
    end
    
    % Z-direction refinement (slices)
    if peak_z > 1 && peak_z < slices
        % Get three points around the peak
        front = score_map(peak_y, peak_x, peak_z - 1);
        center = score_map(peak_y, peak_x, peak_z);
        back = score_map(peak_y, peak_x, peak_z + 1);
        
        % Parabolic interpolation to find sub-pixel maximum
        denom = 2 * (front - 2*center + back);
        if abs(denom) > eps % Avoid division by zero
            dz = (front - back) / denom;
            % Clamp the offset to reasonable bounds
            dz = max(-0.5, min(0.5, dz));
            refined_z = peak_z + dz;
        end
    end
end

function quality_metrics = calculate_quality_metrics_3d(score_map, gradient_count_map, center, total_gradients, contributing_gradients, max_score)
    % Calculate comprehensive quality metrics for 3D radial symmetry detection
    
    [rows, cols, slices] = size(score_map);
    center_y = center(1);
    center_x = center(2);
    center_z = center(3);
    
    % 1. Peak sharpness - how well-defined is the maximum?
    window_size = 2; % Smaller window for 3D to avoid edge effects
    y_start = max(1, center_y - window_size);
    y_end = min(rows, center_y + window_size);
    x_start = max(1, center_x - window_size);
    x_end = min(cols, center_x + window_size);
    z_start = max(1, center_z - window_size);
    z_end = min(slices, center_z + window_size);
    
    surrounding_region = score_map(y_start:y_end, x_start:x_end, z_start:z_end);
    surrounding_region(surrounding_region == max_score) = []; % Remove the peak itself
    
    if ~isempty(surrounding_region) && mean(surrounding_region) > 0
        peak_sharpness = max_score / mean(surrounding_region);
    else
        peak_sharpness = 1;
    end
    
    % 2. Gradient convergence - what fraction of gradients contribute to the detection?
    convergence_ratio = contributing_gradients / max(1, total_gradients);
    
    % 3. Spatial concentration - how concentrated are the contributing gradients around the center?
    gradient_concentration = gradient_count_map(center_y, center_x, center_z) / max(1, max(gradient_count_map(:)));
    
    % 4. Score relative to image intensity
    max_intensity = max(score_map(:));
    relative_score = max_score / max(1, max_intensity);
    
    % 5. Background uniformity - how uniform is the score map background?
    background_scores = score_map(score_map > 0 & score_map < max_score * 0.5);
    if length(background_scores) > 1
        background_uniformity = 1 / (1 + std(background_scores) / mean(background_scores));
    else
        background_uniformity = 0.5;
    end
    
    % Combine metrics into a single quality score (0 to 1)
    % Weight the different components
    weights = [0.3, 0.25, 0.2, 0.15, 0.1]; % peak_sharpness, convergence, concentration, relative_score, uniformity
    
    % Normalize individual metrics to 0-1 range
    normalized_sharpness = min(1, max(0, (peak_sharpness - 1) / (peak_sharpness + 9))); % Sharpness > 1 is good, asymptotic to 1
    normalized_convergence = convergence_ratio; % Already 0-1
    normalized_concentration = gradient_concentration; % Already 0-1
    normalized_relative = relative_score; % Already 0-1, don't scale up
    normalized_uniformity = background_uniformity; % Already 0-1
    
    normalized_quality = weights(1) * normalized_sharpness + ...
                        weights(2) * normalized_convergence + ...
                        weights(3) * normalized_concentration + ...
                        weights(4) * normalized_relative + ...
                        weights(5) * normalized_uniformity;
    
    % Store all metrics
    quality_metrics = struct();
    quality_metrics.normalized_quality = normalized_quality;
    quality_metrics.peak_sharpness = peak_sharpness;
    quality_metrics.convergence_ratio = convergence_ratio;
    quality_metrics.gradient_concentration = gradient_concentration;
    quality_metrics.relative_score = relative_score;
    quality_metrics.background_uniformity = background_uniformity;
end

