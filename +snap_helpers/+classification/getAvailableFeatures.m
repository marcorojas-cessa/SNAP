function [features, featureInfo] = getAvailableFeatures(fittingMethod, has3D, hasPhysicalSpacing)
% getAvailableFeatures - SINGLE SOURCE OF TRUTH for classification features
%
% ============================================================================
% CRITICAL: This is THE ONLY PLACE where classification features are defined
% ============================================================================
%
% This function dynamically discovers what features are available for SVM
% classification based on the current fitting method. It mirrors the logic
% in buildSignalColumnList.m but focuses on numeric features suitable for
% machine learning.
%
% KEY DESIGN PRINCIPLE:
%   Features are COMPLETELY CUSTOMIZABLE - this function returns ALL
%   possible features, and the user selects which to use via the UI.
%
% FEATURE CATEGORIES:
%   'position'  - Spatial coordinates (x, y, z)
%   'intensity' - Brightness measurements (amplitude, background, etc.)
%   'shape'     - Gaussian shape parameters (sigma_x, sigma_y, etc.)
%   'quality'   - Fit quality metrics (r_squared, radial_symmetry_score)
%   'derived'   - Computed from other features (SNR, ratios, etc.)
%   'physical'  - Physical units if spacing available (sigma_x_nm, etc.)
%
% USAGE:
%   [features, featureInfo] = snap_helpers.classification.getAvailableFeatures('3D Gaussian', true, false)
%
% INPUTS:
%   fittingMethod       - String: '1D (X,Y,Z)', '2D (XY) + 1D (Z)', '3D Gaussian',
%                         'Distorted 3D Gaussian', 'Radial Symmetry', etc.
%   has3D               - Logical: whether data is 3D (default: false)
%   hasPhysicalSpacing  - Logical: whether physical spacing is available (default: false)
%
% OUTPUTS:
%   features    - Cell array of feature names (strings)
%   featureInfo - Struct with metadata for each feature:
%                 .category    - Feature category string
%                 .description - Human-readable description
%                 .compute     - (optional) Function handle for derived features
%                                @(data) data.field1 ./ data.field2
%
% ============================================================================

    % Handle defaults
    if nargin < 2, has3D = false; end
    if nargin < 3, hasPhysicalSpacing = false; end
    
    % Initialize outputs
    features = {};
    featureInfo = struct();
    
    % Detect fitting method characteristics (mirrors buildSignalColumnList.m)
    is1DFit = contains(fittingMethod, '1D (X,Y,Z)', 'IgnoreCase', true) || ...
              contains(fittingMethod, '1D x,y,z', 'IgnoreCase', true);
    is2DPlus1DFit = contains(fittingMethod, '2D (XY) + 1D (Z)', 'IgnoreCase', true) || ...
                    contains(fittingMethod, '2D+1D', 'IgnoreCase', true);
    isRadialSymmetry = contains(fittingMethod, 'Radial Symmetry', 'IgnoreCase', true);
    isSkewed = contains(fittingMethod, 'Skewed', 'IgnoreCase', true);
    isDistorted = contains(fittingMethod, 'Distorted', 'IgnoreCase', true) && ~isSkewed;
    
    % ========================================================================
    % POSITION FEATURES (always available but often excluded from classification)
    % ========================================================================
    features{end+1} = 'x';
    featureInfo.x = struct('category', 'position', ...
        'description', 'X coordinate (pixels)', ...
        'source', 'fitted_coords');
    
    features{end+1} = 'y';
    featureInfo.y = struct('category', 'position', ...
        'description', 'Y coordinate (pixels)', ...
        'source', 'fitted_coords');
    
    if has3D
        features{end+1} = 'z';
        featureInfo.z = struct('category', 'position', ...
            'description', 'Z coordinate (slices)', ...
            'source', 'fitted_coords');
    end
    
    % ========================================================================
    % INTENSITY FEATURES (always available with fitting)
    % ========================================================================
    
    % Amplitude (method-dependent)
    if is1DFit
        features{end+1} = 'amplitude_x';
        featureInfo.amplitude_x = struct('category', 'intensity', ...
            'description', 'Fitted amplitude in X dimension');
        
        features{end+1} = 'amplitude_y';
        featureInfo.amplitude_y = struct('category', 'intensity', ...
            'description', 'Fitted amplitude in Y dimension');
        
        features{end+1} = 'amplitude_z';
        featureInfo.amplitude_z = struct('category', 'intensity', ...
            'description', 'Fitted amplitude in Z dimension');
            
    elseif is2DPlus1DFit
        features{end+1} = 'amplitude_xy';
        featureInfo.amplitude_xy = struct('category', 'intensity', ...
            'description', 'Fitted amplitude in XY plane (2D fit)');
        
        features{end+1} = 'amplitude_z';
        featureInfo.amplitude_z = struct('category', 'intensity', ...
            'description', 'Fitted amplitude in Z dimension (1D fit)');
            
    elseif ~isRadialSymmetry
        % Standard Gaussian, Distorted, or Skewed
        features{end+1} = 'amplitude';
        featureInfo.amplitude = struct('category', 'intensity', ...
            'description', 'Peak amplitude from Gaussian fit');
    end
    
    % Integrated intensity (always available with fitting)
    features{end+1} = 'integrated_intensity';
    featureInfo.integrated_intensity = struct('category', 'intensity', ...
        'description', 'Sum of background-corrected intensity in fit window');
    
    % Background (always available with fitting)
    features{end+1} = 'background';
    featureInfo.background = struct('category', 'intensity', ...
        'description', 'Estimated local background level');
    
    % ========================================================================
    % FIT QUALITY FEATURES
    % ========================================================================
    if isRadialSymmetry
        features{end+1} = 'radial_symmetry_score';
        featureInfo.radial_symmetry_score = struct('category', 'quality', ...
            'description', 'Radial symmetry convergence score (0-1)');
    else
        features{end+1} = 'r_squared';
        featureInfo.r_squared = struct('category', 'quality', ...
            'description', 'Goodness of fit (R² value, 0-1)');
    end
    
    % ========================================================================
    % SHAPE FEATURES (Gaussian methods only, not Radial Symmetry)
    % ========================================================================
    if ~isRadialSymmetry
        features{end+1} = 'sigma_x';
        featureInfo.sigma_x = struct('category', 'shape', ...
            'description', 'Gaussian width in X dimension (pixels)');
        
        features{end+1} = 'sigma_y';
        featureInfo.sigma_y = struct('category', 'shape', ...
            'description', 'Gaussian width in Y dimension (pixels)');
        
        if has3D
            features{end+1} = 'sigma_z';
            featureInfo.sigma_z = struct('category', 'shape', ...
                'description', 'Gaussian width in Z dimension (slices)');
        end
    end
    
    % ========================================================================
    % CORRELATION FEATURES (Distorted/Skewed 3D only)
    % ========================================================================
    if has3D && (isDistorted || isSkewed)
        features{end+1} = 'rho_xy';
        featureInfo.rho_xy = struct('category', 'shape', ...
            'description', 'XY correlation coefficient (-1 to 1)');
        
        features{end+1} = 'rho_xz';
        featureInfo.rho_xz = struct('category', 'shape', ...
            'description', 'XZ correlation coefficient (-1 to 1)');
        
        features{end+1} = 'rho_yz';
        featureInfo.rho_yz = struct('category', 'shape', ...
            'description', 'YZ correlation coefficient (-1 to 1)');
    end
    
    % ========================================================================
    % SKEWNESS FEATURES (Skewed method only)
    % ========================================================================
    if isSkewed
        features{end+1} = 'alpha_x';
        featureInfo.alpha_x = struct('category', 'shape', ...
            'description', 'X-axis skewness parameter');
        
        features{end+1} = 'alpha_y';
        featureInfo.alpha_y = struct('category', 'shape', ...
            'description', 'Y-axis skewness parameter');
        
        if has3D
            features{end+1} = 'alpha_z';
            featureInfo.alpha_z = struct('category', 'shape', ...
                'description', 'Z-axis skewness parameter');
        end
    end
    
    % ========================================================================
    % DERIVED FEATURES (computed from raw features)
    % These are ALWAYS available as they're computed on-the-fly
    % ========================================================================
    
    % Signal-to-Noise Ratio
    features{end+1} = 'snr';
    featureInfo.snr = struct('category', 'derived', ...
        'description', 'Signal-to-noise ratio (integrated_intensity / background)', ...
        'compute', @(d) computeSNR(d));
    
    % Amplitude over background (if amplitude available)
    if ~isRadialSymmetry && ~is1DFit && ~is2DPlus1DFit
        features{end+1} = 'amplitude_over_background';
        featureInfo.amplitude_over_background = struct('category', 'derived', ...
            'description', 'Amplitude / Background ratio', ...
            'compute', @(d) computeAmpOverBg(d));
    end
    
    % Sigma-based derived features (if sigmas available)
    if ~isRadialSymmetry
        % XY aspect ratio (ellipticity)
        features{end+1} = 'sigma_xy_ratio';
        featureInfo.sigma_xy_ratio = struct('category', 'derived', ...
            'description', 'sigma_x / sigma_y ratio (ellipticity, 1.0 = circular)', ...
            'compute', @(d) computeSigmaXYRatio(d));
        
        % Total sigma (spot size proxy)
        features{end+1} = 'sigma_sum';
        featureInfo.sigma_sum = struct('category', 'derived', ...
            'description', 'Sum of all sigma values (total spot size)', ...
            'compute', @(d) computeSigmaSum(d, has3D));
        
        % Sigma product (spot area/volume proxy)
        features{end+1} = 'sigma_product';
        featureInfo.sigma_product = struct('category', 'derived', ...
            'description', 'Product of sigma values (spot area/volume proxy)', ...
            'compute', @(d) computeSigmaProduct(d, has3D));
        
        if has3D
            % Z-ratio (how elongated in Z vs XY)
            features{end+1} = 'sigma_z_ratio';
            featureInfo.sigma_z_ratio = struct('category', 'derived', ...
                'description', 'sigma_z / mean(sigma_x, sigma_y) (Z elongation)', ...
                'compute', @(d) computeSigmaZRatio(d));
        end
    end
    
    % ========================================================================
    % PHYSICAL UNIT FEATURES (if physical spacing available)
    % ========================================================================
    if hasPhysicalSpacing && ~isRadialSymmetry
        features{end+1} = 'sigma_x_nm';
        featureInfo.sigma_x_nm = struct('category', 'physical', ...
            'description', 'Gaussian width in X (nanometers)', ...
            'requires_spacing', true);
        
        features{end+1} = 'sigma_y_nm';
        featureInfo.sigma_y_nm = struct('category', 'physical', ...
            'description', 'Gaussian width in Y (nanometers)', ...
            'requires_spacing', true);
        
        if has3D
            features{end+1} = 'sigma_z_nm';
            featureInfo.sigma_z_nm = struct('category', 'physical', ...
                'description', 'Gaussian width in Z (nanometers)', ...
                'requires_spacing', true);
        end
    end
end

% ============================================================================
% HELPER FUNCTIONS FOR DERIVED FEATURES
% ============================================================================

function result = computeSNR(data)
    try
        if istable(data)
            ii = data.integratedIntensity;
            bg = data.background;
        elseif isstruct(data)
            if isfield(data, 'integratedIntensity')
                ii = [data.integratedIntensity]';
            elseif isfield(data, 'integrated_intensity')
                ii = [data.integrated_intensity]';
            else
                ii = NaN;
            end
            bg = [data.background]';
        else
            result = NaN;
            return;
        end
        result = ii ./ max(bg, 1);
    catch
        result = NaN;
    end
end

function result = computeAmpOverBg(data)
    try
        if istable(data)
            amp = data.amplitude;
            bg = data.background;
        elseif isstruct(data)
            amp = [data.amplitude]';
            bg = [data.background]';
        else
            result = NaN;
            return;
        end
        result = amp ./ max(bg, 1);
    catch
        result = NaN;
    end
end

function result = computeSigmaXYRatio(data)
    try
        if istable(data)
            sx = data.sigma_x;
            sy = data.sigma_y;
        elseif isstruct(data)
            sx = [data.sigma_x]';
            sy = [data.sigma_y]';
        else
            result = NaN;
            return;
        end
        result = sx ./ max(sy, 0.01);
    catch
        result = NaN;
    end
end

function result = computeSigmaSum(data, has3D)
    try
        if istable(data)
            sx = data.sigma_x;
            sy = data.sigma_y;
            result = sx + sy;
            if has3D && ismember('sigma_z', data.Properties.VariableNames)
                sz = data.sigma_z;
                valid_z = ~isnan(sz);
                result(valid_z) = result(valid_z) + sz(valid_z);
            end
        elseif isstruct(data)
            sx = [data.sigma_x]';
            sy = [data.sigma_y]';
            result = sx + sy;
            if has3D && isfield(data, 'sigma_z')
                sz = [data.sigma_z]';
                valid_z = ~isnan(sz);
                result(valid_z) = result(valid_z) + sz(valid_z);
            end
        else
            result = NaN;
        end
    catch
        result = NaN;
    end
end

function result = computeSigmaProduct(data, has3D)
    try
        if istable(data)
            sx = data.sigma_x;
            sy = data.sigma_y;
            result = sx .* sy;
            if has3D && ismember('sigma_z', data.Properties.VariableNames)
                sz = data.sigma_z;
                valid_z = ~isnan(sz);
                result(valid_z) = result(valid_z) .* sz(valid_z);
            end
        elseif isstruct(data)
            sx = [data.sigma_x]';
            sy = [data.sigma_y]';
            result = sx .* sy;
            if has3D && isfield(data, 'sigma_z')
                sz = [data.sigma_z]';
                valid_z = ~isnan(sz);
                result(valid_z) = result(valid_z) .* sz(valid_z);
            end
        else
            result = NaN;
        end
    catch
        result = NaN;
    end
end

function result = computeSigmaZRatio(data)
    try
        if istable(data)
            sx = data.sigma_x;
            sy = data.sigma_y;
            sz = data.sigma_z;
        elseif isstruct(data)
            sx = [data.sigma_x]';
            sy = [data.sigma_y]';
            sz = [data.sigma_z]';
        else
            result = NaN;
            return;
        end
        result = sz ./ max((sx + sy)/2, 0.01);
    catch
        result = NaN;
    end
end

