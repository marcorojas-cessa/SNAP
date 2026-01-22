function exportChannelDataStandardized(outputPath, imageName, channelResults, fitParams, options)
% exportChannelDataStandardized - Standardized export for channel signal data
%
% ============================================================================
% SIGNAL DATA EXPORT (Part of Consistency System)
% ============================================================================
%
% This function exports SIGNAL data (spot/maxima detections and Gaussian fits).
% It is independent of nucleus measurements but part of the overall system.
%
% RELATIONSHIP TO OTHER EXPORTS:
%   - exportNucleiDataStandardized()       → nucleus data only
%   - exportChannelDataStandardized()      → signal data only (THIS FUNCTION)
%   - exportNucleiSignalDataStandardized() → combines both (links signals to nuclei)
%
% SIGNAL DATA CONSISTENCY:
%   Signal measurements use computeSignalMeasurements() as single source of truth.
%   The same signal data appears in:
%   1. This export (channel-only)
%   2. Nuclei+signal export (via signalComposition)
%   Both guaranteed identical because they use the same source.
%
% ============================================================================
%
% USAGE:
%   exportChannelDataStandardized(outputPath, imageName, channelResults, fitParams, options)
%
% INPUTS:
%   outputPath      - Base path for output files (without extension)
%   imageName       - Name of the image (e.g., 'image001')
%   channelResults  - Struct with fields:
%                     .maxima_coords (Nx3 array in ARRAY CONVENTION [row,col,slice])
%                     .fit_results (array of fit result structs, optional)
%   fitParams       - Fitting parameters struct containing:
%                     .imageType (required)
%                     .gaussFitMethod (required if fitting performed)
%                     .gaussFitMode (optional, legacy field)
%                     .gaussFitBgCorrMethod, .gaussFitBgCorrWidth, etc. (optional)
%   options         - Optional struct with fields:
%                     .spacing ([dy, dx, dz] in microns, optional)
%                     .is_3d (logical, default: true)
%
% OUTPUTS:
%   Creates:
%     - [outputPath].mat (MATLAB struct)
%     - [outputPath].csv (CSV file)
%     - [outputPath]_receipt.txt (Parameter receipt)
%
% COORDINATE CONVENTION:
%   All coordinates stored/exported in ARRAY CONVENTION: [row, col, slice]
%   This matches MATLAB array indexing: imageData(row, col, slice)
%
% CSV FORMAT (Method-Adaptive):
%   Always included: image_name, signal_id, maxima_x, maxima_y, maxima_z
%   If fitting: fitted_x, fitted_y, fitted_z, amplitude*, integrated_intensity,
%               background, r_squared (non-radial symmetry), radial_symmetry_score (radial symmetry),
%               sigma_x, sigma_y, sigma_z, rho_xy, rho_xz, rho_yz (3D distorted),
%               alpha_x, alpha_y, alpha_z (skewed)
%   *amplitude may be amplitude_x, amplitude_y, amplitude_z for 1D methods
%    or amplitude_xy, amplitude_z for 2D+1D methods
%
% MATLAB EXPORT INCLUDES:
%   - All CSV fields
%   - Background ROI information (bounds, width)
%   - Background formula/coefficients (for polynomial/plane fitting)
%   - Metadata on coordinate convention

    % Parse options
    if nargin < 5 || isempty(options)
        options = struct();
    end
    if ~isfield(options, 'is_3d')
        options.is_3d = true;
    end
    if ~isfield(options, 'spacing')
        options.spacing = [];
    end
    
    % Check if fitting was performed
    hasFitting = isfield(channelResults, 'fit_results') && ~isempty(channelResults.fit_results);
    
    % Determine fitting method characteristics
    if hasFitting
        fitMethod = channelResults.fit_results(1).fitMethod;
        is1DFit = contains(fitMethod, '1D x,y,z', 'IgnoreCase', true);
        is2DPlus1DFit = contains(fitMethod, '2D+1D', 'IgnoreCase', true);
        isRadialSymmetry = contains(fitMethod, 'Radial Symmetry', 'IgnoreCase', true);
        isSkewed = contains(fitMethod, 'Skewed', 'IgnoreCase', true);
        isDistorted = contains(fitMethod, 'Distorted', 'IgnoreCase', true) && ~isSkewed;
    else
        is1DFit = false;
        is2DPlus1DFit = false;
        isRadialSymmetry = false;
        isSkewed = false;
        isDistorted = false;
    end
    
    % Get number of signals
    numSignals = size(channelResults.maxima_coords, 1);
    
    %% ===== MATLAB EXPORT =====
    % Create structured data
    matData = struct();
    matData.metadata.coordinate_convention = 'ARRAY_CONVENTION';
    matData.metadata.coordinate_format = '[row, col, slice]';
    matData.metadata.description = 'Coordinates match MATLAB array indexing: imageData(row, col, slice)';
    matData.metadata.image_name = imageName;
    matData.metadata.num_signals = numSignals;
    matData.metadata.fitting_performed = hasFitting;
    
    % Store parameters
    matData.parameters.image_type = fitParams.imageType;
    if hasFitting
        if isfield(fitParams, 'gaussFitMode')
            matData.parameters.fit_mode = fitParams.gaussFitMode;
        else
            matData.parameters.fit_mode = 'Not specified';
        end
        matData.parameters.fit_method = fitParams.gaussFitMethod;
        if isfield(fitParams, 'gaussFitVoxelWindowSize')
            matData.parameters.fit_window_size = fitParams.gaussFitVoxelWindowSize;
        end
        if isfield(fitParams, 'gaussFitMaxIterations')
            matData.parameters.max_iterations = fitParams.gaussFitMaxIterations;
        end
        if isfield(fitParams, 'gaussFitTolerance')
            matData.parameters.tolerance = fitParams.gaussFitTolerance;
        end
        if isRadialSymmetry && isfield(fitParams, 'gaussFitRadialRadius')
            matData.parameters.radial_radius = fitParams.gaussFitRadialRadius;
        end
    end
    matData.parameters.bg_method = fitParams.gaussFitBgCorrMethod;
    matData.parameters.bg_width = fitParams.gaussFitBgCorrWidth;
    if isfield(fitParams, 'gaussFitPolyDegree')
        matData.parameters.bg_poly_degree = fitParams.gaussFitPolyDegree;
    end
    
    % Pre-allocate signals struct array with all fields (prevents dissimilar structures error)
    matData.signals = struct('signal_id', cell(numSignals, 1), ...
                             'image_name', cell(numSignals, 1), ...
                             'maxima_coords', cell(numSignals, 1), ...
                             'fitted_coords', cell(numSignals, 1), ...
                             'amplitude', cell(numSignals, 1), ...
                             'amplitude_x', cell(numSignals, 1), ...
                             'amplitude_y', cell(numSignals, 1), ...
                             'amplitude_z', cell(numSignals, 1), ...
                             'amplitude_xy', cell(numSignals, 1), ...
                             'integrated_intensity', cell(numSignals, 1), ...
                             'background_value', cell(numSignals, 1), ...
                             'background_method', cell(numSignals, 1), ...
                             'background_coeffs', cell(numSignals, 1), ...
                             'background_degree', cell(numSignals, 1), ...
                             'background_roi_width', cell(numSignals, 1), ...
                             'background_roi_bounds', cell(numSignals, 1), ...
                             'r_squared', cell(numSignals, 1), ...
                             'radial_symmetry_score', cell(numSignals, 1), ...
                             'sigma_x', cell(numSignals, 1), ...
                             'sigma_y', cell(numSignals, 1), ...
                             'sigma_z', cell(numSignals, 1), ...
                             'rho_xy', cell(numSignals, 1), ...
                             'rho_xz', cell(numSignals, 1), ...
                             'rho_yz', cell(numSignals, 1), ...
                             'alpha_x', cell(numSignals, 1), ...
                             'alpha_y', cell(numSignals, 1), ...
                             'alpha_z', cell(numSignals, 1));
    
    for i = 1:numSignals
        % Initialize basic fields
        matData.signals(i).signal_id = i;
        matData.signals(i).image_name = imageName;
        matData.signals(i).maxima_coords = channelResults.maxima_coords(i, :);
        
        % Fitting results (if available)
        if hasFitting && i <= length(channelResults.fit_results)
            fitRes = channelResults.fit_results(i);
            
            % Fitted center (global coordinates, ARRAY CONVENTION)
            if isfield(fitRes, 'globalFitCenter') && ~isempty(fitRes.globalFitCenter)
                matData.signals(i).fitted_coords = fitRes.globalFitCenter;
            end
            
            % Amplitude (method-dependent)
            if is1DFit
                if ~isnan(fitRes.amplitude_x)
                    matData.signals(i).amplitude_x = fitRes.amplitude_x;
                    matData.signals(i).amplitude_y = fitRes.amplitude_y;
                    matData.signals(i).amplitude_z = fitRes.amplitude_z;
                end
            elseif is2DPlus1DFit
                if ~isnan(fitRes.amplitude_xy)
                    matData.signals(i).amplitude_xy = fitRes.amplitude_xy;
                    matData.signals(i).amplitude_z = fitRes.amplitude_z;
                end
            else
                if ~isnan(fitRes.amplitude)
                    matData.signals(i).amplitude = fitRes.amplitude;
                end
            end
            
            % Integrated intensity
            if ~isnan(fitRes.integratedIntensity)
                matData.signals(i).integrated_intensity = fitRes.integratedIntensity;
            end
            
            % Background (adaptive: value or formula)
            if ~isnan(fitRes.background)
                matData.signals(i).background_value = fitRes.background;
            end
            if isfield(fitRes, 'background_method') && ~isempty(fitRes.background_method)
                matData.signals(i).background_method = fitRes.background_method;
            end
            if isfield(fitRes, 'background_coeffs') && ~isempty(fitRes.background_coeffs)
                matData.signals(i).background_coeffs = fitRes.background_coeffs;
            end
            if isfield(fitRes, 'background_degree') && ~isnan(fitRes.background_degree)
                matData.signals(i).background_degree = fitRes.background_degree;
            end
            if isfield(fitRes, 'background_roi_width')
                matData.signals(i).background_roi_width = fitRes.background_roi_width;
            end
            if isfield(fitRes, 'background_roi_bounds') && ~isempty(fitRes.background_roi_bounds)
                matData.signals(i).background_roi_bounds = fitRes.background_roi_bounds;
            end
            
            % R-squared (non-radial symmetry methods only)
            if ~isRadialSymmetry && ~isnan(fitRes.r_squared)
                matData.signals(i).r_squared = fitRes.r_squared;
            end
            
            % Radial symmetry quality score
            if isRadialSymmetry && isfield(fitRes, 'radialSymmetryScore') && ~isnan(fitRes.radialSymmetryScore)
                matData.signals(i).radial_symmetry_score = fitRes.radialSymmetryScore;
            end
            
            % Sigma (for Gaussian methods)
            if ~isRadialSymmetry
                if ~isnan(fitRes.sigma_x)
                    matData.signals(i).sigma_x = fitRes.sigma_x;
                end
                if ~isnan(fitRes.sigma_y)
                    matData.signals(i).sigma_y = fitRes.sigma_y;
                end
                if options.is_3d && ~isnan(fitRes.sigma_z)
                    matData.signals(i).sigma_z = fitRes.sigma_z;
                end
            end
            
            % Rho (for distorted/skewed methods in 3D)
            if options.is_3d && (isDistorted || isSkewed)
                if ~isnan(fitRes.rho_xy)
                    matData.signals(i).rho_xy = fitRes.rho_xy;
                end
                if ~isnan(fitRes.rho_xz)
                    matData.signals(i).rho_xz = fitRes.rho_xz;
                end
                if ~isnan(fitRes.rho_yz)
                    matData.signals(i).rho_yz = fitRes.rho_yz;
                end
            end
            
            % Alpha (for skewed methods only)
            if isSkewed
                if ~isnan(fitRes.alpha_x)
                    matData.signals(i).alpha_x = fitRes.alpha_x;
                end
                if ~isnan(fitRes.alpha_y)
                    matData.signals(i).alpha_y = fitRes.alpha_y;
                end
                if options.is_3d && ~isnan(fitRes.alpha_z)
                    matData.signals(i).alpha_z = fitRes.alpha_z;
                end
            end
        end
    end
    
    %% ===== SAVE FILES =====
    matFilename = [outputPath '.mat'];
    csvFilename = [outputPath '.csv'];
    fid = fopen(csvFilename, 'w');
    
    % Build header using SHARED helper (ZERO REDUNDANCY!)
    % This ensures columns EXACTLY match exportNucleiSignalDataStandardized (expanded CSV)
    % NOTE: Channel export is for ONE channel, so channel_id column not needed (false parameter)
    [signal_header_cols, include_flags] = snap_helpers.buildSignalColumnList(hasFitting, fitParams.gaussFitMethod, options.is_3d, false);
    
    % Conditionally add image_name (only for batch exports)
    include_image_name = ~isempty(imageName);
    if include_image_name
        header_cols = ['image_name', signal_header_cols];
            else
        header_cols = signal_header_cols;
    end
    
    % Write header
    fprintf(fid, '%s\n', strjoin(header_cols, ','));
    
    % Save MAT file
    save(matFilename, '-struct', 'matData');
    
    % Write data rows - Use include_flags from shared helper for consistency
    for i = 1:numSignals
        % Convert signal to standard format expected by writeSignalDataRow
        sig = struct();
        sig.signal_id = i;
        sig.channel_id = []; % Not used in this export (no channel_id column)
        sig.maxima_coords = channelResults.maxima_coords(i, :);
        sig.fitted_coords = [NaN, NaN, NaN];
        sig.amplitude = NaN;
        sig.amplitude_x = NaN;
        sig.amplitude_y = NaN;
        sig.amplitude_z = NaN;
        sig.amplitude_xy = NaN;
        sig.integrated_intensity = NaN;
        sig.background = NaN;
        sig.r_squared = NaN;
        sig.radial_symmetry_score = NaN;
        sig.sigma_x = NaN;
        sig.sigma_y = NaN;
        sig.sigma_z = NaN;
        sig.rho_xy = NaN;
        sig.rho_xz = NaN;
        sig.rho_yz = NaN;
        sig.alpha_x = NaN;
        sig.alpha_y = NaN;
        sig.alpha_z = NaN;
        
        % Populate with fit data if available
        if hasFitting && i <= length(channelResults.fit_results)
            fitRes = channelResults.fit_results(i);
            
            if isfield(fitRes, 'globalFitCenter') && ~isempty(fitRes.globalFitCenter)
                sig.fitted_coords = fitRes.globalFitCenter;
            end
            
            % Copy all fit fields
            if isfield(fitRes, 'amplitude'), sig.amplitude = fitRes.amplitude; end
            if isfield(fitRes, 'amplitude_x'), sig.amplitude_x = fitRes.amplitude_x; end
            if isfield(fitRes, 'amplitude_y'), sig.amplitude_y = fitRes.amplitude_y; end
            if isfield(fitRes, 'amplitude_z'), sig.amplitude_z = fitRes.amplitude_z; end
            if isfield(fitRes, 'amplitude_xy'), sig.amplitude_xy = fitRes.amplitude_xy; end
            if isfield(fitRes, 'integratedIntensity'), sig.integrated_intensity = fitRes.integratedIntensity; end
            if isfield(fitRes, 'background'), sig.background = fitRes.background; end
            if isfield(fitRes, 'r_squared'), sig.r_squared = fitRes.r_squared; end
            if isfield(fitRes, 'radialSymmetryScore'), sig.radial_symmetry_score = fitRes.radialSymmetryScore; end
            if isfield(fitRes, 'sigma_x'), sig.sigma_x = fitRes.sigma_x; end
            if isfield(fitRes, 'sigma_y'), sig.sigma_y = fitRes.sigma_y; end
            if isfield(fitRes, 'sigma_z'), sig.sigma_z = fitRes.sigma_z; end
            if isfield(fitRes, 'rho_xy'), sig.rho_xy = fitRes.rho_xy; end
            if isfield(fitRes, 'rho_xz'), sig.rho_xz = fitRes.rho_xz; end
            if isfield(fitRes, 'rho_yz'), sig.rho_yz = fitRes.rho_yz; end
            if isfield(fitRes, 'alpha_x'), sig.alpha_x = fitRes.alpha_x; end
            if isfield(fitRes, 'alpha_y'), sig.alpha_y = fitRes.alpha_y; end
            if isfield(fitRes, 'alpha_z'), sig.alpha_z = fitRes.alpha_z; end
        end
        
        % Use SHARED helper to write signal data (ZERO REDUNDANCY!)
        % This ensures data order and conditionals EXACTLY match include_flags
        signal_row = snap_helpers.writeSignalDataRow(sig, include_flags);
        
        % Prepend image_name ONLY if included in header
        if include_image_name
            row_data = [{imageName}, signal_row];
        else
            row_data = signal_row;
        end
        
        % Write row
        fprintf(fid, '%s\n', strjoin(row_data, ','));
    end
    
    fclose(fid);
    % CSV saved silently
    
    %% ===== PARAMETER RECEIPT (TXT) =====
    txtFilename = [outputPath '_receipt.txt'];
    fid = fopen(txtFilename, 'w');
    
    fprintf(fid, '========================================\n');
    fprintf(fid, 'SNAP Channel Signal Analysis Parameters\n');
    fprintf(fid, '========================================\n\n');
    fprintf(fid, 'Image: %s\n\n', imageName);
    
    fprintf(fid, '--- Image Type ---\n');
    fprintf(fid, 'Type: %s\n\n', fitParams.imageType);
    
    if hasFitting
        fprintf(fid, '--- Fitting Parameters ---\n');
        if isfield(fitParams, 'gaussFitMode')
            fprintf(fid, 'Fitting Mode: %s\n', fitParams.gaussFitMode);
        end
        fprintf(fid, 'Fitting Method: %s\n', fitParams.gaussFitMethod);
        if isfield(fitParams, 'gaussFitVoxelWindowSize')
            fprintf(fid, 'Window Size: %d pixels/voxels\n', fitParams.gaussFitVoxelWindowSize);
        end
        if isfield(fitParams, 'gaussFitMaxIterations')
            fprintf(fid, 'Max Iterations: %d\n', fitParams.gaussFitMaxIterations);
        end
        if isfield(fitParams, 'gaussFitTolerance')
            fprintf(fid, 'Tolerance: %.2e\n', fitParams.gaussFitTolerance);
        end
        if isRadialSymmetry && isfield(fitParams, 'gaussFitRadialRadius')
            fprintf(fid, 'Radial Radius: %d pixels\n', fitParams.gaussFitRadialRadius);
        end
        fprintf(fid, '\n');
        
        fprintf(fid, '--- Background Correction ---\n');
        fprintf(fid, 'Method: %s\n', fitParams.gaussFitBgCorrMethod);
        fprintf(fid, 'Width: %d pixels/voxels\n', fitParams.gaussFitBgCorrWidth);
        if isfield(fitParams, 'gaussFitPolyDegree')
            fprintf(fid, 'Polynomial Degree: %d\n', fitParams.gaussFitPolyDegree);
        end
        fprintf(fid, '\n');
    else
        fprintf(fid, '--- Fitting ---\n');
        fprintf(fid, 'No fitting performed (maxima detection only)\n\n');
    end
    
    fprintf(fid, '--- Coordinate Convention ---\n');
    fprintf(fid, 'Convention: ARRAY (row, col, slice)\n');
    fprintf(fid, 'Access: imageData(row, col, slice)\n\n');
    
    fprintf(fid, '--- Export Date & Time ---\n');
    fprintf(fid, '%s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    
    fclose(fid);
    % Receipt and MAT saved silently
end

%% Helper function for formatting values
function str = formatValue(val)
    if isnan(val) || isempty(val)
        str = 'NaN';
    else
        str = sprintf('%.4f', val);
    end
end

