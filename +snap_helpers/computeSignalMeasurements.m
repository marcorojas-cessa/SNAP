function [signalData, methodInfo] = computeSignalMeasurements(maximaCoords, fitResults, fitMethod)
% computeSignalMeasurements - SINGLE SOURCE OF TRUTH for signal field extraction
%
% ============================================================================
% CRITICAL: This is THE ONLY PLACE where signal fit data is extracted
% ============================================================================
%
% WHY THIS FUNCTION EXISTS:
%   To ensure signal measurements are IDENTICAL across all export types.
%   
%   Signal data appears in TWO exports:
%   1. Channel export (exportChannelDataStandardized)
%   2. Nuclei+Signal export (exportNucleiSignalDataStandardized)
%
%   Both MUST show identical signal measurements. This function guarantees that.
%
% HOW IT WORKS:
%   - Receives raw fitResults from fitGaussians()
%   - Extracts relevant fields based on fit method (method-adaptive)
%   - Returns standardized struct array with ALL fields initialized
%   - All exports calling this function get IDENTICAL signal data
%
% METHOD-ADAPTIVE EXTRACTION:
%   Different fitting methods produce different parameters:
%   - 1D (X,Y,Z) Gaussian: amplitude_x, amplitude_y, amplitude_z
%   - 2D (XY) + 1D (Z) Gaussian: amplitude_xy, amplitude_z
%   - Radial Symmetry: radial_symmetry_score (instead of r_squared)
%   - Distorted/Skewed: rho_xy, rho_xz, rho_yz (correlation coefficients)
%   - Skewed only: alpha_x, alpha_y, alpha_z (skew parameters)
%
% PRE-ALLOCATION STRATEGY:
%   ALL fields are pre-allocated in struct array to prevent "dissimilar structures"
%   errors in MATLAB. Fields not applicable to the fit method are set to NaN.
%
% USAGE CHAIN:
%   fitGaussians() → fitResults
%        ↓
%   computeSignalMeasurements(coords, fitResults, method) → signalData
%        ↓
%   ├─→ exportChannelDataStandardized() uses signalData
%   └─→ analyzeNucleiSignalComposition() uses signalData
%            ↓
%       exportNucleiSignalDataStandardized() uses embedded signalData
%
%   Result: ALL exports have IDENTICAL signal measurements
%
% ============================================================================
%
% INPUTS:
%   maximaCoords - Nx3 array of maxima coordinates [row, col, slice] ARRAY CONVENTION
%   fitResults   - Array of fit result structs (optional, can be empty)
%   fitMethod    - String describing fit method (used for field detection)
%
% OUTPUTS:
%   signalData   - Struct array with one entry per signal, containing:
%                  .signal_id
%                  .maxima_coords ([row, col, slice])
%                  .fitted_coords ([row, col, slice]) - if fitting
%                  .amplitude (or amplitude_x/y/z/xy for method-specific)
%                  .integrated_intensity
%                  .background
%                  .r_squared OR radial_symmetry_score
%                  .sigma_x/y/z (if Gaussian)
%                  .rho_xy/xz/yz (if distorted/skewed)
%                  .alpha_x/y/z (if skewed)
%   methodInfo   - Struct with boolean flags:
%                  .hasFitting, .is1DFit, .is2DPlus1DFit
%                  .isRadialSymmetry, .isSkewed, .isDistorted

    numSignals = size(maximaCoords, 1);
    
    % Determine fitting method characteristics
    hasFitting = ~isempty(fitResults) && length(fitResults) > 0;
    
    if hasFitting
        is1DFit = contains(fitMethod, '1D x,y,z', 'IgnoreCase', true) || ...
                  contains(fitMethod, '1D (X,Y,Z)', 'IgnoreCase', true);
        is2DPlus1DFit = contains(fitMethod, '2D+1D', 'IgnoreCase', true) || ...
                        contains(fitMethod, '2D (XY) + 1D (Z)', 'IgnoreCase', true);
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
    
    % Store method info
    methodInfo = struct();
    methodInfo.hasFitting = hasFitting;
    methodInfo.is1DFit = is1DFit;
    methodInfo.is2DPlus1DFit = is2DPlus1DFit;
    methodInfo.isRadialSymmetry = isRadialSymmetry;
    methodInfo.isSkewed = isSkewed;
    methodInfo.isDistorted = isDistorted;
    
    % Pre-allocate struct array with all fields to ensure consistency
    % This prevents MATLAB "dissimilar structures" errors
    signalData = struct('signal_id', cell(numSignals, 1), ...
                        'maxima_coords', cell(numSignals, 1), ...
                        'fitted_coords', cell(numSignals, 1), ...
                        'amplitude', cell(numSignals, 1), ...
                        'amplitude_x', cell(numSignals, 1), ...
                        'amplitude_y', cell(numSignals, 1), ...
                        'amplitude_z', cell(numSignals, 1), ...
                        'amplitude_xy', cell(numSignals, 1), ...
                        'integrated_intensity', cell(numSignals, 1), ...
                        'background', cell(numSignals, 1), ...
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
        % Initialize all fields to proper values
        signalData(i).signal_id = i;
        signalData(i).maxima_coords = maximaCoords(i, :);  % [row, col, slice] ARRAY CONVENTION
        signalData(i).fitted_coords = [NaN, NaN, NaN];
        signalData(i).amplitude = NaN;
        signalData(i).amplitude_x = NaN;
        signalData(i).amplitude_y = NaN;
        signalData(i).amplitude_z = NaN;
        signalData(i).amplitude_xy = NaN;
        signalData(i).integrated_intensity = NaN;
        signalData(i).background = NaN;
        signalData(i).r_squared = NaN;
        signalData(i).radial_symmetry_score = NaN;
        signalData(i).sigma_x = NaN;
        signalData(i).sigma_y = NaN;
        signalData(i).sigma_z = NaN;
        signalData(i).rho_xy = NaN;
        signalData(i).rho_xz = NaN;
        signalData(i).rho_yz = NaN;
        signalData(i).alpha_x = NaN;
        signalData(i).alpha_y = NaN;
        signalData(i).alpha_z = NaN;
        
        % Populate with fit data if available
        if hasFitting && i <= length(fitResults)
            fitRes = fitResults(i);
            
            % Fitted center (global coordinates, ARRAY CONVENTION)
            if isfield(fitRes, 'globalFitCenter') && ~isempty(fitRes.globalFitCenter)
                signalData(i).fitted_coords = fitRes.globalFitCenter;
            end
            
            % Amplitude (method-dependent)
            if is1DFit
                if isfield(fitRes, 'amplitude_x') && ~isnan(fitRes.amplitude_x)
                    signalData(i).amplitude_x = fitRes.amplitude_x;
                    signalData(i).amplitude_y = fitRes.amplitude_y;
                    signalData(i).amplitude_z = fitRes.amplitude_z;
                end
            elseif is2DPlus1DFit
                if isfield(fitRes, 'amplitude_xy') && ~isnan(fitRes.amplitude_xy)
                    signalData(i).amplitude_xy = fitRes.amplitude_xy;
                    signalData(i).amplitude_z = fitRes.amplitude_z;
                end
            else
                if isfield(fitRes, 'amplitude') && ~isnan(fitRes.amplitude)
                    signalData(i).amplitude = fitRes.amplitude;
                end
            end
            
            % Integrated intensity
            if isfield(fitRes, 'integratedIntensity') && ~isnan(fitRes.integratedIntensity)
                signalData(i).integrated_intensity = fitRes.integratedIntensity;
            end
            
            % Background
            if isfield(fitRes, 'background') && ~isnan(fitRes.background)
                signalData(i).background = fitRes.background;
            end
            
            % Quality metric (R-squared or radial symmetry score)
            if isRadialSymmetry
                if isfield(fitRes, 'radialSymmetryScore') && ~isnan(fitRes.radialSymmetryScore)
                    signalData(i).radial_symmetry_score = fitRes.radialSymmetryScore;
                end
            else
                if isfield(fitRes, 'r_squared')
                    if isscalar(fitRes.r_squared)
                        signalData(i).r_squared = fitRes.r_squared;
                    else
                        signalData(i).r_squared = mean(fitRes.r_squared(~isnan(fitRes.r_squared)));
                    end
                end
            end
            
            % Sigma (for Gaussian methods)
            if ~isRadialSymmetry
                if isfield(fitRes, 'sigma_x') && ~isnan(fitRes.sigma_x)
                    signalData(i).sigma_x = fitRes.sigma_x;
                end
                if isfield(fitRes, 'sigma_y') && ~isnan(fitRes.sigma_y)
                    signalData(i).sigma_y = fitRes.sigma_y;
                end
                if isfield(fitRes, 'sigma_z') && ~isnan(fitRes.sigma_z)
                    signalData(i).sigma_z = fitRes.sigma_z;
                end
            end
            
            % Rho (for distorted/skewed methods)
            if isDistorted || isSkewed
                if isfield(fitRes, 'rho_xy') && ~isnan(fitRes.rho_xy)
                    signalData(i).rho_xy = fitRes.rho_xy;
                end
                if isfield(fitRes, 'rho_xz') && ~isnan(fitRes.rho_xz)
                    signalData(i).rho_xz = fitRes.rho_xz;
                end
                if isfield(fitRes, 'rho_yz') && ~isnan(fitRes.rho_yz)
                    signalData(i).rho_yz = fitRes.rho_yz;
                end
            end
            
            % Alpha (for skewed methods only)
            if isSkewed
                if isfield(fitRes, 'alpha_x') && ~isnan(fitRes.alpha_x)
                    signalData(i).alpha_x = fitRes.alpha_x;
                end
                if isfield(fitRes, 'alpha_y') && ~isnan(fitRes.alpha_y)
                    signalData(i).alpha_y = fitRes.alpha_y;
                end
                if isfield(fitRes, 'alpha_z') && ~isnan(fitRes.alpha_z)
                    signalData(i).alpha_z = fitRes.alpha_z;
                end
            end
        end
    end
end

