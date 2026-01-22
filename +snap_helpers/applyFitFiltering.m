function [filtered_fits, filter_mask] = applyFitFiltering(fits, channel_idx, handles)
% Applies fit filtering based on user-defined criteria
% Inputs:
%   fits - structure containing fit results
%   channel_idx - channel index for parameter lookup
%   handles - handles structure containing UI parameters
% Outputs:
%   filtered_fits - structure containing filtered fit results
%   filter_mask - logical array indicating which fits passed the filter

    % Check if fit filtering is enabled for this channel
    if ~handles.fitFilterEnabledChecks(channel_idx).Value
        filtered_fits = fits;
        filter_mask = true(1, length(fits));
        return;
    end
    
    % Check fitting method to determine which filters are applicable
    fitting_method = handles.gaussFitMethodDrop(channel_idx).Value;
    is_radial_symmetry = strcmp(fitting_method, 'Radial Symmetry');
    
    % Initialize filter mask (all fits pass by default)
    num_fits = length(fits);
    filter_mask = true(1, num_fits);
    
    % Apply R² filtering (or quality score for radial symmetry)
    if handles.fitFilterRSquaredEnabledChecks(channel_idx).Value
        min_r_squared = handles.fitFilterRSquaredMinInputs(channel_idx).Value;
        max_r_squared = handles.fitFilterRSquaredMaxInputs(channel_idx).Value;
        
        if is_radial_symmetry
            % For radial symmetry, use quality score
            if isfield(fits, 'radialSymmetryQuality')
                quality_values = [fits.radialSymmetryQuality];
                quality_mask = (quality_values >= min_r_squared) & (quality_values <= max_r_squared);
                filter_mask = filter_mask & quality_mask;
            else
                warning('Radial symmetry quality scores not available. Skipping quality filtering.');
            end
        else
            % For Gaussian fits, use R²
            if isfield(fits, 'r_squared')
                % Handle both scalar and array R² values
                r_squared_mask = false(1, num_fits);
                for i = 1:num_fits
                    r2_val = fits(i).r_squared;
                    if isscalar(r2_val)
                        % Single R² value
                        r_squared_mask(i) = (r2_val >= min_r_squared) && (r2_val <= max_r_squared);
                    else
                        % Multiple R² values (e.g., for 1D fits: X, Y, Z)
                        % Use the mean or minimum R² across dimensions
                        mean_r2 = mean(r2_val);
                        r_squared_mask(i) = (mean_r2 >= min_r_squared) && (mean_r2 <= max_r_squared);
                    end
                end
                filter_mask = filter_mask & r_squared_mask;
            else
                warning('R-squared values not available in fit results. Skipping R² filtering.');
            end
        end
    end
    
    % Apply sigma sum filtering (not applicable for radial symmetry with NaN sigmas)
    if handles.fitFilterSigmaSumEnabledChecks(channel_idx).Value && ~is_radial_symmetry
        min_sigma_sum = handles.fitFilterSigmaSumMinInputs(channel_idx).Value;
        max_sigma_sum = handles.fitFilterSigmaSumMaxInputs(channel_idx).Value;
        if isfield(fits, 'sigma_x') && isfield(fits, 'sigma_y')
            sigma_x_vals = [fits.sigma_x];
            sigma_y_vals = [fits.sigma_y];
            sigma_sum = sigma_x_vals + sigma_y_vals;
            if isfield(fits, 'sigma_z')
                sigma_z_vals = [fits.sigma_z];
                % Only add z component if not NaN
                valid_z = ~isnan(sigma_z_vals);
                sigma_sum(valid_z) = sigma_sum(valid_z) + sigma_z_vals(valid_z);
            end
            % Filter out NaN values
            valid_sigma = ~isnan(sigma_sum);
            sigma_sum_mask = valid_sigma & (sigma_sum >= min_sigma_sum) & (sigma_sum <= max_sigma_sum);
            filter_mask = filter_mask & sigma_sum_mask;
        else
            warning('Sigma values not available in fit results. Skipping sigma sum filtering.');
        end
    end
    
    % Apply amplitude filtering (not applicable for radial symmetry)
    if handles.fitFilterAmplitudeEnabledChecks(channel_idx).Value && ~is_radial_symmetry
        min_amplitude = handles.fitFilterAmplitudeMinInputs(channel_idx).Value;
        max_amplitude = handles.fitFilterAmplitudeMaxInputs(channel_idx).Value;
        if isfield(fits, 'amplitude')
            amplitude_values = [fits.amplitude];
            % Filter out NaN values for radial symmetry compatibility
            valid_amplitudes = ~isnan(amplitude_values);
            amplitude_mask = valid_amplitudes & (amplitude_values >= min_amplitude) & (amplitude_values <= max_amplitude);
            filter_mask = filter_mask & amplitude_mask;
        else
            warning('Amplitude values not available in fit results. Skipping amplitude filtering.');
        end
    end
    
    % Apply intensity filtering
    if handles.fitFilterIntensityEnabledChecks(channel_idx).Value
        min_intensity = handles.fitFilterIntensityMinInputs(channel_idx).Value;
        max_intensity = handles.fitFilterIntensityMaxInputs(channel_idx).Value;
        % Check for both 'intensity' and 'integratedIntensity' field names (for compatibility)
        if isfield(fits, 'intensity')
            intensity_values = [fits.intensity];
            intensity_mask = (intensity_values >= min_intensity) & (intensity_values <= max_intensity);
            filter_mask = filter_mask & intensity_mask;
        elseif isfield(fits, 'integratedIntensity')
            intensity_values = [fits.integratedIntensity];
            intensity_mask = (intensity_values >= min_intensity) & (intensity_values <= max_intensity);
            filter_mask = filter_mask & intensity_mask;
        else
            warning('Intensity values not available in fit results. Skipping intensity filtering.');
        end
    end
    
    % Apply the filter mask
    filtered_fits = fits(filter_mask);
    
    % Display filtering results
    num_filtered = sum(filter_mask);
    num_removed = num_fits - num_filtered;
    if num_removed > 0
        fprintf('Fit filtering: %d fits removed, %d fits remaining (%.1f%% kept)\n', ...
                num_removed, num_filtered, 100 * num_filtered / num_fits);
    end
end
