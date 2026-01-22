function row_data = writeSignalDataRow(sig, include_flags)
% writeSignalDataRow - Build signal data row for CSV using shared logic
%
% ============================================================================
% SINGLE SOURCE OF TRUTH for signal data row writing
% ============================================================================
%
% This function writes signal measurements to CSV in the EXACT order
% and with the EXACT conditional logic as buildSignalColumnList().
%
% CRITICAL: This function MUST write fields in the same order and with
%           the same conditionals as buildSignalColumnList() builds headers.
%
% Used by:
%   1. exportChannelDataStandardized
%   2. exportNucleiSignalDataStandardized (expanded CSV)
%
% CONSISTENCY GUARANTEE:
%   - include_flags comes from buildSignalColumnList()
%   - This function follows the EXACT same conditionals
%   - Result: Data columns perfectly match header columns
%
%   Examples:
%     - If buildSignalColumnList includes amplitude_x/y/z, this writes them
%     - If buildSignalColumnList excludes r_squared (radial symmetry), this doesn't write it
%     - If buildSignalColumnList only includes sigma_x/y (2D), this doesn't write sigma_z
%
%   If this function and buildSignalColumnList() ever get out of sync,
%   you'll get CSV parsing errors or misaligned columns. That's why they
%   must be updated together.
%
% ============================================================================
%
% INPUTS:
%   sig           - Signal struct with all measurement fields
%   include_flags - Flags from buildSignalColumnList() indicating which fields to include
%
% OUTPUTS:
%   row_data - Cell array of formatted values matching header order

row_data = {};

% Channel ID (conditional - only if include_flags says to include it)
if include_flags.includeChannelID
    if isfield(sig, 'channel_id') && ~isempty(sig.channel_id)
        row_data{end+1} = sprintf('%d', sig.channel_id);
    else
        row_data{end+1} = 'NaN';  % Should not happen if includeChannelID is true
    end
end

% Signal ID (always included)
row_data{end+1} = sprintf('%d', sig.signal_id);

% Maxima coordinates (always included)
maxima = sig.maxima_coords;
row_data{end+1} = sprintf('%.4f', maxima(1));  % row
row_data{end+1} = sprintf('%.4f', maxima(2));  % col
row_data{end+1} = sprintf('%.4f', maxima(3));  % slice

if ~include_flags.hasFitting
    % No fitting - done
    return;
end

% Fitted coordinates (always if fitting)
fitted = sig.fitted_coords;
row_data{end+1} = formatValue(fitted(1));
row_data{end+1} = formatValue(fitted(2));
row_data{end+1} = formatValue(fitted(3));

% Amplitude (method-dependent) - ONLY write applicable fields
if include_flags.is1DFit
    row_data{end+1} = formatValue(sig.amplitude_x);
    row_data{end+1} = formatValue(sig.amplitude_y);
    row_data{end+1} = formatValue(sig.amplitude_z);
elseif include_flags.is2DPlus1DFit
    row_data{end+1} = formatValue(sig.amplitude_xy);
    row_data{end+1} = formatValue(sig.amplitude_z);
else
    row_data{end+1} = formatValue(sig.amplitude);
end

% Integrated intensity and background (always if fitting)
row_data{end+1} = formatValue(sig.integrated_intensity);
row_data{end+1} = formatValue(sig.background);

% Quality metric (method-dependent) - ONLY write applicable field
if include_flags.isRadialSymmetry
    row_data{end+1} = formatValue(sig.radial_symmetry_score);
else
    row_data{end+1} = formatValue(sig.r_squared);
end

% Sigma (for Gaussian methods, not radial symmetry)
if ~include_flags.isRadialSymmetry
    row_data{end+1} = formatValue(sig.sigma_x);
    row_data{end+1} = formatValue(sig.sigma_y);
    if include_flags.is_3d
        row_data{end+1} = formatValue(sig.sigma_z);
    end
end

% Rho (for distorted/skewed in 3D only)
if include_flags.is_3d && (include_flags.isDistorted || include_flags.isSkewed)
    row_data{end+1} = formatValue(sig.rho_xy);
    row_data{end+1} = formatValue(sig.rho_xz);
    row_data{end+1} = formatValue(sig.rho_yz);
end

% Alpha (for skewed only)
if include_flags.isSkewed
    row_data{end+1} = formatValue(sig.alpha_x);
    row_data{end+1} = formatValue(sig.alpha_y);
    if include_flags.is_3d
        row_data{end+1} = formatValue(sig.alpha_z);
    end
end

end

function str = formatValue(val)
    if isnan(val) || isempty(val)
        str = 'NaN';
    else
        str = sprintf('%.4f', val);
    end
end

