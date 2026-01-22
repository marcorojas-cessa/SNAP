function [header_cols, include_flags] = buildSignalColumnList(hasFitting, fitMethod, is_3d, includeChannelID)
% buildSignalColumnList - Determine which signal columns to include in CSV
%
% ============================================================================
% SINGLE SOURCE OF TRUTH for signal column selection
% ============================================================================
%
% This function encapsulates the EXACT logic for deciding which signal
% columns to include in CSV exports. Used by:
%   1. exportChannelDataStandardized (includeChannelID = false, single channel export)
%   2. exportNucleiSignalDataStandardized (includeChannelID = true, multi-channel export)
%
% WHY THIS EXISTS:
%   Ensures both exports have IDENTICAL column selection logic for signals.
%   Without this, we'd duplicate the method-adaptive logic in two places.
%
%   PROBLEM THIS SOLVES:
%     Before: Each export had its own signal column logic
%     - exportChannelData: Only includes amplitude_x/y/z if 1D fit method
%     - exportNucleiSignal: Includes ALL amplitude variants (x,y,z,xy) with NaNs
%     - User sees different columns (CONFUSING!)
%
%     Now: All exports call this function
%     - Same method → same columns
%     - exportChannelData: only applicable amplitude columns
%     - exportNucleiSignal: only applicable amplitude columns (MATCHES!)
%     - User sees consistent column sets (CLEAR!)
%
% METHOD-ADAPTIVE CONDITIONAL LOGIC:
%   Different fitting methods produce different parameters:
%   
%   If NO fitting:
%     - Only include: channel_id, signal_id, maxima_x/y/z
%   
%   If 1D (X,Y,Z) Gaussian:
%     - Include: amplitude_x, amplitude_y, amplitude_z (NOT amplitude or amplitude_xy)
%   
%   If 2D (XY) + 1D (Z) Gaussian:
%     - Include: amplitude_xy, amplitude_z (NOT amplitude or amplitude_x/y)
%   
%   If Standard Gaussian:
%     - Include: amplitude (NOT amplitude_x/y/z/xy)
%   
%   If Radial Symmetry:
%     - Include: radial_symmetry_score (NOT r_squared)
%     - Exclude: sigma fields (radial symmetry doesn't fit Gaussians)
%   
%   If Distorted/Skewed (3D only):
%     - Include: rho_xy, rho_xz, rho_yz
%   
%   If Skewed only:
%     - Include: alpha_x, alpha_y, alpha_z
%
% ============================================================================
%
% INPUTS:
%   hasFitting      - Boolean indicating if Gaussian fitting was performed
%   fitMethod       - String describing fit method (e.g., '1D (X,Y,Z) Gaussian')
%   is_3d           - Boolean indicating if data is 3D
%   includeChannelID - Boolean indicating if channel_id column should be included
%                      true for multi-channel exports (nuclei+signal)
%                      false for single-channel exports (channel-only)
%
% OUTPUTS:
%   header_cols   - Cell array of column names to include
%   include_flags - Struct with boolean flags for which columns to include

% Handle optional parameter (backward compatibility)
if nargin < 4
    includeChannelID = true;  % Default to including it
end

% Initialize - conditionally include channel_id
if includeChannelID
    header_cols = {'channel_id', 'signal_id', 'maxima_x', 'maxima_y', 'maxima_z'};
else
    header_cols = {'signal_id', 'maxima_x', 'maxima_y', 'maxima_z'};
end

include_flags = struct();
include_flags.hasFitting = hasFitting;
include_flags.includeChannelID = includeChannelID;

if ~hasFitting
    % No fitting - only maxima coordinates
    return;
end

% Determine fitting method characteristics
is1DFit = contains(fitMethod, '1D (X,Y,Z)', 'IgnoreCase', true) || ...
          contains(fitMethod, '1D x,y,z', 'IgnoreCase', true);
is2DPlus1DFit = contains(fitMethod, '2D (XY) + 1D (Z)', 'IgnoreCase', true) || ...
                contains(fitMethod, '2D+1D', 'IgnoreCase', true);
isRadialSymmetry = contains(fitMethod, 'Radial Symmetry', 'IgnoreCase', true);
isSkewed = contains(fitMethod, 'Skewed', 'IgnoreCase', true);
isDistorted = contains(fitMethod, 'Distorted', 'IgnoreCase', true) && ~isSkewed;

% Store flags
include_flags.is1DFit = is1DFit;
include_flags.is2DPlus1DFit = is2DPlus1DFit;
include_flags.isRadialSymmetry = isRadialSymmetry;
include_flags.isSkewed = isSkewed;
include_flags.isDistorted = isDistorted;
include_flags.is_3d = is_3d;

% Fitted coordinates (always include if fitting)
header_cols = [header_cols, {'fitted_x', 'fitted_y', 'fitted_z'}];

% Amplitude columns (method-dependent) - ONLY include applicable ones
if is1DFit
    header_cols = [header_cols, {'amplitude_x', 'amplitude_y', 'amplitude_z'}];
elseif is2DPlus1DFit
    header_cols = [header_cols, {'amplitude_xy', 'amplitude_z'}];
else
    header_cols = [header_cols, {'amplitude'}];
end

% Integrated intensity and background (always include if fitting)
header_cols = [header_cols, {'integrated_intensity', 'background'}];

% Quality metric (method-dependent) - ONLY include applicable one
if isRadialSymmetry
    header_cols = [header_cols, {'radial_symmetry_score'}];
else
    header_cols = [header_cols, {'r_squared'}];
end

% Sigma (for Gaussian methods only, not radial symmetry)
if ~isRadialSymmetry
    if is_3d
        header_cols = [header_cols, {'sigma_x', 'sigma_y', 'sigma_z'}];
    else
        header_cols = [header_cols, {'sigma_x', 'sigma_y'}];
    end
end

% Rho (for distorted/skewed in 3D only)
if is_3d && (isDistorted || isSkewed)
    header_cols = [header_cols, {'rho_xy', 'rho_xz', 'rho_yz'}];
end

% Alpha (for skewed only)
if isSkewed
    if is_3d
        header_cols = [header_cols, {'alpha_x', 'alpha_y', 'alpha_z'}];
    else
        header_cols = [header_cols, {'alpha_x', 'alpha_y'}];
    end
end

end

