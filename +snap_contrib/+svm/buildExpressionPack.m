function pack = buildExpressionPack(varargin)
% buildExpressionPack - Build a per-channel SVM expression pack for SNAP.
%
% This is an example contribution that expands beyond the default
% "all non-position base features" behavior by adding a curated set of
% scientifically motivated custom expressions.
%
% Usage:
%   pack = snap_contrib.svm.buildExpressionPack('ParameterFile', 'params.mat');
%   pack = snap_contrib.svm.buildExpressionPack('ParameterStruct', paramsStruct);
%
% Name-value options:
%   'ParameterFile'                - SNAP parameter MAT file
%   'ParameterStruct'              - already-loaded parameter struct
%   'Channels'                     - channel indices to include (default: all)
%   'IncludeAdvancedStatsFeatures' - include model-selection/statistical
%                                    base features produced by
%                                    snap_contrib.svm.augmentFitResultsWithModelStats
%
% Output:
%   pack - struct with per-channel selectedFeatures/customExpressions.

    p = inputParser;
    p.addParameter('ParameterFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ParameterStruct', struct(), @isstruct);
    p.addParameter('Channels', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
    p.addParameter('IncludeAdvancedStatsFeatures', false, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    params = opts.ParameterStruct;
    parameterSource = 'ParameterStruct';

    if isempty(fieldnames(params))
        parameterFile = strtrim(char(string(opts.ParameterFile)));
        if isempty(parameterFile)
            error('buildExpressionPack:MissingParameters', ...
                'Provide either ParameterStruct or ParameterFile.');
        end
        params = localLoadParameterStruct(parameterFile);
        parameterSource = parameterFile;
    end

    numChannels = localInferNumChannels(params, 1);
    channels = 1:numChannels;
    if ~isempty(opts.Channels)
        channels = unique(round(opts.Channels(:)'));
        channels = channels(channels >= 1 & channels <= numChannels);
        if isempty(channels)
            error('buildExpressionPack:InvalidChannels', ...
                'No valid channels selected. Valid range is 1..%d.', numChannels);
        end
    end

    includeAdvanced = logical(opts.IncludeAdvancedStatsFeatures);
    advancedStatsFeatures = localAdvancedStatsFeatureNames();
    exprLibrary = localExpressionLibrary();

    channelPacks = repmat(struct( ...
        'channelIdx', 1, ...
        'fittingMethod', '', ...
        'has3D', true, ...
        'selectedFeatures', {{}}, ...
        'customExpressions', struct('name', {}, 'expression', {}), ...
        'advancedStatsFeatureNames', {{}}, ...
        'notes', ''), 1, numel(channels));

    for i = 1:numel(channels)
        ch = channels(i);
        [fittingMethod, has3D] = localInferChannelContext(params, ch);
        [allFeatures, featureInfo] = snap_helpers.classification.getAvailableFeatures( ...
            fittingMethod, has3D, false);

        selected = localNonPositionFeatures(allFeatures, featureInfo);
        if includeAdvanced
            selected = unique([selected, advancedStatsFeatures], 'stable');
        end

        availableSet = unique([allFeatures, selected], 'stable');
        customExpressions = localSelectExpressions(exprLibrary, availableSet);

        channelPacks(i).channelIdx = ch;
        channelPacks(i).fittingMethod = fittingMethod;
        channelPacks(i).has3D = logical(has3D);
        channelPacks(i).selectedFeatures = selected;
        channelPacks(i).customExpressions = customExpressions;
        if includeAdvanced
            channelPacks(i).advancedStatsFeatureNames = advancedStatsFeatures;
            channelPacks(i).notes = sprintf([ ...
                'Includes optional model-stat features from ', ...
                'snap_contrib.svm.augmentFitResultsWithModelStats.']);
        else
            channelPacks(i).advancedStatsFeatureNames = {};
            channelPacks(i).notes = 'UI-compatible expression set using canonical SNAP features only.';
        end
    end

    pack = struct();
    pack.name = 'SNAP Literature-Informed SVM Expression Pack';
    pack.version = '1.0.0';
    pack.generatedOn = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    pack.parameterSource = parameterSource;
    pack.numChannelsDetected = numChannels;
    pack.channelsIncluded = channels;
    pack.includeAdvancedStatsFeatures = includeAdvanced;
    pack.channelPacks = channelPacks;
    pack.rationale = [ ...
        'Expressions emphasize dynamic-range stabilization (log transforms), ', ...
        'signal-vs-background normalization, PSF/shape compactness, anisotropy, ', ...
        'and fit-quality coupling.'];
end

function params = localLoadParameterStruct(paramFile)
    if exist(paramFile, 'file') ~= 2
        error('buildExpressionPack:ParameterFileNotFound', ...
            'Parameter file not found: %s', paramFile);
    end

    raw = load(paramFile);
    if isfield(raw, 'batchConfig') && isfield(raw.batchConfig, 'parameters')
        params = raw.batchConfig.parameters;
        return;
    end
    if isfield(raw, 'paramData') && isfield(raw.paramData, 'parameters')
        params = raw.paramData.parameters;
        return;
    end
    if isfield(raw, 'paramData') && isfield(raw.paramData, 'workflowConfig')
        params = raw.paramData.workflowConfig;
        return;
    end
    if isfield(raw, 'workflowConfig')
        params = raw.workflowConfig;
        return;
    end
    if isfield(raw, 'parameters')
        params = raw.parameters;
        return;
    end
    if isfield(raw, 'lastUsed')
        params = raw.lastUsed;
        return;
    end

    error('buildExpressionPack:ParameterFileUnsupported', ...
        ['Could not find a SNAP parameter struct in %s. ', ...
         'Expected one of: batchConfig.parameters, paramData.parameters, ', ...
         'paramData.workflowConfig, workflowConfig, parameters, or lastUsed.'], paramFile);
end

function numChannels = localInferNumChannels(params, fallback)
    if nargin < 2
        fallback = 1;
    end
    numChannels = fallback;
    if ~isstruct(params)
        return;
    end

    explicitN = nan;
    if isfield(params, 'numChannels')
        explicitN = localToScalarNumeric(params.numChannels, nan);
    elseif isfield(params, 'numChan')
        explicitN = localToScalarNumeric(params.numChan, nan);
    elseif isfield(params, 'numChanDrop')
        explicitN = localToScalarNumeric(params.numChanDrop, nan);
    end
    if isfinite(explicitN) && explicitN >= 1
        numChannels = max(1, round(explicitN));
        return;
    end

    hints = {'gaussFitMethod', 'maximaMode', 'preProcMode', 'classifyEnabled'};
    for i = 1:numel(hints)
        f = hints{i};
        if isfield(params, f)
            numChannels = max(numChannels, localCountChannelsFromValue(params.(f)));
        end
    end
    numChannels = max(1, round(numChannels));
end

function n = localCountChannelsFromValue(v)
    n = 1;
    if isempty(v)
        return;
    end
    if iscell(v) || isstring(v)
        n = numel(v);
    elseif isnumeric(v) || islogical(v)
        if isvector(v)
            n = numel(v);
        end
    end
end

function [fittingMethod, has3D] = localInferChannelContext(params, ch)
    fittingMethod = '3D Gaussian';
    has3D = true;

    if isfield(params, 'gaussFitMethod')
        fm = localGetChannelValue(params.gaussFitMethod, ch, '');
        if ~isempty(fm)
            fittingMethod = char(string(fm));
        end
    end

    modeVal = '';
    if isfield(params, 'maximaMode')
        modeVal = char(string(localGetChannelValue(params.maximaMode, ch, '')));
    elseif isfield(params, 'preProcMode')
        modeVal = char(string(localGetChannelValue(params.preProcMode, ch, '')));
    end

    if ~isempty(modeVal)
        has3D = strcmpi(modeVal, '3D');
    end
end

function value = localGetChannelValue(raw, ch, fallback)
    value = fallback;
    if isempty(raw)
        return;
    end
    if iscell(raw)
        if ch <= numel(raw) && ~isempty(raw{ch})
            value = raw{ch};
        elseif ~isempty(raw{1})
            value = raw{1};
        end
    elseif isstring(raw)
        if ch <= numel(raw)
            value = char(raw(ch));
        else
            value = char(raw(1));
        end
    elseif ischar(raw)
        value = raw;
    elseif isnumeric(raw) || islogical(raw)
        if isscalar(raw)
            value = raw;
        elseif ch <= numel(raw)
            value = raw(ch);
        else
            value = raw(1);
        end
    end
end

function value = localToScalarNumeric(raw, fallback)
    value = fallback;
    if iscell(raw)
        if isempty(raw)
            return;
        end
        raw = raw{1};
    end
    if isstring(raw) || ischar(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            value = parsed;
        end
        return;
    end
    if isnumeric(raw) && isscalar(raw) && isfinite(raw)
        value = double(raw);
    end
end

function selected = localNonPositionFeatures(allFeatures, featureInfo)
    selected = {};
    for i = 1:numel(allFeatures)
        f = allFeatures{i};
        if isfield(featureInfo, f) && isfield(featureInfo.(f), 'category') && ...
                strcmpi(featureInfo.(f).category, 'position')
            continue;
        end
        selected{end+1} = f; %#ok<AGROW>
    end
end

function expressions = localSelectExpressions(exprLibrary, availableFeatures)
    if isempty(exprLibrary)
        expressions = struct('name', {}, 'expression', {});
        return;
    end

    availableSet = string(availableFeatures);
    expressions = struct('name', {}, 'expression', {});
    addedNames = strings(0, 1);

    for i = 1:numel(exprLibrary)
        req = string(exprLibrary(i).requiredFeatures);
        if ~all(ismember(req, availableSet))
            continue;
        end

        exprName = string(exprLibrary(i).name);
        if any(exprName == addedNames)
            continue;
        end

        expressions(end+1) = struct( ... %#ok<AGROW>
            'name', char(exprLibrary(i).name), ...
            'expression', char(exprLibrary(i).expression));
        addedNames(end+1, 1) = exprName; %#ok<AGROW>
    end
end

function names = localAdvancedStatsFeatureNames()
    names = { ...
        'aic_gaussian', ...
        'bic_gaussian', ...
        'aic_flat', ...
        'bic_flat', ...
        'delta_aic_flat_minus_gaussian', ...
        'delta_bic_flat_minus_gaussian', ...
        'nll_poisson_gaussian', ...
        'nll_poisson_flat', ...
        'delta_nll_poisson_flat_minus_gaussian', ...
        'lrt_stat_flat_vs_gaussian', ...
        'residual_mean', ...
        'residual_std', ...
        'residual_mad', ...
        'residual_skewness', ...
        'residual_kurtosis_excess', ...
        'residual_over_signal', ...
        'normality_ks_pvalue' ...
    };
end

function lib = localExpressionLibrary()
    lib = struct('name', {}, 'expression', {}, 'requiredFeatures', {}, 'rationale', {});

    lib(end+1) = localExpr( ...
        'log_integrated_intensity', ...
        'log10(integrated_intensity + 1)', ...
        {'integrated_intensity'}, ...
        'Stabilizes heavy-tailed spot intensity distributions.');
    lib(end+1) = localExpr( ...
        'log_background', ...
        'log10(background + 1)', ...
        {'background'}, ...
        'Compresses local background dynamic range.');
    lib(end+1) = localExpr( ...
        'log_snr_proxy', ...
        'log10((integrated_intensity + 1) / (background + 1))', ...
        {'integrated_intensity', 'background'}, ...
        'Background-normalized signal proxy using log scaling.');
    lib(end+1) = localExpr( ...
        'shot_noise_proxy', ...
        '(integrated_intensity - background) / sqrt(background + 1)', ...
        {'integrated_intensity', 'background'}, ...
        'Approximate Poisson-dominated contrast score.');
    lib(end+1) = localExpr( ...
        'background_fraction', ...
        'background / (integrated_intensity + background + 1)', ...
        {'integrated_intensity', 'background'}, ...
        'Fraction of local intensity explained by background.');

    lib(end+1) = localExpr( ...
        'quality_weighted_snr', ...
        'r_squared * log10((integrated_intensity + 1) / (background + 1))', ...
        {'r_squared', 'integrated_intensity', 'background'}, ...
        'Couples fit quality and normalized signal strength.');
    lib(end+1) = localExpr( ...
        'quality_penalty', ...
        '(1 - r_squared) * log10(integrated_intensity + 1)', ...
        {'r_squared', 'integrated_intensity'}, ...
        'Penalizes bright but poorly fit detections.');

    lib(end+1) = localExpr( ...
        'ellipticity_log', ...
        'abs(log((sigma_x + 1e-3) / (sigma_y + 1e-3)))', ...
        {'sigma_x', 'sigma_y'}, ...
        'Shape anisotropy magnitude in XY.');
    lib(end+1) = localExpr( ...
        'compactness_2d', ...
        'integrated_intensity / (sigma_x * sigma_y + 1e-3)', ...
        {'integrated_intensity', 'sigma_x', 'sigma_y'}, ...
        'Intensity concentration per estimated PSF area.');
    lib(end+1) = localExpr( ...
        'size_penalized_snr_2d', ...
        'log10((integrated_intensity + 1) / (background + 1)) / (sigma_x + sigma_y + 1e-3)', ...
        {'integrated_intensity', 'background', 'sigma_x', 'sigma_y'}, ...
        'Normalized signal penalized by broad PSF widths.');
    lib(end+1) = localExpr( ...
        'quality_per_sigma_sum_2d', ...
        'r_squared / (sigma_x + sigma_y + 1e-3)', ...
        {'r_squared', 'sigma_x', 'sigma_y'}, ...
        'Fit quality normalized by spot footprint.');

    lib(end+1) = localExpr( ...
        'compactness_3d', ...
        'integrated_intensity / (sigma_x * sigma_y * sigma_z + 1e-3)', ...
        {'integrated_intensity', 'sigma_x', 'sigma_y', 'sigma_z'}, ...
        'Intensity concentration per estimated PSF volume.');
    lib(end+1) = localExpr( ...
        'z_anisotropy_abs', ...
        'abs(log((sigma_z + 1e-3) / ((sigma_x + sigma_y) / 2 + 1e-3)))', ...
        {'sigma_x', 'sigma_y', 'sigma_z'}, ...
        'Axial-vs-lateral anisotropy magnitude.');
    lib(end+1) = localExpr( ...
        'size_penalized_snr_3d', ...
        'log10((integrated_intensity + 1) / (background + 1)) / (sigma_x + sigma_y + sigma_z + 1e-3)', ...
        {'integrated_intensity', 'background', 'sigma_x', 'sigma_y', 'sigma_z'}, ...
        '3D normalized signal penalized by total width.');

    lib(end+1) = localExpr( ...
        'amp_bg_log', ...
        'log10(amplitude + 1) - log10(background + 1)', ...
        {'amplitude', 'background'}, ...
        'Log contrast of fitted amplitude against local background.');
    lib(end+1) = localExpr( ...
        'amp_to_intensity_ratio', ...
        'amplitude / (integrated_intensity + 1)', ...
        {'amplitude', 'integrated_intensity'}, ...
        'Peak-to-total intensity ratio.');
    lib(end+1) = localExpr( ...
        'amp_quality_coupled', ...
        'r_squared * (amplitude / (background + 1))', ...
        {'r_squared', 'amplitude', 'background'}, ...
        'High when amplitude contrast and fit quality are both strong.');

    lib(end+1) = localExpr( ...
        'axis_amplitude_imbalance', ...
        '(abs(amplitude_x - amplitude_y) + abs(amplitude_x - amplitude_z) + abs(amplitude_y - amplitude_z)) / (amplitude_x + amplitude_y + amplitude_z + 1)', ...
        {'amplitude_x', 'amplitude_y', 'amplitude_z'}, ...
        'Axis-wise 1D amplitude imbalance score.');
    lib(end+1) = localExpr( ...
        'axis_amplitude_bg_ratio', ...
        '((amplitude_x + amplitude_y + amplitude_z) / 3) / (background + 1)', ...
        {'amplitude_x', 'amplitude_y', 'amplitude_z', 'background'}, ...
        'Mean axis amplitude relative to background.');

    lib(end+1) = localExpr( ...
        'xy_z_amplitude_ratio', ...
        'amplitude_xy / (amplitude_z + 1e-3)', ...
        {'amplitude_xy', 'amplitude_z'}, ...
        '2D XY vs 1D Z amplitude balance.');
    lib(end+1) = localExpr( ...
        'xy_z_amplitude_balance_log', ...
        'abs(log((amplitude_xy + 1e-3) / (amplitude_z + 1e-3)))', ...
        {'amplitude_xy', 'amplitude_z'}, ...
        'Magnitude of XY-vs-Z amplitude imbalance.');

    lib(end+1) = localExpr( ...
        'distortion_energy', ...
        'rho_xy^2 + rho_xz^2 + rho_yz^2', ...
        {'rho_xy', 'rho_xz', 'rho_yz'}, ...
        'Correlation distortion energy for distorted Gaussian fits.');
    lib(end+1) = localExpr( ...
        'distortion_penalized_quality', ...
        'r_squared / (1 + rho_xy^2 + rho_xz^2 + rho_yz^2)', ...
        {'r_squared', 'rho_xy', 'rho_xz', 'rho_yz'}, ...
        'Fit quality penalized by strong axis correlations.');

    lib(end+1) = localExpr( ...
        'radial_symmetry_size_norm', ...
        'radial_symmetry_score / (sigma_x + sigma_y + 1e-3)', ...
        {'radial_symmetry_score', 'sigma_x', 'sigma_y'}, ...
        'Radial symmetry quality normalized by spot size.');
end

function entry = localExpr(name, expression, requiredFeatures, rationale)
    entry = struct();
    entry.name = name;
    entry.expression = expression;
    entry.requiredFeatures = requiredFeatures;
    entry.rationale = rationale;
end
