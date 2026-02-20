function [fitResultsOut, summary] = augmentFitResultsWithModelStats(fitResultsIn, varargin)
% augmentFitResultsWithModelStats - Add model-selection/statistical features.
%
% This example contribution computes optional SVM features from each fit's
% local intensity window (if present), including:
%   - Gaussian-vs-flat AIC/BIC deltas
%   - Poisson NLL deltas
%   - Residual distribution descriptors
%
% Intended usage:
%   fitData = snap_contrib.svm.augmentFitResultsWithModelStats(fitData);
%
% Notes:
%   - Works best when fitResults include rawDataWindow from fitGaussians.
%   - Uses a simplified Gaussian reconstruction from fitted parameters.
%   - Keeps behavior non-destructive: missing/invalid cases are NaN.
%
% Name-value options:
%   'ComputeNormalityPValue' - run residual KS normality p-value (default: false)
%   'Verbose'                - print summary (default: false)

    p = inputParser;
    p.addParameter('ComputeNormalityPValue', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('Verbose', false, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    if isempty(fitResultsIn)
        fitResultsOut = fitResultsIn;
        summary = struct('nTotal', 0, 'nComputed', 0, 'nMissingWindow', 0, 'nModelFailures', 0);
        return;
    end
    if ~isstruct(fitResultsIn)
        error('augmentFitResultsWithModelStats:InvalidInput', ...
            'fitResultsIn must be a struct array.');
    end

    fitResultsOut = fitResultsIn;
    n = numel(fitResultsOut);
    nComputed = 0;
    nMissingWindow = 0;
    nModelFailures = 0;

    for i = 1:n
        fitResultsOut(i) = localInitializeFields(fitResultsOut(i));

        rawWindow = [];
        if isfield(fitResultsOut(i), 'rawDataWindow') && ~isempty(fitResultsOut(i).rawDataWindow)
            rawWindow = fitResultsOut(i).rawDataWindow;
        elseif isfield(fitResultsOut(i), 'raw_data_window') && ~isempty(fitResultsOut(i).raw_data_window)
            rawWindow = fitResultsOut(i).raw_data_window;
        end

        if isempty(rawWindow)
            nMissingWindow = nMissingWindow + 1;
            continue;
        end

        rawWindow = double(rawWindow);
        if ~isnumeric(rawWindow) || isempty(rawWindow)
            nMissingWindow = nMissingWindow + 1;
            continue;
        end

        [predictedRaw, modelOk, kGaussian] = localReconstructGaussianModel(rawWindow, fitResultsOut(i));
        if ~modelOk
            nModelFailures = nModelFailures + 1;
            continue;
        end

        metrics = localComputeModelStats(rawWindow, predictedRaw, kGaussian, logical(opts.ComputeNormalityPValue));
        fitResultsOut(i) = localAssignMetrics(fitResultsOut(i), metrics);
        nComputed = nComputed + 1;
    end

    summary = struct();
    summary.nTotal = n;
    summary.nComputed = nComputed;
    summary.nMissingWindow = nMissingWindow;
    summary.nModelFailures = nModelFailures;

    if logical(opts.Verbose)
        fprintf(['augmentFitResultsWithModelStats: computed=%d, missingWindow=%d, ', ...
                 'modelFailures=%d, total=%d\n'], ...
            nComputed, nMissingWindow, nModelFailures, n);
    end
end

function fr = localInitializeFields(fr)
    % Model-selection scores
    fr.aic_gaussian = NaN;
    fr.bic_gaussian = NaN;
    fr.aic_flat = NaN;
    fr.bic_flat = NaN;
    fr.delta_aic_flat_minus_gaussian = NaN;
    fr.delta_bic_flat_minus_gaussian = NaN;
    fr.nll_poisson_gaussian = NaN;
    fr.nll_poisson_flat = NaN;
    fr.delta_nll_poisson_flat_minus_gaussian = NaN;
    fr.lrt_stat_flat_vs_gaussian = NaN;

    % Residual descriptors
    fr.residual_mean = NaN;
    fr.residual_std = NaN;
    fr.residual_mad = NaN;
    fr.residual_skewness = NaN;
    fr.residual_kurtosis_excess = NaN;
    fr.residual_over_signal = NaN;
    fr.normality_ks_pvalue = NaN;
end

function fr = localAssignMetrics(fr, m)
    names = fieldnames(m);
    for j = 1:numel(names)
        fr.(names{j}) = m.(names{j});
    end
end

function [predictedRaw, ok, kGaussian] = localReconstructGaussianModel(rawWindow, fitResult)
    ok = false;
    predictedRaw = [];
    kGaussian = NaN;

    dims = size(rawWindow);
    is3D = ndims(rawWindow) == 3 && dims(3) > 1;

    amp = localGetAmplitude(fitResult);
    cx = localGetFieldNumeric(fitResult, 'center_x');
    cy = localGetFieldNumeric(fitResult, 'center_y');
    cz = localGetFieldNumeric(fitResult, 'center_z');
    sx = localGetFieldNumeric(fitResult, 'sigma_x');
    sy = localGetFieldNumeric(fitResult, 'sigma_y');
    sz = localGetFieldNumeric(fitResult, 'sigma_z');
    bg = localGetFieldNumeric(fitResult, 'background');
    if ~isfinite(bg)
        bg = median(rawWindow(:), 'omitnan');
        if ~isfinite(bg)
            bg = 0;
        end
    end

    if ~(isfinite(amp) && isfinite(cx) && isfinite(cy) && isfinite(sx) && isfinite(sy) && sx > 0 && sy > 0)
        return;
    end

    if is3D && isfinite(cz) && isfinite(sz) && sz > 0
        [COL, ROW, SLICE] = meshgrid(1:dims(2), 1:dims(1), 1:dims(3));
        rho_xy = localGetFieldNumeric(fitResult, 'rho_xy');
        rho_xz = localGetFieldNumeric(fitResult, 'rho_xz');
        rho_yz = localGetFieldNumeric(fitResult, 'rho_yz');

        useDistorted3D = all(isfinite([rho_xy, rho_xz, rho_yz]));
        if useDistorted3D
            Sigma = [sx^2, rho_xy * sx * sy, rho_xz * sx * sz; ...
                     rho_xy * sx * sy, sy^2, rho_yz * sy * sz; ...
                     rho_xz * sx * sz, rho_yz * sy * sz, sz^2];
            [~, cholFlag] = chol(Sigma);
            if cholFlag == 0
                xyz = [ROW(:), COL(:), SLICE(:)];
                mu = [cx, cy, cz];
                centered = xyz - mu;
                gaussianValues = amp .* exp(-0.5 .* sum((centered / Sigma) .* centered, 2));
                gaussianModel = reshape(gaussianValues, dims);
                kGaussian = 11; % A, centers(3), sigmas(3), rho(3), background
            else
                gaussianModel = amp .* exp(-((ROW - cx).^2 ./ (2 * sx^2) + ...
                                             (COL - cy).^2 ./ (2 * sy^2) + ...
                                             (SLICE - cz).^2 ./ (2 * sz^2)));
                kGaussian = 8; % A, centers(3), sigmas(3), background
            end
        else
            gaussianModel = amp .* exp(-((ROW - cx).^2 ./ (2 * sx^2) + ...
                                         (COL - cy).^2 ./ (2 * sy^2) + ...
                                         (SLICE - cz).^2 ./ (2 * sz^2)));
            kGaussian = 8;
        end
        predictedRaw = gaussianModel + bg;
        ok = true;
        return;
    end

    [COL, ROW] = meshgrid(1:dims(2), 1:dims(1));
    theta = localGetFieldNumeric(fitResult, 'rho_xy');
    if isfinite(theta) && abs(theta) <= pi
        xCentered = ROW - cx;
        yCentered = COL - cy;
        xRot = xCentered .* cos(theta) + yCentered .* sin(theta);
        yRot = -xCentered .* sin(theta) + yCentered .* cos(theta);
        gaussianModel = amp .* exp(-(xRot.^2 ./ (2 * sx^2) + yRot.^2 ./ (2 * sy^2)));
    else
        gaussianModel = amp .* exp(-((ROW - cx).^2 ./ (2 * sx^2) + ...
                                     (COL - cy).^2 ./ (2 * sy^2)));
    end
    predictedRaw = gaussianModel + bg;
    kGaussian = 6; % A, centers(2), sigmas(2), background
    ok = true;
end

function metrics = localComputeModelStats(rawWindow, predictedRaw, kGaussian, computeNormalityP)
    obs = double(rawWindow(:));
    pred = double(predictedRaw(:));
    n = max(1, numel(obs));
    epsv = 1e-12;

    residual = obs - pred;
    sseGaussian = sum(residual.^2, 'omitnan');
    if ~isfinite(sseGaussian) || sseGaussian <= 0
        sseGaussian = epsv;
    end

    flatMean = mean(obs, 'omitnan');
    if ~isfinite(flatMean)
        flatMean = 0;
    end
    residualFlat = obs - flatMean;
    sseFlat = sum(residualFlat.^2, 'omitnan');
    if ~isfinite(sseFlat) || sseFlat <= 0
        sseFlat = epsv;
    end

    sigma2Gaussian = max(sseGaussian / n, epsv);
    sigma2Flat = max(sseFlat / n, epsv);
    llGaussian = -0.5 * n * (log(2 * pi * sigma2Gaussian) + 1);
    llFlat = -0.5 * n * (log(2 * pi * sigma2Flat) + 1);

    kFlat = 1;
    aicGaussian = 2 * kGaussian - 2 * llGaussian;
    bicGaussian = kGaussian * log(n) - 2 * llGaussian;
    aicFlat = 2 * kFlat - 2 * llFlat;
    bicFlat = kFlat * log(n) - 2 * llFlat;

    obsCounts = max(obs, 0);
    lambdaGaussian = max(pred, epsv);
    lambdaFlat = mean(obsCounts, 'omitnan');
    if ~isfinite(lambdaFlat) || lambdaFlat <= 0
        lambdaFlat = epsv;
    end
    nllPoissonGaussian = sum(lambdaGaussian - obsCounts .* log(lambdaGaussian) + gammaln(obsCounts + 1), 'omitnan');
    nllPoissonFlat = sum(lambdaFlat - obsCounts .* log(lambdaFlat) + gammaln(obsCounts + 1), 'omitnan');

    residualMean = mean(residual, 'omitnan');
    residualStd = std(residual, 0, 'omitnan');
    residualMed = median(residual, 'omitnan');
    residualMad = median(abs(residual - residualMed), 'omitnan');

    if isfinite(residualStd) && residualStd > 0
        z = (residual - residualMean) ./ residualStd;
        residualSkewness = mean(z.^3, 'omitnan');
        residualKurtExcess = mean(z.^4, 'omitnan') - 3;
    else
        residualSkewness = NaN;
        residualKurtExcess = NaN;
    end

    predMedian = median(pred, 'omitnan');
    if ~isfinite(predMedian)
        predMedian = 0;
    end
    ampScale = max(pred - predMedian, [], 'omitnan');
    if ~isfinite(ampScale) || ampScale <= 0
        ampScale = 1;
    end
    residualOverSignal = residualStd / ampScale;

    normalityP = NaN;
    if computeNormalityP && exist('kstest', 'file') == 2 && isfinite(residualStd) && residualStd > 0
        try
            z = (residual - residualMean) ./ residualStd;
            z = z(isfinite(z));
            if isempty(z)
                normalityP = NaN;
            else
                [~, normalityP] = kstest(z);
            end
        catch
            normalityP = NaN;
        end
    end

    metrics = struct();
    metrics.aic_gaussian = aicGaussian;
    metrics.bic_gaussian = bicGaussian;
    metrics.aic_flat = aicFlat;
    metrics.bic_flat = bicFlat;
    metrics.delta_aic_flat_minus_gaussian = aicFlat - aicGaussian;
    metrics.delta_bic_flat_minus_gaussian = bicFlat - bicGaussian;
    metrics.nll_poisson_gaussian = nllPoissonGaussian;
    metrics.nll_poisson_flat = nllPoissonFlat;
    metrics.delta_nll_poisson_flat_minus_gaussian = nllPoissonFlat - nllPoissonGaussian;
    metrics.lrt_stat_flat_vs_gaussian = 2 * max(0, nllPoissonFlat - nllPoissonGaussian);
    metrics.residual_mean = residualMean;
    metrics.residual_std = residualStd;
    metrics.residual_mad = residualMad;
    metrics.residual_skewness = residualSkewness;
    metrics.residual_kurtosis_excess = residualKurtExcess;
    metrics.residual_over_signal = residualOverSignal;
    metrics.normality_ks_pvalue = normalityP;
end

function amplitude = localGetAmplitude(fitResult)
    amplitude = localGetFieldNumeric(fitResult, 'amplitude');
    if isfinite(amplitude)
        return;
    end

    amplitude = localGetFieldNumeric(fitResult, 'amplitude_xy');
    if isfinite(amplitude)
        return;
    end

    amps = [ ...
        localGetFieldNumeric(fitResult, 'amplitude_x'), ...
        localGetFieldNumeric(fitResult, 'amplitude_y'), ...
        localGetFieldNumeric(fitResult, 'amplitude_z')];
    amps = amps(isfinite(amps));
    if isempty(amps)
        amplitude = NaN;
    else
        amplitude = mean(amps);
    end
end

function value = localGetFieldNumeric(s, fieldName)
    value = NaN;
    if ~isstruct(s) || ~isfield(s, fieldName)
        return;
    end

    raw = s.(fieldName);
    if isempty(raw)
        return;
    end
    if isnumeric(raw) || islogical(raw)
        if isscalar(raw) && isfinite(double(raw))
            value = double(raw);
        end
    elseif ischar(raw) || isstring(raw)
        parsed = str2double(string(raw));
        if isfinite(parsed)
            value = parsed;
        end
    end
end
