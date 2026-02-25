function [packOut, report] = validateExpressionPackAgainstCapabilities(packIn, capabilities, varargin)
% validateExpressionPackAgainstCapabilities - Validate/sanitize pack against channel capabilities.
%
% Usage:
%   [packOut, report] = snap_helpers.classification.validateExpressionPackAgainstCapabilities(pack, caps)
%
% Name-value options:
%   'Mode'                       - 'auto' (default), 'permissive', or 'strict'
%   'AutoGuardUnsafeExpressions' - in permissive mode, try real(...) guard on complex outputs
%   'StressSampleCount'          - synthetic samples for expression stress-check (default: 8)

    p = inputParser;
    p.addParameter('Mode', 'auto', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoGuardUnsafeExpressions', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('StressSampleCount', 8, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 4);
    p.parse(varargin{:});
    opts = p.Results;

    [packOut, normalizeReport] = snap_helpers.classification.normalizeExpressionPack(packIn);
    capList = localNormalizeCapabilities(capabilities);

    if isempty(capList)
        error('validateExpressionPackAgainstCapabilities:MissingCapabilities', ...
            'Capabilities input is empty. Resolve capabilities from parameter context first.');
    end

    mode = lower(strtrim(char(string(opts.Mode))));
    if strcmp(mode, 'auto')
        if isfield(packOut, 'strictModeDefault') && logical(packOut.strictModeDefault)
            mode = 'strict';
        else
            mode = 'permissive';
        end
    end
    if ~ismember(mode, {'strict', 'permissive'})
        error('validateExpressionPackAgainstCapabilities:InvalidMode', ...
            'Mode must be "strict", "permissive", or "auto".');
    end

    isStrict = strcmp(mode, 'strict');
    doAutoGuard = logical(opts.AutoGuardUnsafeExpressions) && ~isStrict;
    nStress = max(4, round(opts.StressSampleCount));

    report = struct();
    report.mode = mode;
    report.success = true;
    report.errors = {};
    report.warnings = normalizeReport.warnings;
    report.channelReports = repmat(localEmptyChannelReport(), 1, numel(packOut.channelPacks));
    report.nDroppedBase = 0;
    report.nDroppedCustom = 0;
    report.nAutoGuarded = 0;

    for i = 1:numel(packOut.channelPacks)
        cp = packOut.channelPacks(i);
        ch = cp.channelIdx;
        chReport = localEmptyChannelReport();
        chReport.channelIdx = ch;

        cap = localFindCapability(capList, ch, i);
        if isempty(cap)
            msg = sprintf('[Channel %d] No capability context available for this channel.', ch);
            chReport.errors{end+1} = msg; %#ok<AGROW>
            report.errors{end+1} = msg; %#ok<AGROW>
            if ~isStrict
                report.warnings{end+1} = sprintf('%s Keeping pack entries unchanged.', msg); %#ok<AGROW>
            end
            report.channelReports(i) = chReport;
            continue;
        end

        available = cap.availableFeatures(:)';
        chReport.availableFeatures = available;
        chReport.fittingMethod = cap.fittingMethod;
        chReport.has3D = logical(cap.has3D);

        reqCaps = cp.requiredCapabilities;
        reqMsgs = localCheckRequiredCapabilities(reqCaps, cap, ch);
        if ~isempty(reqMsgs)
            if isStrict
                for k = 1:numel(reqMsgs)
                    report.errors{end+1} = reqMsgs{k}; %#ok<AGROW>
                    chReport.errors{end+1} = reqMsgs{k}; %#ok<AGROW>
                end
            else
                for k = 1:numel(reqMsgs)
                    report.warnings{end+1} = reqMsgs{k}; %#ok<AGROW>
                    chReport.warnings{end+1} = reqMsgs{k}; %#ok<AGROW>
                end
            end
        end

        selected = cp.selectedFeatures(:)';
        keepBaseMask = ismember(selected, available);
        droppedBase = selected(~keepBaseMask);
        cp.selectedFeatures = selected(keepBaseMask);
        cp.requiredFeatures = unique([cp.requiredFeatures(:)', cp.selectedFeatures], 'stable');

        if ~isempty(droppedBase)
            msg = sprintf('[Channel %d] Dropped incompatible base features: %s', ch, strjoin(droppedBase, ', '));
            chReport.droppedBase = droppedBase;
            report.nDroppedBase = report.nDroppedBase + numel(droppedBase);
            if isStrict
                report.errors{end+1} = msg; %#ok<AGROW>
                chReport.errors{end+1} = msg; %#ok<AGROW>
            else
                report.warnings{end+1} = msg; %#ok<AGROW>
                chReport.warnings{end+1} = msg; %#ok<AGROW>
            end
        end

        keptExpr = struct('name', {}, 'expression', {}, 'requiredFeatures', {}, 'requiredCapabilities', {});
        nextExprIdx = 0;

        for e = 1:numel(cp.customExpressions)
            exprEntry = cp.customExpressions(e);
            exprName = strtrim(char(string(exprEntry.name)));
            exprText = strtrim(char(string(exprEntry.expression)));
            reqFeatures = localNormalizeCellstr(localGetField(exprEntry, 'requiredFeatures', {}));
            if isempty(reqFeatures)
                reqFeatures = localInferExpressionRequiredFeatures(exprText);
            end

            missingReq = setdiff(reqFeatures, available);
            if ~isempty(missingReq)
                msg = sprintf('[Channel %d] Dropped custom expression "%s": missing required feature(s): %s', ...
                    ch, exprName, strjoin(missingReq, ', '));
                chReport.droppedCustom{end+1} = exprName; %#ok<AGROW>
                report.nDroppedCustom = report.nDroppedCustom + 1;
                if isStrict
                    report.errors{end+1} = msg; %#ok<AGROW>
                    chReport.errors{end+1} = msg; %#ok<AGROW>
                else
                    report.warnings{end+1} = msg; %#ok<AGROW>
                    chReport.warnings{end+1} = msg; %#ok<AGROW>
                end
                continue;
            end

            safetyWarnings = localAnalyzeExpressionSafety(exprText);
            for w = 1:numel(safetyWarnings)
                msg = sprintf('[Channel %d] Expression "%s": %s', ch, exprName, safetyWarnings{w});
                report.warnings{end+1} = msg; %#ok<AGROW>
                chReport.warnings{end+1} = msg; %#ok<AGROW>
            end

            evalInfo = localStressEvaluateExpression(exprText, available, nStress);
            if evalInfo.valid
                nextExprIdx = nextExprIdx + 1;
                keptExpr(nextExprIdx) = struct( ... %#ok<AGROW>
                    'name', exprName, ...
                    'expression', exprText, ...
                    'requiredFeatures', {reqFeatures}, ...
                    'requiredCapabilities', localGetField(exprEntry, 'requiredCapabilities', struct()));
                continue;
            end

            guardedApplied = false;
            guardedExpr = exprText;
            if doAutoGuard && evalInfo.hasComplex
                candidateExpr = sprintf('real(%s)', exprText);
                guardedEval = localStressEvaluateExpression(candidateExpr, available, nStress);
                if guardedEval.valid && ~guardedEval.hasComplex
                    guardedApplied = true;
                    guardedExpr = candidateExpr;
                    evalInfo = guardedEval;
                end
            end

            if guardedApplied
                nextExprIdx = nextExprIdx + 1;
                keptExpr(nextExprIdx) = struct( ... %#ok<AGROW>
                    'name', exprName, ...
                    'expression', guardedExpr, ...
                    'requiredFeatures', {reqFeatures}, ...
                    'requiredCapabilities', localGetField(exprEntry, 'requiredCapabilities', struct()));
                chReport.autoGuarded{end+1} = exprName; %#ok<AGROW>
                report.nAutoGuarded = report.nAutoGuarded + 1;
                msg = sprintf('[Channel %d] Auto-guarded expression "%s" with real(...).', ch, exprName);
                report.warnings{end+1} = msg; %#ok<AGROW>
                chReport.warnings{end+1} = msg; %#ok<AGROW>
                continue;
            end

            msg = sprintf('[Channel %d] Dropped custom expression "%s": %s', ch, exprName, evalInfo.reason);
            chReport.droppedCustom{end+1} = exprName; %#ok<AGROW>
            report.nDroppedCustom = report.nDroppedCustom + 1;
            if isStrict
                report.errors{end+1} = msg; %#ok<AGROW>
                chReport.errors{end+1} = msg; %#ok<AGROW>
            else
                report.warnings{end+1} = msg; %#ok<AGROW>
                chReport.warnings{end+1} = msg; %#ok<AGROW>
            end
        end

        cp.customExpressions = keptExpr;
        cp.requiredFeatures = unique([cp.selectedFeatures, localCollectRequiredFeatures(keptExpr)], 'stable');
        packOut.channelPacks(i) = cp;
        report.channelReports(i) = chReport;
    end

    if isStrict && ~isempty(report.errors)
        report.success = false;
    else
        report.success = true;
    end
end

function capList = localNormalizeCapabilities(capabilities)
    capList = struct('channelIdx', {}, 'fittingMethod', {}, 'has3D', {}, 'availableFeatures', {}, 'featureInfo', {});

    if isempty(capabilities) || ~isstruct(capabilities)
        return;
    end

    capList = capabilities(:)';
    for i = 1:numel(capList)
        if ~isfield(capList(i), 'channelIdx') || isempty(capList(i).channelIdx)
            capList(i).channelIdx = i;
        end

        if ~isfield(capList(i), 'fittingMethod') || isempty(capList(i).fittingMethod)
            capList(i).fittingMethod = '3D Gaussian';
        else
            capList(i).fittingMethod = char(string(capList(i).fittingMethod));
        end

        if ~isfield(capList(i), 'has3D') || isempty(capList(i).has3D)
            capList(i).has3D = true;
        else
            capList(i).has3D = logical(capList(i).has3D);
        end

        if ~isfield(capList(i), 'availableFeatures') || isempty(capList(i).availableFeatures)
            [features, info] = snap_helpers.classification.getAvailableFeatures( ...
                capList(i).fittingMethod, capList(i).has3D, false);
            capList(i).availableFeatures = features(:)';
            capList(i).featureInfo = info;
        else
            capList(i).availableFeatures = localNormalizeCellstr(capList(i).availableFeatures);
            if ~isfield(capList(i), 'featureInfo')
                capList(i).featureInfo = struct();
            end
        end
    end
end

function cap = localFindCapability(capList, channelIdx, fallbackIdx)
    cap = struct([]);
    if isempty(capList)
        return;
    end

    idx = find([capList.channelIdx] == channelIdx, 1, 'first');
    if ~isempty(idx)
        cap = capList(idx);
        return;
    end

    if numel(capList) == 1
        cap = capList(1);
        return;
    end

    if nargin >= 3 && fallbackIdx >= 1 && fallbackIdx <= numel(capList)
        cap = capList(fallbackIdx);
    end
end

function messages = localCheckRequiredCapabilities(reqCaps, cap, channelIdx)
    messages = {};
    if ~isstruct(reqCaps)
        return;
    end

    reqMethod = localNormalizeString(localGetField(reqCaps, 'fittingMethod', ''), '');
    if ~isempty(reqMethod)
        if ~strcmpi(reqMethod, cap.fittingMethod)
            messages{end+1} = sprintf([ ...
                '[Channel %d] Pack requires fittingMethod="%s" but parameter context is "%s".'], ...
                channelIdx, reqMethod, cap.fittingMethod); %#ok<AGROW>
        end
    end

    if isfield(reqCaps, 'has3D') && ~isempty(reqCaps.has3D)
        reqHas3D = logical(reqCaps.has3D);
        if reqHas3D ~= logical(cap.has3D)
            messages{end+1} = sprintf([ ...
                '[Channel %d] Pack requires has3D=%d but parameter context is has3D=%d.'], ...
                channelIdx, reqHas3D, logical(cap.has3D)); %#ok<AGROW>
        end
    end
end

function info = localStressEvaluateExpression(expression, availableFeatures, nSamples)
    info = struct();
    info.valid = false;
    info.hasComplex = false;
    info.reason = 'unknown evaluation failure';

    dummy = localBuildDummyFeatureStruct(availableFeatures, nSamples);

    warnState = warning('query', 'all');
    warning('off', 'all');
    cleanup = onCleanup(@() warning(warnState));
    result = snap_helpers.classification.evaluateExpression(expression, dummy, availableFeatures);
    clear cleanup;

    if isempty(result)
        info.reason = 'expression returned empty output';
        return;
    end

    result = result(:);
    if numel(result) ~= nSamples
        info.reason = sprintf('expression returned %d value(s), expected %d', numel(result), nSamples);
        return;
    end

    imagPart = imag(result);
    complexMask = isfinite(imagPart) & (abs(imagPart) > 1e-12);
    info.hasComplex = any(complexMask);
    if info.hasComplex
        info.reason = 'expression produced complex values under stress test';
        return;
    end

    realResult = real(result);
    finiteMask = isfinite(realResult);
    finiteCount = sum(finiteMask);

    if finiteCount == 0
        info.reason = 'expression produced no finite values under stress test';
        return;
    end

    info.valid = true;
    info.reason = '';
end

function warnings = localAnalyzeExpressionSafety(expression)
    warnings = {};
    expr = char(string(expression));
    exprLower = lower(expr);

    if ~isempty(regexp(exprLower, '(?<![a-z0-9_])log\s*\(', 'once')) || ...
            ~isempty(regexp(exprLower, '(?<![a-z0-9_])log10\s*\(', 'once')) || ...
            ~isempty(regexp(exprLower, '(?<![a-z0-9_])log2\s*\(', 'once'))
        if isempty(strfind(exprLower, 'max(')) %#ok<STREMP>
            warnings{end+1} = 'log/log10/log2 appears without an explicit max(...) guard.'; %#ok<AGROW>
        end
    end

    if ~isempty(regexp(exprLower, '(?<![a-z0-9_])sqrt\s*\(', 'once'))
        if isempty(strfind(exprLower, 'max(')) && isempty(strfind(exprLower, 'abs(')) %#ok<STREMP>
            warnings{end+1} = 'sqrt(...) appears without max(...) or abs(...) guard.'; %#ok<AGROW>
        end
    end

    if contains(expr, '/') && ~contains(exprLower, '+1') && ~contains(exprLower, '+ 1')
        warnings{end+1} = 'expression includes division; consider denominator guard (e.g., +1e-6 or +1).'; %#ok<AGROW>
    end
end

function data = localBuildDummyFeatureStruct(featureNames, n)
    data = repmat(struct(), n, 1);
    baseT = linspace(0, 1, n)';

    for i = 1:n
        for f = 1:numel(featureNames)
            name = featureNames{f};
            lowerName = lower(name);

            if contains(lowerName, 'rho_')
                value = -0.3 + 0.6 * baseT(i);
            elseif contains(lowerName, 'r_squared') || contains(lowerName, 'rsquared')
                value = 0.7 + 0.29 * baseT(i);
            elseif contains(lowerName, 'sigma')
                value = 0.8 + 1.4 * baseT(i);
            elseif contains(lowerName, 'background')
                value = 10 + 5 * baseT(i);
            elseif contains(lowerName, 'amplitude')
                value = 20 + 80 * baseT(i);
            elseif contains(lowerName, 'intensity')
                value = 100 + 400 * baseT(i);
            elseif strcmp(lowerName, 'z')
                value = 5 + 10 * baseT(i);
            else
                value = 1 + 4 * baseT(i);
            end

            data(i).(name) = value;
        end
    end
end

function req = localCollectRequiredFeatures(customExpr)
    req = {};
    for i = 1:numel(customExpr)
        if isfield(customExpr(i), 'requiredFeatures')
            req = [req, localNormalizeCellstr(customExpr(i).requiredFeatures)]; %#ok<AGROW>
        end
    end
    req = unique(req, 'stable');
end

function req = localInferExpressionRequiredFeatures(expr)
    req = {};
    if isempty(expr)
        return;
    end

    exprText = char(expr);
    % Strip scientific-notation exponent fragments (e.g., 1e-6) so "e"
    % is not interpreted as a required feature token.
    exprText = regexprep(exprText, '(?<=[0-9])[eEdD][+\-]?[0-9]+', ' ');

    tokens = regexp(exprText, '[A-Za-z_][A-Za-z0-9_]*', 'match');
    if isempty(tokens)
        return;
    end

    reserved = lower({ ...
        'log','log10','log2','sqrt','abs','exp','sin','cos','tan', ...
        'min','max','mean','std','median','sum','prod','sign','floor','ceil','round', ...
        'real','imag','conj','pi','inf','nan','true','false','workspace'});

    keep = true(size(tokens));
    for i = 1:numel(tokens)
        tok = lower(tokens{i});
        if strcmp(tok, 'e') || strcmp(tok, 'd')
            keep(i) = false;
        elseif ismember(tok, reserved)
            keep(i) = false;
        end
    end

    req = unique(tokens(keep), 'stable');
end

function out = localNormalizeCellstr(v)
    out = {};
    if isempty(v)
        return;
    end

    if ischar(v)
        out = {strtrim(v)};
        return;
    end

    if isstring(v)
        out = cellstr(v(:));
    elseif iscell(v)
        out = cell(size(v));
        for i = 1:numel(v)
            out{i} = char(string(v{i}));
        end
    else
        out = {char(string(v))};
    end

    out = out(:)';
    out = cellfun(@(s) strtrim(char(string(s))), out, 'UniformOutput', false);
    out = out(~cellfun(@isempty, out));
    out = unique(out, 'stable');
end

function value = localGetField(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function out = localNormalizeString(v, defaultValue)
    if nargin < 2
        defaultValue = '';
    end
    if isempty(v)
        out = defaultValue;
        return;
    end
    out = strtrim(char(string(v)));
    if isempty(out)
        out = defaultValue;
    end
end

function report = localEmptyChannelReport()
    report = struct( ...
        'channelIdx', NaN, ...
        'fittingMethod', '', ...
        'has3D', true, ...
        'availableFeatures', {{}}, ...
        'droppedBase', {{}}, ...
        'droppedCustom', {{}}, ...
        'autoGuarded', {{}}, ...
        'warnings', {{}}, ...
        'errors', {{}});
end
