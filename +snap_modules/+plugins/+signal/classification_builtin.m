function spec = classification_builtin()
% Built-in SVM classification plugin for batch/signal pipelines.

    spec = struct();
    spec.id = 'builtin.signal.classification';
    spec.stage = 'classification';
    spec.displayName = 'SNAP SVM Classification';
    spec.priority = 100;
    spec.version = '1.0.0';
    spec.source = 'builtin';
    spec.supportsFcn = [];
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = isEnabled(~, context)
    tf = false;
    if ~isfield(context, 'enableClassification') || ~logical(context.enableClassification)
        return;
    end
    if ~isfield(context, 'params') || ~isstruct(context.params)
        return;
    end

    params = context.params;
    ch = context.channelIdx;
    tf = getCellLogical(params, 'classifyEnabled', ch, false);
end

function state = runPlugin(state, context)
    if isempty(state.fitResults)
        return;
    end

    params = context.params;
    ch = context.channelIdx;

    if ~hasCellEntry(params, 'classifiers', ch)
        snap_modules.signal.emitProgress(context, '[classification] Channel %d has no loaded classifier. Skipping.', ch);
        return;
    end

    classifier = params.classifiers{ch};
    if isempty(classifier)
        snap_modules.signal.emitProgress(context, '[classification] Channel %d classifier entry is empty. Skipping.', ch);
        return;
    end

    if ~hasCellEntry(params, 'classifierFeatures', ch) || ~hasCellEntry(params, 'classifierFeatureInfo', ch)
        snap_modules.signal.emitProgress(context, '[classification] Channel %d classifier metadata missing. Skipping.', ch);
        return;
    end

    features = params.classifierFeatures{ch};
    featureInfo = params.classifierFeatureInfo{ch};

    customExpr = struct('name', {}, 'expression', {});
    if hasCellEntry(params, 'classifierCustomExpressions', ch) && ~isempty(params.classifierCustomExpressions{ch})
        customExpr = params.classifierCustomExpressions{ch};
    end

    normParams = struct('mu', [], 'sigma', [], 'standardized', false);
    if hasCellEntry(params, 'classifierNormParams', ch) && ~isempty(params.classifierNormParams{ch})
        normParams = params.classifierNormParams{ch};
    end

    [X, featureNames, validMask] = snap_helpers.classification.buildFeatureMatrix( ...
        state.fitResults, features, featureInfo, customExpr);

    [predictions, scores, confidence, classLabels] = snap_helpers.classification.applyClassifier( ...
        classifier, X, featureNames, featureNames, normParams);

    state.classification.predictions = predictions;
    state.classification.scores = scores;
    state.classification.confidence = confidence;
    state.classification.classLabels = classLabels;
    state.classification.validMask = validMask;

    keepMask = true(size(predictions));
    filterNoise = getCellLogical(params, 'classifyFilterNoise', ch, true);
    if filterNoise
        keepMask = (predictions == 1) | ~validMask;
        state.classification.keepMask = keepMask;

        nBefore = numel(state.fitResults);
        state.fitResults = state.fitResults(keepMask);

        if ~isempty(state.maximaCoords) && size(state.maximaCoords, 1) == numel(keepMask)
            state.maximaCoords = state.maximaCoords(keepMask, :);
        end

        nAfter = numel(state.fitResults);
        snap_modules.signal.emitProgress(context, ...
            '[classification] Channel %d kept %d/%d fit(s) after noise filtering.', ...
            ch, nAfter, nBefore);
    else
        snap_modules.signal.emitProgress(context, ...
            '[classification] Channel %d classified %d fit(s) (no filtering applied).', ...
            ch, numel(predictions));
    end
end

function tf = hasCellEntry(params, fieldName, idx)
    tf = false;
    if ~isstruct(params) || ~isfield(params, fieldName)
        return;
    end
    raw = params.(fieldName);
    if ~iscell(raw) || numel(raw) < idx
        return;
    end
    tf = true;
end

function value = getCellLogical(params, fieldName, idx, defaultValue)
    value = defaultValue;
    if ~isstruct(params) || ~isfield(params, fieldName)
        return;
    end
    raw = params.(fieldName);
    try
        if iscell(raw)
            if numel(raw) < idx || isempty(raw{idx})
                return;
            end
            value = logical(raw{idx});
        elseif isnumeric(raw) || islogical(raw)
            if numel(raw) < idx
                return;
            end
            value = logical(raw(idx));
        end
    catch
        value = defaultValue;
    end
end
