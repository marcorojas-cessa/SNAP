function [predictions, scores, confidence, classLabels] = applyClassifier(model, X, featureNames, modelFeatureNames, normParams)
% applyClassifier - Apply trained SVM classifier to new data with z-score normalization
%
% ============================================================================
% APPLY TRAINED MODEL WITH CONSISTENT NORMALIZATION
% ============================================================================
%
% This function applies a trained SVM model to classify new spots.
% If normalization parameters are provided, features are z-score normalized
% using the SAME parameters from training (critical for consistency).
%
% USAGE:
%   [pred, scores, conf] = snap_helpers.classification.applyClassifier(model, X)
%   [pred, scores, conf, labels] = snap_helpers.classification.applyClassifier(model, X, featureNames, modelFeatureNames, normParams)
%
% INPUTS:
%   model             - Trained SVM model (from trainClassifier)
%   X                 - [N x F] feature matrix to classify
%   featureNames      - (optional) Cell array of feature names for columns of X
%   modelFeatureNames - (optional) Cell array of feature names the model expects
%   normParams        - (optional) Struct with normalization parameters:
%                       .mu    - [1 x F] feature means from training
%                       .sigma - [1 x F] feature stds from training
%                       .standardized - logical, whether to apply normalization
%
% OUTPUTS:
%   predictions - [N x 1] predicted class labels (0 = Noise, 1 = Real)
%   scores      - [N x 2] SVM decision scores for each class
%   confidence  - [N x 1] confidence values (abs distance from decision boundary)
%   classLabels - Cell array: {'Noise', 'Real'} (class name mapping)
%
% ============================================================================

    % Initialize outputs
    classLabels = {'Noise', 'Real'};
    
    % Handle empty inputs
    if isempty(model)
        predictions = [];
        scores = [];
        confidence = [];
        warning('Empty model provided');
        return;
    end
    
    if isempty(X)
        predictions = [];
        scores = [];
        confidence = [];
        return;
    end
    
    % Handle optional normParams
    if nargin < 5
        normParams = struct();
    end
    
    % Reorder features if names provided
    if nargin >= 4 && ~isempty(featureNames) && ~isempty(modelFeatureNames)
        [X, alignmentOk] = alignFeatures(X, featureNames, modelFeatureNames);
        if ~alignmentOk
            predictions = nan(size(X, 1), 1);
            scores = nan(size(X, 1), 2);
            confidence = zeros(size(X, 1), 1);
            warning('Feature alignment failed');
            return;
        end
    end
    
    % Handle NaN values - mark as invalid
    validMask = all(~isnan(X), 2);
    nSamples = size(X, 1);
    
    % Initialize outputs
    predictions = nan(nSamples, 1);
    scores = nan(nSamples, 2);
    confidence = zeros(nSamples, 1);
    
    if ~any(validMask)
        warning('All samples have NaN features');
        return;
    end
    
    % Extract valid samples
    X_valid = X(validMask, :);
    
    % ========================================================================
    % Z-SCORE NORMALIZATION using training parameters
    % ========================================================================
    if isfield(normParams, 'standardized') && normParams.standardized && ...
       isfield(normParams, 'mu') && isfield(normParams, 'sigma') && ...
       ~isempty(normParams.mu) && ~isempty(normParams.sigma)
        
        % Verify dimensions match
        if numel(normParams.mu) == size(X_valid, 2) && numel(normParams.sigma) == size(X_valid, 2)
            X_valid = (X_valid - normParams.mu) ./ normParams.sigma;
        else
            warning('Normalization parameter dimension mismatch. Skipping normalization.');
        end
    end
    
    % Apply model to valid samples
    try
        [pred, score] = predict(model, X_valid);
        
        % Store results for valid samples
        predictions(validMask) = pred;
        if size(score, 2) >= 2
            scores(validMask, :) = score;
            % Confidence = distance from decision boundary
            confidence(validMask) = abs(score(:, 2) - score(:, 1)) / 2;
        else
            scores(validMask, 1) = score;
            confidence(validMask) = abs(score);
        end
        
    catch ME
        warning('Prediction failed: %s', ME.message);
        predictions = nan(nSamples, 1);
        scores = nan(nSamples, 2);
        confidence = zeros(nSamples, 1);
    end
end

% ============================================================================
% HELPER FUNCTION: Align feature columns to model's expected order
% ============================================================================
function [X_aligned, success] = alignFeatures(X, featureNames, modelFeatureNames)
    success = true;
    
    [isPresent, idx] = ismember(modelFeatureNames, featureNames);
    
    if ~all(isPresent)
        missingFeatures = modelFeatureNames(~isPresent);
        warning('Missing features for classification: %s', strjoin(missingFeatures, ', '));
        success = false;
        X_aligned = X;
        return;
    end
    
    X_aligned = X(:, idx);
end
