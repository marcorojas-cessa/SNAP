function [model, trainStats, normParams] = trainClassifier(X, labels, options)
% trainClassifier - Train SVM classifier with z-score normalization and configurable hyperparameters
%
% ============================================================================
% ROBUST SVM TRAINING WITH Z-SCORE NORMALIZATION
% ============================================================================
%
% This function trains a Support Vector Machine classifier for spot
% classification (Real vs Noise). Features are ALWAYS z-score normalized.
%
% KEY FEATURES:
%   - Z-score normalization (mean=0, std=1) for all features
%   - Automatic class balancing (handles imbalanced datasets)
%   - Configurable kernel and hyperparameters
%   - Cross-validation support
%   - Comprehensive training statistics
%
% USAGE:
%   [model, stats, norm] = snap_helpers.classification.trainClassifier(X, labels)
%   [model, stats, norm] = snap_helpers.classification.trainClassifier(X, labels, options)
%
% INPUTS:
%   X       - [N x F] feature matrix (N samples, F features)
%   labels  - [N x 1] class labels: 1 = Real spot, 0 = Noise
%   options - (optional) Struct with training parameters:
%
%     KERNEL OPTIONS:
%       .kernelFunction  - 'rbf' (default), 'linear', 'polynomial', 'gaussian'
%       .kernelScale     - 'auto' (default), or numeric value for RBF gamma
%       .polynomialOrder - 2-5 for polynomial kernel (default: 3)
%
%     REGULARIZATION:
%       .boxConstraint   - C parameter: higher = less regularization (default: 1)
%                          Range: 0.001 to 1000 typical
%
%     NORMALIZATION:
%       .standardize     - true (default): z-score normalize features
%                          ALWAYS RECOMMENDED for SVM
%
%     VALIDATION:
%       .crossValidate   - true (default) or false
%       .kFold           - Number of CV folds (default: 5)
%
%     OTHER:
%       .verbose         - true (default) or false
%       .solver          - 'SMO' (default), 'ISDA', 'L1QP'
%       .cacheSize       - Cache size in MB (default: 1000)
%
% OUTPUTS:
%   model      - Trained SVM model (compact form for storage efficiency)
%   trainStats - Struct with training statistics
%   normParams - Struct with normalization parameters:
%                .mu  - [1 x F] feature means
%                .sigma - [1 x F] feature standard deviations
%                These are needed to normalize new data before prediction
%
% ============================================================================

    % Default options
    if nargin < 3
        options = struct();
    end
    
    defaultOptions = struct();
    % Kernel options
    defaultOptions.kernelFunction = 'rbf';
    defaultOptions.kernelScale = 'auto';
    defaultOptions.polynomialOrder = 3;
    % Regularization
    defaultOptions.boxConstraint = 1;
    % Normalization - ALWAYS ON BY DEFAULT
    defaultOptions.standardize = true;
    % Validation
    defaultOptions.crossValidate = true;
    defaultOptions.kFold = 5;
    % Other
    defaultOptions.verbose = true;
    defaultOptions.solver = 'SMO';
    defaultOptions.cacheSize = 1000;
    
    % Merge with provided options
    optFields = fieldnames(defaultOptions);
    for i = 1:numel(optFields)
        if ~isfield(options, optFields{i})
            options.(optFields{i}) = defaultOptions.(optFields{i});
        end
    end
    
    % Initialize outputs
    trainStats = struct();
    trainStats.success = false;
    trainStats.error = '';
    trainStats.options = options;
    trainStats.trainTime = datetime('now');
    
    normParams = struct();
    normParams.mu = [];
    normParams.sigma = [];
    normParams.standardized = options.standardize;
    
    % Validate inputs
    if isempty(X) || isempty(labels)
        trainStats.error = 'Empty training data';
        model = [];
        return;
    end
    
    % Ensure labels is a column vector
    labels = labels(:);
    
    % Convert labels to numeric if necessary
    [labels, classNames] = normalizeLabels(labels);
    
    % Validate dimensions
    if size(X, 1) ~= numel(labels)
        trainStats.error = sprintf('Dimension mismatch: X has %d rows but labels has %d elements', ...
            size(X, 1), numel(labels));
        model = [];
        return;
    end
    
    % Remove samples with NaN features
    validMask = all(~isnan(X), 2) & ~isnan(labels);
    X = X(validMask, :);
    labels = labels(validMask);
    
    if isempty(X)
        trainStats.error = 'No valid samples after removing NaN values';
        model = [];
        return;
    end
    
    % ========================================================================
    % Z-SCORE NORMALIZATION
    % ========================================================================
    % Compute normalization parameters from training data
    normParams.mu = mean(X, 1, 'omitnan');
    normParams.sigma = std(X, 0, 1, 'omitnan');
    
    % Handle zero/near-zero variance features
    zeroVarIdx = normParams.sigma < 1e-10;
    if any(zeroVarIdx)
        normParams.sigma(zeroVarIdx) = 1;  % Prevent division by zero
        if options.verbose
            fprintf('  Warning: %d features have zero variance\n', sum(zeroVarIdx));
        end
    end
    
    % Apply z-score normalization if enabled
    if options.standardize
        X_normalized = (X - normParams.mu) ./ normParams.sigma;
    else
        X_normalized = X;
    end
    
    % Store basic statistics
    trainStats.nSamples = size(X, 1);
    trainStats.nFeatures = size(X, 2);
    trainStats.classNames = classNames;
    trainStats.featureMeans = normParams.mu;
    trainStats.featureStds = normParams.sigma;
    
    % Compute class statistics
    uniqueClasses = unique(labels);
    nClasses = numel(uniqueClasses);
    classCounts = zeros(nClasses, 1);
    for i = 1:nClasses
        classCounts(i) = sum(labels == uniqueClasses(i));
    end
    trainStats.classCounts = classCounts;
    
    % Validate minimum samples
    minSamplesPerClass = min(classCounts);
    if minSamplesPerClass < 2
        trainStats.error = sprintf('Insufficient samples: need at least 2 per class, have %d', minSamplesPerClass);
        model = [];
        return;
    end
    
    % Compute class weights (inverse frequency for balanced training)
    classWeights = trainStats.nSamples ./ (nClasses * classCounts);
    trainStats.classWeights = classWeights;
    
    % Create sample weight vector
    sampleWeights = ones(trainStats.nSamples, 1);
    for i = 1:nClasses
        sampleWeights(labels == uniqueClasses(i)) = classWeights(i);
    end
    
    if options.verbose
        fprintf('Training SVM classifier:\n');
        fprintf('  Samples: %d (%d Real, %d Noise)\n', trainStats.nSamples, classCounts(end), classCounts(1));
        fprintf('  Features: %d\n', trainStats.nFeatures);
        fprintf('  Z-score normalization: %s\n', string(options.standardize));
        fprintf('  Kernel: %s\n', options.kernelFunction);
        if strcmp(options.kernelFunction, 'polynomial')
            fprintf('  Polynomial order: %d\n', options.polynomialOrder);
        end
        if isnumeric(options.kernelScale)
            fprintf('  Kernel scale: %.4f\n', options.kernelScale);
        else
            fprintf('  Kernel scale: %s\n', options.kernelScale);
        end
        fprintf('  Box constraint (C): %.4f\n', options.boxConstraint);
        fprintf('  Class weights: [%.2f, %.2f]\n', classWeights(1), classWeights(end));
    end
    
    % Build SVM arguments
    svmArgs = {...
        'KernelFunction', options.kernelFunction, ...
        'BoxConstraint', options.boxConstraint, ...
        'Standardize', false, ...  % We handle normalization ourselves
        'Weights', sampleWeights, ...
        'Solver', options.solver, ...
        'CacheSize', options.cacheSize};
    
    % Add kernel-specific parameters
    if strcmp(options.kernelFunction, 'polynomial')
        svmArgs = [svmArgs, {'PolynomialOrder', options.polynomialOrder}];
    end
    
    if ~strcmp(options.kernelScale, 'auto') && isnumeric(options.kernelScale)
        svmArgs = [svmArgs, {'KernelScale', options.kernelScale}];
    elseif strcmp(options.kernelScale, 'auto')
        svmArgs = [svmArgs, {'KernelScale', 'auto'}];
    end
    
    % Train SVM
    try
        if options.crossValidate && trainStats.nSamples >= options.kFold * 2
            % Cross-validated training
            if options.verbose
                fprintf('  Performing %d-fold cross-validation...\n', options.kFold);
            end
            
            cvModel = fitcsvm(X_normalized, labels, svmArgs{:}, ...
                'CrossVal', 'on', 'KFold', options.kFold);
            
            % Get cross-validation loss
            trainStats.cvLoss = kfoldLoss(cvModel);
            trainStats.cvAccuracy = 1 - trainStats.cvLoss;
            trainStats.crossValidated = true;
            
            % Get CV predictions for confusion matrix
            cvPredictions = kfoldPredict(cvModel);
            trainStats.confusionMatrix = confusionmat(labels, cvPredictions);
            
            % Train final model on all data
            finalModel = fitcsvm(X_normalized, labels, svmArgs{:});
            
            if options.verbose
                fprintf('  Cross-validation accuracy: %.1f%%\n', trainStats.cvAccuracy * 100);
            end
            
        else
            % Direct training (no CV)
            if options.verbose
                fprintf('  Training without cross-validation...\n');
            end
            
            finalModel = fitcsvm(X_normalized, labels, svmArgs{:});
            
            % Training accuracy
            [predLabels, ~] = predict(finalModel, X_normalized);
            trainStats.trainAccuracy = mean(predLabels == labels);
            trainStats.crossValidated = false;
            
            % Confusion matrix
            trainStats.confusionMatrix = confusionmat(labels, predLabels);
            
            if options.verbose
                fprintf('  Training accuracy: %.1f%%\n', trainStats.trainAccuracy * 100);
            end
        end
        
        % Compact model for storage
        model = compact(finalModel);
        
        % Calculate precision, recall, F1 for each class
        cm = trainStats.confusionMatrix;
        if size(cm, 1) >= 2
            % For binary classification: row=actual, col=predicted
            TN = cm(1,1); FP = cm(1,2);
            FN = cm(2,1); TP = cm(2,2);
            
            trainStats.precision_real = TP / max(TP + FP, 1);
            trainStats.recall_real = TP / max(TP + FN, 1);
            trainStats.f1_real = 2 * trainStats.precision_real * trainStats.recall_real / ...
                max(trainStats.precision_real + trainStats.recall_real, 0.001);
            
            trainStats.precision_noise = TN / max(TN + FN, 1);
            trainStats.recall_noise = TN / max(TN + FP, 1);
            trainStats.f1_noise = 2 * trainStats.precision_noise * trainStats.recall_noise / ...
                max(trainStats.precision_noise + trainStats.recall_noise, 0.001);
        end
        
        trainStats.success = true;
        
        if options.verbose
            fprintf('  Training complete!\n');
            if isfield(trainStats, 'f1_real')
                fprintf('  Real spots - Precision: %.1f%%, Recall: %.1f%%, F1: %.2f\n', ...
                    trainStats.precision_real*100, trainStats.recall_real*100, trainStats.f1_real);
                fprintf('  Noise - Precision: %.1f%%, Recall: %.1f%%, F1: %.2f\n', ...
                    trainStats.precision_noise*100, trainStats.recall_noise*100, trainStats.f1_noise);
            end
        end
        
    catch ME
        model = [];
        trainStats.success = false;
        trainStats.error = ME.message;
        
        if options.verbose
            fprintf('  ERROR: %s\n', ME.message);
        end
    end
end

% ============================================================================
% HELPER FUNCTION: Normalize labels to numeric format
% ============================================================================
function [numericLabels, classNames] = normalizeLabels(labels)
    classNames = {'Noise', 'Real'};
    
    if isnumeric(labels)
        if all(labels == 0 | labels == 1)
            numericLabels = labels;
        elseif all(labels == -1 | labels == 1)
            numericLabels = (labels + 1) / 2;
        else
            uLabels = unique(labels);
            numericLabels = double(labels == uLabels(end));
        end
        
    elseif iscategorical(labels)
        cats = categories(labels);
        classNames = cats;
        numericLabels = double(labels == cats{end});
        
    elseif iscell(labels) || isstring(labels)
        labels = string(labels);
        uLabels = unique(labels);
        classNames = cellstr(uLabels);
        
        realIdx = find(contains(lower(uLabels), 'real') | contains(lower(uLabels), 'true') | ...
                      contains(lower(uLabels), 'positive') | contains(lower(uLabels), 'spot'));
        if isempty(realIdx)
            realIdx = numel(uLabels);
        end
        
        numericLabels = double(labels == uLabels(realIdx));
        
    else
        error('Unsupported label type: %s', class(labels));
    end
    
    numericLabels = numericLabels(:);
end
