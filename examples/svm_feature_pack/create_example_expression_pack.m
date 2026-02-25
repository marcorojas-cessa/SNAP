function pack = create_example_expression_pack(parameterFile, outputPath, varargin)
% create_example_expression_pack - Build and save the example SVM pack.
%
% Example:
%   pack = create_example_expression_pack('svm_parameters.mat', ...
%       'examples/svm_feature_pack/snap_svm_expression_pack.mat');
%
% Optional name-value options:
%   'SelectionMode'                - 'focused' (default) or 'all_non_position'
%   'IncludeAdvancedStatsFeatures' - true (default) to include AIC/BIC/residual features
%   'LintMode'                     - 'auto' (default), 'strict', or 'permissive'
%                                    Note: lint mode validates compatibility only; it does not train SVM.

    if nargin < 1 || isempty(parameterFile)
        parameterFile = fullfile(pwd, 'svm_parameters.mat');
    end
    if nargin < 2 || isempty(outputPath)
        outputPath = fullfile(pwd, 'examples', 'svm_feature_pack', ...
            'snap_svm_expression_pack.mat');
    end

    p = inputParser;
    p.addParameter('IncludeAdvancedStatsFeatures', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('SelectionMode', 'focused', @(x) ischar(x) || isstring(x));
    p.addParameter('LintMode', 'auto', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;
    includeAdvanced = logical(opts.IncludeAdvancedStatsFeatures);
    lintModeInput = lower(strtrim(char(string(opts.LintMode))));
    if ~ismember(lintModeInput, {'auto', 'strict', 'permissive'})
        error('create_example_expression_pack:InvalidLintMode', ...
            'LintMode must be ''auto'', ''strict'', or ''permissive''.');
    end

    pack = snap_contrib.svm.buildExpressionPack( ...
        'ParameterFile', parameterFile, ...
        'IncludeAdvancedStatsFeatures', includeAdvanced, ...
        'SelectionMode', opts.SelectionMode);

    snap_contrib.svm.saveExpressionPack(pack, outputPath);

    % Advanced model-stats features are context-dependent and can be pruned
    % by capability checks when unavailable, so lint permissively by default
    % in that mode. Non-advanced mode remains strict.
    lintMode = lintModeInput;
    if strcmp(lintMode, 'auto')
        lintMode = 'strict';
        if includeAdvanced
            lintMode = 'permissive';
        end
    end
    autoGuard = strcmp(lintMode, 'permissive');

    lintReport = snap_contrib.svm.lintExpressionPack(outputPath, ...
        'ParameterFile', parameterFile, ...
        'Mode', lintMode, ...
        'AutoGuardUnsafeExpressions', autoGuard, ...
        'Verbose', false);
    if ~lintReport.success
        error('create_example_expression_pack:LintFailed', ...
            'Generated expression pack failed compatibility lint checks.');
    end

    fprintf('Saved example expression pack: %s\n', outputPath);
    fprintf('Channels in pack: %s\n', mat2str(pack.channelsIncluded));
    fprintf('Selection mode: %s\n', char(string(opts.SelectionMode)));
    fprintf('Include advanced model-stat features: %s\n', string(includeAdvanced));
    fprintf('Lint mode (validation only): %s\n', lintMode);
    for i = 1:numel(pack.channelPacks)
        chPack = pack.channelPacks(i);
        fprintf('  Ch%d: base features=%d, custom expressions=%d\n', ...
            chPack.channelIdx, numel(chPack.selectedFeatures), numel(chPack.customExpressions));
    end
end
