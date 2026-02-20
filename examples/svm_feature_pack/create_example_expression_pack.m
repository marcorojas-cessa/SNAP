function pack = create_example_expression_pack(parameterFile, outputPath)
% create_example_expression_pack - Build and save the example SVM pack.
%
% Example:
%   pack = create_example_expression_pack('svm_parameters.mat', ...
%       'examples/svm_feature_pack/snap_svm_expression_pack.mat');

    if nargin < 1 || isempty(parameterFile)
        parameterFile = fullfile(pwd, 'svm_parameters.mat');
    end
    if nargin < 2 || isempty(outputPath)
        outputPath = fullfile(pwd, 'examples', 'svm_feature_pack', ...
            'snap_svm_expression_pack.mat');
    end

    pack = snap_contrib.svm.buildExpressionPack( ...
        'ParameterFile', parameterFile, ...
        'IncludeAdvancedStatsFeatures', false);

    snap_contrib.svm.saveExpressionPack(pack, outputPath);

    fprintf('Saved example expression pack: %s\n', outputPath);
    fprintf('Channels in pack: %s\n', mat2str(pack.channelsIncluded));
    for i = 1:numel(pack.channelPacks)
        chPack = pack.channelPacks(i);
        fprintf('  Ch%d: base features=%d, custom expressions=%d\n', ...
            chPack.channelIdx, numel(chPack.selectedFeatures), numel(chPack.customExpressions));
    end
end
