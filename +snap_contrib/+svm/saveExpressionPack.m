function saveExpressionPack(pack, outputPath)
% saveExpressionPack - Save an SVM expression pack to MAT and CSV.
%
% Usage:
%   pack = snap_contrib.svm.buildExpressionPack('ParameterFile', 'params.mat');
%   snap_contrib.svm.saveExpressionPack(pack, 'svm_expression_pack.mat');

    if nargin < 2 || isempty(outputPath)
        error('saveExpressionPack:MissingOutputPath', 'Provide outputPath.');
    end
    if ~isstruct(pack) || ~isfield(pack, 'channelPacks')
        error('saveExpressionPack:InvalidPack', ...
            'Input pack must be a struct from snap_contrib.svm.buildExpressionPack.');
    end

    outputPath = char(string(outputPath));
    [outDir, ~, outExt] = fileparts(outputPath);
    if isempty(outDir)
        outDir = pwd;
        outputPath = fullfile(outDir, outputPath);
    end
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    if isempty(outExt)
        outputPath = [outputPath '.mat'];
    end

    expressionPack = pack; %#ok<NASGU>
    save(outputPath, 'expressionPack');

    [baseDir, baseName, ~] = fileparts(outputPath);
    csvPath = fullfile(baseDir, [baseName '_expressions.csv']);
    fid = fopen(csvPath, 'w');
    if fid == -1
        warning('saveExpressionPack:CsvWriteFailed', ...
            'Could not open CSV path for writing: %s', csvPath);
        return;
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'channel_idx,fitting_method,has3d,expression_name,expression\n');
    for i = 1:numel(pack.channelPacks)
        chPack = pack.channelPacks(i);
        for j = 1:numel(chPack.customExpressions)
            exprName = localCsvEscape(chPack.customExpressions(j).name);
            exprText = localCsvEscape(chPack.customExpressions(j).expression);
            fitMethod = localCsvEscape(chPack.fittingMethod);
            fprintf(fid, '%d,%s,%d,%s,%s\n', ...
                chPack.channelIdx, fitMethod, logical(chPack.has3D), exprName, exprText);
        end
    end
end

function out = localCsvEscape(in)
    txt = char(string(in));
    txt = strrep(txt, '"', '""');
    out = ['"' txt '"'];
end
