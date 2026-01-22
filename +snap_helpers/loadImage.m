function imageData = loadImage(filePath)
%loadImage a utility to load image data from a file.
%   Supports standard image formats and .mat files.

[~, ~, ext] = fileparts(filePath);

if strcmpi(ext, '.mat')
    % For .mat files, assume the image data is in a variable
    % named 'imageData', 'img', or the first variable in the file.
    vars = load(filePath);
    if isfield(vars, 'imageData')
        imageData = vars.imageData;
    elseif isfield(vars, 'img')
        imageData = vars.img;
    else
        f = fieldnames(vars);
        imageData = vars.(f{1});
    end
else
    % For other formats, use imread
    imageData = imread(filePath);
end

end
