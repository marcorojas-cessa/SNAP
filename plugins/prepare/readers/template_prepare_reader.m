function spec = template_prepare_reader()
% template_prepare_reader - External SNAP_prepare reader provider template.
%
% To activate:
%   1) Copy/rename this file and function.
%   2) Implement canRead/readFcn.
%   3) Place in plugins/prepare/readers or SNAP_PREPARE_PROVIDER_PATH.

    spec = struct();
    spec.id = 'external.example.prepare_reader';
    spec.displayName = 'Example Prepare Reader';
    spec.priority = 200; % higher than built-in if needed
    spec.version = '1.0.0';
    spec.source = 'external-template';
    spec.canReadFcn = @canRead;
    spec.readFcn = @readLibrary;
end

function tf = canRead(filePath)
    tf = false;
    if ~(ischar(filePath) || isstring(filePath))
        return;
    end
    filePath = char(string(filePath));
    [~, ~, ext] = fileparts(filePath);
    tf = strcmpi(ext, '.yourformat');
end

function out = readLibrary(filePath, progressCb)
    % Replace with format-specific logic.
    if ~isempty(progressCb)
        progressCb(sprintf('Example reader invoked for %s', filePath));
    end

    out = struct();
    out.imageData = {};
    out.tableData = {};
    out.numChannels = 0;
end
