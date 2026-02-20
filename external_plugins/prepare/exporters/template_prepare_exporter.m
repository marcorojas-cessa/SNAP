function spec = template_prepare_exporter()
% template_prepare_exporter - External SNAP_prepare exporter template.
%
% To activate:
%   1) Copy/rename this file and function.
%   2) Implement canExport/exportFcn.
%   3) Place in external_plugins/prepare/exporters or SNAP_EXTERNAL_PREPARE_EXPORTER_PATH.

    spec = struct();
    spec.id = 'external.example.prepare_exporter';
    spec.displayName = 'Example Prepare Exporter';
    spec.priority = 200; % higher than built-in if needed
    spec.version = '1.0.0';
    spec.source = 'external-template';
    spec.canExportFcn = @canExport;
    spec.exportFcn = @exportChannels;
end

function tf = canExport(imageInfo, mapping, outputPath)
    tf = isstruct(imageInfo) && isstruct(mapping) && ...
         (ischar(outputPath) || isstring(outputPath));
end

function exportChannels(imageInfo, mapping, outputPath, progressCb)
    %#ok<INUSD>
    % Replace with custom export behavior.
    if ~isempty(progressCb)
        progressCb('Example exporter ran (no files written).');
    end
end
