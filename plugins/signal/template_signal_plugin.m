function spec = template_signal_plugin()
% template_signal_plugin - External signal plugin template for SNAP.
%
% To activate this plugin:
%   1) Copy this file and rename function/file to a unique name.
%   2) Update id/stage/priority and implement the run callback.
%   3) Place it in plugins/signal or a folder listed in SNAP_SIGNAL_PLUGIN_PATH.
%
% Stage options (built-in order):
%   signal_processing, maxima_detection, gaussian_fitting,
%   fit_filtering, classification, nuclei_filtering
%
% Priority: highest value wins per stage.

    spec = struct();
    spec.id = 'external.example.signal_plugin';
    spec.stage = 'maxima_detection';
    spec.displayName = 'Example Signal Plugin';
    spec.priority = 200; % override built-ins by using >100
    spec.version = '1.0.0';
    spec.source = 'external-template';
    spec.supportsFcn = @supports;
    spec.isEnabledFcn = @isEnabled;
    spec.run = @runPlugin;
end

function tf = supports(~, ~)
    tf = true;
end

function tf = isEnabled(~, context)
    tf = isfield(context, 'mode') && strcmpi(context.mode, 'batch');
end

function state = runPlugin(state, context)
    % Example: keep built-in maxima by doing nothing.
    % Replace this with your algorithm implementation.
    snap_modules.signal.emitProgress(context, '[template] External plugin executed.');
end
