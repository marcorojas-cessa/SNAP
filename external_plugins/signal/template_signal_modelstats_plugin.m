function spec = template_signal_modelstats_plugin()
% template_signal_modelstats_plugin - Example plugin adding fit model stats.
%
% This template demonstrates how to extend SNAP's modular signal pipeline
% without editing core app files. Because the function name contains
% "template", SNAP's external plugin discovery skips it by default.
%
% To activate:
%   1) Copy this file and rename both file + function.
%   2) Keep stage='gaussian_fitting' and set priority > 100.
%   3) Place the renamed file in external_plugins/signal or SNAP_EXTERNAL_SIGNAL_PLUGIN_PATH.

    spec = struct();
    spec.id = 'external.example.signal_modelstats';
    spec.stage = 'gaussian_fitting';
    spec.displayName = 'Example Gaussian Fit + Model Stats';
    spec.priority = 200; % override built-in gaussian fitting plugin
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
    % Example: enable only for training mode by default.
    tf = isfield(context, 'mode') && strcmpi(char(string(context.mode)), 'train');
end

function state = runPlugin(state, context)
    % First execute the canonical built-in Gaussian fitting behavior.
    builtin = snap_modules.plugins.signal.gaussian_fitting_builtin();
    state = builtin.run(state, context);

    % Then add optional model-selection/statistical fit features.
    if isfield(state, 'fitResults') && ~isempty(state.fitResults)
        [state.fitResults, summary] = snap_contrib.svm.augmentFitResultsWithModelStats( ...
            state.fitResults, 'ComputeNormalityPValue', false, 'Verbose', false);
        snap_modules.signal.emitProgress(context, ...
            '[gaussian_fitting] Model-stats extension computed for %d/%d fit(s).', ...
            summary.nComputed, summary.nTotal);
    end
end
