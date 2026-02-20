function emitProgress(context, msgFmt, varargin)
% emitProgress - Safe progress emission helper for signal pipeline plugins.

    if nargin < 2
        return;
    end

    if ~isstruct(context) || ~isfield(context, 'progressCallback') || ...
            isempty(context.progressCallback) || ~isa(context.progressCallback, 'function_handle')
        return;
    end

    try
        context.progressCallback(sprintf(msgFmt, varargin{:}));
    catch
        % Progress reporting must never fail the pipeline.
    end
end
