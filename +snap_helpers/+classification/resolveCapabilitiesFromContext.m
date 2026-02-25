function capability = resolveCapabilitiesFromContext(fittingMethod, has3D, varargin)
% resolveCapabilitiesFromContext - Resolve feature capabilities for one channel context.
%
% Usage:
%   cap = snap_helpers.classification.resolveCapabilitiesFromContext('3D Gaussian', true)
%
% Name-value options:
%   'ChannelIndex'      - channel index metadata (default: 1)
%   'HasPhysicalSpacing'- expose physical-unit features (default: false)

    p = inputParser;
    p.addParameter('ChannelIndex', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
    p.addParameter('HasPhysicalSpacing', false, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    fm = strtrim(char(string(fittingMethod)));
    if isempty(fm)
        fm = '3D Gaussian';
    end

    if isempty(has3D)
        h3d = true;
    else
        h3d = logical(has3D);
    end

    hasPhysicalSpacing = logical(opts.HasPhysicalSpacing);
    [availableFeatures, featureInfo] = snap_helpers.classification.getAvailableFeatures( ...
        fm, h3d, hasPhysicalSpacing);

    capability = struct();
    capability.channelIdx = max(1, round(opts.ChannelIndex));
    capability.fittingMethod = fm;
    capability.has3D = h3d;
    capability.hasPhysicalSpacing = hasPhysicalSpacing;
    capability.availableFeatures = availableFeatures(:)';
    capability.featureInfo = featureInfo;
    capability.requiredCapabilities = struct( ...
        'fittingMethod', fm, ...
        'has3D', h3d, ...
        'hasPhysicalSpacing', hasPhysicalSpacing);
end
