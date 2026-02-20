function fitParams = getFitParamsFromHandles(handles, channelIdx)
% getFitParamsFromHandles - Build fit parameter struct from handle-like data.

    fitParams = struct();
    fitParams.gaussFitVoxelWindowSize = getHandleValue(handles, 'gaussFitVoxelWindowSlider', channelIdx, 7);
    fitParams.gaussFitBgCorrMethod = getHandleValue(handles, 'gaussFitBgCorrMethodDrop', channelIdx, 'Mean Surrounding Subtraction');
    fitParams.gaussFitBgCorrWidth = getHandleValue(handles, 'gaussFitBgCorrWidthEdit', channelIdx, 2);
    fitParams.gaussFitPolyDegree = getHandleValue(handles, 'gaussFitPolyDegreeEdit', channelIdx, 2);
    fitParams.gaussFitMethod = getHandleValue(handles, 'gaussFitMethodDrop', channelIdx, '1D (X,Y,Z) Gaussian');
    fitParams.gaussFitMaxIterations = getHandleValue(handles, 'gaussFitMaxIterationsEdit', channelIdx, 200);
    fitParams.gaussFitTolerance = getHandleValue(handles, 'gaussFitToleranceEdit', channelIdx, 1e-6);
    fitParams.gaussFitRadialRadius = getHandleValue(handles, 'gaussFitRadialRadiusEdit', channelIdx, 3);
    fitParams.gaussFitPlotCheck = false;
    fitParams.xySpacing = getHandleValue(handles, 'xySpacingInputs', channelIdx, 1);
    fitParams.zSpacing = getHandleValue(handles, 'zSpacingInputs', channelIdx, 1);
end

function value = getHandleValue(handles, fieldName, idx, defaultValue)
    value = defaultValue;
    if ~isstruct(handles) || ~isfield(handles, fieldName)
        return;
    end

    raw = handles.(fieldName);
    try
        if numel(raw) < idx
            return;
        end
        entry = raw(idx);
        if isstruct(entry) && isfield(entry, 'Value')
            value = entry.Value;
        elseif isobject(entry) && isprop(entry, 'Value')
            value = entry.Value;
        end
    catch
        value = defaultValue;
    end
end
