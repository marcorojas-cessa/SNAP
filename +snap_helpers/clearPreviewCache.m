function clearPreviewCache(inputArg)
% Clears the preview cache when data changes to ensure consistency

    handles = [];
    fig_handle = [];

    if isstruct(inputArg)
        handles = inputArg;
        if isfield(handles, 'fig')
            fig_handle = handles.fig;
        end
    elseif isa(inputArg, 'matlab.ui.Figure')
        fig_handle = inputArg;
        handles = guidata(fig_handle);
    else
        return;
    end

    if isempty(handles) || ~isstruct(handles)
        return;
    end

    if isfield(handles, 'previewCache')
        handles = rmfield(handles, 'previewCache');
        if ~isempty(fig_handle) && isvalid(fig_handle)
            guidata(fig_handle, handles);
        end
    end
end
