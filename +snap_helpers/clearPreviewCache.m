function clearPreviewCache(handles)
% Clears the preview cache when data changes to ensure consistency
    
    if isfield(handles, 'previewCache')
        handles = rmfield(handles, 'previewCache');
        guidata(handles.fig, handles);
    end
end
