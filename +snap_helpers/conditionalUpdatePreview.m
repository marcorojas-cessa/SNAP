function conditionalUpdatePreview(fig, content_type)
% Wrapper that calls redrawPreview for all previews showing a specific content type.
    handles = guidata(fig);
    
    for i = 1:5
        preview_content = handles.previewContentDrops(i).Value;
        if contains(preview_content, content_type, 'IgnoreCase', true)
            snap_helpers.redrawPreview(fig, i);
        end
    end
end