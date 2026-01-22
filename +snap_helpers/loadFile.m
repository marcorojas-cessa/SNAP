function loadFile(fig, fileType, varargin)
% Opens a file dialog, updates path, saves handles, and enables/disables 3D controls.
    handles = guidata(fig);
    
    [filename, filepath] = uigetfile({'*.tif;*.tiff;*.ome.tif', 'TIFF Files'}, 'Select Image File');
    if isequal(filename, 0), 
        return; 
    end
    fullPath = fullfile(filepath, filename);
    
    channelIdx = 0; % Default to 0, indicating not a standard channel
    
    switch fileType
        case 'dic'
            handles.dicPathText.Value = fullPath;
        case 'nuc'
            handles.nucPathText.Value = fullPath;
            try
                info = imfinfo(fullPath);
                is3D = numel(info) > 1;
                handles.nucZSpacingInput.Enable = matlab.lang.OnOffSwitchState(is3D);
                if ~is3D
                    handles.nucPreprocessModeDrop.Items = {'2D (Slice-by-slice)', 'On Z-Projection'};
                    handles.nucPreprocessModeDrop.Value = '2D (Slice-by-slice)';
                    handles.nucSegModeDrop.Items = {'2D (Slice-by-slice)', 'On Z-Projection'};
                    handles.nucSegModeDrop.Value = '2D (Slice-by-slice)';
                else
                    handles.nucPreprocessModeDrop.Items = {'3D', '2D (Slice-by-slice)', 'On Z-Projection'};
                    handles.nucSegModeDrop.Items = {'3D', '2D (Slice-by-slice)', 'On Z-Projection'};
                end
            catch ME
                warning('Could not read image info for Nuclei: %s', ME.message);
            end
        case 'channel'
            if ~isempty(varargin)
                channelIdx = varargin{1};
                handles.channelPathTexts(channelIdx).Value = fullPath;
            else
                error('Channel index not provided for channel file type.');
            end
        otherwise
            error('Unknown file type: %s', fileType);
    end
    
    % Update controls based on image dimension if a channel was loaded
    if channelIdx > 0
        try
            info = imfinfo(fullPath);
            is3D = numel(info) > 1;
            handles.zSpacingInputs(channelIdx).Enable = is3D;
            
            if is3D
                set(handles.preprocessModeDrops(channelIdx), 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'});
                set(handles.bgCorrModeDrops(channelIdx), 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'});
                set(handles.maximaModeDrops(channelIdx), 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'});
            else
                set(handles.preprocessModeDrops(channelIdx), 'Items', {'2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', '2D (Slice-by-slice)');
                set(handles.bgCorrModeDrops(channelIdx), 'Items', {'2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', '2D (Slice-by-slice)');
                set(handles.maximaModeDrops(channelIdx), 'Items', {'2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', '2D (Slice-by-slice)');
            end
        catch ME
            warning('Could not read image info for %s: %s', fullPath, ME.message);
        end
    end

    % Clear preview cache when new data is loaded
    snap_helpers.clearPreviewCache(handles);
    
    % Save the updated handles structure
    guidata(fig, handles);
    % No updateControls call here - let it happen when user clicks "Update Previews"
end