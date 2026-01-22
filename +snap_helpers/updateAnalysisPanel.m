function updateAnalysisPanel(handles)
    % Updates the analysis panel with nuclei signal composition data
    
    try
        % Check if analysis panel exists
        if ~isfield(handles, 'analysisTable') || ~isvalid(handles.analysisTable)
            return;
        end
        
        % Get analysis data from cache
        analysis_data = [];
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'analysis')
            analysis_data = handles.previewCache.analysis;
        end
        
        % Determine active channels (default to empty if no data)
        active_channels = [];
        if ~isempty(analysis_data) && isfield(analysis_data, 'active_channels')
            active_channels = analysis_data.active_channels;
        end
        
        % Build column names dynamically based on active channels
        % Fixed columns: Nucleus, Area/Vol, Circularity/Sphericity, Solidity
        column_names = {'Nucleus', 'Area/Vol', 'Circularity/Sphericity', 'Solidity'};
        col_width = {60, 80, 120, 70}; % Widths for fixed columns
        
        % Add columns only for active channels
        for ch = active_channels
            column_names{end+1} = sprintf('# Ch%d Signals', ch);
            col_width{end+1} = 80;
        end
        
        % Show empty table if no data
        if isempty(analysis_data) || ~isfield(analysis_data, 'nuclei_data') || isempty(analysis_data.nuclei_data) || analysis_data.num_nuclei == 0
            handles.analysisTable.Data = {};
            handles.analysisTable.ColumnName = column_names;
            handles.analysisTable.ColumnWidth = col_width;
            handles.analysisInfoLabel.Text = 'No nuclei data available. Enable segmentation and update previews to populate table.';
            return;
        end
        
        nuclei_data = analysis_data.nuclei_data;
        
        % ====================================================================
        % SHAPE METRICS FROM SHARED SOURCE (Consistency Fix)
        % ====================================================================
        % analyzeNucleiSignalComposition() no longer computes shape metrics
        % to avoid discrepancies. Get them from computeNucleusMeasurements()
        % which is the single source of truth for ALL exports.
        % ====================================================================
        
        % Get nucleus labels and mask from cache
        shape_metrics = [];
        if isfield(handles, 'previewCache') && isfield(handles.previewCache, 'nuclei')
            nucleus_labels = handles.previewCache.nuclei.labels;
            nuclei_mask = handles.previewCache.nuclei.mask;
            
            if ~isempty(nucleus_labels) && ~isempty(nuclei_mask)
                % Compute shape metrics using SINGLE SOURCE OF TRUTH
                % These will EXACTLY match export values
                measureOptions = struct();
                measureOptions.is_3d = (ndims(nuclei_mask) == 3 && size(nuclei_mask, 3) > 1);
                measureOptions.xy_spacing = handles.nucXYSpacingInput.Value;
                measureOptions.z_spacing = handles.nucZSpacingInput.Value;
                measureOptions.image_name = '';
                
                shape_metrics = snap_helpers.computeNucleusMeasurements(nucleus_labels, nuclei_mask, measureOptions);
            end
        end
        
        % Determine if 2D or 3D for area/volume display
        is_3d = ~isempty(shape_metrics) && shape_metrics.is_3d;
        if isempty(shape_metrics)
        is_3d = isfield(nuclei_data(1), 'volume_voxels');
        end
        
        % Build table data with dynamic number of columns
        num_columns = length(column_names);
        table_data = cell(length(nuclei_data), num_columns);
        
        for i = 1:length(nuclei_data)
            nucleus = nuclei_data(i);
            
            % Column 1: Nucleus ID (ensure integer display)
            table_data{i, 1} = int32(round(nucleus.nucleus_id));
            
            % Column 2: Area/Volume - Use shape_metrics if available for consistency
            if ~isempty(shape_metrics) && i <= length(shape_metrics.nucleus_ids)
                if shape_metrics.is_3d
                    table_data{i, 2} = int32(round(shape_metrics.volume_voxels(i)));
                else
                    table_data{i, 2} = int32(round(shape_metrics.area_pixels(i)));
                end
            elseif is_3d && isfield(nucleus, 'volume_voxels')
                table_data{i, 2} = int32(round(nucleus.volume_voxels));
            elseif ~is_3d && isfield(nucleus, 'area_pixels')
                table_data{i, 2} = int32(round(nucleus.area_pixels));
            else
                table_data{i, 2} = NaN;
            end
            
            % Column 3: Circularity/Sphericity - Use shape_metrics for consistency with exports
            if ~isempty(shape_metrics) && i <= length(shape_metrics.nucleus_ids)
                if shape_metrics.is_3d
                    table_data{i, 3} = round(shape_metrics.sphericity(i), 3);
                else
                    table_data{i, 3} = round(shape_metrics.circularity(i), 3);
                end
            else
                table_data{i, 3} = NaN;
            end
            
            % Column 4: Solidity - Use shape_metrics for consistency with exports
            if ~isempty(shape_metrics) && i <= length(shape_metrics.nucleus_ids)
                table_data{i, 4} = round(shape_metrics.solidity(i), 3);
            else
                table_data{i, 4} = NaN;
            end
            
            % Columns 5+: Signal counts for active channels only
            col_idx = 5;
            for ch = active_channels
                signal_count = 0;
                if isfield(nucleus.channels, sprintf('channel_%d', ch))
                    ch_field = sprintf('channel_%d', ch);
                    if isfield(nucleus.channels.(ch_field), 'signal_count')
                        signal_count = nucleus.channels.(ch_field).signal_count;
                    end
                end
                table_data{i, col_idx} = int32(round(signal_count));
                col_idx = col_idx + 1;
            end
        end
        
        % Update table
        handles.analysisTable.ColumnName = column_names;
        handles.analysisTable.Data = table_data;
        handles.analysisTable.ColumnWidth = col_width;
        
        % Update info label
        handles.analysisInfoLabel.Text = sprintf('%d nuclei analyzed with shape metrics and signal counts', ...
            length(nuclei_data));
        
        % Store data for export
        handles.analysisData = analysis_data;
        
    catch ME
        warning('Failed to update analysis panel: %s', ME.message);
        if isfield(handles, 'analysisTable') && isvalid(handles.analysisTable)
            handles.analysisTable.Data = {};
            handles.analysisInfoLabel.Text = 'Error updating analysis display.';
        end
    end
end
