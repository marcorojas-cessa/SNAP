function preprocessed_nuc = preprocessNuclei(handles)
% Applies the selected pre-processing steps to the raw nuclei image.

    if ~isfield(handles, 'rawNuclei') || isempty(handles.rawNuclei)
        preprocessed_nuc = [];
        return;
    end
    
    % Check if nuclei preprocessing is enabled
    if ~handles.nucPreprocEnabledCheck.Value
        % Return raw image if preprocessing is disabled
        preprocessed_nuc = double(handles.rawNuclei);
        return;
    end
    
    raw_img = double(handles.rawNuclei);
    
    % --- Pre-processing ---
    mode = handles.nucPreprocessModeDrop.Value;
    is_scaled = handles.nucPreprocessScaleCheck.Value;
    xy_spacing = handles.nucXYSpacingInput.Value;
    z_spacing = handles.nucZSpacingInput.Value;
    
    img_after_preprocess = [];
    
    if strcmp(mode, 'On Z-Projection')
        proj_type = handles.nucPreprocessProjectionDrop.Value;
        img_2d = projectZ(raw_img, proj_type);
        img_after_preprocess = apply_nuc_2d_preproc(img_2d, handles);
    else % 3D or 2D Slice-by-slice
        img_3d = raw_img;
        is_3d_mode = strcmp(mode, '3D');
        
        % Get the selected preprocessing method
        method = handles.nucPreprocMethodDrop.Value;
        
        % Apply preprocessing based on method
        switch method
            case 'Gaussian'
                sigma_val = handles.nucPreprocParam1Inputs.Value;
                if sigma_val > 0
                    if is_3d_mode
                        if is_scaled
                            sigma_pixels_xy = sigma_val / xy_spacing;
                            sigma_pixels_z = sigma_val / z_spacing;
                            sigma_vec = [sigma_pixels_xy, sigma_pixels_xy, sigma_pixels_z];
                            img_3d = imgaussfilt3(img_3d, sigma_vec);
                        else
                            img_3d = imgaussfilt3(img_3d, sigma_val);
                        end
                    else % Slice-by-slice (pixels)
                        for z = 1:size(img_3d, 3)
                            img_3d(:,:,z) = imgaussfilt(img_3d(:,:,z), sigma_val);
                        end
                    end
                end
        
            case 'Median'
                med_val = handles.nucPreprocParam1Inputs.Value;
                if med_val > 0
                    if is_3d_mode
                        if is_scaled
                            neigh_xy_f = med_val / xy_spacing;
                            neigh_z_f = med_val / z_spacing;
                            neigh_xy = max(1, 2*round(neigh_xy_f/2) + 1);
                            neigh_z = max(1, 2*round(neigh_z_f/2) + 1);
                            neigh_vec = [neigh_xy, neigh_xy, neigh_z];
                            img_3d = medfilt3(img_3d, neigh_vec);
                        else
                            med_val_odd = max(1, 2*round(med_val/2) + 1);
                            img_3d = medfilt3(img_3d, [med_val_odd med_val_odd med_val_odd]);
                        end
                    else % Slice-by-slice (pixels)
                        for z = 1:size(img_3d, 3)
                            med_val_odd = max(1, 2*round(med_val/2) + 1);
                            img_3d(:,:,z) = medfilt2(img_3d(:,:,z), [med_val_odd med_val_odd]);
                        end
                    end
                end
                
            case 'Wavelet Denoising'
                % Only available in 2D modes
                if ~is_3d_mode
                    wname = handles.nucWaveletNameDrop.Value;
                    level = handles.nucPreprocParam2Inputs.Value;
                    rule = handles.nucPreprocParam3Inputs.Value;
                    method_str = handles.nucPreprocParam4Inputs.Value;
                    if strcmp(method_str, 'soft'), sorh = 's'; else, sorh = 'h'; end
                    
                    if level > 0 && floor(level) == level
                        for z = 1:size(img_3d, 3)
                            slice = img_3d(:,:,z);
                            try
                                % Use wdencmp for 2D image denoising
                                thr = thselect(slice, rule);
                                if isscalar(thr) && isfinite(thr)
                                    denoised_slice = wdencmp('gbl', slice, wname, level, thr, sorh, 1);
                                    if isnumeric(denoised_slice) && ~isempty(denoised_slice) && isreal(denoised_slice) && all(isfinite(denoised_slice(:)))
                                        img_3d(:,:,z) = double(denoised_slice);
                                    end
                                end
                            catch ME
                                warning('Wavelet denoising failed for nuclei Z-frame %d: %s. Using original slice.', z, ME.message);
                            end
                        end
                    end
                end
                
            case 'None'
                % No preprocessing applied
                % No break needed - fall through to end
        end
        
        img_after_preprocess = img_3d;
    end
    
    preprocessed_nuc = img_after_preprocess;
end

% --- Nested Helper Functions ---

function img_out = apply_nuc_2d_preproc(img_in, h)
    img_out = img_in;
    
    % Get the selected preprocessing method
    method = h.nucPreprocMethodDrop.Value;
    
    % Apply preprocessing based on method
    switch method
        case 'Gaussian'
            sigma_val = h.nucPreprocParam1Inputs.Value;
            if sigma_val > 0
                img_out = imgaussfilt(img_out, sigma_val);
            end
            
        case 'Median'
            med_val = h.nucPreprocParam1Inputs.Value;
            if med_val > 0
                med_val_odd = max(1, 2*round(med_val/2) + 1);
                img_out = medfilt2(img_out, [med_val_odd med_val_odd]);
            end
            
        case 'Wavelet Denoising'
            wname = h.nucWaveletNameDrop.Value;
            level = h.nucPreprocParam2Inputs.Value;
            rule = h.nucPreprocParam3Inputs.Value;
            method_str = h.nucPreprocParam4Inputs.Value;
            if strcmp(method_str, 'soft'), sorh = 's'; else, sorh = 'h'; end
            
            if level > 0 && floor(level) == level
                try
                    % Use wdencmp for 2D image denoising
                    thr = thselect(img_out, rule);
                    if isscalar(thr) && isfinite(thr)
                        denoised_img = wdencmp('gbl', img_out, wname, level, thr, sorh, 1);
                        if isnumeric(denoised_img) && ~isempty(denoised_img) && isreal(denoised_img) && all(isfinite(denoised_img(:)))
                            img_out = double(denoised_img);
                        end
                    end
                catch ME
                    warning('Wavelet denoising failed for nuclei 2D projection: %s. Using original image.', ME.message);
                end
            end
            
        case 'None'
            % No preprocessing applied
            % No break needed - fall through to end
    end
end

function proj = projectZ(img, type)
    if ndims(img) < 3, proj = img; return; end
    if strcmp(type, 'Max'), proj = max(img, [], 3);
    elseif strcmp(type, 'Min'), proj = min(img, [], 3);
    elseif strcmp(type, 'Median'), proj = median(img, 3);
    elseif strcmp(type, 'Mean'), proj = mean(img, 3);
    else, proj = max(img, [], 3);
    end
end
