function preprocessed_nuc = preprocessNucleiWithBgCorr(handles)
% Applies pre-processing and background correction to the raw nuclei image.

    if ~isfield(handles, 'rawNuclei') || isempty(handles.rawNuclei)
        preprocessed_nuc = [];
        return;
    end
    
    raw_img = double(handles.rawNuclei);
    
    img_after_deconv = raw_img; % Initialize with raw image
    
    % --- Step 0: Deconvolution (Applied First, Before Everything) ---
    if isfield(handles, 'nucDeconvEnabledCheck') && getScalar(handles.nucDeconvEnabledCheck.Value)
        method = handles.nucDeconvMethodDrop.Value;
        psf_source = handles.nucDeconvPSFSourceDrop.Value;
        
        % Get or generate PSF
        if strcmp(psf_source, 'Load File')
            psf_file = handles.nucDeconvPSFPathText.Value;
            if ~isempty(psf_file) && isfile(psf_file)
                psf_raw = double(tiffreadVolume(psf_file));
                % Crop/pad PSF to match image dimensions
                psf = cropPSFToImage(psf_raw, size(raw_img));
                psf = psf / sum(psf(:)); % Normalize
            else
                warning('Nuclei PSF file not found. Skipping deconvolution.');
                psf = [];
            end
        else % 'Generate'
            psf_sigma_xy = handles.nucDeconvPSFSigmaXYInput.Value;
            psf_sigma_z = handles.nucDeconvPSFSigmaZInput.Value;
            psf_size_xy = handles.nucDeconvPSFSizeXYInput.Value;
            psf_size_z = handles.nucDeconvPSFSizeZInput.Value;
            psf = generate3DPSF(psf_sigma_xy, psf_sigma_z, psf_size_xy, psf_size_z);
        end
        
        % Apply deconvolution if PSF is available
        if ~isempty(psf)
            switch method
                case 'Lucy-Richardson'
                    iterations = handles.nucDeconvLRIterationsInput.Value;
                    damping = handles.nucDeconvLRDampingInput.Value;
                    if damping > 0
                        img_after_deconv = deconvlucy(raw_img, psf, iterations, damping);
                    else
                        img_after_deconv = deconvlucy(raw_img, psf, iterations);
                    end
                    
                case 'Wiener'
                    nsr = handles.nucDeconvWienerNSRInput.Value;
                    img_after_deconv = deconvwnr(raw_img, psf, nsr);
                    
                case 'Blind'
                    iterations = handles.nucDeconvBlindIterationsInput.Value;
                    under_relax = handles.nucDeconvBlindUnderRelaxInput.Value;
                    if under_relax > 0
                        [img_after_deconv, ~] = deconvblind(raw_img, psf, iterations, under_relax);
                    else
                        [img_after_deconv, ~] = deconvblind(raw_img, psf, iterations);
                    end
            end
        end
    end
    
    % Check for abort after deconvolution
    if isfield(handles, 'fig') && isfield(handles, 'abortRequested')
        latest_handles = guidata(handles.fig);
        if latest_handles.abortRequested
            preprocessed_nuc = [];
            return;
        end
    end
    
    % --- Step 1: Pre-processing ---
    img_after_preprocess = img_after_deconv;
    
    if getScalar(handles.nucPreprocEnabledCheck.Value)
        mode = handles.nucPreprocessModeDrop.Value;
        is_scaled = getScalar(handles.nucPreprocessScaleCheck.Value);
        xy_spacing = handles.nucXYSpacingInput.Value;
        z_spacing = handles.nucZSpacingInput.Value;
        
        if strcmp(mode, 'On Z-Projection')
            proj_type = handles.nucPreprocessProjectionDrop.Value;
            img_2d = projectZ(raw_img, proj_type);
            img_after_preprocess = apply_nuc_2d_preproc(img_2d, handles);
        else % 3D or 2D Slice-by-slice
            img_3d = raw_img;
            is_3d_mode = strcmp(mode, '3D');
            method = handles.nucPreprocMethodDrop.Value;

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
                                % Check for abort every 10 slices
                                if mod(z, 10) == 0 && isfield(handles, 'fig') && isfield(handles, 'abortRequested')
                                    latest_handles = guidata(handles.fig);
                                    if latest_handles.abortRequested
                                        preprocessed_nuc = [];
                                        return;
                                    end
                                end
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
                case 'Non-Local Means'
                    % Get NLM parameters
                    filter_strength = handles.nucNlmFilterStrengthInput.Value;
                    search_window = handles.nucNlmSearchWindowInput.Value;
                    comparison_window = handles.nucNlmComparisonWindowInput.Value;
                    
                    if filter_strength > 0
                        if is_3d_mode
                            % Check if imnlmfilt3 is available (newer MATLAB versions)
                            if exist('imnlmfilt3', 'file') == 2
                                % 3D Non-Local Means (process entire volume)
                                img_3d = imnlmfilt3(img_3d, 'DegreeOfSmoothing', filter_strength, ...
                                                   'SearchWindowSize', search_window, ...
                                                   'ComparisonWindowSize', comparison_window);
                            else
                                % Fallback: Use slice-by-slice 2D processing for 3D mode
                                fprintf('Warning: imnlmfilt3 not available. Using slice-by-slice processing for 3D Non-Local Means.\n');
                                for z = 1:size(img_3d, 3)
                                    img_3d(:,:,z) = imnlmfilt(img_3d(:,:,z), 'DegreeOfSmoothing', filter_strength, ...
                                                              'SearchWindowSize', search_window, ...
                                                              'ComparisonWindowSize', comparison_window);
                                end
                            end
                        else
                            % Slice-by-slice Non-Local Means
                            for z = 1:size(img_3d, 3)
                                img_3d(:,:,z) = imnlmfilt(img_3d(:,:,z), 'DegreeOfSmoothing', filter_strength, ...
                                                          'SearchWindowSize', search_window, ...
                                                          'ComparisonWindowSize', comparison_window);
                            end
                        end
                    end
                case 'Wavelet Denoising'
                    if ~is_3d_mode
                        wname = handles.nucWaveletNameDrop.Value;
                        level = handles.nucPreprocParam2Inputs.Value;
                        rule = handles.nucPreprocParam3Inputs.Value;
                        method_str = handles.nucPreprocParam4Inputs.Value;
                        if strcmp(method_str, 'soft'), sorh = 's'; else, sorh = 'h'; end
                        
                        if level > 0 && floor(level) == level
                            for z = 1:size(img_3d, 3)
                                % Check for abort every 10 slices
                                if mod(z, 10) == 0 && isfield(handles, 'fig') && isfield(handles, 'abortRequested')
                                    latest_handles = guidata(handles.fig);
                                    if latest_handles.abortRequested
                                        preprocessed_nuc = [];
                                        return;
                                    end
                                end
                                
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
            end
            img_after_preprocess = img_3d;
        end
        
        % Enforce non-negativity after all pre-processing if checked
        if getScalar(handles.nucPreprocClipChecks.Value)
            img_after_preprocess(img_after_preprocess < 0) = 0;
        end
    end
    
    % Check for abort after preprocessing
    if isfield(handles, 'fig') && isfield(handles, 'abortRequested')
        latest_handles = guidata(handles.fig);
        if latest_handles.abortRequested
            preprocessed_nuc = [];
            return;
        end
    end
    
    % --- Step 2: Background Correction ---
    preprocessed_nuc = img_after_preprocess;
    
    if getScalar(handles.nucBgCorrEnabledCheck.Value)
        bg_mode = handles.nucBgCorrModeDrop.Value;
        is_scaled_bg = getScalar(handles.nucBgCorrScaleCheck.Value);
        xy_spacing = handles.nucXYSpacingInput.Value;
        z_spacing = handles.nucZSpacingInput.Value;

        if strcmp(bg_mode, 'On Z-Projection')
            proj_type = handles.nucBgCorrProjectionDrop.Value;
            img_2d = projectZ(img_after_preprocess, proj_type);
            preprocessed_nuc = apply_nuc_2d_bg_corr(img_2d, handles);
        else % 3D or 2D Slice-by-slice
            img_3d = img_after_preprocess;
            is_3d_mode_bg = strcmp(bg_mode, '3D');
            
            method = handles.nucBgMethodDrop.Value;
            param_val = handles.nucBgParamInput.Value;
            if param_val > 0 && ~strcmp(method, 'None')
                if strcmp(method, 'Gaussian')
                    sigma_val = param_val;
                    if is_3d_mode_bg
                        if is_scaled_bg
                            sigma_pixels_xy = sigma_val / xy_spacing;
                            sigma_pixels_z = sigma_val / z_spacing;
                            sigma_vec = [sigma_pixels_xy, sigma_pixels_xy, sigma_pixels_z];
                            if any(sigma_vec > 0)
                                bg = imgaussfilt3(img_3d, sigma_vec);
                                img_3d = img_3d - bg;
                            end
                        else % Unscaled 3D
                            bg = imgaussfilt3(img_3d, sigma_val);
                            img_3d = img_3d - bg;
                        end
                    else % Slice-by-slice
                        if is_scaled_bg
                            sigma_pixels_2d = sigma_val / xy_spacing;
                            for z = 1:size(img_3d, 3)
                                slice = img_3d(:,:,z);
                                bg_slice = imgaussfilt(slice, sigma_pixels_2d);
                                processed_slice = slice - bg_slice;
                                img_3d(:,:,z) = processed_slice;
                            end
                        else
                            for z = 1:size(img_3d, 3)
                                slice = img_3d(:,:,z);
                                bg_slice = imgaussfilt(slice, sigma_val);
                                processed_slice = slice - bg_slice;
                                img_3d(:,:,z) = processed_slice;
                            end
                        end
                    end
                else % Rolling-ball or Top-hat
                     if is_3d_mode_bg && is_scaled_bg
                        % True 3D Anisotropic Filtering using offsetstrel
                        radius_xy_pixels = round(param_val / xy_spacing);
                        height_z_pixels = round(param_val / z_spacing); % This is the Z-radius
                        
                        if radius_xy_pixels > 0 && height_z_pixels > 0
                            se = offsetstrel('ball', radius_xy_pixels, height_z_pixels);
                            if strcmp(method, 'Rolling-ball'), img_3d = img_3d - imopen(img_3d, se);
                            elseif strcmp(method, 'Top-hat'), img_3d = imtophat(img_3d, se);
                            end
                        end
                     else % Voxel-based or Slice-by-slice
                        radius = round(param_val);
                        if radius > 0
                            se = strel('disk', radius);
                            for z = 1:size(img_3d, 3)
                                slice = img_3d(:,:,z);
                                if strcmp(method, 'Rolling-ball'), slice = slice - imopen(slice, se);
                                elseif strcmp(method, 'Top-hat'), slice = imtophat(slice, se);
                                end
                                img_3d(:,:,z) = slice;
                            end
                        end
                     end
                end
            end
            preprocessed_nuc = img_3d;
        end
        
        % Apply final clipping if checked
        if getScalar(handles.nucBgCorrClipChecks.Value)
            preprocessed_nuc(preprocessed_nuc < 0) = 0;
        end
    end
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
    end
end

function img_out = apply_nuc_2d_bg_corr(img_in, h)
    img_out = img_in;
    method = h.nucBgMethodDrop.Value;
    param_val = h.nucBgParamInput.Value;
    if param_val > 0 && ~strcmp(method, 'None')
        if strcmp(method, 'Gaussian')
            img_out = img_in - imgaussfilt(img_in, param_val);
        else
            se = strel('disk', round(param_val));
            if strcmp(method, 'Rolling-ball'), img_out = img_in - imopen(img_in, se);
            elseif strcmp(method, 'Top-hat'), img_out = imtophat(img_in, se);
            end
        end
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

function psf = generate3DPSF(sigma_xy, sigma_z, size_xy, size_z)
    % Generate a 3D Gaussian PSF for deconvolution
    [X, Y, Z] = meshgrid(-size_xy:size_xy, -size_xy:size_xy, -size_z:size_z);
    psf = exp(-(X.^2 + Y.^2)/(2*sigma_xy^2) - Z.^2/(2*sigma_z^2));
    psf = psf / sum(psf(:)); % Normalize
end

function psf_out = cropPSFToImage(psf_in, img_size)
    % Crop or pad PSF to match image dimensions (centered on middle frame)
    psf_size = size(psf_in);
    if length(psf_size) == 2, psf_size(3) = 1; end
    if length(img_size) == 2, img_size(3) = 1; end
    
    psf_center = ceil(psf_size / 2);
    img_center = ceil(img_size / 2);
    psf_out = zeros(img_size);
    
    for dim = 1:3
        if psf_size(dim) > img_size(dim)
            start_idx(dim) = psf_center(dim) - floor(img_size(dim)/2);
            end_idx(dim) = start_idx(dim) + img_size(dim) - 1;
        else
            start_idx(dim) = 1;
            end_idx(dim) = psf_size(dim);
        end
    end
    
    if all(psf_size >= img_size)
        psf_out = psf_in(start_idx(1):end_idx(1), start_idx(2):end_idx(2), start_idx(3):end_idx(3));
    else
        out_start = max(1, img_center - psf_center + 1);
        out_end = min(img_size, out_start + psf_size - 1);
        psf_start = max(1, psf_center - img_center + 1);
        psf_end = min(psf_size, psf_start + (out_end - out_start));
        psf_out(out_start(1):out_end(1), out_start(2):out_end(2), out_start(3):out_end(3)) = ...
            psf_in(psf_start(1):psf_end(1), psf_start(2):psf_end(2), psf_start(3):psf_end(3));
    end
    
    if sum(psf_out(:)) > 0
        psf_out = psf_out / sum(psf_out(:));
    end
end

%% Helper Functions for SNAP_batch compatibility

function val = getScalar(input)
    % Ensure value is a proper scalar (handle cells, arrays, etc.)
    if iscell(input)
        val = input{1};
    elseif numel(input) > 1
        val = input(1);
    else
        val = input;
    end
    % Ensure truly scalar
    if ~isscalar(val)
        val = val(1);
    end
end
