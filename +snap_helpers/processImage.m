function processed_img = processImage(handles, ch, statusLabel)
% Processes a single channel image using the selected mode for each step.
    
    if nargin < 3, statusLabel = []; end % Optional status label handle

    raw_img = double(handles.rawChannel{ch});
    
    % --- Get Spacing Parameters (needed by multiple steps) ---
    xy_spacing = handles.xySpacingInputs(ch).Value;
    z_spacing = handles.zSpacingInputs(ch).Value;
    
    img_after_deconv = raw_img; % Initialize with raw image

    % --- Step 0: Deconvolution (Applied First, Before Everything) ---
    if handles.deconvEnabledChecks(ch).Value
        if ~isempty(statusLabel)
            statusLabel.Text = ['Status: Deconvolving Ch. ' num2str(ch) '...'];
            drawnow;
        end
        
        % Get deconvolution parameters
        method = handles.deconvMethodDrops(ch).Value;
        psf_source = handles.deconvPSFSourceDrops(ch).Value;
        
        % Get or generate PSF
        if strcmp(psf_source, 'Load File')
            psf_file = handles.deconvPSFPathTexts(ch).Value;
            if isempty(psf_file) || ~isfile(psf_file)
                warning('PSF file not found. Skipping deconvolution for channel %d.', ch);
                img_after_deconv = raw_img;
            else
                psf_raw = double(tiffreadVolume(psf_file));
                % Crop/pad PSF to match image dimensions (centered on middle frame)
                psf = cropPSFToImage(psf_raw, size(raw_img));
                psf = psf / sum(psf(:)); % Normalize
            end
        else % 'Generate'
            psf_sigma_xy = handles.deconvPSFSigmaXYInputs(ch).Value;
            psf_sigma_z = handles.deconvPSFSigmaZInputs(ch).Value;
            psf_size_xy = handles.deconvPSFSizeXYInputs(ch).Value;
            psf_size_z = handles.deconvPSFSizeZInputs(ch).Value;
            psf = generate3DPSF(psf_sigma_xy, psf_sigma_z, psf_size_xy, psf_size_z);
        end
        
        % Apply deconvolution based on selected method
        if exist('psf', 'var')
            switch method
                case 'Lucy-Richardson'
                    iterations = handles.deconvLRIterationsInputs(ch).Value;
                    damping = handles.deconvLRDampingInputs(ch).Value;
                    if damping > 0
                        img_after_deconv = deconvlucy(raw_img, psf, iterations, damping);
                    else
                        img_after_deconv = deconvlucy(raw_img, psf, iterations);
                    end
                    
                case 'Wiener'
                    nsr = handles.deconvWienerNSRInputs(ch).Value;
                    img_after_deconv = deconvwnr(raw_img, psf, nsr);
                    
                case 'Blind'
                    iterations = handles.deconvBlindIterationsInputs(ch).Value;
                    under_relax = handles.deconvBlindUnderRelaxInputs(ch).Value;
                    if under_relax > 0
                        [img_after_deconv, ~] = deconvblind(raw_img, psf, iterations, under_relax);
                    else
                        [img_after_deconv, ~] = deconvblind(raw_img, psf, iterations);
                    end
            end
            
            if ~isempty(statusLabel)
                statusLabel.Text = ['Status: Deconvolution complete for Ch. ' num2str(ch)];
                drawnow;
            end
        end
    end
    
    % Check for abort after deconvolution
    if isfield(handles, 'fig') && isfield(handles, 'abortRequested')
        latest_handles = guidata(handles.fig);
        if latest_handles.abortRequested
            processed_img = [];
            return;
        end
    end
    
    img_after_preprocess = img_after_deconv; % Use deconvolved image for next steps

    % --- Step 1: Pre-processing ---
    if handles.preprocEnabledChecks(ch).Value
        if ~isempty(statusLabel)
            statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) '...'];
            drawnow;
        end
        mode = handles.preprocessModeDrops(ch).Value;
        is_scaled = handles.preprocessScaleChecks(ch).Value;
        
        % Apply pre-processing based on the selected mode
        img_after_preprocess = raw_img; % Start with the raw image

        if strcmp(mode, 'On Z-Projection')
            proj_type = handles.preprocessProjectionDrops(ch).Value;
            img_2d = projectZ(raw_img, proj_type);
            img_after_preprocess = apply_2d_preproc(img_2d, handles, ch);
        else % 3D or 2D Slice-by-slice
            img_3d = raw_img;
            is_3d_mode = strcmp(mode, '3D');
            method = handles.preprocMethodDrops(ch).Value;

            switch method
                case 'Gaussian'
                    sigma_val = handles.gaussInputs(ch).Value;
                    if sigma_val > 0
                        if is_3d_mode
                            if is_scaled
                                sigma_pixels_xy = sigma_val / xy_spacing;
                                sigma_pixels_z = sigma_val / z_spacing;
                                sigma_vec = [sigma_pixels_xy, sigma_pixels_xy, sigma_pixels_z];
                                if any(sigma_vec > 0), img_3d = imgaussfilt3(img_3d, sigma_vec); end
                            else % Unscaled 3D (voxels)
                                img_3d = imgaussfilt3(img_3d, sigma_val);
                            end
                        else % Slice-by-slice (pixels)
                            if is_scaled
                                sigma_pixels_2d = sigma_val / xy_spacing;
                                for z = 1:size(img_3d, 3)
                                    if ~isempty(statusLabel)
                                        statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                        drawnow;
                                    end
                                    img_3d(:,:,z) = imgaussfilt(img_3d(:,:,z), sigma_pixels_2d);
                                end
                            else
                                for z = 1:size(img_3d, 3)
                                    if ~isempty(statusLabel)
                                        statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                        drawnow;
                                    end
                                    img_3d(:,:,z) = imgaussfilt(img_3d(:,:,z), sigma_val);
                                end
                            end
                        end
                    end
                case 'Median'
                    med_val = handles.medianInputs(ch).Value;
                    if med_val > 0
                        if is_3d_mode
                            if is_scaled
                                neigh_xy = max(1, 2*round( (med_val / xy_spacing) / 2) + 1);
                                neigh_z = max(1, 2*round( (med_val / z_spacing) / 2) + 1);
                                neigh_vec = [neigh_xy, neigh_xy, neigh_z];
                                img_3d = medfilt3(img_3d, neigh_vec);
                            else % Unscaled 3D (voxels)
                                med_val_odd = max(1, 2*round(med_val/2) + 1);
                                img_3d = medfilt3(img_3d, [med_val_odd med_val_odd med_val_odd]);
                            end
                        else % Slice-by-slice (pixels)
                            med_val_odd = max(1, 2*round(med_val/2) + 1);
                            for z = 1:size(img_3d, 3)
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                    drawnow;
                                end
                                img_3d(:,:,z) = medfilt2(img_3d(:,:,z), [med_val_odd med_val_odd]);
                            end
                        end
                    end
                case 'Non-Local Means'
                    % Get NLM parameters
                    filter_strength = handles.nlmFilterStrengthInputs(ch).Value;
                    search_window = handles.nlmSearchWindowInputs(ch).Value;
                    comparison_window = handles.nlmComparisonWindowInputs(ch).Value;
                    
                    if filter_strength > 0
                        if is_3d_mode
                            % Check if imnlmfilt3 is available (newer MATLAB versions)
                            if exist('imnlmfilt3', 'file') == 2
                                % 3D Non-Local Means (process entire volume)
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - 3D Non-Local Means...'];
                                    drawnow;
                                end
                                img_3d = imnlmfilt3(img_3d, 'DegreeOfSmoothing', filter_strength, ...
                                                   'SearchWindowSize', search_window, ...
                                                   'ComparisonWindowSize', comparison_window);
                            else
                                % Fallback: Use slice-by-slice 2D processing for 3D mode
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - 3D Non-Local Means (slice-by-slice fallback)...'];
                                    drawnow;
                                end
                                fprintf('Warning: imnlmfilt3 not available. Using slice-by-slice processing for 3D Non-Local Means on Channel %d.\n', ch);
                                for z = 1:size(img_3d, 3)
                                    if ~isempty(statusLabel)
                                        statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) ' (NLM fallback)...'];
                                        drawnow;
                                    end
                                    img_3d(:,:,z) = imnlmfilt(img_3d(:,:,z), 'DegreeOfSmoothing', filter_strength, ...
                                                              'SearchWindowSize', search_window, ...
                                                              'ComparisonWindowSize', comparison_window);
                                end
                            end
                        else
                            % Slice-by-slice Non-Local Means
                            for z = 1:size(img_3d, 3)
                                % Check for abort every 10 slices
                                if mod(z, 10) == 0 && isfield(handles, 'fig') && isfield(handles, 'abortRequested')
                                    latest_handles = guidata(handles.fig);
                                    if latest_handles.abortRequested
                                        processed_img = [];
                                        return;
                                    end
                                end
                                
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) ' (NLM)...'];
                                    drawnow;
                                end
                                img_3d(:,:,z) = imnlmfilt(img_3d(:,:,z), 'DegreeOfSmoothing', filter_strength, ...
                                                          'SearchWindowSize', search_window, ...
                                                          'ComparisonWindowSize', comparison_window);
                            end
                        end
                    end
                case 'Wavelet Denoising'
                    if ~is_3d_mode
                         for z = 1:size(img_3d, 3)
                            % Check for abort every 10 slices
                            if mod(z, 10) == 0 && isfield(handles, 'fig') && isfield(handles, 'abortRequested')
                                latest_handles = guidata(handles.fig);
                                if latest_handles.abortRequested
                                    processed_img = [];
                                    return;
                                end
                            end
                            
                            if ~isempty(statusLabel)
                                statusLabel.Text = ['Status: Pre-processing Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                drawnow;
                            end
                            slice = img_3d(:,:,z);
                            
                            wname = handles.waveletNameDrops(ch).Value;
                            level = handles.waveletLevelInputs(ch).Value;
                            rule = handles.waveletThresholdRuleDrops(ch).Value;
                            sorh_str = handles.waveletThresholdMethodDrops(ch).Value;
                            if strcmp(sorh_str, 'soft'), sorh = 's'; else, sorh = 'h'; end
                            
                            if level > 0 && floor(level) == level
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
                                    warning('Wavelet denoising failed for Channel %d, Z-frame %d: %s. Using original slice.', ch, z, ME.message);
                                end
                            end
                        end
                    end
            end
            img_after_preprocess = img_3d;
        end
        
        % Enforce non-negativity after all pre-processing if checked
        if handles.preprocClipChecks(ch).Value
            img_after_preprocess(img_after_preprocess < 0) = 0;
        end
    else
        img_after_preprocess = raw_img; % Pass raw if preproc is disabled
    end
    
    % Check for abort after preprocessing
    if isfield(handles, 'fig') && isfield(handles, 'abortRequested')
        latest_handles = guidata(handles.fig);
        if latest_handles.abortRequested
            processed_img = [];
            return;
        end
    end

    processed_img = img_after_preprocess; % Pass through by default

    % --- Step 2: Background Correction ---
    if handles.bgCorrEnabledChecks(ch).Value
        if ~isempty(statusLabel)
            statusLabel.Text = ['Status: BG Correcting Ch. ' num2str(ch) '...'];
            drawnow;
        end
        bg_mode = handles.bgCorrModeDrops(ch).Value;
        is_scaled_bg = handles.bgCorrScaleChecks(ch).Value;

        if strcmp(bg_mode, 'On Z-Projection')
            proj_type = handles.bgCorrProjectionDrops(ch).Value;
            img_2d = projectZ(img_after_preprocess, proj_type);
            processed_img = apply2DBG(img_2d, handles, ch);
        else % 3D or 2D Slice-by-slice
            img_3d = img_after_preprocess;
            is_3d_mode_bg = strcmp(bg_mode, '3D');
            
            method = handles.bgMethodDrops(ch).Value;
            param_val = handles.bgParamInputs(ch).Value;
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
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: BG Correcting Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                    drawnow;
                                end
                                slice = img_3d(:,:,z);
                                bg_slice = imgaussfilt(slice, sigma_pixels_2d);
                                processed_slice = slice - bg_slice;
                                img_3d(:,:,z) = processed_slice;
                            end
                        else
                            for z = 1:size(img_3d, 3)
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: BG Correcting Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                    drawnow;
                                end
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
                                if ~isempty(statusLabel)
                                    statusLabel.Text = ['Status: BG Correcting Ch. ' num2str(ch) ' - Z-frame ' num2str(z) '/' num2str(size(img_3d, 3)) '...'];
                                    drawnow;
                                end
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
        end
        processed_img = img_3d;
        
        % Apply final clipping if checked
        if handles.bgCorrClipChecks(ch).Value
            processed_img(processed_img < 0) = 0;
        end
    else
        processed_img = img_after_preprocess; % Pass through if BG correction is disabled
    end

    % --- Nested Helper Functions ---
    
    function img_out = apply_2d_preproc(img_in, h, c)
        img_out = img_in;
        method = h.preprocMethodDrops(c).Value;

        switch method
            case 'Gaussian'
                if h.gaussInputs(c).Value > 0
                    img_out = imgaussfilt(img_out, h.gaussInputs(c).Value);
                end
            case 'Median'
                if h.medianInputs(c).Value > 0
                    med_val = h.medianInputs(c).Value;
                    med_val_odd = max(1, 2*round(med_val/2) + 1);
                    img_out = medfilt2(img_out, [med_val_odd med_val_odd]);
                end
            case 'Wavelet Denoising'
                 wname = h.waveletNameDrops(c).Value;
                 level = h.waveletLevelInputs(c).Value;
                 rule = h.waveletThresholdRuleDrops(c).Value;
                 sorh_str = h.waveletThresholdMethodDrops(c).Value;
                 if strcmp(sorh_str, 'soft'), sorh = 's'; else, sorh = 'h'; end
                 
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
                        warning('Wavelet denoising failed for 2D projection: %s. Using original image.', ME.message);
                    end
                 end
        end
    end

    function img_out = apply2DBG(img_in, h, c)
        img_out = img_in;
        method = h.bgMethodDrops(c).Value;
        param_val = h.bgParamInputs(c).Value;
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
        % sigma_xy: PSF sigma in XY plane (pixels)
        % sigma_z: PSF sigma in Z direction (slices)
        % size_xy: PSF half-size in XY (pixels)
        % size_z: PSF half-size in Z (slices)
        
        [X, Y, Z] = meshgrid(-size_xy:size_xy, -size_xy:size_xy, -size_z:size_z);
        psf = exp(-(X.^2 + Y.^2)/(2*sigma_xy^2) - Z.^2/(2*sigma_z^2));
        psf = psf / sum(psf(:)); % Normalize
    end
    
    function psf_out = cropPSFToImage(psf_in, img_size)
        % Crop or pad PSF to match image dimensions
        % Both PSF and image are assumed to have their focus at the center Z-frame
        % psf_in: Input PSF (may be larger or smaller than image)
        % img_size: Size of the image [rows, cols, slices]
        
        psf_size = size(psf_in);
        
        % Handle 2D case
        if length(psf_size) == 2
            psf_size(3) = 1;
        end
        if length(img_size) == 2
            img_size(3) = 1;
        end
        
        % Find center indices for both PSF and image
        psf_center = ceil(psf_size / 2);
        img_center = ceil(img_size / 2);
        
        % Calculate how much to crop/pad in each dimension
        % We want to keep the PSF centered on its middle frame
        psf_out = zeros(img_size);
        
        for dim = 1:3
            if psf_size(dim) > img_size(dim)
                % PSF is larger - need to crop centered on middle frame
                start_idx(dim) = psf_center(dim) - floor(img_size(dim)/2);
                end_idx(dim) = start_idx(dim) + img_size(dim) - 1;
            else
                % PSF is smaller or equal - will pad with zeros
                start_idx(dim) = 1;
                end_idx(dim) = psf_size(dim);
            end
        end
        
        % Extract the cropped PSF region
        if all(psf_size >= img_size)
            % Simple crop case
            psf_out = psf_in(start_idx(1):end_idx(1), start_idx(2):end_idx(2), start_idx(3):end_idx(3));
        else
            % Need to handle padding - place PSF centered in output
            out_start = max(1, img_center - psf_center + 1);
            out_end = min(img_size, out_start + psf_size - 1);
            
            psf_start = max(1, psf_center - img_center + 1);
            psf_end = min(psf_size, psf_start + (out_end - out_start));
            
            psf_out(out_start(1):out_end(1), out_start(2):out_end(2), out_start(3):out_end(3)) = ...
                psf_in(psf_start(1):psf_end(1), psf_start(2):psf_end(2), psf_start(3):psf_end(3));
        end
        
        % Ensure PSF is normalized
        if sum(psf_out(:)) > 0
            psf_out = psf_out / sum(psf_out(:));
        end
    end
end