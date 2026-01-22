function handles = updatePreview(handles)
% Processes image, finds local maxima, and displays the result.

ch = str2double(handles.previewChanDrop.Value);

if isempty(handles.rawChannel{ch})
    cla(handles.previewAx);
    title(handles.previewAx, ['Channel ' num2str(ch) ' (No Image)']);
    return;
end

img = double(handles.rawChannel{ch});

% --- Step 1: Pre-processing (Smoothing, Wavelet, BG Correction) ---
if handles.smoothGaussianChecks(ch).Value && handles.smoothGaussianInputs(ch).Value > 0
    sigma = handles.smoothGaussianInputs(ch).Value;
    if ndims(img) == 3, img = imgaussfilt3(img, sigma); else, img = imgaussfilt(img, sigma); end
end
if handles.smoothMedianChecks(ch).Value && handles.smoothMedianInputs(ch).Value > 0
    filterSize = repmat(handles.smoothMedianInputs(ch).Value, 1, ndims(img));
    if ndims(img) == 3, img = medfilt3(img, filterSize); else, img = medfilt2(img, filterSize(1:2)); end
end
if handles.waveletDenoiseChecks(ch).Value
    level = handles.waveletLevelInputs(ch).Value;
    if level > 0 && isfinite(level) && level == round(level)
        try
            wname = handles.waveletNameDrops(ch).Value;
            thrRule = handles.waveletThresholdRuleDrops(ch).Value;
            thrMethod = 's';
            if strcmp(handles.waveletThresholdMethodDrops(ch).Value, 'hard'), thrMethod = 'h'; end
            
            if ndims(img) == 3
                for z = 1:size(img,3)
                    % Note: No status label available in updatePreview, so no z-frame progress shown
                    currentSlice = img(:,:,z);
                    maxLvl = min(floor(log2(size(currentSlice))));
                    if level > maxLvl, continue; end
                    thr = thselect(currentSlice, thrRule);
                    if ~isscalar(thr) || ~isfinite(thr), continue; end
                    denoisedSlice = wdencmp('gbl', currentSlice, wname, level, thr, thrMethod, 1);
                    isValid = true;
                    if ~isa(denoisedSlice,'numeric') || isempty(denoisedSlice) || ~isreal(denoisedSlice) || ~all(isfinite(denoisedSlice(:))) || ~isequal(size(denoisedSlice), size(currentSlice)), isValid = false; end
                    if isValid, img(:,:,z) = double(denoisedSlice); end
                end
            else
                maxLvl = min(floor(log2(size(img))));
                if level <= maxLvl
                    thr = thselect(img, thrRule);
                    if isscalar(thr) && isfinite(thr)
                        denoisedImg = wdencmp('gbl', img, wname, level, thr, thrMethod, 1);
                        isValid = true;
                        if ~isa(denoisedImg,'numeric') || isempty(denoisedImg) || ~isreal(denoisedImg) || ~all(isfinite(denoisedImg(:))) || ~isequal(size(denoisedImg), size(img)), isValid = false; end
                        if isValid, img = double(denoisedImg); end
                    end
                end
            end
        catch ME, warning('Wavelet denoising failed for Channel %d: %s.', ch, ME.message); end
    end
end
method = handles.bgMethodDrop(ch).Value;
param = handles.bgParamInput(ch).Value;
if param > 0 && ~strcmp(method, 'None')
    if strcmp(method, 'Gaussian')
         if ndims(img) == 3, bg = imgaussfilt3(img,param); else, bg = imgaussfilt(img,param); end
         img = img - bg;
    else
        if ndims(img) == 3
            if strcmp(method, 'Rolling-ball'), se = strel('ball',param,param); else, se = strel('sphere',param); end
        else, se = strel('disk',param); end
        if strcmp(method, 'Rolling-ball'), img = img - imopen(img,se); else, img = imtophat(img,se); end
    end
end

% --- Max projection for consistent processing ---
if ndims(img) == 3, procProj = max(img,[],3); else, procProj = img; end

% --- Step 2: Local Maxima Detection ---
maxima_method = handles.maximaMethodDrops(ch).Value;
    neighborhood_size = handles.maximaNeighborhoodInputs(ch).Value;
bw_maxima = false(size(procProj)); % Initialize an empty mask

try
    switch maxima_method
        case 'Simple Regional'
            bw_maxima = findLocalMaximaNeighborhood(procProj, neighborhood_size);
            
        case 'Extended Maxima'
            h = handles.hmaxInputs(ch).Value;
            if h > 0
                bw_maxima = findExtendedMaximaNeighborhood(procProj, neighborhood_size, h);
            end
            
        case 'Laplacian of Gaussian'
            sigma = handles.logSigmaInputs(ch).Value;
            thresh = handles.logThreshInputs(ch).Value;
            if sigma > 0
                filt_size = 2 * ceil(3*sigma) + 1;
                log_filter = fspecial('log', filt_size, sigma);
                img_log = imfilter(procProj, log_filter, 'replicate', 'same');
                bw_maxima = findLocalMaximaNeighborhood(img_log, neighborhood_size) & (img_log > thresh);
            end
    end
catch ME
    warning('Local maxima detection failed for Channel %d: %s.', ch, ME.message);
end

% --- Step 3: Display Results ---
imagesc(handles.previewAx, procProj);
axis(handles.previewAx,'image'); colormap(handles.previewAx,'gray');
title(handles.previewAx,['Preview - Channel ' num2str(ch)]);

% Overlay the found maxima only if any were found
hold(handles.previewAx, 'on');
[rows, cols] = find(bw_maxima);
if ~isempty(rows)
    plot(handles.previewAx, cols, rows, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
end
hold(handles.previewAx, 'off');

% --- Helper Functions for Neighborhood-Based Maxima Detection ---

function bw = findLocalMaximaNeighborhood(img, neighborhood_size)
    % Find local maxima using a neighborhood window approach
    % neighborhood_size: radius in pixels (e.g., 2 means 5x5 window)
    
    % Ensure neighborhood_size is integer to avoid colon operator errors
    neighborhood_size = round(neighborhood_size);
    
    if neighborhood_size <= 0
        bw = false(size(img));
        return;
    end
    
    [rows, cols] = size(img);
    bw = false(rows, cols);
    
    % Create neighborhood window
    window_size = 2 * neighborhood_size + 1;
    
    for r = 1:rows
        for c = 1:cols
            % Define neighborhood boundaries
            r_start = max(1, r - neighborhood_size);
            r_end = min(rows, r + neighborhood_size);
            c_start = max(1, c - neighborhood_size);
            c_end = min(cols, c + neighborhood_size);
            
            % Extract neighborhood
            neighborhood = img(r_start:r_end, c_start:c_end);
            center_value = img(r, c);
            
            % Check if center is maximum in neighborhood
            if center_value == max(neighborhood(:)) && center_value > 0
                % Additional check: ensure no ties (only one maximum)
                max_positions = find(neighborhood == center_value);
                if length(max_positions) == 1 || (length(max_positions) > 1 && max_positions(1) == sub2ind(size(neighborhood), neighborhood_size + 1, neighborhood_size + 1))
                    bw(r, c) = true;
                end
            end
        end
    end
end

function bw = findExtendedMaximaNeighborhood(img, neighborhood_size, h_threshold)
    % Find extended maxima using a neighborhood window approach with height threshold
    % neighborhood_size: radius in pixels (e.g., 2 means 5x5 window)
    % h_threshold: minimum height difference from surrounding
    
    if neighborhood_size <= 0 || h_threshold <= 0
        bw = false(size(img));
        return;
    end
    
    [rows, cols] = size(img);
    bw = false(rows, cols);
    
    % Create neighborhood window
    window_size = 2 * neighborhood_size + 1;
    
    for r = 1:rows
        for c = 1:cols
            % Define neighborhood boundaries
            r_start = max(1, r - neighborhood_size);
            r_end = min(rows, r + neighborhood_size);
            c_start = max(1, c - neighborhood_size);
            c_end = min(cols, c + neighborhood_size);
            
            % Extract neighborhood
            neighborhood = img(r_start:r_end, c_start:c_end);
            center_value = img(r, c);
            
            % Check if center is maximum in neighborhood
            if center_value == max(neighborhood(:)) && center_value > 0
                % Check height threshold: center must be h_threshold above surrounding
                surrounding_values = neighborhood(:);
                surrounding_values(sub2ind(size(neighborhood), neighborhood_size + 1, neighborhood_size + 1)) = [];
                
                if center_value >= max(surrounding_values) + h_threshold
                    % Additional check: ensure no ties (only one maximum)
                    max_positions = find(neighborhood == center_value);
                    if length(max_positions) == 1 || (length(max_positions) > 1 && max_positions(1) == sub2ind(size(neighborhood), neighborhood_size + 1, neighborhood_size + 1))
                        bw(r, c) = true;
                    end
                end
            end
        end
    end
end