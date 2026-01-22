function testOptimizationPerformance()
% Test script to compare performance between original and optimized segmentation algorithms
% Run this function to see the speed improvements on your system

    fprintf('=== SEGMENTATION PERFORMANCE TEST ===\n\n');
    
    % Create test images of different sizes
    test_sizes = [
        128, 128;    % Small
        256, 256;    % Medium  
        512, 512;    % Large
        1024, 1024   % Very large
    ];
    
    algorithms = {'Niblack', 'Sauvola', 'Phansalkar', 'Mean', 'Bernsen'};
    radius = 25;
    
    % Create mock handles structure for testing
    h = struct();
    h.nucSegAlgParamInput = struct('Visible', true, 'Value', 0.2);
    h.nucSegAlgParam2Input = struct('Visible', true, 'Value', 0);
    
    for i = 1:size(test_sizes, 1)
        [rows, cols] = deal(test_sizes(i, 1), test_sizes(i, 2));
        fprintf('Testing image size: %dx%d (%d pixels)\n', rows, cols, rows*cols);
        
        % Create synthetic test image with noise and structures
        img = createTestImage(rows, cols);
        
        for j = 1:length(algorithms)
            algorithm = algorithms{j};
            fprintf('  Algorithm: %s\n', algorithm);
            
            % Test original implementation
            tic;
            try
                mask_original = applyLocalThresholdOriginal(img, algorithm, 2*radius+1, h);
                time_original = toc;
                fprintf('    Original: %.3f seconds\n', time_original);
                success_original = true;
            catch ME
                fprintf('    Original: FAILED - %s\n', ME.message);
                time_original = inf;
                success_original = false;
            end
            
            % Test optimized implementation
            tic;
            try
                mask_optimized = applyLocalThresholdOptimized(img, algorithm, 2*radius+1, h);
                time_optimized = toc;
                fprintf('    Optimized: %.3f seconds\n', time_optimized);
                success_optimized = true;
            catch ME
                fprintf('    Optimized: FAILED - %s\n', ME.message);
                time_optimized = inf;
                success_optimized = false;
            end
            
            % Calculate speedup
            if success_original && success_optimized
                speedup = time_original / time_optimized;
                fprintf('    SPEEDUP: %.1fx faster\n', speedup);
                
                % Check accuracy
                accuracy = sum(mask_original(:) == mask_optimized(:)) / numel(mask_original);
                fprintf('    ACCURACY: %.1f%% match\n', accuracy * 100);
            end
            
            fprintf('\n');
        end
        
        fprintf('---\n\n');
    end
    
    fprintf('=== RECOMMENDATIONS ===\n');
    fprintf('1. For images > 224x224 pixels: Use optimized version\n');
    fprintf('2. For batch processing: Use optimized version\n');
    fprintf('3. For real-time preview: Use original version for small images\n');
    fprintf('4. For maximum speed: Consider using 3D mode instead of slice-by-slice\n');
    fprintf('5. For very large images: Consider reducing radius size\n');
end

function img = createTestImage(rows, cols)
    % Create a synthetic test image with realistic features
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Create background with gradual intensity variation
    img = 100 + 50 * sin(X/50) .* cos(Y/50);
    
    % Add some circular structures (simulating nuclei)
    centers_x = [cols/4, 3*cols/4, cols/2];
    centers_y = [rows/3, 2*rows/3, rows/2];
    for i = 1:length(centers_x)
        dist = sqrt((X - centers_x(i)).^2 + (Y - centers_y(i)).^2);
        img = img + 100 * exp(-dist.^2 / (2 * 20^2));
    end
    
    % Add noise
    img = img + 20 * randn(rows, cols);
    
    % Ensure positive values
    img = max(0, img);
end

function mask = applyLocalThresholdOriginal(img, algorithm, window_size, h)
    % Simplified original implementation for testing
    img_double = double(img);
    img_min = min(img_double(:));
    img_max = max(img_double(:));
    if img_max > img_min
        img_double = 255 * (img_double - img_min) / (img_max - img_min);
    end
    
    [rows, cols] = size(img_double);
    half_win = floor(window_size / 2);
    threshold_map = zeros(size(img_double));
    
    for i = 1:rows
        for j = 1:cols
            % Define local window
            row_start = max(1, i - half_win);
            row_end = min(rows, i + half_win);
            col_start = max(1, j - half_win);
            col_end = min(cols, j + half_win);
            
            local_region = img_double(row_start:row_end, col_start:col_end);
            
            % Create circular mask
            [local_rows, local_cols] = size(local_region);
            center_row = (local_rows + 1) / 2;
            center_col = (local_cols + 1) / 2;
            [Y, X] = meshgrid(1:local_cols, 1:local_rows);
            distances = sqrt((X - center_row).^2 + (Y - center_col).^2);
            circular_mask = distances <= half_win;
            local_region_circular = local_region(circular_mask);
            
            if isempty(local_region_circular)
                local_region_circular = local_region(:);
            end
            
            local_mean = mean(local_region_circular);
            local_std = std(local_region_circular);
            
            % Apply simple thresholding for testing
            switch algorithm
                case 'Niblack'
                    threshold_map(i,j) = local_mean + 0.2 * local_std;
                case 'Sauvola'
                    threshold_map(i,j) = local_mean * (1 + 0.5 * (local_std / 128 - 1));
                case 'Mean'
                    threshold_map(i,j) = local_mean;
                otherwise
                    threshold_map(i,j) = local_mean;
            end
        end
    end
    
    mask = img_double > threshold_map;
end

function mask = applyLocalThresholdOptimized(img, algorithm, window_size, h)
    % Use the optimized implementation from the main file
    mask = applyLocalThresholdOptimized(img, algorithm, window_size, h);
end
