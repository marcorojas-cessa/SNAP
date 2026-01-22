function performanceComparison(img, method, radius, algorithm, h)
% Performance comparison between original and optimized segmentation algorithms
% This function helps you choose the best optimization level for your data

    fprintf('\n=== PERFORMANCE COMPARISON ===\n');
    fprintf('Testing %s algorithm with radius %d on image size %dx%d\n', ...
        algorithm, radius, size(img, 1), size(img, 2));
    
    % Test original implementation
    fprintf('\n1. Testing ORIGINAL implementation...\n');
    tic;
    try
        mask_original = snap_helpers.segmentNuclei(img, method, radius, algorithm, h);
        time_original = toc;
        fprintf('   Original: %.3f seconds\n', time_original);
        success_original = true;
    catch ME
        fprintf('   Original: FAILED - %s\n', ME.message);
        time_original = inf;
        success_original = false;
        mask_original = [];
    end
    
    % Test optimized implementation
    fprintf('\n2. Testing OPTIMIZED implementation...\n');
    tic;
    try
        mask_optimized = snap_helpers.segmentNuclei_optimized(img, method, radius, algorithm, h);
        time_optimized = toc;
        fprintf('   Optimized: %.3f seconds\n', time_optimized);
        success_optimized = true;
    catch ME
        fprintf('   Optimized: FAILED - %s\n', ME.message);
        time_optimized = inf;
        success_optimized = false;
        mask_optimized = [];
    end
    
    % Test ultra-fast implementation
    fprintf('\n3. Testing ULTRA-FAST implementation...\n');
    tic;
    try
        mask_ultrafast = snap_helpers.segmentNuclei_ultrafast(img, method, radius, algorithm, h);
        time_ultrafast = toc;
        fprintf('   Ultra-fast: %.3f seconds\n', time_ultrafast);
        success_ultrafast = true;
    catch ME
        fprintf('   Ultra-fast: FAILED - %s\n', ME.message);
        time_ultrafast = inf;
        success_ultrafast = false;
        mask_ultrafast = [];
    end
    
    % Calculate speedup factors
    fprintf('\n=== PERFORMANCE RESULTS ===\n');
    if success_original && success_optimized
        speedup_opt = time_original / time_optimized;
        fprintf('Optimized vs Original: %.1fx speedup\n', speedup_opt);
    end
    
    if success_original && success_ultrafast
        speedup_ultra = time_original / time_ultrafast;
        fprintf('Ultra-fast vs Original: %.1fx speedup\n', speedup_ultra);
    end
    
    if success_optimized && success_ultrafast
        speedup_ultra_vs_opt = time_optimized / time_ultrafast;
        fprintf('Ultra-fast vs Optimized: %.1fx speedup\n', speedup_ultra_vs_opt);
    end
    
    % Check accuracy (if both original and optimized succeeded)
    if success_original && (success_optimized || success_ultrafast)
        fprintf('\n=== ACCURACY COMPARISON ===\n');
        
        if success_optimized
            accuracy_opt = compareMasks(mask_original, mask_optimized);
            fprintf('Optimized accuracy: %.2f%% match with original\n', accuracy_opt * 100);
        end
        
        if success_ultrafast
            accuracy_ultra = compareMasks(mask_original, mask_ultrafast);
            fprintf('Ultra-fast accuracy: %.2f%% match with original\n', accuracy_ultra * 100);
        end
    end
    
    % Recommendations
    fprintf('\n=== RECOMMENDATIONS ===\n');
    if success_ultrafast && time_ultrafast < time_optimized
        fprintf('RECOMMENDED: Use ULTRA-FAST implementation\n');
        fprintf('  - Best performance\n');
        fprintf('  - Suitable for large images and batch processing\n');
    elseif success_optimized && time_optimized < time_original
        fprintf('RECOMMENDED: Use OPTIMIZED implementation\n');
        fprintf('  - Good balance of speed and accuracy\n');
        fprintf('  - More stable than ultra-fast version\n');
    else
        fprintf('RECOMMENDED: Use ORIGINAL implementation\n');
        fprintf('  - Most reliable\n');
        fprintf('  - Use when accuracy is critical\n');
    end
    
    fprintf('\n=== OPTIMIZATION TIPS ===\n');
    fprintf('For even better performance:\n');
    fprintf('1. Use 3D mode instead of slice-by-slice when possible\n');
    fprintf('2. Reduce radius size if image quality allows\n');
    fprintf('3. Use parallel processing for multiple channels\n');
    fprintf('4. Consider GPU acceleration for very large images\n');
    fprintf('5. Process only regions of interest when possible\n');
end

function accuracy = compareMasks(mask1, mask2)
    % Compare two binary masks and return accuracy percentage
    if isempty(mask1) || isempty(mask2) || ~isequal(size(mask1), size(mask2))
        accuracy = 0;
        return;
    end
    
    % Calculate overlap
    overlap = sum(mask1(:) & mask2(:));
    total = sum(mask1(:) | mask2(:));
    
    if total == 0
        accuracy = 1; % Both masks are empty
    else
        accuracy = overlap / total;
    end
end
