function SNAP_classify()
% SNAP_classify - Standalone SVM training app for SNAP spot classification
%
% ============================================================================
% COMPREHENSIVE SPOT LABELING AND CLASSIFIER TRAINING APPLICATION
% ============================================================================
%
% This application allows users to:
%   1. Load exported fit results from SNAP or SNAP_batch
%   2. Load associated images for visual inspection
%   3. Select which features to use for classification
%   4. Manually label spots as Real or Noise
%   5. Train and evaluate SVM classifiers
%   6. Save classifiers for use in SNAP or SNAP_batch
%
% WORKFLOW:
%   1. Run SNAP_classify()
%   2. Load exported .mat file with fit results
%   3. (Optional) Load associated image for visualization
%   4. Select features via "Select Features" button
%   5. Label spots using keyboard shortcuts or buttons:
%      - R = Real spot
%      - N = Noise
%      - S = Skip
%      - Left/Right arrows = Navigate
%      - Ctrl+Z = Undo
%   6. Train SVM when you have 5+ labels per class
%   7. Export classifier for SNAP/SNAP_batch integration
%
% USAGE:
%   SNAP_classify()
%
% KEYBOARD SHORTCUTS:
%   R         - Mark current spot as REAL
%   N         - Mark current spot as NOISE
%   S         - Skip current spot
%   Left/Right- Navigate between spots
%   Ctrl+Z    - Undo last label
%   Space     - Toggle live prediction overlay
%
% ============================================================================

    % Create main figure
    fig = uifigure('Name', 'SNAP Classify - SVM Training', ...
        'Position', [50 50 1400 900], ...
        'CloseRequestFcn', @onClose);
    
    % Initialize application state
    state = struct();
    state.fitData = [];
    state.selectedFeatures = {};
    state.customExpressions = struct('name', {}, 'expression', {});  % Custom expression definitions
    state.featureInfo = struct();
    state.labeledReal = [];
    state.labeledNoise = [];
    state.skippedSpots = [];
    state.currentIdx = 1;
    state.classifier = [];
    state.trainStats = [];
    state.normParams = struct('mu', [], 'sigma', [], 'standardized', true);  % Z-score normalization params
    state.imageData = [];
    state.dicImage = [];
    state.fittingMethod = '';
    state.has3D = false;
    state.undoStack = {};
    state.modified = false;
    state.showPredictions = false;
    state.filepath = '';
    state.cropData = [];
    
    % Store handles
    handles = struct();
    
    % Create UI
    createClassifyUI();
    
    % Store state in figure
    fig.UserData.state = state;
    fig.UserData.handles = handles;
    
    % ========================================================================
    % UI CREATION
    % ========================================================================
    function createClassifyUI()
        % Main grid layout
        mainGrid = uigridlayout(fig, [1 3]);
        mainGrid.ColumnWidth = {380, '1x', 380};
        mainGrid.Padding = [10 10 10 10];
        mainGrid.ColumnSpacing = 10;
        
        % === LEFT PANEL ===
        leftPanel = uipanel(mainGrid, 'Title', 'Data & Features');
        leftGrid = uigridlayout(leftPanel, [10 1]);
        leftGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x', 'fit', 'fit', 'fit', 'fit', 'fit'};
        leftGrid.Padding = [10 10 10 10];
        leftGrid.RowSpacing = 8;
        
        % Load data section
        loadPanel = uipanel(leftGrid, 'Title', 'Load Data');
        loadGrid = uigridlayout(loadPanel, [3 2]);
        loadGrid.ColumnWidth = {130, '1x'};
        loadGrid.RowHeight = {28, 28, 28};
        
        uibutton(loadGrid, 'Text', 'Load SNAP Export...', 'ButtonPushedFcn', @(~,~) loadSNAPExport());
        handles.dataStatusLabel = uilabel(loadGrid, 'Text', 'No data loaded', 'FontColor', [0.5 0.5 0.5]);
        
        uibutton(loadGrid, 'Text', 'Load Image...', 'ButtonPushedFcn', @(~,~) loadImageData());
        handles.imageStatusLabel = uilabel(loadGrid, 'Text', 'No image loaded', 'FontColor', [0.5 0.5 0.5]);
        
        uibutton(loadGrid, 'Text', 'Load DIC...', 'ButtonPushedFcn', @(~,~) loadDICImage());
        handles.dicStatusLabel = uilabel(loadGrid, 'Text', 'No DIC loaded', 'FontColor', [0.5 0.5 0.5]);
        
        % Fitting method display
        methodPanel = uipanel(leftGrid, 'Title', 'Detected Settings');
        methodGrid = uigridlayout(methodPanel, [2 2]);
        methodGrid.ColumnWidth = {'fit', '1x'};
        methodGrid.RowHeight = {22, 22};
        
        uilabel(methodGrid, 'Text', 'Fit Method:');
        handles.fittingMethodLabel = uilabel(methodGrid, 'Text', '-', 'FontColor', [0.3 0.3 0.3], 'FontWeight', 'bold');
        
        uilabel(methodGrid, 'Text', '3D Data:');
        handles.is3DLabel = uilabel(methodGrid, 'Text', '-', 'FontColor', [0.3 0.3 0.3], 'FontWeight', 'bold');
        
        % Feature selection
        featurePanel = uipanel(leftGrid, 'Title', 'Classification Features');
        featureGrid = uigridlayout(featurePanel, [2 1]);
        featureGrid.RowHeight = {30, 'fit'};
        
        handles.selectFeaturesBtn = uibutton(featureGrid, 'Text', 'Select Features...', ...
            'ButtonPushedFcn', @(~,~) selectFeatures(), 'Enable', 'off');
        handles.featureCountLabel = uilabel(featureGrid, 'Text', '0 features selected', 'FontColor', [0.5 0.5 0.5]);
        
        % Feature list
        handles.featureListArea = uitextarea(leftGrid, 'Editable', 'off', ...
            'Value', {'No features selected'}, 'FontSize', 10);
        
        % Statistics
        statsPanel = uipanel(leftGrid, 'Title', 'Labeling Statistics');
        statsGrid = uigridlayout(statsPanel, [3 2]);
        statsGrid.ColumnWidth = {'fit', '1x'};
        statsGrid.RowHeight = {22, 22, 22};
        
        uilabel(statsGrid, 'Text', 'Real spots:');
        handles.realCountLabel = uilabel(statsGrid, 'Text', '0', 'FontColor', [0.2 0.7 0.3], 'FontWeight', 'bold');
        
        uilabel(statsGrid, 'Text', 'Noise:');
        handles.noiseCountLabel = uilabel(statsGrid, 'Text', '0', 'FontColor', [0.8 0.2 0.2], 'FontWeight', 'bold');
        
        uilabel(statsGrid, 'Text', 'Skipped:');
        handles.skippedCountLabel = uilabel(statsGrid, 'Text', '0', 'FontColor', [0.5 0.5 0.5]);
        
        % Progress
        handles.progressLabel = uilabel(leftGrid, 'Text', 'Spot 0 / 0', ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        % Navigation
        navGrid = uigridlayout(leftGrid, [1 4]);
        navGrid.ColumnWidth = {'1x', '1x', '1x', '1x'};
        navGrid.Padding = [0 0 0 0];
        
        uibutton(navGrid, 'Text', '<<', 'ButtonPushedFcn', @(~,~) navigateSpot(-10));
        uibutton(navGrid, 'Text', '<', 'ButtonPushedFcn', @(~,~) navigateSpot(-1));
        uibutton(navGrid, 'Text', '>', 'ButtonPushedFcn', @(~,~) navigateSpot(1));
        uibutton(navGrid, 'Text', '>>', 'ButtonPushedFcn', @(~,~) navigateSpot(10));
        
        % Jump to spot
        jumpGrid = uigridlayout(leftGrid, [1 2]);
        jumpGrid.ColumnWidth = {'fit', '1x'};
        jumpGrid.Padding = [0 0 0 0];
        
        uilabel(jumpGrid, 'Text', 'Go to spot:');
        handles.jumpSpinner = uispinner(jumpGrid, 'Limits', [1 1], 'Value', 1, ...
            'ValueChangedFcn', @(src,~) jumpToSpot(src.Value));
        
        % === CENTER PANEL ===
        centerPanel = uipanel(mainGrid, 'Title', 'Spot Visualization');
        centerGrid = uigridlayout(centerPanel, [3 2]);
        centerGrid.RowHeight = {'1x', '1x', 30};
        centerGrid.ColumnWidth = {'1x', '1x'};
        centerGrid.Padding = [5 5 5 5];
        
        handles.ax3D = uiaxes(centerGrid);
        title(handles.ax3D, 'Z-Stack Crop');
        handles.ax3D.Layout.Row = 1;
        handles.ax3D.Layout.Column = 1;
        
        handles.axMaxProj = uiaxes(centerGrid);
        title(handles.axMaxProj, 'Max Projection (Crop)');
        handles.axMaxProj.Layout.Row = 1;
        handles.axMaxProj.Layout.Column = 2;
        
        handles.axContext = uiaxes(centerGrid);
        title(handles.axContext, 'Full Image Context');
        handles.axContext.Layout.Row = 2;
        handles.axContext.Layout.Column = [1 2];
        
        sliderGrid = uigridlayout(centerGrid, [1 3]);
        sliderGrid.ColumnWidth = {'fit', '1x', 'fit'};
        sliderGrid.Layout.Row = 3;
        sliderGrid.Layout.Column = [1 2];
        sliderGrid.Padding = [0 0 0 0];
        
        uilabel(sliderGrid, 'Text', 'Z:');
        handles.zSlider = uislider(sliderGrid, 'Limits', [1 10], 'Value', 1, ...
            'ValueChangedFcn', @(~,~) updateZDisplay());
        handles.zLabel = uilabel(sliderGrid, 'Text', '1/1');
        
        % === RIGHT PANEL ===
        rightPanel = uipanel(mainGrid, 'Title', 'Labeling & Training');
        rightGrid = uigridlayout(rightPanel, [12 1]);
        rightGrid.RowHeight = {70, 70, 35, 'fit', 'fit', '1x', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
        rightGrid.Padding = [10 10 10 10];
        rightGrid.RowSpacing = 8;
        
        handles.realBtn = uibutton(rightGrid, 'Text', 'REAL SPOT (R)', ...
            'FontSize', 18, 'FontWeight', 'bold', ...
            'BackgroundColor', [0.2 0.7 0.3], 'FontColor', 'white', ...
            'ButtonPushedFcn', @(~,~) labelSpot(1));
        
        handles.noiseBtn = uibutton(rightGrid, 'Text', 'NOISE (N)', ...
            'FontSize', 18, 'FontWeight', 'bold', ...
            'BackgroundColor', [0.8 0.2 0.2], 'FontColor', 'white', ...
            'ButtonPushedFcn', @(~,~) labelSpot(0));
        
        handles.skipBtn = uibutton(rightGrid, 'Text', 'Skip (S)', 'ButtonPushedFcn', @(~,~) skipSpot());
        
        handles.spotInfoArea = uitextarea(rightGrid, 'Editable', 'off', ...
            'Value', {'Spot Information:', '', 'Load data to begin.'}, ...
            'FontName', 'Consolas', 'FontSize', 10);
        
        handles.predictionLabel = uilabel(rightGrid, 'Text', 'Prediction: -', ...
            'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        
        trainPanel = uipanel(rightGrid, 'Title', 'SVM Training Parameters');
        trainGrid = uigridlayout(trainPanel, [9 2]);
        trainGrid.ColumnWidth = {'fit', '1x'};
        trainGrid.RowHeight = {22, 22, 22, 22, 22, 22, 22, 22, 35};
        trainGrid.RowSpacing = 4;
        trainGrid.Padding = [8 8 8 8];
        
        % Row 1: Kernel type
        uilabel(trainGrid, 'Text', 'Kernel:');
        handles.kernelDropdown = uidropdown(trainGrid, ...
            'Items', {'RBF (Gaussian)', 'Linear', 'Polynomial'}, 'Value', 'RBF (Gaussian)', ...
            'Tooltip', 'Kernel function: RBF is most common, Linear is faster, Polynomial for complex boundaries');
        
        % Row 2: Box Constraint (C parameter)
        uilabel(trainGrid, 'Text', 'Box Constraint (C):');
        handles.boxConstraintEdit = uieditfield(trainGrid, 'numeric', 'Value', 1.0, ...
            'Limits', [0.001 10000], 'Tooltip', 'Regularization: Higher = less regularization, risk of overfitting. Try 0.1, 1, 10, 100');
        
        % Row 3: Kernel Scale (gamma for RBF)
        uilabel(trainGrid, 'Text', 'Kernel Scale:');
        kernelScaleGrid = uigridlayout(trainGrid, [1 2]);
        kernelScaleGrid.ColumnWidth = {'fit', '1x'};
        kernelScaleGrid.Padding = [0 0 0 0];
        handles.kernelScaleAutoCheck = uicheckbox(kernelScaleGrid, 'Text', 'Auto', 'Value', true, ...
            'ValueChangedFcn', @(src,~) toggleKernelScaleAuto(src));
        handles.kernelScaleEdit = uieditfield(kernelScaleGrid, 'numeric', 'Value', 1.0, ...
            'Limits', [0.0001 1000], 'Enable', 'off', 'Tooltip', 'RBF gamma: smaller = wider influence, larger = tighter fit');
        
        % Row 4: Polynomial order (only for polynomial kernel)
        uilabel(trainGrid, 'Text', 'Polynomial Order:');
        handles.polyOrderSpinner = uispinner(trainGrid, 'Value', 3, 'Limits', [2 5], 'Step', 1, ...
            'Tooltip', 'Order for polynomial kernel (2-5)', 'Enable', 'off');
        
        % Row 5: Z-score normalization
        uilabel(trainGrid, 'Text', 'Z-Score Normalize:');
        handles.standardizeCheck = uicheckbox(trainGrid, 'Text', 'Always recommended', 'Value', true, ...
            'Tooltip', 'Normalize features to mean=0, std=1. Critical for SVM performance');
        
        % Row 6: Cross-validation
        uilabel(trainGrid, 'Text', 'Cross-Validation:');
        cvGrid = uigridlayout(trainGrid, [1 2]);
        cvGrid.ColumnWidth = {'fit', '1x'};
        cvGrid.Padding = [0 0 0 0];
        handles.cvCheck = uicheckbox(cvGrid, 'Value', true, 'Text', '');
        handles.kFoldSpinner = uispinner(cvGrid, 'Value', 5, 'Limits', [2 10], 'Step', 1, ...
            'Tooltip', 'Number of cross-validation folds');
        
        % Row 7: Status
        uilabel(trainGrid, 'Text', 'Status:');
        handles.trainStatusLabel = uilabel(trainGrid, 'Text', 'Not trained', 'FontColor', [0.5 0.5 0.5]);
        
        % Row 8: Accuracy
        uilabel(trainGrid, 'Text', 'Accuracy:');
        handles.accuracyLabel = uilabel(trainGrid, 'Text', '-', 'FontWeight', 'bold');
        
        % Row 9: Train button
        handles.trainBtn = uibutton(trainGrid, 'Text', 'Train SVM', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) trainSVM(), 'Enable', 'off');
        handles.trainBtn.Layout.Column = [1 2];
        
        % Callback for kernel dropdown to enable/disable polynomial order
        handles.kernelDropdown.ValueChangedFcn = @(src,~) updateKernelOptions(src);
        
        handles.liveUpdateBtn = uibutton(rightGrid, 'Text', 'Toggle Live Predictions (Space)', ...
            'ButtonPushedFcn', @(~,~) toggleLivePredictions());
        
        handles.undoBtn = uibutton(rightGrid, 'Text', 'Undo (Ctrl+Z)', ...
            'ButtonPushedFcn', @(~,~) undoLabel(), 'Enable', 'off');
        
        progressBtnGrid = uigridlayout(rightGrid, [1 2]);
        progressBtnGrid.ColumnWidth = {'1x', '1x'};
        progressBtnGrid.Padding = [0 0 0 0];
        
        uibutton(progressBtnGrid, 'Text', 'Save Progress', 'ButtonPushedFcn', @(~,~) saveProgress());
        uibutton(progressBtnGrid, 'Text', 'Load Progress', 'ButtonPushedFcn', @(~,~) loadProgress());
        
        handles.exportBtn = uibutton(rightGrid, 'Text', 'Export Classifier for SNAP...', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) exportClassifier(), 'Enable', 'off');
        
        % Keyboard shortcuts
        fig.KeyPressFcn = @handleKeyPress;
    end
    
    % ========================================================================
    % DATA LOADING
    % ========================================================================
    
    function loadSNAPExport()
        [file, path] = uigetfile({'*.mat', 'SNAP Export Files (*.mat)'}, 'Select SNAP Export File');
        if file == 0, return; end
        
        filepath = fullfile(path, file);
        
        try
            data = load(filepath);
            fitData = findFitResults(data);
            
            if isempty(fitData)
                uialert(fig, 'Could not find fit results in file.', 'Load Error');
                return;
            end
            
            state.fitData = fitData(:);
            state.filepath = filepath;
            
            % Try to get fitting method from metadata/parameters first
            state.fittingMethod = '';
            if isfield(data, 'parameters') && isfield(data.parameters, 'fit_method')
                state.fittingMethod = data.parameters.fit_method;
                fprintf('Fitting method from parameters: %s\n', state.fittingMethod);
            elseif isfield(data, 'metadata') && isfield(data.metadata, 'fit_method')
                state.fittingMethod = data.metadata.fit_method;
                fprintf('Fitting method from metadata: %s\n', state.fittingMethod);
            end
            
            % Fall back to detection from data if not in metadata
            if isempty(state.fittingMethod)
                state.fittingMethod = detectFittingMethod(fitData);
            end
            
            state.has3D = detectIs3D(fitData);
            
            [~, state.featureInfo] = snap_helpers.classification.getAvailableFeatures(...
                state.fittingMethod, state.has3D, false);
            
            handles.dataStatusLabel.Text = sprintf('%d spots loaded', numel(fitData));
            handles.dataStatusLabel.FontColor = [0.2 0.6 0.2];
            handles.fittingMethodLabel.Text = state.fittingMethod;
            handles.is3DLabel.Text = string(state.has3D);
            handles.jumpSpinner.Limits = [1 max(1, numel(fitData))];
            handles.selectFeaturesBtn.Enable = 'on';
            
            state.labeledReal = [];
            state.labeledNoise = [];
            state.skippedSpots = [];
            state.currentIdx = 1;
            state.classifier = [];
            state.undoStack = {};
            state.modified = false;
            
            fig.UserData.state = state;
            updateDisplay();
            updateStats();
            
            fprintf('Loaded %d spots from: %s\n', numel(fitData), file);
            fprintf('Fitting method: %s, 3D: %s\n', state.fittingMethod, string(state.has3D));
            
            % Remind user to load image if not already loaded
            if isempty(state.imageData)
                fprintf('TIP: Load the channel image for visualization (Browse button next to "Channel Image")\n');
            end
            
        catch ME
            uialert(fig, sprintf('Failed to load data: %s', ME.message), 'Load Error');
        end
    end
    
    function fitData = findFitResults(data)
        fitData = [];
        
        % Priority 1: 'signals' field (from exportChannelDataStandardized)
        % This is the standard SNAP/SNAP_batch export format
        if isfield(data, 'signals')
            fitData = data.signals;
            fprintf('Found signal data in "signals" field (standard SNAP export)\n');
        % Priority 2: Legacy 'fitResults' field
        elseif isfield(data, 'fitResults')
            fitData = data.fitResults;
            fprintf('Found fit data in "fitResults" field\n');
        elseif isfield(data, 'fit_results')
            fitData = data.fit_results;
            fprintf('Found fit data in "fit_results" field\n');
        % Priority 3: Nested in exportData
        elseif isfield(data, 'exportData')
            if isfield(data.exportData, 'signals')
                fitData = data.exportData.signals;
                fprintf('Found signal data in "exportData.signals" field\n');
            elseif isfield(data.exportData, 'fitResults')
                fitData = data.exportData.fitResults;
                fprintf('Found fit data in "exportData.fitResults" field\n');
            end
        % Priority 4: Other common names
        elseif isfield(data, 'spots')
            fitData = data.spots;
            fprintf('Found data in "spots" field\n');
        else
            % Search all fields for struct arrays with signal/fit-like fields
            fnames = fieldnames(data);
            fprintf('Available top-level fields: %s\n', strjoin(fnames, ', '));
            
            for i = 1:numel(fnames)
                candidate = data.(fnames{i});
                if isstruct(candidate) && numel(candidate) >= 1
                    % Check if it looks like signal/fit data
                    if isfield(candidate, 'fitted_coords') || isfield(candidate, 'maxima_coords') || ...
                       isfield(candidate, 'amplitude') || isfield(candidate, 'sigma_x') || ...
                       isfield(candidate, 'globalFitCenter') || isfield(candidate, 'signal_id')
                        fitData = candidate;
                        fprintf('Found signal/fit data in field: %s\n', fnames{i});
                        break;
                    end
                end
            end
        end
        
        % Convert table to struct if needed
        if istable(fitData)
            fitData = table2struct(fitData);
        end
        
        % Check if we found anything
        if isempty(fitData)
            fprintf('ERROR: Could not find signal/fit data in file.\n');
            fprintf('Expected fields: "signals" (SNAP export) or "fitResults" (legacy)\n');
            return;
        end
        
        fprintf('Loaded %d signals/spots\n', numel(fitData));
        
        % Normalize the fit data structure
        fitData = normalizeFitData(fitData);
    end
    
    function fitData = normalizeFitData(fitData)
        % Ensure consistent field names for position data
        % Report what we found first
        if ~isempty(fitData)
            fnames = fieldnames(fitData);
            fprintf('Fit data fields: %s\n', strjoin(fnames, ', '));
        end
        
        for i = 1:numel(fitData)
            spot = fitData(i);
            
            % Create globalFitCenter if it doesn't exist
            if ~isfield(spot, 'globalFitCenter') || isempty(spot.globalFitCenter)
                pos = [];
                
                % Priority 1: fitted_coords array (from exportChannelDataStandardized)
                if isfield(spot, 'fitted_coords') && ~isempty(spot.fitted_coords)
                    coords = spot.fitted_coords;
                    if numel(coords) >= 2
                        % fitted_coords is in [row, col, slice] = [y, x, z] format
                        pos = coords(:)';
                    end
                % Priority 2: maxima_coords array (fallback if no fitting)
                elseif isfield(spot, 'maxima_coords') && ~isempty(spot.maxima_coords)
                    coords = spot.maxima_coords;
                    if numel(coords) >= 2
                        % maxima_coords is in [row, col, slice] = [y, x, z] format
                        pos = coords(:)';
                    end
                % Priority 3: Individual coordinate fields
                elseif isfield(spot, 'fitted_y') && isfield(spot, 'fitted_x')
                    y = spot.fitted_y;
                    x = spot.fitted_x;
                    if isfield(spot, 'fitted_z')
                        z = spot.fitted_z;
                    else
                        z = 1;
                    end
                    pos = [y, x, z];
                elseif isfield(spot, 'y') && isfield(spot, 'x')
                    y = spot.y;
                    x = spot.x;
                    if isfield(spot, 'z')
                        z = spot.z;
                    else
                        z = 1;
                    end
                    pos = [y, x, z];
                elseif isfield(spot, 'maxima_y') && isfield(spot, 'maxima_x')
                    y = spot.maxima_y;
                    x = spot.maxima_x;
                    if isfield(spot, 'maxima_z')
                        z = spot.maxima_z;
                    else
                        z = 1;
                    end
                    pos = [y, x, z];
                elseif isfield(spot, 'row') && isfield(spot, 'col')
                    y = spot.row;
                    x = spot.col;
                    if isfield(spot, 'slice')
                        z = spot.slice;
                    else
                        z = 1;
                    end
                    pos = [y, x, z];
                elseif isfield(spot, 'center') && ~isempty(spot.center)
                    pos = spot.center(:)';
                elseif isfield(spot, 'position') && ~isempty(spot.position)
                    pos = spot.position(:)';
                end
                
                % Ensure pos has at least 3 elements [y, x, z]
                if ~isempty(pos)
                    if numel(pos) == 2
                        pos = [pos, 1];
                    end
                    fitData(i).globalFitCenter = pos(1:min(3, numel(pos)));
                end
            end
            
            % Normalize r_squared field name
            if ~isfield(spot, 'r_squared') && isfield(spot, 'rsquared')
                fitData(i).r_squared = spot.rsquared;
            end
        end
        
        % Report position info
        if ~isempty(fitData) && isfield(fitData, 'globalFitCenter') && ~isempty(fitData(1).globalFitCenter)
            pos = fitData(1).globalFitCenter;
            fprintf('First spot position (y,x,z): [%.1f, %.1f, %.1f]\n', pos(1), pos(2), pos(min(3,numel(pos))));
        else
            fprintf('WARNING: Could not extract position data from spots.\n');
        end
    end
    
    function method = detectFittingMethod(fitData)
        if isfield(fitData, 'fitMethod') && ~isempty(fitData(1).fitMethod)
            method = fitData(1).fitMethod;
        elseif isfield(fitData, 'fit_method') && ~isempty(fitData(1).fit_method)
            method = fitData(1).fit_method;
        elseif isfield(fitData, 'rho_xy') && any(~isnan([fitData.rho_xy]))
            if isfield(fitData, 'alpha_x') && any(~isnan([fitData.alpha_x]))
                method = 'Skewed 3D Gaussian';
            else
                method = 'Distorted 3D Gaussian';
            end
        elseif isfield(fitData, 'radial_symmetry_score') && any(~isnan([fitData.radial_symmetry_score]))
            method = 'Radial Symmetry';
        elseif isfield(fitData, 'radialSymmetryScore') && any(~isnan([fitData.radialSymmetryScore]))
            method = 'Radial Symmetry';
        elseif isfield(fitData, 'amplitude_xy') && any(~isnan([fitData.amplitude_xy]))
            method = '2D (XY) + 1D (Z)';
        elseif isfield(fitData, 'amplitude_x') && any(~isnan([fitData.amplitude_x]))
            method = '1D (X,Y,Z)';
        elseif isfield(fitData, 'sigma_z') && any(~isnan([fitData.sigma_z]))
            method = '3D Gaussian';
        elseif isfield(fitData, 'sigma_x') && isfield(fitData, 'sigma_y')
            % Has sigma_x and sigma_y but no sigma_z - could be 2D or 3D
            if isfield(fitData, 'sigma_z') && any(~isnan([fitData.sigma_z]))
                method = '3D Gaussian';
            else
                method = '2D Gaussian';
            end
        else
            method = 'Unknown';
        end
        
        fprintf('Detected fitting method: %s\n', method);
    end
    
    function is3D = detectIs3D(fitData)
        is3D = false;
        if isfield(fitData, 'sigma_z') && any(~isnan([fitData.sigma_z]))
            is3D = true;
        elseif isfield(fitData, 'globalFitCenter') && ~isempty(fitData(1).globalFitCenter)
            coords = fitData(1).globalFitCenter;
            if numel(coords) >= 3 && ~isnan(coords(3)) && coords(3) > 1
                is3D = true;
            end
        elseif isfield(fitData, 'fitted_coords') && ~isempty(fitData(1).fitted_coords)
            coords = fitData(1).fitted_coords;
            if numel(coords) >= 3 && ~isnan(coords(3)) && coords(3) > 1
                is3D = true;
            end
        elseif isfield(fitData, 'maxima_coords') && ~isempty(fitData(1).maxima_coords)
            coords = fitData(1).maxima_coords;
            if numel(coords) >= 3 && ~isnan(coords(3)) && coords(3) > 1
                is3D = true;
            end
        end
    end
    
    function loadImageData()
        [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF Images'; '*.*', 'All Files'}, 'Select Channel Image');
        if file == 0, return; end
        
        try
            imgPath = fullfile(path, file);
            info = imfinfo(imgPath);
            
            % Read first frame to get data type
            firstFrame = imread(imgPath, 1);
            
            if numel(info) > 1
                % Multi-page TIFF - allocate with correct type
                state.imageData = zeros(info(1).Height, info(1).Width, numel(info), class(firstFrame));
                state.imageData(:,:,1) = firstFrame;
                for z = 2:numel(info)
                    state.imageData(:,:,z) = imread(imgPath, z);
                end
            else
                state.imageData = firstFrame;
            end
            
            % Update UI status
            handles.imageStatusLabel.Text = sprintf('%dx%dx%d loaded', size(state.imageData, 1), size(state.imageData, 2), size(state.imageData, 3));
            handles.imageStatusLabel.FontColor = [0.2 0.6 0.2];
            
            % Update Z slider
            nSlices = size(state.imageData, 3);
            if nSlices > 1
                handles.zSlider.Limits = [1 nSlices];
                handles.zSlider.Value = ceil(nSlices / 2);
                handles.zSlider.Enable = 'on';
            else
                handles.zSlider.Limits = [1 1];
                handles.zSlider.Value = 1;
            end
            handles.zLabel.Text = sprintf('%d/%d', round(handles.zSlider.Value), nSlices);
            
            % Sync state and refresh display
            fig.UserData.state = state;
            
            % Force refresh of visualization
            if ~isempty(state.fitData)
                updateDisplay();
            else
                % Even without fit data, show the image in context
                fullMaxProj = max(double(state.imageData), [], 3);
                lims = computeContrastLimits(fullMaxProj);
                imagesc(handles.axContext, fullMaxProj, lims);
                colormap(handles.axContext, 'gray');
                axis(handles.axContext, 'image');
                title(handles.axContext, 'Image loaded - Now load spot data');
            end
            
            fprintf('Image loaded: %s (%dx%dx%d, %s)\n', file, ...
                size(state.imageData, 1), size(state.imageData, 2), size(state.imageData, 3), class(state.imageData));
            
        catch ME
            uialert(fig, sprintf('Failed to load image: %s', ME.message), 'Load Error');
            fprintf('Image load error: %s\n', ME.message);
        end
    end
    
    function loadDICImage()
        [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF Images'; '*.*', 'All Files'}, 'Select DIC Image');
        if file == 0, return; end
        
        try
            state.dicImage = imread(fullfile(path, file));
            handles.dicStatusLabel.Text = sprintf('%dx%d loaded', size(state.dicImage, 1), size(state.dicImage, 2));
            handles.dicStatusLabel.FontColor = [0.2 0.6 0.2];
            fig.UserData.state = state;
            
            % Force refresh
            if ~isempty(state.fitData)
                updateDisplay();
            end
            
            fprintf('DIC image loaded: %s (%dx%d)\n', file, size(state.dicImage, 1), size(state.dicImage, 2));
        catch ME
            uialert(fig, sprintf('Failed to load DIC: %s', ME.message), 'Load Error');
        end
    end
    
    % ========================================================================
    % FEATURE SELECTION
    % ========================================================================
    
    function selectFeatures()
        if isempty(state.fitData)
            uialert(fig, 'Load data first', 'No Data');
            return;
        end
        
        [selected, customExpr, cancelled] = snap_helpers.classification.featureSelectionUI(...
            state.fittingMethod, state.has3D, false, state.selectedFeatures, state.customExpressions);
        
        if ~cancelled
            state.selectedFeatures = selected;
            state.customExpressions = customExpr;
            
            % Count total features
            nBase = numel(selected);
            nCustom = numel(customExpr);
            nTotal = nBase + nCustom;
            
            handles.featureCountLabel.Text = sprintf('%d base + %d custom = %d features', nBase, nCustom, nTotal);
            handles.featureCountLabel.FontColor = [0.2 0.6 0.2];
            
            % Update feature list display
            featureList = selected;
            for i = 1:numel(customExpr)
                featureList{end+1} = sprintf('[EXPR] %s = %s', customExpr(i).name, customExpr(i).expression);
            end
            handles.featureListArea.Value = featureList;
            
            fig.UserData.state = state;
            checkTrainingEligibility();
        end
    end
    
    % ========================================================================
    % LABELING
    % ========================================================================
    
    function labelSpot(label)
        if isempty(state.fitData), return; end
        
        pushUndo();
        idx = state.currentIdx;
        
        state.labeledReal(state.labeledReal == idx) = [];
        state.labeledNoise(state.labeledNoise == idx) = [];
        state.skippedSpots(state.skippedSpots == idx) = [];
        
        if label == 1
            state.labeledReal(end+1) = idx;
        else
            state.labeledNoise(end+1) = idx;
        end
        
        state.modified = true;
        fig.UserData.state = state;
        
        advanceToNextUnlabeled();
        updateDisplay();
        updateStats();
        checkTrainingEligibility();
    end
    
    function skipSpot()
        if isempty(state.fitData), return; end
        
        pushUndo();
        idx = state.currentIdx;
        if ~ismember(idx, state.skippedSpots)
            state.skippedSpots(end+1) = idx;
        end
        
        fig.UserData.state = state;
        advanceToNextUnlabeled();
        updateDisplay();
        updateStats();
    end
    
    function advanceToNextUnlabeled()
        if isempty(state.fitData), return; end
        
        nSpots = numel(state.fitData);
        labeled = [state.labeledReal, state.labeledNoise, state.skippedSpots];
        
        for i = 1:nSpots
            nextIdx = mod(state.currentIdx + i - 1, nSpots) + 1;
            if ~ismember(nextIdx, labeled)
                state.currentIdx = nextIdx;
                fig.UserData.state = state;
                return;
            end
        end
        
        uialert(fig, 'All spots have been labeled or skipped!', 'Complete', 'Icon', 'success');
    end
    
    function navigateSpot(delta)
        if isempty(state.fitData), return; end
        nSpots = numel(state.fitData);
        state.currentIdx = max(1, min(nSpots, state.currentIdx + delta));
        fig.UserData.state = state;
        updateDisplay();
    end
    
    function jumpToSpot(idx)
        if isempty(state.fitData), return; end
        state.currentIdx = max(1, min(numel(state.fitData), idx));
        fig.UserData.state = state;
        updateDisplay();
    end
    
    function pushUndo()
        undoState = struct();
        undoState.labeledReal = state.labeledReal;
        undoState.labeledNoise = state.labeledNoise;
        undoState.skippedSpots = state.skippedSpots;
        undoState.currentIdx = state.currentIdx;
        state.undoStack{end+1} = undoState;
        
        if numel(state.undoStack) > 100
            state.undoStack = state.undoStack(end-99:end);
        end
        
        handles.undoBtn.Enable = 'on';
        fig.UserData.state = state;
    end
    
    function undoLabel()
        if isempty(state.undoStack), return; end
        
        undoState = state.undoStack{end};
        state.undoStack(end) = [];
        
        state.labeledReal = undoState.labeledReal;
        state.labeledNoise = undoState.labeledNoise;
        state.skippedSpots = undoState.skippedSpots;
        state.currentIdx = undoState.currentIdx;
        
        if isempty(state.undoStack)
            handles.undoBtn.Enable = 'off';
        end
        
        fig.UserData.state = state;
        updateDisplay();
        updateStats();
        checkTrainingEligibility();
    end
    
    % ========================================================================
    % DISPLAY
    % ========================================================================
    
    function updateDisplay()
        if isempty(state.fitData), return; end
        
        nSpots = numel(state.fitData);
        idx = state.currentIdx;
        spot = state.fitData(idx);
        
        handles.progressLabel.Text = sprintf('Spot %d / %d', idx, nSpots);
        handles.jumpSpinner.Value = idx;
        
        % Spot info
        infoLines = {sprintf('=== Spot #%d ===', idx), ''};
        
        % Try to display position from various field names
        pos = [];
        if isfield(spot, 'globalFitCenter') && ~isempty(spot.globalFitCenter)
            pos = spot.globalFitCenter;
        elseif isfield(spot, 'fitted_coords') && ~isempty(spot.fitted_coords)
            pos = spot.fitted_coords(:)';
        elseif isfield(spot, 'maxima_coords') && ~isempty(spot.maxima_coords)
            pos = spot.maxima_coords(:)';
        elseif isfield(spot, 'fitted_y') && isfield(spot, 'fitted_x')
            if isfield(spot, 'fitted_z')
                pos = [spot.fitted_y, spot.fitted_x, spot.fitted_z];
            else
                pos = [spot.fitted_y, spot.fitted_x];
            end
        elseif isfield(spot, 'maxima_y') && isfield(spot, 'maxima_x')
            if isfield(spot, 'maxima_z')
                pos = [spot.maxima_y, spot.maxima_x, spot.maxima_z];
            else
                pos = [spot.maxima_y, spot.maxima_x];
            end
        end
        
        if ~isempty(pos) && any(~isnan(pos))
            if numel(pos) >= 3 && ~isnan(pos(3))
                infoLines{end+1} = sprintf('Position: [%.1f, %.1f, %.1f]', pos(2), pos(1), pos(3));
            else
                infoLines{end+1} = sprintf('Position: [%.1f, %.1f]', pos(2), pos(1));
            end
        else
            infoLines{end+1} = 'Position: (not found)';
        end
        
        % Extended list of metrics to show
        metricsToShow = {'amplitude', 'sigma_x', 'sigma_y', 'sigma_z', 'r_squared', 'rsquared', ...
            'integratedIntensity', 'integrated_intensity', 'background', 'radialSymmetryScore', ...
            'radial_symmetry_score', 'fitted_x', 'fitted_y', 'fitted_z', 'maxima_x', 'maxima_y', 'maxima_z'};
        shownFields = {};
        for i = 1:numel(metricsToShow)
            fname = metricsToShow{i};
            if isfield(spot, fname) && ~ismember(fname, shownFields)
                val = spot.(fname);
                if isnumeric(val) && isscalar(val) && ~isnan(val)
                    infoLines{end+1} = sprintf('%s: %.4g', fname, val);
                    shownFields{end+1} = fname;
                elseif isnumeric(val) && numel(val) <= 3 && any(~isnan(val))
                    infoLines{end+1} = sprintf('%s: [%s]', fname, num2str(val, '%.3f '));
                    shownFields{end+1} = fname;
                end
            end
        end
        
        % If no metrics shown, list available fields
        if numel(shownFields) == 0
            infoLines{end+1} = '';
            infoLines{end+1} = 'Available fields:';
            fnames = fieldnames(spot);
            for i = 1:min(10, numel(fnames))
                infoLines{end+1} = sprintf('  %s', fnames{i});
            end
            if numel(fnames) > 10
                infoLines{end+1} = sprintf('  ... and %d more', numel(fnames) - 10);
            end
        end
        
        handles.spotInfoArea.Value = infoLines;
        
        updatePredictionDisplay(spot);
        visualizeSpot(spot);
        updateLabelButtonHighlights(idx);
    end
    
    function updatePredictionDisplay(spot)
        hasFeatures = ~isempty(state.selectedFeatures) || ~isempty(state.customExpressions);
        if ~isempty(state.classifier) && hasFeatures
            try
                [X, featureNames, ~] = snap_helpers.classification.buildFeatureMatrix(...
                    spot, state.selectedFeatures, state.featureInfo, state.customExpressions);
                
                if all(~isnan(X))
                    % Apply classifier with normalization parameters
                    [pred, ~, conf] = snap_helpers.classification.applyClassifier(...
                        state.classifier, X, featureNames, featureNames, state.normParams);
                    
                    if pred == 1
                        handles.predictionLabel.Text = sprintf('Prediction: REAL (%.0f%% conf)', conf*100);
                        handles.predictionLabel.FontColor = [0.2 0.7 0.3];
                    else
                        handles.predictionLabel.Text = sprintf('Prediction: NOISE (%.0f%% conf)', conf*100);
                        handles.predictionLabel.FontColor = [0.8 0.2 0.2];
                    end
                else
                    handles.predictionLabel.Text = 'Prediction: (missing features)';
                    handles.predictionLabel.FontColor = [0.5 0.5 0.5];
                end
            catch
                handles.predictionLabel.Text = 'Prediction: (error)';
                handles.predictionLabel.FontColor = [0.5 0.5 0.5];
            end
        else
            handles.predictionLabel.Text = 'Prediction: (train classifier first)';
            handles.predictionLabel.FontColor = [0.5 0.5 0.5];
        end
    end
    
    function visualizeSpot(spot)
        % Try to get spot position from various possible field names
        pos = [];
        
        % Check different possible position field names
        % Priority 1: globalFitCenter (our normalized field)
        if isfield(spot, 'globalFitCenter') && ~isempty(spot.globalFitCenter)
            pos = spot.globalFitCenter;
        % Priority 2: fitted_coords array (from exportChannelDataStandardized)
        elseif isfield(spot, 'fitted_coords') && ~isempty(spot.fitted_coords)
            pos = spot.fitted_coords(:)';
        % Priority 3: maxima_coords array
        elseif isfield(spot, 'maxima_coords') && ~isempty(spot.maxima_coords)
            pos = spot.maxima_coords(:)';
        % Priority 4: Other common formats
        elseif isfield(spot, 'fitted_position') && ~isempty(spot.fitted_position)
            pos = spot.fitted_position(:)';
        elseif isfield(spot, 'position') && ~isempty(spot.position)
            pos = spot.position(:)';
        elseif isfield(spot, 'center') && ~isempty(spot.center)
            pos = spot.center(:)';
        elseif isfield(spot, 'fitCenter') && ~isempty(spot.fitCenter)
            pos = spot.fitCenter(:)';
        % Priority 5: Individual coordinate fields
        elseif isfield(spot, 'fitted_y') && isfield(spot, 'fitted_x')
            if isfield(spot, 'fitted_z')
                pos = [spot.fitted_y, spot.fitted_x, spot.fitted_z];
            else
                pos = [spot.fitted_y, spot.fitted_x, 1];
            end
        elseif isfield(spot, 'y') && isfield(spot, 'x')
            if isfield(spot, 'z')
                pos = [spot.y, spot.x, spot.z];
            else
                pos = [spot.y, spot.x, 1];
            end
        elseif isfield(spot, 'maxima_y') && isfield(spot, 'maxima_x')
            if isfield(spot, 'maxima_z')
                pos = [spot.maxima_y, spot.maxima_x, spot.maxima_z];
            else
                pos = [spot.maxima_y, spot.maxima_x, 1];
            end
        elseif isfield(spot, 'row') && isfield(spot, 'col')
            if isfield(spot, 'slice')
                pos = [spot.row, spot.col, spot.slice];
            else
                pos = [spot.row, spot.col, 1];
            end
        end
        
        % Validate position
        if isempty(pos) || all(isnan(pos(1:min(2,numel(pos)))))
            % No valid position - show debug info
            cla(handles.ax3D);
            cla(handles.axMaxProj);
            cla(handles.axContext);
            
            % Show what fields ARE available
            fnames = fieldnames(spot);
            title(handles.ax3D, 'No position found');
            title(handles.axMaxProj, sprintf('Fields: %s', strjoin(fnames(1:min(5,numel(fnames))), ', ')));
            title(handles.axContext, 'Check console for field names');
            
            % Print debug info to console
            fprintf('WARNING: Spot has no recognized position field.\n');
            fprintf('Available fields: %s\n', strjoin(fnames, ', '));
            return;
        end
        
        % Ensure pos has at least 2 elements
        if numel(pos) < 2
            fprintf('WARNING: Position has fewer than 2 coordinates: %s\n', num2str(pos));
            return;
        end
        
        % Extract coordinates (pos is in [row, col, slice] = [y, x, z] format)
        x = round(pos(2));
        y = round(pos(1));
        if numel(pos) >= 3 && ~isnan(pos(3))
            z = round(pos(3));
        else
            z = 1;
        end
        
        windowSize = 15;
        
        if ~isempty(state.imageData)
            img = double(state.imageData);  % Convert to double for proper scaling
            [h, w, d] = size(img);
            
            x = max(1, min(w, x));
            y = max(1, min(h, y));
            z = max(1, min(d, z));
            
            x1 = max(1, x - windowSize);
            x2 = min(w, x + windowSize);
            y1 = max(1, y - windowSize);
            y2 = min(h, y + windowSize);
            z1 = max(1, z - windowSize);
            z2 = min(d, z + windowSize);
            
            crop3D = img(y1:y2, x1:x2, z1:z2);
            state.cropData = struct('crop', crop3D, 'centerX', x - x1 + 1, 'centerY', y - y1 + 1, 'z1', z1);
            fig.UserData.state = state;
            
            nCropSlices = size(crop3D, 3);
            handles.zSlider.Limits = [1 max(1, nCropSlices)];
            handles.zSlider.Value = min(handles.zSlider.Value, nCropSlices);
            handles.zLabel.Text = sprintf('%d/%d', round(handles.zSlider.Value), nCropSlices);
            
            currentZ = max(1, min(nCropSlices, round(handles.zSlider.Value)));
            
            % Calculate contrast limits for crop (percentile-based)
            cropSlice = crop3D(:,:,currentZ);
            cropLims = computeContrastLimits(cropSlice);
            
            imagesc(handles.ax3D, cropSlice, cropLims);
            colormap(handles.ax3D, 'hot');
            axis(handles.ax3D, 'image');
            title(handles.ax3D, sprintf('Z-slice (z=%d)', z1 + currentZ - 1));
            hold(handles.ax3D, 'on');
            plot(handles.ax3D, x - x1 + 1, y - y1 + 1, 'g+', 'MarkerSize', 20, 'LineWidth', 2);
            hold(handles.ax3D, 'off');
            
            maxProj = max(crop3D, [], 3);
            maxProjLims = computeContrastLimits(maxProj);
            imagesc(handles.axMaxProj, maxProj, maxProjLims);
            colormap(handles.axMaxProj, 'hot');
            axis(handles.axMaxProj, 'image');
            title(handles.axMaxProj, 'Max Projection (Crop)');
            hold(handles.axMaxProj, 'on');
            plot(handles.axMaxProj, x - x1 + 1, y - y1 + 1, 'g+', 'MarkerSize', 20, 'LineWidth', 2);
            hold(handles.axMaxProj, 'off');
            
            % Context image - compute limits once for full image
            fullMaxProj = max(img, [], 3);
            contextLims = computeContrastLimits(fullMaxProj);
            imagesc(handles.axContext, fullMaxProj, contextLims);
            colormap(handles.axContext, 'gray');
            axis(handles.axContext, 'image');
            title(handles.axContext, 'Full Image - Spot Location');
            hold(handles.axContext, 'on');
            plot(handles.axContext, x, y, 'go', 'MarkerSize', 12, 'LineWidth', 2);
            
            for i = 1:numel(state.labeledReal)
                rpos = state.fitData(state.labeledReal(i)).globalFitCenter;
                plot(handles.axContext, rpos(2), rpos(1), 'g.', 'MarkerSize', 10);
            end
            for i = 1:numel(state.labeledNoise)
                npos = state.fitData(state.labeledNoise(i)).globalFitCenter;
                plot(handles.axContext, npos(2), npos(1), 'r.', 'MarkerSize', 10);
            end
            
            if state.showPredictions && ~isempty(state.classifier)
                showLivePredictions(handles.axContext);
            end
            hold(handles.axContext, 'off');
            
        elseif ~isempty(state.dicImage)
            dicImg = double(state.dicImage);
            dicLims = computeContrastLimits(dicImg);
            imagesc(handles.axContext, dicImg, dicLims);
            colormap(handles.axContext, 'gray');
            axis(handles.axContext, 'image');
            title(handles.axContext, 'DIC Image');
            hold(handles.axContext, 'on');
            plot(handles.axContext, x, y, 'go', 'MarkerSize', 12, 'LineWidth', 2);
            hold(handles.axContext, 'off');
            
            % Clear other axes when only DIC is available
            cla(handles.ax3D);
            cla(handles.axMaxProj);
            title(handles.ax3D, 'Load channel image for Z-view');
            title(handles.axMaxProj, 'Load channel image');
        else
            % No image loaded - show placeholder with spot info
            cla(handles.ax3D);
            cla(handles.axMaxProj);
            cla(handles.axContext);
            title(handles.ax3D, 'No image loaded');
            title(handles.axMaxProj, 'No image loaded');
            title(handles.axContext, sprintf('Spot at [%.1f, %.1f] - Load image for visualization', x, y));
        end
    end
    
    function lims = computeContrastLimits(data)
        % Compute display limits using percentiles for good contrast
        data = double(data(:));
        data = data(~isnan(data) & ~isinf(data));
        if isempty(data)
            lims = [0 1];
            return;
        end
        
        % Use 1st and 99.5th percentile for contrast
        lowVal = prctile(data, 1);
        highVal = prctile(data, 99.5);
        
        % Ensure valid limits
        if highVal <= lowVal
            highVal = lowVal + 1;
        end
        
        lims = [lowVal, highVal];
    end
    
    function updateZDisplay()
        if ~isempty(state.cropData)
            crop = double(state.cropData.crop);
            centerX = state.cropData.centerX;
            centerY = state.cropData.centerY;
            
            currentZ = max(1, min(size(crop, 3), round(handles.zSlider.Value)));
            
            cropSlice = crop(:,:,currentZ);
            lims = computeContrastLimits(cropSlice);
            imagesc(handles.ax3D, cropSlice, lims);
            colormap(handles.ax3D, 'hot');
            axis(handles.ax3D, 'image');
            hold(handles.ax3D, 'on');
            plot(handles.ax3D, centerX, centerY, 'g+', 'MarkerSize', 20, 'LineWidth', 2);
            hold(handles.ax3D, 'off');
            
            handles.zLabel.Text = sprintf('%d/%d', currentZ, size(crop, 3));
        end
    end
    
    function updateLabelButtonHighlights(idx)
        if ismember(idx, state.labeledReal)
            handles.realBtn.BackgroundColor = [0.1 0.9 0.2];
        else
            handles.realBtn.BackgroundColor = [0.2 0.7 0.3];
        end
        
        if ismember(idx, state.labeledNoise)
            handles.noiseBtn.BackgroundColor = [0.9 0.1 0.1];
        else
            handles.noiseBtn.BackgroundColor = [0.8 0.2 0.2];
        end
    end
    
    function updateStats()
        handles.realCountLabel.Text = num2str(numel(state.labeledReal));
        handles.noiseCountLabel.Text = num2str(numel(state.labeledNoise));
        handles.skippedCountLabel.Text = num2str(numel(state.skippedSpots));
    end
    
    function checkTrainingEligibility()
        nReal = numel(state.labeledReal);
        nNoise = numel(state.labeledNoise);
        nBaseFeatures = numel(state.selectedFeatures);
        nCustomFeatures = numel(state.customExpressions);
        nTotalFeatures = nBaseFeatures + nCustomFeatures;
        
        if nReal >= 5 && nNoise >= 5 && nTotalFeatures >= 1
            handles.trainBtn.Enable = 'on';
        else
            handles.trainBtn.Enable = 'off';
        end
    end
    
    function toggleLivePredictions()
        state.showPredictions = ~state.showPredictions;
        fig.UserData.state = state;
        
        if state.showPredictions
            handles.liveUpdateBtn.BackgroundColor = [0.3 0.5 0.7];
        else
            handles.liveUpdateBtn.BackgroundColor = [0.94 0.94 0.94];
        end
        
        updateDisplay();
    end
    
    function showLivePredictions(ax)
        hasFeatures = ~isempty(state.selectedFeatures) || ~isempty(state.customExpressions);
        if isempty(state.classifier) || ~hasFeatures
            return;
        end
        
        labeled = [state.labeledReal, state.labeledNoise, state.skippedSpots];
        unlabeled = setdiff(1:numel(state.fitData), labeled);
        
        if isempty(unlabeled), return; end
        
        try
            [X, featureNames, validMask] = snap_helpers.classification.buildFeatureMatrix(...
                state.fitData(unlabeled), state.selectedFeatures, state.featureInfo, state.customExpressions);
            
            % Apply classifier with normalization parameters
            [predictions, ~, ~] = snap_helpers.classification.applyClassifier(...
                state.classifier, X, featureNames, featureNames, state.normParams);
            
            for i = 1:numel(unlabeled)
                if validMask(i)
                    pos = state.fitData(unlabeled(i)).globalFitCenter;
                    if predictions(i) == 1
                        plot(ax, pos(2), pos(1), 'g^', 'MarkerSize', 6);
                    else
                        plot(ax, pos(2), pos(1), 'rv', 'MarkerSize', 6);
                    end
                end
            end
        catch
            % Silently fail
        end
    end
    
    % ========================================================================
    % TRAINING
    % ========================================================================
    
    function trainSVM()
        if isempty(state.selectedFeatures) && isempty(state.customExpressions)
            uialert(fig, 'Please select features first', 'No Features');
            return;
        end
        
        if numel(state.labeledReal) < 5 || numel(state.labeledNoise) < 5
            uialert(fig, 'Need at least 5 labels per class', 'Insufficient Labels');
            return;
        end
        
        realData = state.fitData(state.labeledReal);
        noiseData = state.fitData(state.labeledNoise);
        
        % Pass custom expressions to feature matrix builder
        [Xreal, ~, validReal] = snap_helpers.classification.buildFeatureMatrix(...
            realData, state.selectedFeatures, state.featureInfo, state.customExpressions);
        [Xnoise, ~, validNoise] = snap_helpers.classification.buildFeatureMatrix(...
            noiseData, state.selectedFeatures, state.featureInfo, state.customExpressions);
        
        X = [Xreal(validReal, :); Xnoise(validNoise, :)];
        labels = [ones(sum(validReal), 1); zeros(sum(validNoise), 1)];
        
        if size(X, 1) < 10
            uialert(fig, 'Too few valid samples', 'Insufficient Data');
            return;
        end
        
        % Build options from UI controls
        options = struct();
        
        % Kernel type
        kernelMap = containers.Map({'RBF (Gaussian)', 'Linear', 'Polynomial'}, {'rbf', 'linear', 'polynomial'});
        options.kernelFunction = kernelMap(handles.kernelDropdown.Value);
        
        % Box constraint (C parameter)
        options.boxConstraint = handles.boxConstraintEdit.Value;
        
        % Kernel scale (gamma)
        if handles.kernelScaleAutoCheck.Value
            options.kernelScale = 'auto';
        else
            options.kernelScale = handles.kernelScaleEdit.Value;
        end
        
        % Polynomial order
        options.polynomialOrder = handles.polyOrderSpinner.Value;
        
        % Z-score normalization (always recommended)
        options.standardize = handles.standardizeCheck.Value;
        
        % Cross-validation
        options.crossValidate = handles.cvCheck.Value;
        options.kFold = handles.kFoldSpinner.Value;
        
        options.verbose = true;
        
        handles.trainStatusLabel.Text = 'Training...';
        handles.trainStatusLabel.FontColor = [0.8 0.6 0.2];
        drawnow;
        
        % Train classifier and get normalization parameters
        [state.classifier, state.trainStats, state.normParams] = ...
            snap_helpers.classification.trainClassifier(X, labels, options);
        
        if state.trainStats.success
            handles.trainStatusLabel.Text = 'Trained!';
            handles.trainStatusLabel.FontColor = [0.2 0.7 0.3];
            
            if state.trainStats.crossValidated
                handles.accuracyLabel.Text = sprintf('%.1f%% (CV)', state.trainStats.cvAccuracy * 100);
            else
                handles.accuracyLabel.Text = sprintf('%.1f%% (train)', state.trainStats.trainAccuracy * 100);
            end
            handles.accuracyLabel.FontColor = [0.2 0.6 0.2];
            handles.exportBtn.Enable = 'on';
        else
            handles.trainStatusLabel.Text = 'Failed';
            handles.trainStatusLabel.FontColor = [0.8 0.2 0.2];
            handles.accuracyLabel.Text = '-';
            uialert(fig, sprintf('Training failed: %s', state.trainStats.error), 'Training Error');
        end
        
        fig.UserData.state = state;
        updateDisplay();
    end
    
    % ========================================================================
    % SAVE/LOAD
    % ========================================================================
    
    function saveProgress()
        [file, path] = uiputfile('*.mat', 'Save Progress', 'classification_progress.mat');
        if file == 0, return; end
        
        saveData = struct();
        saveData.labeledReal = state.labeledReal;
        saveData.labeledNoise = state.labeledNoise;
        saveData.skippedSpots = state.skippedSpots;
        saveData.selectedFeatures = state.selectedFeatures;
        saveData.customExpressions = state.customExpressions;
        saveData.classifier = state.classifier;
        saveData.trainStats = state.trainStats;
        saveData.normParams = state.normParams;  % Save normalization parameters
        saveData.fittingMethod = state.fittingMethod;
        saveData.sourceFile = state.filepath;
        saveData.timestamp = datetime('now');
        
        save(fullfile(path, file), '-struct', 'saveData');
        state.modified = false;
        fig.UserData.state = state;
        
        uialert(fig, 'Progress saved!', 'Saved', 'Icon', 'success');
    end
    
    function loadProgress()
        [file, path] = uigetfile('*.mat', 'Load Progress');
        if file == 0, return; end
        
        try
            data = load(fullfile(path, file));
            
            if isfield(data, 'labeledReal'), state.labeledReal = data.labeledReal; end
            if isfield(data, 'labeledNoise'), state.labeledNoise = data.labeledNoise; end
            if isfield(data, 'skippedSpots'), state.skippedSpots = data.skippedSpots; end
            if isfield(data, 'selectedFeatures')
                state.selectedFeatures = data.selectedFeatures;
            end
            if isfield(data, 'customExpressions')
                state.customExpressions = data.customExpressions;
            end
            
            % Update feature display
            nBase = numel(state.selectedFeatures);
            nCustom = numel(state.customExpressions);
            handles.featureCountLabel.Text = sprintf('%d base + %d custom', nBase, nCustom);
            featureList = state.selectedFeatures;
            for ii = 1:nCustom
                featureList{end+1} = sprintf('[EXPR] %s', state.customExpressions(ii).name);
            end
            handles.featureListArea.Value = featureList;
            
            if isfield(data, 'classifier')
                state.classifier = data.classifier;
                if ~isempty(data.classifier)
                    handles.trainStatusLabel.Text = 'Loaded';
                    handles.trainStatusLabel.FontColor = [0.2 0.7 0.3];
                    handles.exportBtn.Enable = 'on';
                end
            end
            if isfield(data, 'trainStats')
                state.trainStats = data.trainStats;
                if isfield(data.trainStats, 'cvAccuracy')
                    handles.accuracyLabel.Text = sprintf('%.1f%%', data.trainStats.cvAccuracy * 100);
                end
            end
            if isfield(data, 'normParams')
                state.normParams = data.normParams;
            end
            
            fig.UserData.state = state;
            updateDisplay();
            updateStats();
            checkTrainingEligibility();
            
            uialert(fig, 'Progress loaded!', 'Loaded', 'Icon', 'success');
        catch ME
            uialert(fig, sprintf('Failed to load: %s', ME.message), 'Load Error');
        end
    end
    
    function exportClassifier()
        if isempty(state.classifier)
            uialert(fig, 'No classifier to export', 'Export Error');
            return;
        end
        
        [file, path] = uiputfile('*.mat', 'Export Classifier for SNAP', 'snap_classifier.mat');
        if file == 0, return; end
        
        metadata = struct();
        metadata.nRealLabeled = numel(state.labeledReal);
        metadata.nNoiseLabeled = numel(state.labeledNoise);
        metadata.sourceFile = state.filepath;
        
        % Save with custom expressions and normalization parameters
        success = snap_helpers.classification.saveClassifier(...
            fullfile(path, file), state.classifier, state.selectedFeatures, ...
            state.featureInfo, state.trainStats, state.fittingMethod, metadata, ...
            state.customExpressions, state.normParams);
        
        if success
            accuracy = 0;
            if isfield(state.trainStats, 'cvAccuracy')
                accuracy = state.trainStats.cvAccuracy * 100;
            elseif isfield(state.trainStats, 'trainAccuracy')
                accuracy = state.trainStats.trainAccuracy * 100;
            end
            
            % Build feature summary
            nBase = numel(state.selectedFeatures);
            nCustom = numel(state.customExpressions);
            if nCustom > 0
                featureSummary = sprintf('%d base + %d custom features', nBase, nCustom);
            else
                featureSummary = strjoin(state.selectedFeatures(1:min(3,end)), ', ');
                if nBase > 3
                    featureSummary = [featureSummary, '...'];
                end
            end
            
            uialert(fig, sprintf('Classifier exported!\n\nFeatures: %s\nAccuracy: %.1f%%\n\nUse this file in SNAP or SNAP_batch.', ...
                featureSummary, accuracy), ...
                'Export Complete', 'Icon', 'success');
        end
    end
    
    % ========================================================================
    % KEYBOARD
    % ========================================================================
    
    function handleKeyPress(~, event)
        switch lower(event.Key)
            case 'r'
                labelSpot(1);
            case 'n'
                labelSpot(0);
            case 's'
                skipSpot();
            case 'space'
                toggleLivePredictions();
            case 'z'
                if any(strcmp(event.Modifier, 'control')) || any(strcmp(event.Modifier, 'command'))
                    undoLabel();
                end
            case 'leftarrow'
                navigateSpot(-1);
            case 'rightarrow'
                navigateSpot(1);
            case 'uparrow'
                navigateSpot(-10);
            case 'downarrow'
                navigateSpot(10);
        end
    end
    
    % ========================================================================
    % UI HELPER CALLBACKS FOR SVM PARAMETERS
    % ========================================================================
    
    function toggleKernelScaleAuto(src)
        if src.Value
            handles.kernelScaleEdit.Enable = 'off';
        else
            handles.kernelScaleEdit.Enable = 'on';
        end
    end
    
    function updateKernelOptions(src)
        % Enable/disable polynomial order based on kernel type
        if strcmp(src.Value, 'Polynomial')
            handles.polyOrderSpinner.Enable = 'on';
        else
            handles.polyOrderSpinner.Enable = 'off';
        end
    end
    
    function onClose(~, ~)
        if state.modified
            answer = questdlg('Save progress before closing?', 'Unsaved Changes', ...
                'Save', 'Discard', 'Cancel', 'Save');
            switch answer
                case 'Save'
                    saveProgress();
                case 'Cancel'
                    return;
            end
        end
        delete(fig);
    end
end


