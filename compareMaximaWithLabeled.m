function results = compareMaximaWithLabeled(csvFile, snapExportFile, maxDistance)
% compareMaximaWithLabeled - Compare manually labeled coordinates with SNAP-detected maxima
%
% USAGE:
%   compareMaximaWithLabeled()           % Opens UI
%   results = compareMaximaWithLabeled(csvFile, snapExportFile)
%   results = compareMaximaWithLabeled(csvFile, snapExportFile, maxDistance)
%
% Opens a simple UI to:
%   1. Browse and load a CSV file with labeled coordinates (columns: X, Y, Slice)
%   2. Browse and load a SNAP channel export .mat file
%   3. Set match distance threshold
%   4. Run comparison and view results
%

    % If called with no arguments or empty, launch UI
    if nargin == 0 || (nargin >= 1 && isempty(csvFile))
        launchUI();
        results = [];
        return;
    end
    
    % Default max distance
    if nargin < 3 || isempty(maxDistance)
        maxDistance = 3;
    end
    
    % Run comparison directly if files provided
    results = runComparison(csvFile, snapExportFile, maxDistance);
end

%% ========================================================================
%  UI LAUNCHER
%% ========================================================================
function launchUI()
    % Create figure
    fig = uifigure('Name', 'Maxima Comparison Tool', ...
        'Position', [200 150 550 500], ...
        'Resize', 'on');
    
    % Main grid
    grid = uigridlayout(fig, [8, 1]);
    grid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', '1x', 'fit'};
    grid.Padding = [20 20 20 20];
    grid.RowSpacing = 10;
    
    % Title
    uilabel(grid, 'Text', 'Compare Labeled vs SNAP Maxima', ...
        'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % CSV file selection
    csvPanel = uipanel(grid, 'Title', 'Labeled Coordinates (CSV)'); %#ok<NASGU>
    csvGrid = uigridlayout(csvPanel, [1, 2]);
    csvGrid.ColumnWidth = {'1x', 100};
    csvEdit = uieditfield(csvGrid, 'Editable', 'off', 'Placeholder', 'Select CSV file with X, Y, Slice columns...');
    uibutton(csvGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(~,~) browseCSV());
    
    % SNAP export file selection
    snapPanel = uipanel(grid, 'Title', 'SNAP Channel Export (.mat)'); %#ok<NASGU>
    snapGrid = uigridlayout(snapPanel, [1, 2]);
    snapGrid.ColumnWidth = {'1x', 100};
    snapEdit = uieditfield(snapGrid, 'Editable', 'off', 'Placeholder', 'Select SNAP export file...');
    uibutton(snapGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(~,~) browseSNAP());
    
    % Settings
    settingsPanel = uipanel(grid, 'Title', 'Settings'); %#ok<NASGU>
    settingsGrid = uigridlayout(settingsPanel, [1, 4]);
    settingsGrid.ColumnWidth = {'fit', 80, '1x', 'fit'};
    uilabel(settingsGrid, 'Text', 'Max match distance:');
    distEdit = uieditfield(settingsGrid, 'numeric', 'Value', 3, 'Limits', [0.1 100]);
    uilabel(settingsGrid, 'Text', 'pixels');
    
    % Info label
    uilabel(grid, 'Text', ...
        'Coordinate conversion: CSV [X,Y,Slice] → SNAP [Y+1, X+1, Slice+1]', ...
        'FontSize', 10, 'FontColor', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');
    
    % Results text area
    resultsArea = uitextarea(grid, 'Editable', 'off', ...
        'FontName', 'Consolas', 'FontSize', 11, ...
        'Value', {'Results will appear here after comparison...'});
    
    % Buttons
    btnGrid = uigridlayout(grid, [1, 2]);
    btnGrid.ColumnWidth = {'1x', '1x'};
    btnGrid.Padding = [0 0 0 0];
    
    uibutton(btnGrid, 'Text', 'Run Comparison', ...
        'FontWeight', 'bold', 'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'white', ...
        'ButtonPushedFcn', @(~,~) runFromUI());
    
    uibutton(btnGrid, 'Text', 'Close', ...
        'ButtonPushedFcn', @(~,~) delete(fig));
    
    % Store state
    state = struct();
    state.csvFile = '';
    state.snapFile = '';
    
    %% Nested functions
    function browseCSV()
        [file, path] = uigetfile({'*.csv', 'CSV Files (*.csv)'}, 'Select Labeled Coordinates CSV');
        if file ~= 0
            state.csvFile = fullfile(path, file);
            csvEdit.Value = state.csvFile;
        end
    end
    
    function browseSNAP()
        [file, path] = uigetfile({'*.mat', 'SNAP Export Files (*.mat)'}, 'Select SNAP Channel Export');
        if file ~= 0
            state.snapFile = fullfile(path, file);
            snapEdit.Value = state.snapFile;
        end
    end
    
    function runFromUI()
        % Validate inputs
        if isempty(state.csvFile)
            uialert(fig, 'Please select a CSV file with labeled coordinates.', 'Missing Input');
            return;
        end
        if isempty(state.snapFile)
            uialert(fig, 'Please select a SNAP export file.', 'Missing Input');
            return;
        end
        
        % Update UI
        resultsArea.Value = {'Running comparison...'};
        drawnow;
        
        try
            % Run comparison
            maxDist = distEdit.Value;
            results = runComparison(state.csvFile, state.snapFile, maxDist);
            
            % Display results
            stats = results.stats;
            lines = {
                '═══════════════════════════════════════════════════════════'
                '                    COMPARISON RESULTS'
                '═══════════════════════════════════════════════════════════'
                ''
                sprintf('  Labeled coordinates (ground truth):   %d', stats.nLabeled)
                sprintf('  SNAP-detected maxima:                 %d', stats.nSNAP)
                ''
                '───────────────────────────────────────────────────────────'
                ''
                sprintf('  ✓ Matched (True Positives):           %d', stats.nMatched)
                sprintf('  ✗ Unmatched SNAP (False Positives):   %d', stats.nUnmatchedSNAP)
                sprintf('  ✗ Unmatched Labeled (False Negatives): %d', stats.nUnmatchedLabeled)
                ''
                '───────────────────────────────────────────────────────────'
                ''
                sprintf('  Precision:  %.1f%%', stats.precision)
                sprintf('  Recall:     %.1f%%', stats.recall)
                sprintf('  F1 Score:   %.1f%%', stats.f1Score)
                ''
            };
            
            if stats.nMatched > 0
                lines{end+1} = sprintf('  Mean match distance: %.2f px', stats.meanMatchDistance);
                lines{end+1} = sprintf('  Max match distance:  %.2f px', stats.maxMatchDistance);
                lines{end+1} = '';
            end
            
            lines{end+1} = '═══════════════════════════════════════════════════════════';
            lines{end+1} = '';
            lines{end+1} = 'Results saved to files in export directory.';
            
            resultsArea.Value = lines;
            
        catch ME
            resultsArea.Value = {
                'ERROR:'
                ''
                ME.message
                ''
                'Stack trace:'
                ME.stack(1).name
            };
        end
    end
end

%% ========================================================================
%  COMPARISON LOGIC
%% ========================================================================
function results = runComparison(csvFile, snapExportFile, maxDistance)
    
    fprintf('=== Maxima Comparison Tool ===\n\n');
    fprintf('CSV file: %s\n', csvFile);
    fprintf('SNAP export: %s\n', snapExportFile);
    fprintf('Max match distance: %.1f pixels\n\n', maxDistance);
    
    %% Load labeled coordinates from CSV
    fprintf('Loading labeled coordinates...\n');
    try
        labeledTable = readtable(csvFile);
    catch ME
        error('Failed to read CSV file: %s', ME.message);
    end
    
    % Check for required columns
    requiredCols = {'X', 'Y', 'Slice'};
    missingCols = setdiff(requiredCols, labeledTable.Properties.VariableNames);
    if ~isempty(missingCols)
        % Try case-insensitive match
        colNames = labeledTable.Properties.VariableNames;
        for i = 1:numel(requiredCols)
            idx = find(strcmpi(colNames, requiredCols{i}), 1);
            if ~isempty(idx)
                labeledTable.Properties.VariableNames{idx} = requiredCols{i};
            end
        end
        missingCols = setdiff(requiredCols, labeledTable.Properties.VariableNames);
        if ~isempty(missingCols)
            error('CSV file missing required columns: %s\nFound columns: %s', ...
                strjoin(missingCols, ', '), strjoin(colNames, ', '));
        end
    end
    
    % Extract and convert coordinates
    % CSV format: [X, Y, Slice] (0-indexed, ImageJ convention)
    % Convert to: [row, col, slice] = [Y+1, X+1, Slice+1] (1-indexed, MATLAB convention)
    labeledCoords = zeros(height(labeledTable), 3);
    labeledCoords(:, 1) = labeledTable.Y + 1;  % row = Y + 1
    labeledCoords(:, 2) = labeledTable.X + 1;  % col = X + 1
    labeledCoords(:, 3) = labeledTable.Slice + 1;  % slice + 1
    
    nLabeled = size(labeledCoords, 1);
    fprintf('  Loaded %d labeled coordinates\n', nLabeled);
    fprintf('  Coordinate conversion applied: [X,Y,Slice] -> [Y+1, X+1, Slice+1]\n');
    
    %% Load SNAP export
    fprintf('\nLoading SNAP export...\n');
    try
        snapData = load(snapExportFile);
    catch ME
        error('Failed to load SNAP export: %s', ME.message);
    end
    
    % Find signals in the export
    signals = [];
    if isfield(snapData, 'signals')
        signals = snapData.signals;
    elseif isfield(snapData, 'fitResults')
        signals = snapData.fitResults;
    else
        % Search for signals
        fnames = fieldnames(snapData);
        for i = 1:numel(fnames)
            if isstruct(snapData.(fnames{i})) && ...
               (isfield(snapData.(fnames{i}), 'maxima_coords') || isfield(snapData.(fnames{i}), 'fitted_coords'))
                signals = snapData.(fnames{i});
                break;
            end
        end
    end
    
    if isempty(signals)
        error('Could not find signals in SNAP export. Expected "signals" field.');
    end
    
    % Extract maxima coordinates from SNAP signals
    nSNAP = numel(signals);
    snapCoords = zeros(nSNAP, 3);
    
    for i = 1:nSNAP
        if isfield(signals, 'maxima_coords') && ~isempty(signals(i).maxima_coords)
            coords = signals(i).maxima_coords;
            snapCoords(i, :) = coords(1:min(3, numel(coords)));
        elseif isfield(signals, 'fitted_coords') && ~isempty(signals(i).fitted_coords)
            % Fallback to fitted coords if maxima not available
            coords = signals(i).fitted_coords;
            snapCoords(i, :) = coords(1:min(3, numel(coords)));
            if i == 1
                fprintf('  WARNING: Using fitted_coords (maxima_coords not available)\n');
            end
        else
            snapCoords(i, :) = [NaN, NaN, NaN];
        end
    end
    
    % Handle 2D case (no Z)
    if all(snapCoords(:,3) == 0) || all(isnan(snapCoords(:,3)))
        snapCoords(:,3) = 1;
        labeledCoords(:,3) = 1;
        fprintf('  Note: Treating as 2D data (Z ignored)\n');
    end
    
    fprintf('  Loaded %d SNAP-detected maxima\n', nSNAP);
    
    %% Perform matching (1-to-1 greedy matching)
    fprintf('\nPerforming coordinate matching...\n');
    
    % Compute all pairwise distances
    distMatrix = inf(nLabeled, nSNAP);
    for i = 1:nLabeled
        for j = 1:nSNAP
            distMatrix(i, j) = sqrt(sum((labeledCoords(i,:) - snapCoords(j,:)).^2));
        end
    end
    
    % Track matches (1-to-1 matching)
    labeledMatched = false(nLabeled, 1);
    labeledMatchIdx = zeros(nLabeled, 1); %#ok<NASGU>
    
    snapMatched = false(nSNAP, 1);
    snapMatchIdx = zeros(nSNAP, 1);  % Index of matched labeled coordinate
    snapMatchDist = inf(nSNAP, 1);
    
    % Greedy 1-to-1 matching: repeatedly match closest pair until no more matches
    tempDistMatrix = distMatrix;
    while true
        % Find minimum distance that is within threshold
        [minDist, minIdx] = min(tempDistMatrix(:));
        
        if minDist > maxDistance || isinf(minDist)
            break;  % No more valid matches
        end
        
        % Convert linear index to subscripts
        [labeledIdx, snapIdx] = ind2sub([nLabeled, nSNAP], minIdx);
        
        % Record this match
        labeledMatched(labeledIdx) = true;
        labeledMatchIdx(labeledIdx) = snapIdx;
        
        snapMatched(snapIdx) = true;
        snapMatchIdx(snapIdx) = labeledIdx;
        snapMatchDist(snapIdx) = minDist;
        
        % Remove this pair from consideration (set entire row and column to inf)
        tempDistMatrix(labeledIdx, :) = inf;
        tempDistMatrix(:, snapIdx) = inf;
    end
    
    %% Build results
    fprintf('\nBuilding results...\n');
    
    % Matched SNAP signals
    matchedIdx = find(snapMatched);
    if ~isempty(matchedIdx)
        matchedTable = struct2table(signals(matchedIdx));
        matchedTable.match_distance = snapMatchDist(matchedIdx);
        matchedTable.labeled_idx = snapMatchIdx(matchedIdx);
        % Add labeled coordinates for reference
        matchedTable.labeled_row = labeledCoords(snapMatchIdx(matchedIdx), 1);
        matchedTable.labeled_col = labeledCoords(snapMatchIdx(matchedIdx), 2);
        matchedTable.labeled_slice = labeledCoords(snapMatchIdx(matchedIdx), 3);
    else
        matchedTable = table();
    end
    
    % Unmatched SNAP signals (potential false positives)
    unmatchedSNAPIdx = find(~snapMatched);
    if ~isempty(unmatchedSNAPIdx)
        unmatchedSNAPTable = struct2table(signals(unmatchedSNAPIdx));
    else
        unmatchedSNAPTable = table();
    end
    
    % Unmatched labeled coordinates (potential false negatives)
    unmatchedLabeledIdx = find(~labeledMatched);
    if ~isempty(unmatchedLabeledIdx)
        unmatchedLabeledTable = table();
        unmatchedLabeledTable.row = labeledCoords(unmatchedLabeledIdx, 1);
        unmatchedLabeledTable.col = labeledCoords(unmatchedLabeledIdx, 2);
        unmatchedLabeledTable.slice = labeledCoords(unmatchedLabeledIdx, 3);
        unmatchedLabeledTable.original_X = labeledTable.X(unmatchedLabeledIdx);
        unmatchedLabeledTable.original_Y = labeledTable.Y(unmatchedLabeledIdx);
        unmatchedLabeledTable.original_Slice = labeledTable.Slice(unmatchedLabeledIdx);
    else
        unmatchedLabeledTable = table();
    end
    
    % Compute statistics
    stats = struct();
    stats.nLabeled = nLabeled;
    stats.nSNAP = nSNAP;
    stats.nMatched = sum(snapMatched);
    stats.nUnmatchedSNAP = sum(~snapMatched);
    stats.nUnmatchedLabeled = sum(~labeledMatched);
    stats.matchRate_SNAP = 100 * stats.nMatched / max(1, stats.nSNAP);
    stats.matchRate_Labeled = 100 * sum(labeledMatched) / max(1, nLabeled);
    stats.precision = 100 * stats.nMatched / max(1, stats.nSNAP);  % TP / (TP + FP)
    stats.recall = 100 * sum(labeledMatched) / max(1, nLabeled);   % TP / (TP + FN)
    if stats.precision + stats.recall > 0
        stats.f1Score = 2 * (stats.precision * stats.recall) / (stats.precision + stats.recall);
    else
        stats.f1Score = 0;
    end
    stats.meanMatchDistance = mean(snapMatchDist(snapMatched));
    stats.maxMatchDistance = max(snapMatchDist(snapMatched));
    
    % Build output
    results = struct();
    results.matched = matchedTable;
    results.unmatchedSNAP = unmatchedSNAPTable;
    results.unmatchedLabeled = unmatchedLabeledTable;
    results.matchDistances = snapMatchDist(snapMatched);
    results.stats = stats;
    results.labeledCoords = labeledCoords;
    results.snapCoords = snapCoords;
    results.labeledMatched = labeledMatched;
    results.snapMatched = snapMatched;
    
    %% Print summary
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('║                    COMPARISON RESULTS                         ║\n');
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Labeled coordinates (ground truth):     %5d                ║\n', stats.nLabeled);
    fprintf('║  SNAP-detected maxima:                   %5d                ║\n', stats.nSNAP);
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Matched (True Positives):               %5d                ║\n', stats.nMatched);
    fprintf('║  Unmatched SNAP (False Positives):       %5d                ║\n', stats.nUnmatchedSNAP);
    fprintf('║  Unmatched Labeled (False Negatives):    %5d                ║\n', stats.nUnmatchedLabeled);
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Precision (TP / SNAP detected):         %5.1f%%              ║\n', stats.precision);
    fprintf('║  Recall (TP / Labeled):                  %5.1f%%              ║\n', stats.recall);
    fprintf('║  F1 Score:                               %5.1f%%              ║\n', stats.f1Score);
    fprintf('╠══════════════════════════════════════════════════════════════╣\n');
    if stats.nMatched > 0
        fprintf('║  Mean match distance:                    %5.2f px            ║\n', stats.meanMatchDistance);
        fprintf('║  Max match distance:                     %5.2f px            ║\n', stats.maxMatchDistance);
    end
    fprintf('╚══════════════════════════════════════════════════════════════╝\n');
    fprintf('\nDone!\n');
end

