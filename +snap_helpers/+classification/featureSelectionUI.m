function [selectedFeatures, customExpressions, cancelled] = featureSelectionUI(fittingMethod, has3D, hasPhysicalSpacing, previousSelection, previousExpressions)
% featureSelectionUI - Interactive GUI for selecting classification features
%
% ============================================================================
% USER-FRIENDLY FEATURE SELECTION WITH CUSTOM EXPRESSION SUPPORT
% ============================================================================
%
% This function presents a GUI that allows users to:
%   1. Select which base features to use for SVM classification
%   2. Create CUSTOM EXPRESSIONS combining multiple features
%
% CUSTOM EXPRESSIONS:
%   Users can define mathematical expressions like:
%     'amplitude / background'
%     'log(integrated_intensity)'
%     'sigma_x / sigma_y'
%     'sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2)'
%
% USAGE:
%   [selected, expressions, cancelled] = snap_helpers.classification.featureSelectionUI(method, is3D)
%
% INPUTS:
%   fittingMethod       - Current fitting method string
%   has3D               - Logical: whether data is 3D
%   hasPhysicalSpacing  - Logical: whether physical spacing is available
%   previousSelection   - Cell array of previously selected features
%   previousExpressions - Struct array of previously defined expressions
%
% OUTPUTS:
%   selectedFeatures  - Cell array of selected feature names
%   customExpressions - Struct array with fields: .name, .expression
%   cancelled         - Logical: true if user cancelled
%
% ============================================================================

    % Handle defaults
    if nargin < 3, hasPhysicalSpacing = false; end
    if nargin < 4, previousSelection = {}; end
    if nargin < 5, previousExpressions = struct('name', {}, 'expression', {}); end
    
    % Get available features
    [features, featureInfo] = snap_helpers.classification.getAvailableFeatures(fittingMethod, has3D, hasPhysicalSpacing);
    
    % Initialize outputs
    cancelled = true;
    selectedFeatures = {};
    customExpressions = previousExpressions;
    
    % Category display order and labels
    categories = {'intensity', 'shape', 'quality', 'derived', 'position', 'physical'};
    categoryLabels = struct();
    categoryLabels.intensity = 'Intensity Features';
    categoryLabels.shape = 'Shape Features';
    categoryLabels.quality = 'Fit Quality';
    categoryLabels.derived = 'Derived Features';
    categoryLabels.position = 'Position (usually exclude)';
    categoryLabels.physical = 'Physical Units';
    
    categoryColors = struct();
    categoryColors.intensity = [0.2 0.6 0.3];
    categoryColors.shape = [0.3 0.4 0.7];
    categoryColors.quality = [0.7 0.5 0.2];
    categoryColors.derived = [0.5 0.3 0.6];
    categoryColors.position = [0.5 0.5 0.5];
    categoryColors.physical = [0.2 0.5 0.5];
    
    % Create figure
    fig = uifigure('Name', 'Select Classification Features', ...
        'Position', [100 50 900 800], ...
        'CloseRequestFcn', @onCancel);
    
    % Main layout - two columns
    mainGrid = uigridlayout(fig, [1 2]);
    mainGrid.ColumnWidth = {'1x', '1x'};
    mainGrid.Padding = [10 10 10 10];
    mainGrid.ColumnSpacing = 15;
    
    % === LEFT COLUMN: Base Features ===
    leftPanel = uipanel(mainGrid, 'Title', 'Base Features', 'FontWeight', 'bold');
    leftGrid = uigridlayout(leftPanel, [4 1]);
    leftGrid.RowHeight = {'fit', '1x', 'fit', 'fit'};
    leftGrid.Padding = [10 10 10 10];
    leftGrid.RowSpacing = 10;
    
    % Header info
    uilabel(leftGrid, 'Text', sprintf('Method: %s | 3D: %s', fittingMethod, string(has3D)), ...
        'FontSize', 10, 'FontColor', [0.4 0.4 0.4]);
    
    % Scrollable feature panel
    scrollPanel = uipanel(leftGrid, 'Scrollable', 'on', 'BorderType', 'line');
    
    % Group features by category
    featuresByCategory = struct();
    for c = 1:numel(categories)
        featuresByCategory.(categories{c}) = {};
    end
    
    for f = 1:numel(features)
        fname = features{f};
        if isfield(featureInfo, fname) && isfield(featureInfo.(fname), 'category')
            cat = featureInfo.(fname).category;
            if isfield(featuresByCategory, cat)
                featuresByCategory.(cat){end+1} = fname;
            end
        end
    end
    
    % Count non-empty categories
    nonEmptyCategories = {};
    for c = 1:numel(categories)
        if ~isempty(featuresByCategory.(categories{c}))
            nonEmptyCategories{end+1} = categories{c};
        end
    end
    
    scrollGrid = uigridlayout(scrollPanel, [numel(nonEmptyCategories) 1]);
    scrollGrid.RowHeight = repmat({'fit'}, 1, numel(nonEmptyCategories));
    scrollGrid.Padding = [5 5 5 5];
    scrollGrid.RowSpacing = 8;
    
    % Store checkboxes for later access
    checkboxes = struct();
    
    % Create category panels
    for c = 1:numel(nonEmptyCategories)
        cat = nonEmptyCategories{c};
        catFeatures = featuresByCategory.(cat);
        
        if isempty(catFeatures), continue; end
        
        catColor = categoryColors.(cat);
        catPanel = uipanel(scrollGrid, 'Title', categoryLabels.(cat), ...
            'FontWeight', 'bold', 'ForegroundColor', catColor);
        
        catGrid = uigridlayout(catPanel, [numel(catFeatures) 2]);
        catGrid.ColumnWidth = {180, '1x'};
        catGrid.RowHeight = repmat({22}, 1, numel(catFeatures));
        catGrid.Padding = [5 5 5 5];
        
        for i = 1:numel(catFeatures)
            fname = catFeatures{i};
            
            desc = '';
            if isfield(featureInfo, fname) && isfield(featureInfo.(fname), 'description')
                desc = featureInfo.(fname).description;
            end
            
            isSelected = false;
            if ~isempty(previousSelection)
                isSelected = ismember(fname, previousSelection);
            else
                isSelected = ~strcmp(cat, 'position');
            end
            
            cb = uicheckbox(catGrid, 'Text', fname, 'Value', isSelected, 'FontName', 'Consolas');
            checkboxes.(matlab.lang.makeValidName(fname)) = cb;
            
            uilabel(catGrid, 'Text', desc, 'FontSize', 9, 'FontColor', [0.5 0.5 0.5]);
        end
    end
    
    % Statistics panel
    statsGrid = uigridlayout(leftGrid, [1 2]);
    statsGrid.ColumnWidth = {'1x', 'fit'};
    statsGrid.Padding = [0 0 0 0];
    
    statsLabel = uilabel(statsGrid, 'Text', '', 'FontSize', 12);
    updateStats();
    
    % Add callback to update stats
    fnames = fieldnames(checkboxes);
    for i = 1:numel(fnames)
        checkboxes.(fnames{i}).ValueChangedFcn = @(~,~) updateStats();
    end
    
    % Quick select buttons
    quickBtnGrid = uigridlayout(statsGrid, [1 3]);
    quickBtnGrid.ColumnWidth = {'fit', 'fit', 'fit'};
    quickBtnGrid.Padding = [0 0 0 0];
    
    uibutton(quickBtnGrid, 'Text', 'All', 'ButtonPushedFcn', @(~,~) selectAll(true));
    uibutton(quickBtnGrid, 'Text', 'None', 'ButtonPushedFcn', @(~,~) selectAll(false));
    uibutton(quickBtnGrid, 'Text', 'Recommended', 'ButtonPushedFcn', @(~,~) selectRecommended());
    
    % OK/Cancel for left column
    leftBtnGrid = uigridlayout(leftGrid, [1 2]);
    leftBtnGrid.ColumnWidth = {'1x', '1x'};
    leftBtnGrid.Padding = [0 0 0 0];
    
    uibutton(leftBtnGrid, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) onCancel());
    uibutton(leftBtnGrid, 'Text', 'OK', 'ButtonPushedFcn', @(~,~) onOK(), 'FontWeight', 'bold');
    
    % === RIGHT COLUMN: Custom Expressions ===
    rightPanel = uipanel(mainGrid, 'Title', 'Custom Expressions (Optional)', 'FontWeight', 'bold');
    rightGrid = uigridlayout(rightPanel, [6 1]);
    rightGrid.RowHeight = {'fit', 'fit', '1x', 'fit', 'fit', 'fit'};
    rightGrid.Padding = [10 10 10 10];
    rightGrid.RowSpacing = 10;
    
    % Instructions
    instrText = uitextarea(rightGrid, 'Editable', 'off', 'FontSize', 10, ...
        'Value', {
            'Create custom features using mathematical expressions.'
            ''
            'Available operations:'
            '  +, -, *, /, ^, (), log, log10, sqrt, abs, exp'
            ''
            'Examples:'
            '  amplitude / background'
            '  log(integrated_intensity / background)'
            '  sqrt(sigma_x^2 + sigma_y^2)'
            '  (amplitude - background) / sqrt(background)'
        });
    instrText.BackgroundColor = [0.95 0.95 0.98];
    
    % Available variables reference
    varRefPanel = uipanel(rightGrid, 'Title', 'Available Variables');
    varRefGrid = uigridlayout(varRefPanel, [1 1]);
    varRefGrid.Padding = [5 5 5 5];
    
    varList = uitextarea(varRefGrid, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 9);
    varList.Value = features;
    
    % Expression list
    exprPanel = uipanel(rightGrid, 'Title', 'Defined Expressions');
    exprGrid = uigridlayout(exprPanel, [1 1]);
    exprGrid.Padding = [5 5 5 5];
    
    exprTable = uitable(exprGrid, 'ColumnName', {'Name', 'Expression'}, ...
        'ColumnWidth', {100, 'auto'}, 'RowName', {}, ...
        'CellSelectionCallback', @onExprSelect);
    updateExprTable();
    
    % Expression input
    inputPanel = uipanel(rightGrid, 'Title', 'Add/Edit Expression');
    inputGrid = uigridlayout(inputPanel, [3 2]);
    inputGrid.RowHeight = {'fit', 'fit', 'fit'};
    inputGrid.ColumnWidth = {80, '1x'};
    inputGrid.Padding = [10 10 10 10];
    inputGrid.RowSpacing = 8;
    
    uilabel(inputGrid, 'Text', 'Name:');
    exprNameEdit = uieditfield(inputGrid, 'Value', '', ...
        'Placeholder', 'e.g., brightness_ratio', 'FontName', 'Consolas');
    
    uilabel(inputGrid, 'Text', 'Expression:');
    exprEdit = uieditfield(inputGrid, 'Value', '', ...
        'Placeholder', 'e.g., amplitude / background', 'FontName', 'Consolas');
    
    uilabel(inputGrid, 'Text', '');
    validateLabel = uilabel(inputGrid, 'Text', '', 'FontSize', 10);
    
    % Expression buttons
    exprBtnGrid = uigridlayout(rightGrid, [1 4]);
    exprBtnGrid.ColumnWidth = {'fit', 'fit', 'fit', '1x'};
    exprBtnGrid.Padding = [0 0 0 0];
    
    uibutton(exprBtnGrid, 'Text', 'Add', 'ButtonPushedFcn', @(~,~) addExpression());
    uibutton(exprBtnGrid, 'Text', 'Update', 'ButtonPushedFcn', @(~,~) updateExpression());
    uibutton(exprBtnGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(~,~) removeExpression());
    uibutton(exprBtnGrid, 'Text', 'Validate', 'ButtonPushedFcn', @(~,~) validateExpression());
    
    % Preset expressions
    presetPanel = uipanel(rightGrid, 'Title', 'Quick Presets');
    presetGrid = uigridlayout(presetPanel, [2 3]);
    presetGrid.RowHeight = {'fit', 'fit'};
    presetGrid.ColumnWidth = {'1x', '1x', '1x'};
    presetGrid.Padding = [5 5 5 5];
    
    uibutton(presetGrid, 'Text', 'SNR', 'ButtonPushedFcn', @(~,~) addPreset('snr', 'integrated_intensity / background'));
    uibutton(presetGrid, 'Text', 'Amp/BG', 'ButtonPushedFcn', @(~,~) addPreset('amp_over_bg', 'amplitude / background'));
    uibutton(presetGrid, 'Text', 'Log SNR', 'ButtonPushedFcn', @(~,~) addPreset('log_snr', 'log(integrated_intensity / background + 1)'));
    uibutton(presetGrid, 'Text', 'Ellipticity', 'ButtonPushedFcn', @(~,~) addPreset('ellipticity', 'sigma_x / sigma_y'));
    uibutton(presetGrid, 'Text', 'Spot Size', 'ButtonPushedFcn', @(~,~) addPreset('spot_size', 'sqrt(sigma_x^2 + sigma_y^2)'));
    uibutton(presetGrid, 'Text', 'Norm Amp', 'ButtonPushedFcn', @(~,~) addPreset('norm_amp', '(amplitude - background) / sqrt(background)'));
    
    selectedExprIdx = [];
    
    % Wait for user
    uiwait(fig);
    
    % === Nested functions ===
    
    function updateStats()
        count = 0;
        cbNames = fieldnames(checkboxes);
        for ii = 1:numel(cbNames)
            if checkboxes.(cbNames{ii}).Value
                count = count + 1;
            end
        end
        nExpr = numel(customExpressions);
        statsLabel.Text = sprintf('%d base + %d custom = %d features', count, nExpr, count + nExpr);
        
        if count + nExpr == 0
            statsLabel.FontColor = [0.8 0.2 0.2];
        elseif count + nExpr < 3
            statsLabel.FontColor = [0.8 0.6 0.2];
        else
            statsLabel.FontColor = [0.2 0.6 0.2];
        end
    end
    
    function selectAll(val)
        cbNames = fieldnames(checkboxes);
        for ii = 1:numel(cbNames)
            checkboxes.(cbNames{ii}).Value = val;
        end
        updateStats();
    end
    
    function selectRecommended()
        recommended = {'amplitude', 'integrated_intensity', 'background', ...
            'r_squared', 'radial_symmetry_score', 'snr', ...
            'sigma_x', 'sigma_y', 'sigma_z', 'sigma_xy_ratio', 'sigma_sum', ...
            'amplitude_over_background', 'amplitude_x', 'amplitude_y', 'amplitude_z', ...
            'amplitude_xy'};
        
        cbNames = fieldnames(checkboxes);
        for ii = 1:numel(cbNames)
            originalName = checkboxes.(cbNames{ii}).Text;
            checkboxes.(cbNames{ii}).Value = ismember(originalName, recommended);
        end
        updateStats();
    end
    
    function updateExprTable()
        if isempty(customExpressions)
            exprTable.Data = {};
        else
            data = cell(numel(customExpressions), 2);
            for ii = 1:numel(customExpressions)
                data{ii, 1} = customExpressions(ii).name;
                data{ii, 2} = customExpressions(ii).expression;
            end
            exprTable.Data = data;
        end
        updateStats();
    end
    
    function onExprSelect(~, event)
        if ~isempty(event.Indices)
            selectedExprIdx = event.Indices(1, 1);
            if selectedExprIdx <= numel(customExpressions)
                exprNameEdit.Value = customExpressions(selectedExprIdx).name;
                exprEdit.Value = customExpressions(selectedExprIdx).expression;
            end
        end
    end
    
    function addExpression()
        name = strtrim(exprNameEdit.Value);
        expr = strtrim(exprEdit.Value);
        
        if isempty(name) || isempty(expr)
            validateLabel.Text = '⚠ Name and expression required';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        
        % Check for duplicate name
        for ii = 1:numel(customExpressions)
            if strcmp(customExpressions(ii).name, name)
                validateLabel.Text = '⚠ Name already exists';
                validateLabel.FontColor = [0.8 0.2 0.2];
                return;
            end
        end
        
        % Validate expression
        if ~validateExpressionSyntax(expr)
            return;
        end
        
        % Add expression
        newExpr = struct('name', name, 'expression', expr);
        if isempty(customExpressions)
            customExpressions = newExpr;
        else
            customExpressions(end+1) = newExpr;
        end
        
        exprNameEdit.Value = '';
        exprEdit.Value = '';
        validateLabel.Text = '✓ Expression added';
        validateLabel.FontColor = [0.2 0.6 0.2];
        updateExprTable();
    end
    
    function updateExpression()
        if isempty(selectedExprIdx) || selectedExprIdx > numel(customExpressions)
            validateLabel.Text = '⚠ Select an expression first';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        
        name = strtrim(exprNameEdit.Value);
        expr = strtrim(exprEdit.Value);
        
        if isempty(name) || isempty(expr)
            validateLabel.Text = '⚠ Name and expression required';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        
        if ~validateExpressionSyntax(expr)
            return;
        end
        
        customExpressions(selectedExprIdx).name = name;
        customExpressions(selectedExprIdx).expression = expr;
        
        validateLabel.Text = '✓ Expression updated';
        validateLabel.FontColor = [0.2 0.6 0.2];
        updateExprTable();
    end
    
    function removeExpression()
        if isempty(selectedExprIdx) || selectedExprIdx > numel(customExpressions)
            validateLabel.Text = '⚠ Select an expression first';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        
        customExpressions(selectedExprIdx) = [];
        selectedExprIdx = [];
        exprNameEdit.Value = '';
        exprEdit.Value = '';
        validateLabel.Text = '✓ Expression removed';
        validateLabel.FontColor = [0.5 0.5 0.5];
        updateExprTable();
    end
    
    function validateExpression()
        expr = strtrim(exprEdit.Value);
        if isempty(expr)
            validateLabel.Text = '⚠ Enter an expression';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        validateExpressionSyntax(expr);
    end
    
    function valid = validateExpressionSyntax(expr)
        valid = false;
        
        % Check for dangerous patterns
        dangerousPatterns = {'eval', 'feval', 'system', 'delete', 'rmdir', 'save', 'load', ...
            'fopen', 'fwrite', 'fprintf', 'sprintf', 'assignin', 'evalin'};
        for ii = 1:numel(dangerousPatterns)
            if contains(lower(expr), dangerousPatterns{ii})
                validateLabel.Text = sprintf('⚠ Forbidden: %s', dangerousPatterns{ii});
                validateLabel.FontColor = [0.8 0.2 0.2];
                return;
            end
        end
        
        % Check parentheses balance
        if sum(expr == '(') ~= sum(expr == ')')
            validateLabel.Text = '⚠ Unbalanced parentheses';
            validateLabel.FontColor = [0.8 0.2 0.2];
            return;
        end
        
        % Try to parse with dummy values
        testWorkspace = struct();
        for ii = 1:numel(features)
            safeName = matlab.lang.makeValidName(features{ii});
            testWorkspace.(safeName) = rand(10, 1);
        end
        % Add common derived names
        testWorkspace.snr = rand(10, 1);
        testWorkspace.amp_bg = rand(10, 1);
        testWorkspace.sigma_xy_ratio = rand(10, 1);
        testWorkspace.sigma_xy_product = rand(10, 1);
        testWorkspace.sigma_xy_sum = rand(10, 1);
        
        testExpr = expr;
        allNames = fieldnames(testWorkspace);
        [~, sortIdx] = sort(cellfun(@length, allNames), 'descend');
        allNames = allNames(sortIdx);
        
        for ii = 1:numel(allNames)
            fname = allNames{ii};
            pattern = ['(?<![a-zA-Z0-9_])' fname '(?![a-zA-Z0-9_])'];
            testExpr = regexprep(testExpr, pattern, ['testWorkspace.' fname]);
        end
        
        try
            result = eval(testExpr);
            if ~isnumeric(result)
                validateLabel.Text = '⚠ Expression must return numeric';
                validateLabel.FontColor = [0.8 0.2 0.2];
                return;
            end
            validateLabel.Text = '✓ Valid expression';
            validateLabel.FontColor = [0.2 0.6 0.2];
            valid = true;
        catch ME
            validateLabel.Text = sprintf('⚠ %s', ME.message);
            validateLabel.FontColor = [0.8 0.2 0.2];
        end
    end
    
    function addPreset(name, expr)
        exprNameEdit.Value = name;
        exprEdit.Value = expr;
        validateExpression();
    end
    
    function onOK()
        cancelled = false;
        selectedFeatures = {};
        cbNames = fieldnames(checkboxes);
        for ii = 1:numel(cbNames)
            if checkboxes.(cbNames{ii}).Value
                selectedFeatures{end+1} = checkboxes.(cbNames{ii}).Text;
            end
        end
        delete(fig);
    end
    
    function onCancel()
        cancelled = true;
        selectedFeatures = {};
        customExpressions = struct('name', {}, 'expression', {});
        delete(fig);
    end
end
