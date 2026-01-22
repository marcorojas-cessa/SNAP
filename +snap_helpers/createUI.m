function handles = createUI(Nmax, lastUsed)
% Creates the user interface using a grid layout for responsiveness.
% Fixed: Corrected label creation syntax to avoid GridLayout indexing errors

fig = uifigure('Name', 'SNAP - Spot & Nuclei Analysis Pipeline', 'Position', [50 50 1900 1250]);
handles.fig = fig;
handles.Nmax = Nmax;
handles.lastUsed = lastUsed; % Store default parameters for reference

% Main grid layout - 3 columns: controls, previews, analysis
mainGrid = uigridlayout(fig, [1, 3]);
mainGrid.ColumnWidth = {420, '1x', 300};
mainGrid.Padding = [10 10 10 10];

% --- Function Handles for Robust UI Creation ---
labelFcn = @uilabel;
editFcn = @uieditfield;
buttonFcn = @uibutton;
dropdownFcn = @uidropdown;
checkboxFcn = @uicheckbox;
sliderFcn = @uislider;
axesFcn = @uiaxes;
panelFcn = @uipanel;

% --- Left Panel ---
leftPanel = panelFcn(mainGrid, 'Title', 'File Loading & Global Controls');
leftPanel.Layout.Row = 1;
leftPanel.Layout.Column = 1;

leftGrid = uigridlayout(leftPanel, [4, 1]);
leftGrid.RowHeight = {'fit', 'fit', 'fit', '1x'};
leftGrid.ColumnWidth = {'1x'};
leftGrid.Padding = [5 5 5 5];

% Load Parameter File
paramPanel = panelFcn(leftGrid, 'BorderType', 'none');
paramPanel.Layout.Row = 1;
paramGrid = uigridlayout(paramPanel, [1, 3]);
paramGrid.ColumnWidth = {'fit', '1x', 'fit'};
labelFcn(paramGrid, 'Text', 'Load Parameters:');
handles.paramLoadPathText = editFcn(paramGrid, 'Value', '', 'Editable', 'off', ...
    'Tooltip', 'Load saved parameter file (.mat)');
handles.paramLoadButton = buttonFcn(paramGrid, 'Text', 'Browse', ...
    'Tooltip', 'Browse for a saved parameter file to load');

% DIC Image
dicPanel = panelFcn(leftGrid, 'BorderType', 'none');
dicPanel.Layout.Row = 2;
dicGrid = uigridlayout(dicPanel, [1, 3]);
dicGrid.ColumnWidth = {'fit', '1x', 'fit'};
labelFcn(dicGrid, 'Text', 'DIC Image:');
handles.dicPathText = editFcn(dicGrid, 'Value', lastUsed.dicFilePath, 'Editable', 'off');
handles.dicBrowseButton = buttonFcn(dicGrid, 'Text', 'Browse');

% Number of Channels
numChanPanel = panelFcn(leftGrid, 'BorderType', 'none');
numChanPanel.Layout.Row = 3;
numChanGrid = uigridlayout(numChanPanel, [1, 2]);
numChanGrid.ColumnWidth = {'fit', '1x'};
labelFcn(numChanGrid, 'Text', 'Number of Fluorescent Channels:');
handles.numChanDrop = dropdownFcn(numChanGrid, 'Items', arrayfun(@num2str, 1:Nmax, 'UniformOutput', false), 'Value', num2str(lastUsed.numChannels));


% --- Pre-allocate Handle Arrays ---
handles.channelTabs = gobjects(1, Nmax);
handles.tabMainGrids = gobjects(1, Nmax);
handles.upButtons = gobjects(1, Nmax);
handles.downButtons = gobjects(1, Nmax);
handles.channelPathTexts = gobjects(1, Nmax);
handles.channelBrowseButtons = gobjects(1, Nmax);
handles.xySpacingInputs = gobjects(1, Nmax);
handles.zSpacingInputs = gobjects(1, Nmax);
% Deconvolution handles
handles.deconvPanels = gobjects(1, Nmax);
handles.deconvEnabledChecks = gobjects(1, Nmax);
handles.deconvMethodDrops = gobjects(1, Nmax);
handles.deconvParam1Labels = gobjects(1, Nmax);
handles.deconvParam2Labels = gobjects(1, Nmax);
handles.deconvLRIterationsInputs = gobjects(1, Nmax);
handles.deconvLRDampingInputs = gobjects(1, Nmax);
handles.deconvWienerNSRInputs = gobjects(1, Nmax);
handles.deconvBlindIterationsInputs = gobjects(1, Nmax);
handles.deconvBlindUnderRelaxInputs = gobjects(1, Nmax);
handles.deconvPSFSourceDrops = gobjects(1, Nmax);
handles.deconvPSFPathTexts = gobjects(1, Nmax);
handles.deconvPSFBrowseButtons = gobjects(1, Nmax);
handles.deconvPSFSigmaXYInputs = gobjects(1, Nmax);
handles.deconvPSFSigmaZInputs = gobjects(1, Nmax);
handles.deconvPSFSizeXYInputs = gobjects(1, Nmax);
handles.deconvPSFSizeZInputs = gobjects(1, Nmax);
handles.preprocessModeDrops = gobjects(1, Nmax);
handles.preprocessScaleChecks = gobjects(1, Nmax);
handles.preprocessProjectionDrops = gobjects(1, Nmax);
handles.preprocMethodDrops = gobjects(1, Nmax); % New
handles.preprocParam1Labels = gobjects(1, Nmax); % Dynamic parameter labels
handles.preprocParam1Inputs = gobjects(1, Nmax); % Dynamic parameter inputs
handles.preprocParam2Labels = gobjects(1, Nmax); % Dynamic secondary labels
handles.preprocParam2Inputs = gobjects(1, Nmax); % Dynamic secondary inputs
handles.preprocParam3Labels = gobjects(1, Nmax); % Wavelet threshold rule labels
handles.preprocParam3Inputs = gobjects(1, Nmax); % Wavelet threshold rule dropdowns
handles.preprocParam4Labels = gobjects(1, Nmax); % Wavelet threshold method labels
handles.preprocParam4Inputs = gobjects(1, Nmax); % Wavelet threshold method dropdowns
handles.gaussInputs = gobjects(1, Nmax); % Alias to param1
handles.medianInputs = gobjects(1, Nmax); % Alias to param1
handles.waveletNameDrops = gobjects(1, Nmax);
handles.waveletLevelInputs = gobjects(1, Nmax); % Alias to param2
handles.waveletThresholdRuleDrops = gobjects(1, Nmax);
handles.waveletThresholdMethodDrops = gobjects(1, Nmax);
handles.nlmFilterStrengthInputs = gobjects(1, Nmax); % Non-local means filter strength
handles.nlmSearchWindowInputs = gobjects(1, Nmax); % Non-local means search window
handles.nlmComparisonWindowInputs = gobjects(1, Nmax); % Non-local means comparison window
handles.preprocClipChecks = gobjects(1, Nmax); % Section-wide
handles.bgCorrModeDrops = gobjects(1, Nmax);
handles.bgCorrScaleChecks = gobjects(1, Nmax);
handles.bgCorrProjectionDrops = gobjects(1, Nmax);
handles.bgMethodDrops = gobjects(1, Nmax);
handles.bgParamInputs = gobjects(1, Nmax);
handles.bgCorrClipChecks = gobjects(1, Nmax); % Section-wide
handles.maximaModeDrops = gobjects(1, Nmax);
handles.maximaProjectionDrops = gobjects(1, Nmax);
handles.maximaMethodDrops = gobjects(1, Nmax);
handles.hMaxPanel = gobjects(1, Nmax);
handles.hMaxInputs = gobjects(1, Nmax);
handles.logPanel = gobjects(1, Nmax);
handles.logSigmaInputs = gobjects(1, Nmax);
handles.logThresholdInputs = gobjects(1, Nmax);
    handles.maximaNeighborhoodInputs = gobjects(1, Nmax);
    handles.maximaScaleChecks = gobjects(1, Nmax);
    handles.showMaximaChecks = gobjects(1, Nmax);
handles.maximaColorDrops = gobjects(1, Nmax);
handles.gaussFitVoxelWindowSlider = gobjects(1, Nmax);
handles.gaussFitBgCorrMethodDrop = gobjects(1, Nmax);
handles.gaussFitBgCorrWidthEdit = gobjects(1, Nmax);
handles.gaussFitPolyDegreeEdit = gobjects(1, Nmax);
handles.gaussFitMethodDrop = gobjects(1, Nmax);
handles.gaussFitPlotCheck = gobjects(1, Nmax);
handles.gaussFitMaxIterationsEdit = gobjects(1, Nmax);
handles.gaussFitToleranceEdit = gobjects(1, Nmax);
handles.gaussFitRadialRadiusEdit = gobjects(1, Nmax);
handles.gaussFitInitGuessPanel = gobjects(1, Nmax);
% Removed initial guess parameter handles - no longer needed

% --- New Handle Allocations for Enable/Disable Feature ---
% Nuclei
handles.nucPreprocPanel = gobjects(1);
handles.nucPreprocEnabledCheck = gobjects(1);
handles.nucPreprocClipChecks = gobjects(1);
handles.nucPreprocessModeDrop = gobjects(1);
handles.nucPreprocessScaleCheck = gobjects(1);
handles.nucPreprocessProjectionDrop = gobjects(1);
handles.nucPreprocMethodDrop = gobjects(1);
handles.nucPreprocParam1Labels = gobjects(1);
handles.nucPreprocParam1Inputs = gobjects(1);
handles.nucPreprocParam2Labels = gobjects(1);
handles.nucPreprocParam2Inputs = gobjects(1);
handles.nucPreprocParam3Labels = gobjects(1);
handles.nucPreprocParam3Inputs = gobjects(1);
handles.nucPreprocParam4Labels = gobjects(1);
handles.nucPreprocParam4Inputs = gobjects(1);
handles.nucWaveletNameDrop = gobjects(1);
handles.nucNlmFilterStrengthInput = gobjects(1);
handles.nucNlmSearchWindowInput = gobjects(1);
handles.nucNlmComparisonWindowInput = gobjects(1);
handles.segPanel = gobjects(1);
handles.nucSegEnabledCheck = gobjects(1);
% Channels
handles.preprocPanels = gobjects(1, Nmax);
handles.preprocEnabledChecks = gobjects(1, Nmax);
handles.bgPanels = gobjects(1, Nmax);
handles.bgCorrEnabledChecks = gobjects(1, Nmax);
handles.maximaPanels = gobjects(1, Nmax);
handles.maximaEnabledChecks = gobjects(1, Nmax);
handles.gaussFitPanels = gobjects(1, Nmax);
handles.gaussFitEnabledChecks = gobjects(1, Nmax);
handles.fitFilterPanels = gobjects(1, Nmax);
handles.fitFilterEnabledChecks = gobjects(1, Nmax);
handles.fitFilterRSquaredEnabledChecks = gobjects(1, Nmax);
handles.fitFilterRSquaredMinInputs = gobjects(1, Nmax);
handles.fitFilterRSquaredMaxInputs = gobjects(1, Nmax);
handles.fitFilterSigmaSumEnabledChecks = gobjects(1, Nmax);
handles.fitFilterSigmaSumMinInputs = gobjects(1, Nmax);
handles.fitFilterSigmaSumMaxInputs = gobjects(1, Nmax);
handles.fitFilterAmplitudeEnabledChecks = gobjects(1, Nmax);
handles.fitFilterAmplitudeMinInputs = gobjects(1, Nmax);
handles.fitFilterAmplitudeMaxInputs = gobjects(1, Nmax);
handles.fitFilterIntensityEnabledChecks = gobjects(1, Nmax);
handles.fitFilterIntensityMinInputs = gobjects(1, Nmax);
handles.fitFilterIntensityMaxInputs = gobjects(1, Nmax);

% === CLASSIFICATION HANDLES ===
handles.classifyPanels = gobjects(1, Nmax);
handles.classifyEnabledChecks = gobjects(1, Nmax);
handles.classifyStatusLabels = gobjects(1, Nmax);
handles.classifyFeatureCountLabels = gobjects(1, Nmax);
handles.classifySelectFeaturesButtons = gobjects(1, Nmax);
handles.classifyLoadButtons = gobjects(1, Nmax);
handles.classifyTrainButtons = gobjects(1, Nmax);
handles.classifyApplyButtons = gobjects(1, Nmax);
handles.classifyFilterNoiseChecks = gobjects(1, Nmax);
% Store classifier data per channel
handles.classifiers = cell(1, Nmax);
handles.classifierFeatures = cell(1, Nmax);
handles.classifierFeatureInfo = cell(1, Nmax);
handles.classifierTrainStats = cell(1, Nmax);

% Tab Group for Channels
tabGroup = uitabgroup(leftGrid);
tabGroup.Layout.Row = 4;
tabGroup.Layout.Column = 1;

% --- Nuclei Tab ---
nucTab = uitab(tabGroup, 'Title', 'Nuclei');
nucTabOuterGrid = uigridlayout(nucTab, [2, 1]);
nucTabOuterGrid.RowHeight = {'fit', '1x'};
nucTabOuterGrid.Padding = [5 5 5 5];

% File Path (Nuclei) with navigation buttons
nucPathPanel = panelFcn(nucTabOuterGrid, 'BorderType', 'none');
nucPathPanel.Layout.Row = 1;
nucPathGrid = uigridlayout(nucPathPanel, [1, 5]);
nucPathGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit'};
labelFcn(nucPathGrid, 'Text', 'Nuclei Path:');
handles.nucPathText = editFcn(nucPathGrid, 'Value', lastUsed.nucFilePath, 'Editable', 'off');
handles.nucBrowseButton = buttonFcn(nucPathGrid, 'Text', 'Browse');
handles.nucUpButton = buttonFcn(nucPathGrid, 'Text', '↑', 'FontSize', 16);
handles.nucDownButton = buttonFcn(nucPathGrid, 'Text', '↓', 'FontSize', 16);

% Set explicit column positions for navigation buttons
handles.nucUpButton.Layout.Column = 4;
handles.nucDownButton.Layout.Column = 5;

% --- Host panel for the main scrolling content ---
nucContentHostPanel = panelFcn(nucTabOuterGrid, 'BorderType', 'none');
nucContentHostPanel.Layout.Row = 2;

% --- Main Content Grid (this is what will "scroll") ---
nucMainContentGrid = uigridlayout(nucContentHostPanel, [7, 1]);
handles.nucMainContentGrid = nucMainContentGrid;

% Set initial navigation state for nuclei
nucInitialState = lastUsed.nucNavPanelIndex;
% Grid rows: 1=Spacing, 2=Deconv, 3=Preproc, 4=BG, 5=Segmentation, 6=Filter, 7=Spacer
if nucInitialState == 0 % Stage 0 - show panels 1-5 (Spacing, Deconv, Preproc, BG, Segmentation)
    nucMainContentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 0, '1x'};
elseif nucInitialState == 1 % Stage 1 - show panels 2-6 (Deconv, Preproc, BG, Segmentation, Filter)
    nucMainContentGrid.RowHeight = {0, 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
elseif nucInitialState == 2 % Stage 2 - show panels 3-6 (Preproc, BG, Segmentation, Filter)
    nucMainContentGrid.RowHeight = {0, 0, 'fit', 'fit', 'fit', 'fit', '1x'};
elseif nucInitialState == 3 % Stage 3 - show panels 4-6 (BG, Segmentation, Filter)
    nucMainContentGrid.RowHeight = {0, 0, 0, 'fit', 'fit', 'fit', '1x'};
else % Default fallback - Stage 0
    nucMainContentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 0, '1x'};
end

% Image Spacing (Nuclei)
nucSpacingPanel = panelFcn(nucMainContentGrid, 'Title', 'Image Spacing');
nucSpacingPanel.Layout.Row = 1;
nucSpacingGrid = uigridlayout(nucSpacingPanel, [1, 6]);
labelFcn(nucSpacingGrid, 'Text', 'XY Pixel:');
handles.nucXYSpacingInput = editFcn(nucSpacingGrid, 'numeric', 'Value', lastUsed.nucXYSpacing, 'Tooltip', 'Pixel width in XY (µm).');
labelFcn(nucSpacingGrid, 'Text', '(µm)');
labelFcn(nucSpacingGrid, 'Text', 'Z-Frame:');
handles.nucZSpacingInput = editFcn(nucSpacingGrid, 'numeric', 'Value', lastUsed.nucZSpacing, 'Tooltip', 'Slice depth in Z (µm).', 'Enable', 'off');
labelFcn(nucSpacingGrid, 'Text', '(µm)');

% --- Deconvolution Panel (Nuclei) ---
nucDeconvPanel = panelFcn(nucMainContentGrid, 'Title', '');
nucDeconvPanel.Layout.Row = 2;
handles.nucDeconvPanel = nucDeconvPanel;
nucDeconvOuterGrid = uigridlayout(nucDeconvPanel, [2, 1]);
nucDeconvOuterGrid.RowHeight = {'fit', '1x'};
nucDeconvOuterGrid.Padding = [2 2 2 2];
% Title-like row
titleGrid = uigridlayout(nucDeconvOuterGrid, [1, 2]);
titleGrid.ColumnWidth = {'1x', 'fit'};
titleGrid.Padding = [0 0 0 0];
labelFcn(titleGrid, 'Text', 'Deconvolution', 'FontWeight', 'bold');
handles.nucDeconvEnabledCheck = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.nucDeconvEnabled);
% Content grid - 6 rows
nucDeconvGrid = uigridlayout(nucDeconvOuterGrid, [6, 4]);
nucDeconvGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
nucDeconvGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
nucDeconvGrid.Padding = [5 5 5 5];
nucDeconvGrid.RowSpacing = 5;
nucDeconvGrid.ColumnSpacing = 10;

% Row 1: Method dropdown
labelFcn(nucDeconvGrid, 'Text', 'Method:');
handles.nucDeconvMethodDrop = dropdownFcn(nucDeconvGrid, 'Items', {'Lucy-Richardson', 'Wiener', 'Blind'}, 'Value', lastUsed.nucDeconvMethod, 'Tooltip', 'Deconvolution algorithm.');
handles.nucDeconvMethodDrop.Layout.Row = 1;
handles.nucDeconvMethodDrop.Layout.Column = [2 4];

% Row 2: Dynamic parameters
handles.nucDeconvParam1Label = labelFcn(nucDeconvGrid, 'Text', 'Iterations:', 'Tag', 'nuc_deconv_param');
handles.nucDeconvParam1Label.Layout.Row = 2;
handles.nucDeconvParam1Label.Layout.Column = 1;
handles.nucDeconvLRIterationsInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvLRIterations, 'Tag', 'nuc_deconv_param', 'Tooltip', 'Number of iterations.');
handles.nucDeconvLRIterationsInput.Layout.Row = 2;
handles.nucDeconvLRIterationsInput.Layout.Column = 2;

handles.nucDeconvParam2Label = labelFcn(nucDeconvGrid, 'Text', 'Damping:', 'Tag', 'nuc_deconv_param');
handles.nucDeconvParam2Label.Layout.Row = 2;
handles.nucDeconvParam2Label.Layout.Column = 3;
handles.nucDeconvLRDampingInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvLRDamping, 'Tag', 'nuc_deconv_param', 'Tooltip', 'Damping/under-relaxation parameter.');
handles.nucDeconvLRDampingInput.Layout.Row = 2;
handles.nucDeconvLRDampingInput.Layout.Column = 4;

% Shared controls
handles.nucDeconvWienerNSRInput = handles.nucDeconvLRIterationsInput;
handles.nucDeconvBlindIterationsInput = handles.nucDeconvLRIterationsInput;
handles.nucDeconvBlindUnderRelaxInput = handles.nucDeconvLRDampingInput;

% Row 3: PSF Source
labelFcn(nucDeconvGrid, 'Text', 'PSF Source:');
handles.nucDeconvPSFSourceDrop = dropdownFcn(nucDeconvGrid, 'Items', {'Generate', 'Load File'}, 'Value', lastUsed.nucDeconvPSFSource, 'Tooltip', 'Generate Gaussian PSF or load experimental PSF image.');
handles.nucDeconvPSFSourceDrop.Layout.Row = 3;
handles.nucDeconvPSFSourceDrop.Layout.Column = [2 4];

% Row 4: PSF File Path
labelFcn(nucDeconvGrid, 'Text', 'PSF File:', 'Tag', 'nuc_deconv_psf_file', 'Visible', 'off');
handles.nucDeconvPSFPathText = editFcn(nucDeconvGrid, 'Value', lastUsed.nucDeconvPSFFilePath, 'Editable', 'off', 'Tag', 'nuc_deconv_psf_file', 'Visible', 'off');
handles.nucDeconvPSFPathText.Layout.Row = 4;
handles.nucDeconvPSFPathText.Layout.Column = 2;
handles.nucDeconvPSFBrowseButton = buttonFcn(nucDeconvGrid, 'Text', 'Browse', 'Tag', 'nuc_deconv_psf_file', 'Visible', 'off');
handles.nucDeconvPSFBrowseButton.Layout.Row = 4;
handles.nucDeconvPSFBrowseButton.Layout.Column = [3 4];

% Row 5: PSF Generation parameters
labelFcn(nucDeconvGrid, 'Text', 'PSF Sigma XY:', 'Tag', 'nuc_deconv_psf_gen');
handles.nucDeconvPSFSigmaXYInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvPSFSigmaXY, 'Tag', 'nuc_deconv_psf_gen', 'Tooltip', 'PSF sigma in XY plane (pixels).');
handles.nucDeconvPSFSigmaXYInput.Layout.Row = 5;
handles.nucDeconvPSFSigmaXYInput.Layout.Column = 2;
labelFcn(nucDeconvGrid, 'Text', 'Sigma Z:', 'Tag', 'nuc_deconv_psf_gen');
handles.nucDeconvPSFSigmaZInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvPSFSigmaZ, 'Tag', 'nuc_deconv_psf_gen', 'Tooltip', 'PSF sigma in Z direction (slices).');
handles.nucDeconvPSFSigmaZInput.Layout.Row = 5;
handles.nucDeconvPSFSigmaZInput.Layout.Column = 4;

% Row 6: PSF Size parameters
labelFcn(nucDeconvGrid, 'Text', 'PSF Size XY:', 'Tag', 'nuc_deconv_psf_gen');
handles.nucDeconvPSFSizeXYInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvPSFSizeXY, 'Tag', 'nuc_deconv_psf_gen', 'Tooltip', 'PSF half-size in XY (pixels).');
handles.nucDeconvPSFSizeXYInput.Layout.Row = 6;
handles.nucDeconvPSFSizeXYInput.Layout.Column = 2;
labelFcn(nucDeconvGrid, 'Text', 'Size Z:', 'Tag', 'nuc_deconv_psf_gen');
handles.nucDeconvPSFSizeZInput = editFcn(nucDeconvGrid, 'numeric', 'Value', lastUsed.nucDeconvPSFSizeZ, 'Tag', 'nuc_deconv_psf_gen', 'Tooltip', 'PSF half-size in Z (slices).');
handles.nucDeconvPSFSizeZInput.Layout.Row = 6;
handles.nucDeconvPSFSizeZInput.Layout.Column = 4;

% --- Pre-processing Panel (Nuclei) ---
nucPreprocPanel = panelFcn(nucMainContentGrid, 'Title', '');
nucPreprocPanel.Layout.Row = 3;
handles.nucPreprocPanel = nucPreprocPanel; % Store handle to the panel
nucPreprocOuterGrid = uigridlayout(nucPreprocPanel, [2, 1]);
nucPreprocOuterGrid.RowHeight = {'fit', '1x'};
nucPreprocOuterGrid.Padding = [2 2 2 2];
% Title-like row
titleGrid = uigridlayout(nucPreprocOuterGrid, [1, 2]);
titleGrid.ColumnWidth = {'1x', 'fit'};
titleGrid.Padding = [0 0 0 0];
labelFcn(titleGrid, 'Text', 'Pre-processing', 'FontWeight', 'bold');
handles.nucPreprocEnabledCheck = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.nucPreprocEnabled);
% Content grid
nucPreprocGrid = uigridlayout(nucPreprocOuterGrid, [6, 4]); % 6 rows, 4 columns for wavelet parameters
nucPreprocGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
nucPreprocGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
nucPreprocGrid.Padding = [5 5 5 5];
nucPreprocGrid.RowSpacing = 5;
nucPreprocGrid.ColumnSpacing = 10;

% Row 1: Clip checkbox (spans all columns)
handles.nucPreprocClipChecks = checkboxFcn(nucPreprocGrid, 'Text', 'Clip Negative Values at Zero', 'Value', lastUsed.nucPreprocClipAtZero);
handles.nucPreprocClipChecks.Layout.Row = 1;
handles.nucPreprocClipChecks.Layout.Column = [1 4];

% Row 2: Mode controls
labelFcn(nucPreprocGrid, 'Text', 'Mode:');
handles.nucPreprocessModeDrop = dropdownFcn(nucPreprocGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.nucPreProcMode, 'Tooltip', 'Choose how to apply pre-processing to nuclei.');
handles.nucPreprocessModeDrop.Layout.Column = 2;
handles.nucPreprocessModeDrop.Layout.Row = 2;
handles.nucPreprocessScaleCheck = checkboxFcn(nucPreprocGrid, 'Text', 'Scale', 'Value', lastUsed.nucPreProcScale, 'Enable', 'off', 'Tooltip', 'Use physical units (µm) for 3D operations.');
handles.nucPreprocessScaleCheck.Layout.Column = 3;
handles.nucPreprocessScaleCheck.Layout.Row = 2;
handles.nucPreprocessProjectionDrop = dropdownFcn(nucPreprocGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.nucPreProcProjection, 'Tooltip', 'Z-projection type.');
handles.nucPreprocessProjectionDrop.Layout.Column = 4;
handles.nucPreprocessProjectionDrop.Layout.Row = 2;

% Row 3: Method controls
labelFcn(nucPreprocGrid, 'Text', 'Method:');
handles.nucPreprocMethodDrop = dropdownFcn(nucPreprocGrid, 'Items', {'None', 'Gaussian', 'Median', 'Non-Local Means', 'Wavelet Denoising'}, 'Value', lastUsed.nucPreProcMethod, 'Tooltip', 'Algorithm for pre-processing.');
handles.nucPreprocMethodDrop.Layout.Column = [2 4];
handles.nucPreprocMethodDrop.Layout.Row = 3;

% Row 4: Dynamic parameter controls (content changes based on selected method)
% Primary parameter label and input (used by all methods)
handles.nucPreprocParam1Labels = labelFcn(nucPreprocGrid, 'Text', 'Parameter:', 'Tag', 'nuc_preproc_param_ctrl');
handles.nucPreprocParam1Labels.Layout.Row = 4;
handles.nucPreprocParam1Labels.Layout.Column = 1;

% Primary parameter input (Gaussian sigma, Median neighborhood, or Wavelet name)
handles.nucPreprocParam1Inputs = editFcn(nucPreprocGrid, 'numeric', 'Value', lastUsed.nucSmoothGaussianValue, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Parameter value for the selected method.');
handles.nucPreprocParam1Inputs.Layout.Row = 4;
handles.nucPreprocParam1Inputs.Layout.Column = 2;

% Units/secondary label
handles.nucPreprocParam2Labels = labelFcn(nucPreprocGrid, 'Text', '(units)', 'Tag', 'nuc_preproc_param_ctrl');
handles.nucPreprocParam2Labels.Layout.Row = 4;
handles.nucPreprocParam2Labels.Layout.Column = 3;

% Secondary parameter input (only used for wavelet level)
handles.nucPreprocParam2Inputs = editFcn(nucPreprocGrid, 'numeric', 'Value', lastUsed.nucWaveletLevel, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Secondary parameter for the selected method.', 'Visible', 'off');
handles.nucPreprocParam2Inputs.Layout.Row = 4;
handles.nucPreprocParam2Inputs.Layout.Column = 4;

% Row 5: Wavelet-specific parameters (Threshold Rule)
handles.nucPreprocParam3Labels = labelFcn(nucPreprocGrid, 'Text', 'Rule:', 'Tag', 'nuc_preproc_param_ctrl', 'Visible', 'off');
handles.nucPreprocParam3Labels.Layout.Row = 5;
handles.nucPreprocParam3Labels.Layout.Column = 1;
handles.nucPreprocParam3Inputs = dropdownFcn(nucPreprocGrid, 'Items', {'sqtwolog','rigrsure','heursure','minimaxi'}, 'Value', lastUsed.nucWaveletThresholdRule, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Threshold selection rule.', 'Visible', 'off');
handles.nucPreprocParam3Inputs.Layout.Row = 5;
handles.nucPreprocParam3Inputs.Layout.Column = 2;

% Row 6: Wavelet threshold method
handles.nucPreprocParam4Labels = labelFcn(nucPreprocGrid, 'Text', 'Method:', 'Tag', 'nuc_preproc_param_ctrl', 'Visible', 'off');
handles.nucPreprocParam4Labels.Layout.Row = 6;
handles.nucPreprocParam4Labels.Layout.Column = 1;
handles.nucPreprocParam4Inputs = dropdownFcn(nucPreprocGrid, 'Items', {'soft', 'hard'}, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Thresholding method.', 'Visible', 'off');
handles.nucPreprocParam4Inputs.Layout.Row = 6;
handles.nucPreprocParam4Inputs.Layout.Column = 2;

% Store original parameter handles for backward compatibility
handles.nucGaussInput = handles.nucPreprocParam1Inputs;
handles.nucMedianInput = handles.nucPreprocParam1Inputs;
handles.nucWaveletNameDrop = dropdownFcn(nucPreprocGrid, 'Items', {'haar','db2','db3','db4','sym2','sym3','sym4','coif1'}, 'Value', lastUsed.nucWaveletName, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Type of wavelet.', 'Visible', 'off');
handles.nucWaveletNameDrop.Layout.Row = 4;
handles.nucWaveletNameDrop.Layout.Column = 2;
handles.nucWaveletLevelInput = handles.nucPreprocParam2Inputs;
handles.nucWaveletThresholdRuleDrop = handles.nucPreprocParam3Inputs;
handles.nucWaveletThresholdMethodDrop = handles.nucPreprocParam4Inputs;

% Non-Local Means specific parameter controls (hidden by default)
handles.nucNlmFilterStrengthInput = editFcn(nucPreprocGrid, 'numeric', 'Value', lastUsed.nucNlmFilterStrength, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Filter strength for non-local means denoising.', 'Visible', 'off');
handles.nucNlmFilterStrengthInput.Layout.Row = 4;
handles.nucNlmFilterStrengthInput.Layout.Column = 2;

handles.nucNlmSearchWindowInput = editFcn(nucPreprocGrid, 'numeric', 'Value', lastUsed.nucNlmSearchWindow, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Search window size for non-local means.', 'Visible', 'off');
handles.nucNlmSearchWindowInput.Layout.Row = 5;
handles.nucNlmSearchWindowInput.Layout.Column = 2;

handles.nucNlmComparisonWindowInput = editFcn(nucPreprocGrid, 'numeric', 'Value', lastUsed.nucNlmComparisonWindow, 'Tag', 'nuc_preproc_param_ctrl', 'Tooltip', 'Comparison window size for non-local means.', 'Visible', 'off');
handles.nucNlmComparisonWindowInput.Layout.Row = 6;
handles.nucNlmComparisonWindowInput.Layout.Column = 2;


% --- Background Correction Panel (Nuclei) ---
nucBgPanel = panelFcn(nucMainContentGrid, 'Title', '');
nucBgPanel.Layout.Row = 4;
handles.nucBgPanel = nucBgPanel; % Store handle to the panel
nucBgOuterGrid = uigridlayout(nucBgPanel, [2, 1]);
nucBgOuterGrid.RowHeight = {'fit', '1x'};
nucBgOuterGrid.Padding = [2 2 2 2];
% Title-like row
titleGrid = uigridlayout(nucBgOuterGrid, [1, 2]);
titleGrid.ColumnWidth = {'1x', 'fit'};
titleGrid.Padding = [0 0 0 0];
labelFcn(titleGrid, 'Text', 'Background Correction', 'FontWeight', 'bold');
handles.nucBgCorrEnabledCheck = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.nucBgCorrEnabled);
% Content grid
nucBgGrid = uigridlayout(nucBgOuterGrid, [4, 4]); % 4 rows, 4 columns
nucBgGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
nucBgGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
nucBgGrid.Padding = [5 5 5 5];
nucBgGrid.RowSpacing = 5;
nucBgGrid.ColumnSpacing = 10;

% Row 1: Clip checkbox (spans all columns)
handles.nucBgCorrClipChecks = checkboxFcn(nucBgGrid, 'Text', 'Clip Negative Values at Zero', 'Value', lastUsed.nucBgCorrClipAtZero);
handles.nucBgCorrClipChecks.Layout.Row = 1;
handles.nucBgCorrClipChecks.Layout.Column = [1 4];

% Row 2: Mode controls
labelFcn(nucBgGrid, 'Text', 'Mode:');
handles.nucBgCorrModeDrop = dropdownFcn(nucBgGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.nucBgCorrMode, 'Tooltip', 'Choose how to apply background correction to nuclei.');
handles.nucBgCorrModeDrop.Layout.Column = 2;
handles.nucBgCorrModeDrop.Layout.Row = 2;
handles.nucBgCorrScaleCheck = checkboxFcn(nucBgGrid, 'Text', 'Scale', 'Value', lastUsed.nucBgCorrScale, 'Enable', 'off', 'Tooltip', 'Use physical units (µm) for 3D operations.');
handles.nucBgCorrScaleCheck.Layout.Column = 3;
handles.nucBgCorrScaleCheck.Layout.Row = 2;
handles.nucBgCorrProjectionDrop = dropdownFcn(nucBgGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.nucBgCorrProjection, 'Tooltip', 'Z-projection type.');
handles.nucBgCorrProjectionDrop.Layout.Column = 4;
handles.nucBgCorrProjectionDrop.Layout.Row = 2;

% Row 3: Method controls
labelFcn(nucBgGrid, 'Text', 'Method:');
handles.nucBgMethodDrop = dropdownFcn(nucBgGrid, 'Items', {'None', 'Gaussian', 'Rolling-ball', 'Top-hat'}, 'Value', lastUsed.nucBgMethod, 'Tooltip', 'Algorithm for background subtraction.');
handles.nucBgMethodDrop.Layout.Column = [2 4];
handles.nucBgMethodDrop.Layout.Row = 3;

% Row 4: Parameter controls
labelFcn(nucBgGrid, 'Text', 'Parameter:');
handles.nucBgParamInput = editFcn(nucBgGrid, 'numeric', 'Value', lastUsed.nucBgParam, 'Tooltip', 'Gaussian: Sigma. Rolling-ball/Top-hat: Radius. Units are µm if "Scale" is checked (3D), otherwise pixels/voxels.');
handles.nucBgParamInput.Layout.Row = 4;
handles.nucBgParamInput.Layout.Column = 2;
unitsLabel = labelFcn(nucBgGrid, 'Text', '(units)');
unitsLabel.Layout.Row = 4;
unitsLabel.Layout.Column = 3;

% --- Segmentation Panel ---
segPanel = panelFcn(nucMainContentGrid, 'Title', '');
segPanel.Layout.Row = 5;
handles.segPanel = segPanel; % Store handle to the panel
segOuterGrid = uigridlayout(segPanel, [2, 1]);
segOuterGrid.RowHeight = {'fit', '1x'};
segOuterGrid.Padding = [2 2 2 2];
% Title-like row
titleGrid = uigridlayout(segOuterGrid, [1, 2]);
titleGrid.ColumnWidth = {'1x', 'fit'};
titleGrid.Padding = [0 0 0 0];
labelFcn(titleGrid, 'Text', 'Segmentation', 'FontWeight', 'bold');
handles.nucSegEnabledCheck = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.nucSegEnabled);
% Content grid
segGrid = uigridlayout(segOuterGrid);
segGrid.RowHeight = {'fit', 'fit', 'fit'};
segGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};

% Segmentation Mode (Row 1)
labelFcn(segGrid, 'Text', 'Mode:');
handles.nucSegModeDrop = dropdownFcn(segGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.nucSegMode, 'Tooltip', 'Choose how to apply segmentation.');
handles.nucSegModeDrop.Layout.Column = [2 3];
handles.nucSegProjectionDrop = dropdownFcn(segGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.nucSegProjection, 'Tooltip', 'Z-projection type for segmentation.');


% Nuclei Segmentation Methods Panel (Row 2)
handles.nucSegMethodsPanel = panelFcn(segGrid, 'BorderType', 'none', 'Visible', 'on');
handles.nucSegMethodsPanel.Layout.Row = 2;
handles.nucSegMethodsPanel.Layout.Column = [1 4];
methodsGrid = uigridlayout(handles.nucSegMethodsPanel, [6, 4]); % Clean 4-column layout like preprocessing
methodsGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
methodsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};

% Row 1: Main method selection
labelFcn(methodsGrid, 'Text', 'Method:', 'Tag', 'nuc_seg_control');
handles.nucSegMainMethodDrop = dropdownFcn(methodsGrid, 'Items', {'Absolute', 'Mean', 'Median', 'Auto Local Threshold'}, 'Value', lastUsed.nucSegMainMethod, 'Tooltip', 'Segmentation method: Absolute=direct threshold, Mean/Median=statistical, Auto Local=adaptive algorithms', 'Tag', 'nuc_seg_control');
handles.nucSegMainMethodDrop.Layout.Column = [2 4];

% Row 2: Sub-method selection (visible for Mean/Median and Auto Local)
handles.nucSegSubMethodLabel = labelFcn(methodsGrid, 'Text', 'Type:', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegSubMethodLabel.Layout.Row = 2;
handles.nucSegSubMethodDrop = dropdownFcn(methodsGrid, 'Items', {'Std Multiplier', 'Absolute Offset'}, 'Value', lastUsed.nucSegSubMethod, 'Tooltip', 'Std Multiplier: adaptive to image noise. Absolute Offset: fixed intensity difference.', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegSubMethodDrop.Layout.Row = 2;
handles.nucSegSubMethodDrop.Layout.Column = [2 4];

% Row 3: Dynamic parameter controls (content changes based on selected method)
% Primary parameter label and input (used by all methods)
handles.nucSegParam1Label = labelFcn(methodsGrid, 'Text', 'Parameter:', 'Tag', 'nuc_seg_control');
handles.nucSegParam1Label.Layout.Row = 3;
handles.nucSegParam1Label.Layout.Column = 1;

% Primary parameter input (Threshold, Std×, Offset, or Algorithm)
handles.nucSegParam1Input = editFcn(methodsGrid, 'numeric', 'Value', lastUsed.nucSegAbsoluteThreshold, 'Tooltip', 'Parameter value for the selected method.', 'Tag', 'nuc_seg_control');
handles.nucSegParam1Input.Layout.Row = 3;
handles.nucSegParam1Input.Layout.Column = 2;

% Units/secondary label
handles.nucSegParam2Label = labelFcn(methodsGrid, 'Text', '(units)', 'Tag', 'nuc_seg_control');
handles.nucSegParam2Label.Layout.Row = 3;
handles.nucSegParam2Label.Layout.Column = 3;

% Secondary parameter input removed - now handled by algorithm-specific parameters
% (This was redundant with the new dual-parameter system)

% Row 4: Auto Local Threshold algorithm selection (only visible for Auto Local Threshold)
handles.nucSegAlgorithmLabel = labelFcn(methodsGrid, 'Text', 'Algorithm:', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgorithmLabel.Layout.Row = 4;
handles.nucSegAlgorithmLabel.Layout.Column = 1;
handles.nucSegAlgorithmDrop = dropdownFcn(methodsGrid, 'Items', {'Bernsen', 'Contrast', 'Mean', 'Median', 'MidGrey', 'Niblack', 'Otsu', 'Phansalkar', 'Sauvola'}, 'Value', lastUsed.nucSegLocalAlgorithm, 'Tooltip', 'Local thresholding algorithm. Each adapts differently to local image properties.', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgorithmDrop.Layout.Row = 4;
handles.nucSegAlgorithmDrop.Layout.Column = [2 4];

% Row 5: Algorithm-specific parameters (only visible for Auto Local Threshold)
handles.nucSegAlgParamLabel = labelFcn(methodsGrid, 'Text', 'Parameter:', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParamLabel.Layout.Row = 5;
handles.nucSegAlgParamLabel.Layout.Column = 1;
handles.nucSegAlgParamInput = editFcn(methodsGrid, 'numeric', 'Value', lastUsed.nucSegAlgParam1, 'Tooltip', 'Algorithm-specific parameter value.', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParamInput.Layout.Row = 5;
handles.nucSegAlgParamInput.Layout.Column = 2;
handles.nucSegAlgParamDefaultCheck = checkboxFcn(methodsGrid, 'Text', 'Default', 'Value', lastUsed.nucSegAlgParam1Default, 'Tooltip', 'Use ImageJ default value', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParamDefaultCheck.Layout.Row = 5;
handles.nucSegAlgParamDefaultCheck.Layout.Column = 3;

% Row 6: Second algorithm-specific parameter (for dual-parameter algorithms)
handles.nucSegAlgParam2Label = labelFcn(methodsGrid, 'Text', 'Parameter 2:', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParam2Label.Layout.Row = 6;
handles.nucSegAlgParam2Label.Layout.Column = 1;
handles.nucSegAlgParam2Input = editFcn(methodsGrid, 'numeric', 'Value', lastUsed.nucSegAlgParam2, 'Tooltip', 'Second algorithm-specific parameter value.', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParam2Input.Layout.Row = 6;
handles.nucSegAlgParam2Input.Layout.Column = 2;
handles.nucSegAlgParam2DefaultCheck = checkboxFcn(methodsGrid, 'Text', 'Default', 'Value', lastUsed.nucSegAlgParam2Default, 'Tooltip', 'Use ImageJ default value', 'Tag', 'nuc_seg_control', 'Visible', 'off');
handles.nucSegAlgParam2DefaultCheck.Layout.Row = 6;
handles.nucSegAlgParam2DefaultCheck.Layout.Column = 3;

% Store original parameter handles for backward compatibility
handles.nucSegAbsoluteInput = handles.nucSegParam1Input;
handles.nucSegStdMultiplierInput = handles.nucSegParam1Input;
handles.nucSegOffsetInput = handles.nucSegParam1Input;
handles.nucSegLocalAlgorithmDrop = handles.nucSegAlgorithmDrop;
% Note: nucSegLocalRadiusInput removed - radius now handled by algorithm-specific parameters

% Note: Old Adaptive Thresholding Panel removed - replaced by Auto Local Threshold methods

% Show Segmentation Checkbox (Row 3)
handles.nucShowSegCheck = checkboxFcn(segGrid, 'Text', 'Show Segmentation Overlay', 'Value', lastUsed.nucShowSeg);
handles.nucShowSegCheck.Layout.Row = 3;
handles.nucShowSegCheck.Layout.Column = [1 4];



% --- Filter Panel ---
nucFilterPanel = panelFcn(nucMainContentGrid, 'Title', '');
nucFilterPanel.Layout.Row = 6;
handles.nucFilterPanel = nucFilterPanel; % Store handle to the panel
nucFilterOuterGrid = uigridlayout(nucFilterPanel, [2, 1]);
nucFilterOuterGrid.RowHeight = {'fit', '1x'};
nucFilterOuterGrid.Padding = [2 2 2 2];
% Title-like row
titleGrid = uigridlayout(nucFilterOuterGrid, [1, 2]);
titleGrid.ColumnWidth = {'1x', 'fit'};
titleGrid.Padding = [0 0 0 0];
labelFcn(titleGrid, 'Text', 'Filter', 'FontWeight', 'bold');
handles.nucFilterEnabledCheck = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.nucFilterEnabled);
% Content grid
nucFilterGrid = uigridlayout(nucFilterOuterGrid, [6, 4]);
nucFilterGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
nucFilterGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
nucFilterGrid.Padding = [5 5 5 5];
nucFilterGrid.RowSpacing = 5;
nucFilterGrid.ColumnSpacing = 10;

% Row 1: Size filter controls
handles.nucFilterSizeEnabledCheck = checkboxFcn(nucFilterGrid, 'Text', 'Min Size:', 'Value', lastUsed.nucFilterSizeEnabled, 'Tooltip', 'Enable filtering by minimum size threshold');
handles.nucFilterSizeEnabledCheck.Layout.Row = 1;
handles.nucFilterSizeEnabledCheck.Layout.Column = 1;
handles.nucFilterMinSizeInput = editFcn(nucFilterGrid, 'numeric', 'Value', lastUsed.nucFilterMinSize, 'Tooltip', 'Objects smaller than this will be removed');
handles.nucFilterMinSizeInput.Layout.Row = 1;
handles.nucFilterMinSizeInput.Layout.Column = 2;
labelFcn(nucFilterGrid, 'Text', 'Unit:');
handles.nucFilterSizeUnitDrop = dropdownFcn(nucFilterGrid, 'Items', {'pixels', 'voxels', 'microns^2', 'microns^3'}, 'Value', lastUsed.nucFilterSizeUnit, 'Tooltip', 'Unit for size measurement');

% Row 2: Circularity filter controls
handles.nucFilterCircularityEnabledCheck = checkboxFcn(nucFilterGrid, 'Text', 'Min Circularity:', 'Value', lastUsed.nucFilterCircularityEnabled, 'Tooltip', 'Enable filtering by shape circularity/sphericity');
handles.nucFilterCircularityEnabledCheck.Layout.Row = 2;
handles.nucFilterCircularityEnabledCheck.Layout.Column = 1;
handles.nucFilterMinCircularityInput = editFcn(nucFilterGrid, 'numeric', 'Value', lastUsed.nucFilterMinCircularity, 'Tooltip', 'Minimum circularity (2D) or sphericity (3D). Range: 0.0-1.0, where 1.0 = perfect circle/sphere');
handles.nucFilterMinCircularityInput.Layout.Row = 2;
handles.nucFilterMinCircularityInput.Layout.Column = 2;
circularityRangeLabel = labelFcn(nucFilterGrid, 'Text', '(0.0-1.0)', 'Tooltip', 'Circularity range: 0.0 (very elongated) to 1.0 (perfect circle/sphere)');
circularityRangeLabel.Layout.Row = 2;
circularityRangeLabel.Layout.Column = [3 4];

% Row 3: Solidity filter controls
handles.nucFilterSolidityEnabledCheck = checkboxFcn(nucFilterGrid, 'Text', 'Min Solidity:', 'Value', lastUsed.nucFilterSolidityEnabled, 'Tooltip', 'Enable filtering by solidity (measures convexity)');
handles.nucFilterSolidityEnabledCheck.Layout.Row = 3;
handles.nucFilterSolidityEnabledCheck.Layout.Column = 1;
handles.nucFilterMinSolidityInput = editFcn(nucFilterGrid, 'numeric', 'Value', lastUsed.nucFilterMinSolidity, 'Tooltip', 'Minimum solidity. Range: 0.0-1.0, where 1.0 = perfectly convex, lower values = has invaginations/concavities');
handles.nucFilterMinSolidityInput.Layout.Row = 3;
handles.nucFilterMinSolidityInput.Layout.Column = 2;
solidityRangeLabel = labelFcn(nucFilterGrid, 'Text', '(0.0-1.0)', 'Tooltip', 'Solidity range: 1.0 (smooth/convex) to lower values (has invaginations)');
solidityRangeLabel.Layout.Row = 3;
solidityRangeLabel.Layout.Column = [3 4];

% Row 4: Exclude Edge Nuclei checkbox
handles.nucExcludeEdgesCheck = checkboxFcn(nucFilterGrid, 'Text', 'Exclude Edge Nuclei', 'Value', lastUsed.nucExcludeEdges, 'Tooltip', 'Remove nuclei that touch the edges of the image.');
handles.nucExcludeEdgesCheck.Layout.Row = 4;
handles.nucExcludeEdgesCheck.Layout.Column = [1 4];

% Row 5: Filter Local Maxima by Nuclei checkbox and mode
handles.nucInclusionExclusionEnabledCheck = checkboxFcn(nucFilterGrid, 'Text', 'Filter Local Maxima by Nuclei', 'Value', lastUsed.nucInclusionExclusionEnabled);
handles.nucInclusionExclusionEnabledCheck.Layout.Row = 5;
handles.nucInclusionExclusionEnabledCheck.Layout.Column = [1 2];
handles.nucInclusionExclusionModeDrop = dropdownFcn(nucFilterGrid, 'Items', {'Include Inside Nuclei', 'Exclude Inside Nuclei'}, 'Value', lastUsed.nucInclusionExclusionMode, 'Tooltip', 'Include: Keep only maxima inside nuclei. Exclude: Remove maxima inside nuclei.');
handles.nucInclusionExclusionModeDrop.Layout.Row = 5;
handles.nucInclusionExclusionModeDrop.Layout.Column = [3 4];

% Row 6: Apply to channels
handles.nucInclusionExclusionApplyLabel = labelFcn(nucFilterGrid, 'Text', 'Apply to Channels:', 'Tag', 'nuc_inclusion_control');
handles.nucInclusionExclusionApplyLabel.Layout.Row = 6;
handles.nucInclusionExclusionApplyLabel.Layout.Column = 1;
handles.nucInclusionExclusionApplyDrop = dropdownFcn(nucFilterGrid, 'Items', {'All Channels', 'Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Channel 5'}, 'Value', lastUsed.nucInclusionExclusionApplyTo, 'Tooltip', 'Which channels to apply the nuclei filtering to.', 'Tag', 'nuc_inclusion_control');
handles.nucInclusionExclusionApplyDrop.Layout.Row = 6;
handles.nucInclusionExclusionApplyDrop.Layout.Column = [2 4];

% Spacer panel to fill remaining vertical space
spacer = panelFcn(nucMainContentGrid, 'BorderType', 'none');
spacer.Layout.Row = 7;

for k = 1:Nmax
    tab = uitab(tabGroup, 'Title', ['Channel ' num2str(k)]);
    handles.channelTabs(k) = tab;
    
    % This grid holds the top controls and the main content grid
    outerGrid = uigridlayout(tab, [2, 1]);
    outerGrid.RowHeight = {'fit', '1x'};

    % --- Top Controls: File Path, Browse, and Navigation Arrows ---
    topControlsPanel = panelFcn(outerGrid, 'BorderType', 'none');
    topControlsPanel.Layout.Row = 1;
    topControlsGrid = uigridlayout(topControlsPanel, [1, 5]);
    topControlsGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit'};
    labelFcn(topControlsGrid, 'Text', ['Path ' num2str(k) ':']);
    handles.channelPathTexts(k) = editFcn(topControlsGrid, 'Value', lastUsed.channelFilePaths{k}, 'Editable', 'off');
    handles.channelBrowseButtons(k) = buttonFcn(topControlsGrid, 'Text', 'Browse');
    handles.upButtons(k) = buttonFcn(topControlsGrid, 'Text', '↑', 'FontSize', 16);
    handles.downButtons(k) = buttonFcn(topControlsGrid, 'Text', '↓', 'FontSize', 16);

    % --- Host panel for the main scrolling content ---
    contentHostPanel = panelFcn(outerGrid, 'BorderType', 'none');
    contentHostPanel.Layout.Row = 2;

    % --- Main Content Grid (this is what will "scroll") ---
    mainContentGrid = uigridlayout(contentHostPanel, [9, 1]);
    handles.tabMainGrids(k) = mainContentGrid;

    % --- Panels to be scrolled ---
    % 1. Image Spacing
    spacingPanel = panelFcn(mainContentGrid, 'Title', 'Image Spacing');
    spacingPanel.Layout.Row = 1;
    spacingGrid = uigridlayout(spacingPanel, [1, 6]);
    labelFcn(spacingGrid, 'Text', 'XY Pixel:');
    handles.xySpacingInputs(k) = editFcn(spacingGrid, 'numeric', 'Value', lastUsed.xySpacing{k}, 'Tooltip', 'Width of a pixel in the XY plane, in micrometers (µm).');
    labelFcn(spacingGrid, 'Text', '(µm)');
    labelFcn(spacingGrid, 'Text', 'Z-Frame:');
    handles.zSpacingInputs(k) = editFcn(spacingGrid, 'numeric', 'Value', lastUsed.zSpacing{k}, 'Enable', 'off', 'Tooltip', 'Depth of a voxel in the Z direction (distance between slices), in micrometers (µm).');
    labelFcn(spacingGrid, 'Text', '(µm)');

    % 2. Deconvolution
    deconvPanel = panelFcn(mainContentGrid, 'Title', '');
    deconvPanel.Layout.Row = 2;
    handles.deconvPanels(k) = deconvPanel;
    deconvOuterGrid = uigridlayout(deconvPanel, [2, 1]);
    deconvOuterGrid.RowHeight = {'fit', '1x'};
    deconvOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(deconvOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Deconvolution', 'FontWeight', 'bold');
    handles.deconvEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.deconvEnabled{k});
    % Content grid - 6 rows
    deconvGrid = uigridlayout(deconvOuterGrid, [6, 4]);
    deconvGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    deconvGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
    deconvGrid.Padding = [5 5 5 5];
    deconvGrid.RowSpacing = 5;
    deconvGrid.ColumnSpacing = 10;
    
    % Row 1: Method dropdown (spans columns 2-4)
    labelFcn(deconvGrid, 'Text', 'Method:');
    handles.deconvMethodDrops(k) = dropdownFcn(deconvGrid, 'Items', {'Lucy-Richardson', 'Wiener', 'Blind'}, 'Value', lastUsed.deconvMethod{k}, 'Tooltip', 'Deconvolution algorithm.');
    handles.deconvMethodDrops(k).Layout.Row = 1;
    handles.deconvMethodDrops(k).Layout.Column = [2 4];
    
    % Row 2: Dynamic Parameter 1 (Iterations for LR/Blind, NSR for Wiener)
    handles.deconvParam1Labels(k) = labelFcn(deconvGrid, 'Text', 'Iterations:', 'Tag', ['deconv_param_' num2str(k)]);
    handles.deconvParam1Labels(k).Layout.Row = 2;
    handles.deconvParam1Labels(k).Layout.Column = 1;
    handles.deconvLRIterationsInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvLRIterations{k}, 'Tag', ['deconv_param_' num2str(k)], 'Tooltip', 'Number of iterations.');
    handles.deconvLRIterationsInputs(k).Layout.Row = 2;
    handles.deconvLRIterationsInputs(k).Layout.Column = 2;
    
    % Row 2 Col 3-4: Dynamic Parameter 2 (Damping for LR, Under-relax for Blind)
    handles.deconvParam2Labels(k) = labelFcn(deconvGrid, 'Text', 'Damping:', 'Tag', ['deconv_param_' num2str(k)]);
    handles.deconvParam2Labels(k).Layout.Row = 2;
    handles.deconvParam2Labels(k).Layout.Column = 3;
    handles.deconvLRDampingInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvLRDamping{k}, 'Tag', ['deconv_param_' num2str(k)], 'Tooltip', 'Damping/under-relaxation parameter.');
    handles.deconvLRDampingInputs(k).Layout.Row = 2;
    handles.deconvLRDampingInputs(k).Layout.Column = 4;
    
    % Shared controls for all methods (reused with different meanings)
    handles.deconvWienerNSRInputs(k) = handles.deconvLRIterationsInputs(k);
    handles.deconvBlindIterationsInputs(k) = handles.deconvLRIterationsInputs(k);
    handles.deconvBlindUnderRelaxInputs(k) = handles.deconvLRDampingInputs(k);
    
    % Row 3: PSF Source (Generate or Load File)
    labelFcn(deconvGrid, 'Text', 'PSF Source:');
    handles.deconvPSFSourceDrops(k) = dropdownFcn(deconvGrid, 'Items', {'Generate', 'Load File'}, 'Value', lastUsed.deconvPSFSource{k}, 'Tooltip', 'Generate Gaussian PSF or load experimental PSF image.');
    handles.deconvPSFSourceDrops(k).Layout.Row = 3;
    handles.deconvPSFSourceDrops(k).Layout.Column = [2 4];
    
    % Row 4: PSF File Path (for Load File option)
    labelFcn(deconvGrid, 'Text', 'PSF File:', 'Tag', ['deconv_psf_file_' num2str(k)], 'Visible', 'off');
    handles.deconvPSFPathTexts(k) = editFcn(deconvGrid, 'Value', lastUsed.deconvPSFFilePath{k}, 'Editable', 'off', 'Tag', ['deconv_psf_file_' num2str(k)], 'Visible', 'off');
    handles.deconvPSFPathTexts(k).Layout.Row = 4;
    handles.deconvPSFPathTexts(k).Layout.Column = 2;
    handles.deconvPSFBrowseButtons(k) = buttonFcn(deconvGrid, 'Text', 'Browse', 'Tag', ['deconv_psf_file_' num2str(k)], 'Visible', 'off');
    handles.deconvPSFBrowseButtons(k).Layout.Row = 4;
    handles.deconvPSFBrowseButtons(k).Layout.Column = [3 4];
    
    % Row 5: PSF Generation parameters (Sigma XY and Sigma Z)
    labelFcn(deconvGrid, 'Text', 'PSF Sigma XY:', 'Tag', ['deconv_psf_gen_' num2str(k)]);
    handles.deconvPSFSigmaXYInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvPSFSigmaXY{k}, 'Tag', ['deconv_psf_gen_' num2str(k)], 'Tooltip', 'PSF sigma in XY plane (pixels).');
    handles.deconvPSFSigmaXYInputs(k).Layout.Row = 5;
    handles.deconvPSFSigmaXYInputs(k).Layout.Column = 2;
    labelFcn(deconvGrid, 'Text', 'Sigma Z:', 'Tag', ['deconv_psf_gen_' num2str(k)]);
    handles.deconvPSFSigmaZInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvPSFSigmaZ{k}, 'Tag', ['deconv_psf_gen_' num2str(k)], 'Tooltip', 'PSF sigma in Z direction (slices).');
    handles.deconvPSFSigmaZInputs(k).Layout.Row = 5;
    handles.deconvPSFSigmaZInputs(k).Layout.Column = 4;
    
    % Row 6: PSF Size parameters (Size XY and Size Z)
    labelFcn(deconvGrid, 'Text', 'PSF Size XY:', 'Tag', ['deconv_psf_gen_' num2str(k)]);
    handles.deconvPSFSizeXYInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvPSFSizeXY{k}, 'Tag', ['deconv_psf_gen_' num2str(k)], 'Tooltip', 'PSF half-size in XY (pixels).');
    handles.deconvPSFSizeXYInputs(k).Layout.Row = 6;
    handles.deconvPSFSizeXYInputs(k).Layout.Column = 2;
    labelFcn(deconvGrid, 'Text', 'Size Z:', 'Tag', ['deconv_psf_gen_' num2str(k)]);
    handles.deconvPSFSizeZInputs(k) = editFcn(deconvGrid, 'numeric', 'Value', lastUsed.deconvPSFSizeZ{k}, 'Tag', ['deconv_psf_gen_' num2str(k)], 'Tooltip', 'PSF half-size in Z (slices).');
    handles.deconvPSFSizeZInputs(k).Layout.Row = 6;
    handles.deconvPSFSizeZInputs(k).Layout.Column = 4;

    % 3. Pre-processing
    preprocPanel = panelFcn(mainContentGrid, 'Title', '');
    preprocPanel.Layout.Row = 3;
    handles.preprocPanels(k) = preprocPanel;
    preprocOuterGrid = uigridlayout(preprocPanel, [2, 1]);
    preprocOuterGrid.RowHeight = {'fit', '1x'};
    preprocOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(preprocOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Pre-processing', 'FontWeight', 'bold');
    handles.preprocEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.preprocEnabled{k});
    % (contents of preprocPanel...)
    preprocGrid = uigridlayout(preprocOuterGrid, [6, 4]); % 6 rows, 4 columns for wavelet parameters
    preprocGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    preprocGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
    preprocGrid.Padding = [5 5 5 5];
    preprocGrid.RowSpacing = 5;
    preprocGrid.ColumnSpacing = 10;
    
    % Row 1: Clip checkbox (spans all columns)
    handles.preprocClipChecks(k) = checkboxFcn(preprocGrid, 'Text', 'Clip Negative Values at Zero', 'Value', lastUsed.preprocClipAtZero{k});
    handles.preprocClipChecks(k).Layout.Row = 1;
    handles.preprocClipChecks(k).Layout.Column = [1 4];
    
    % Row 2: Mode controls
    labelFcn(preprocGrid, 'Text', 'Mode:');
    handles.preprocessModeDrops(k) = dropdownFcn(preprocGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.preProcMode{k}, 'Tooltip', 'Choose how to apply pre-processing.');
    handles.preprocessScaleChecks(k) = checkboxFcn(preprocGrid, 'Text', 'Scale', 'Value', lastUsed.preProcScale{k}, 'Enable', 'off', 'Tooltip', 'Use physical units (µm) for 3D operations.');
    handles.preprocessProjectionDrops(k) = dropdownFcn(preprocGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.preProcProjection{k}, 'Tooltip', 'Z-projection type.');
    
    % Row 3: Method controls
    labelFcn(preprocGrid, 'Text', 'Method:');
    handles.preprocMethodDrops(k) = dropdownFcn(preprocGrid, 'Items', {'None', 'Gaussian', 'Median', 'Non-Local Means', 'Wavelet Denoising'}, 'Value', lastUsed.preProcMethod{k}, 'Tooltip', 'Algorithm for pre-processing.');
    handles.preprocMethodDrops(k).Layout.Column = [2 4];
    handles.preprocMethodDrops(k).Layout.Row = 3;
    
    % Row 4: Dynamic parameter controls (content changes based on selected method)
    % Primary parameter label and input (used by all methods)
    handles.preprocParam1Labels(k) = labelFcn(preprocGrid, 'Text', 'Parameter:', 'Tag', ['preproc_param_ctrl_' num2str(k)]);
    handles.preprocParam1Labels(k).Layout.Row = 4;
    handles.preprocParam1Labels(k).Layout.Column = 1;
    
    % Primary parameter input (Gaussian sigma, Median neighborhood, or Wavelet name)
    handles.preprocParam1Inputs(k) = editFcn(preprocGrid, 'numeric', 'Value', lastUsed.smoothGaussianValues{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Parameter value for the selected method.');
    handles.preprocParam1Inputs(k).Layout.Row = 4;
    handles.preprocParam1Inputs(k).Layout.Column = 2;
    
    % Units/secondary label
    handles.preprocParam2Labels(k) = labelFcn(preprocGrid, 'Text', '(units)', 'Tag', ['preproc_param_ctrl_' num2str(k)]);
    handles.preprocParam2Labels(k).Layout.Row = 4;
    handles.preprocParam2Labels(k).Layout.Column = 3;
    
    % Secondary parameter input (only used for wavelet level)
    handles.preprocParam2Inputs(k) = editFcn(preprocGrid, 'numeric', 'Value', lastUsed.waveletLevel{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Secondary parameter for the selected method.', 'Visible', 'off');
    handles.preprocParam2Inputs(k).Layout.Row = 4;
    handles.preprocParam2Inputs(k).Layout.Column = 4;
    
    % Row 5: Wavelet-specific parameters (Threshold Rule)
    handles.preprocParam3Labels(k) = labelFcn(preprocGrid, 'Text', 'Rule:', 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Visible', 'off');
    handles.preprocParam3Labels(k).Layout.Row = 5;
    handles.preprocParam3Labels(k).Layout.Column = 1;
    handles.preprocParam3Inputs(k) = dropdownFcn(preprocGrid, 'Items', {'sqtwolog','rigrsure','heursure','minimaxi'}, 'Value', lastUsed.waveletThresholdRule{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Threshold selection rule.', 'Visible', 'off');
    handles.preprocParam3Inputs(k).Layout.Row = 5;
    handles.preprocParam3Inputs(k).Layout.Column = 2;
    
    % Row 6: Wavelet threshold method
    handles.preprocParam4Labels(k) = labelFcn(preprocGrid, 'Text', 'Method:', 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Visible', 'off');
    handles.preprocParam4Labels(k).Layout.Row = 6;
    handles.preprocParam4Labels(k).Layout.Column = 1;
    handles.preprocParam4Inputs(k) = dropdownFcn(preprocGrid, 'Items', {'soft', 'hard'}, 'Value', lastUsed.waveletThresholdMethod{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Thresholding method.', 'Visible', 'off');
    handles.preprocParam4Inputs(k).Layout.Row = 6;
    handles.preprocParam4Inputs(k).Layout.Column = 2;
    
    % Store original parameter handles for backward compatibility
    handles.gaussInputs(k) = handles.preprocParam1Inputs(k);
    handles.medianInputs(k) = handles.preprocParam1Inputs(k);
    handles.waveletNameDrops(k) = dropdownFcn(preprocGrid, 'Items', {'haar','db2','db3','db4','sym2','sym3','sym4','coif1'}, 'Value', lastUsed.waveletName{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Type of wavelet.', 'Visible', 'off');
    handles.waveletNameDrops(k).Layout.Row = 4;
    handles.waveletNameDrops(k).Layout.Column = 2;
    handles.waveletLevelInputs(k) = handles.preprocParam2Inputs(k);
    handles.waveletThresholdRuleDrops(k) = handles.preprocParam3Inputs(k);
    handles.waveletThresholdMethodDrops(k) = handles.preprocParam4Inputs(k);
    
    % Non-Local Means specific parameter controls (hidden by default)
    handles.nlmFilterStrengthInputs(k) = editFcn(preprocGrid, 'numeric', 'Value', lastUsed.nlmFilterStrength{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Filter strength for non-local means denoising.', 'Visible', 'off');
    handles.nlmFilterStrengthInputs(k).Layout.Row = 4;
    handles.nlmFilterStrengthInputs(k).Layout.Column = 2;
    
    handles.nlmSearchWindowInputs(k) = editFcn(preprocGrid, 'numeric', 'Value', lastUsed.nlmSearchWindow{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Search window size for non-local means.', 'Visible', 'off');
    handles.nlmSearchWindowInputs(k).Layout.Row = 5;
    handles.nlmSearchWindowInputs(k).Layout.Column = 2;
    
    handles.nlmComparisonWindowInputs(k) = editFcn(preprocGrid, 'numeric', 'Value', lastUsed.nlmComparisonWindow{k}, 'Tag', ['preproc_param_ctrl_' num2str(k)], 'Tooltip', 'Comparison window size for non-local means.', 'Visible', 'off');
    handles.nlmComparisonWindowInputs(k).Layout.Row = 6;
    handles.nlmComparisonWindowInputs(k).Layout.Column = 2;
    

    % 4. Background Correction
    bgPanel = panelFcn(mainContentGrid, 'Title', '');
    bgPanel.Layout.Row = 4;
    handles.bgPanels(k) = bgPanel;
    bgOuterGrid = uigridlayout(bgPanel, [2, 1]);
    bgOuterGrid.RowHeight = {'fit', '1x'};
    bgOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(bgOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Background Correction', 'FontWeight', 'bold');
    handles.bgCorrEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.bgCorrEnabled{k});
    % (contents of bgPanel...)
    bgGrid = uigridlayout(bgOuterGrid, [4, 4]); % 4 rows now
    bgGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    bgGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
    bgGrid.Padding = [5 5 5 5];
    bgGrid.RowSpacing = 5;
    bgGrid.ColumnSpacing = 10;
    
    % Row 1: Clip checkbox (spans all columns)
    handles.bgCorrClipChecks(k) = checkboxFcn(bgGrid, 'Text', 'Clip Negative Values at Zero', 'Value', lastUsed.bgCorrClipAtZero{k});
    handles.bgCorrClipChecks(k).Layout.Row = 1;
    handles.bgCorrClipChecks(k).Layout.Column = [1 4];
    
    % Row 2: Mode controls
    labelFcn(bgGrid, 'Text', 'Mode:');
    handles.bgCorrModeDrops(k) = dropdownFcn(bgGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.bgCorrMode{k}, 'Tooltip', 'Choose how to apply background correction.');
    handles.bgCorrScaleChecks(k) = checkboxFcn(bgGrid, 'Text', 'Scale', 'Value', lastUsed.bgCorrScale{k}, 'Enable', 'off', 'Tooltip', 'If checked, 3D Gaussian method uses physical units (µm).');
    handles.bgCorrProjectionDrops(k) = dropdownFcn(bgGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.bgCorrProjection{k}, 'Tooltip', 'Type of Z-projection to perform if mode is On Z-Projection.');
    
    % Row 3: Method controls
    labelFcn(bgGrid, 'Text', 'Method:');
    handles.bgMethodDrops(k) = dropdownFcn(bgGrid, 'Items', {'None', 'Gaussian', 'Rolling-ball', 'Top-hat'}, 'Value', lastUsed.bgMethod{k}, 'Tooltip', 'Algorithm for background subtraction.');
    handles.bgMethodDrops(k).Layout.Column = [2 4];
    handles.bgMethodDrops(k).Layout.Row = 3;
    
    % Row 4: Parameter controls
    labelFcn(bgGrid, 'Text', 'Parameter:');
    handles.bgParamInputs(k) = editFcn(bgGrid, 'numeric', 'Value', lastUsed.bgParam{k}, 'Tooltip', 'Gaussian: Sigma. Rolling-ball/Top-hat: Radius. Units are µm if "Scale" is checked (3D), otherwise pixels/voxels.');
    handles.bgParamInputs(k).Layout.Row = 4;
    handles.bgParamInputs(k).Layout.Column = 2;
    unitsLabel = labelFcn(bgGrid, 'Text', '(units)');
    unitsLabel.Layout.Row = 4;
    unitsLabel.Layout.Column = 3;

    % 5. Local Maxima Detection
    maximaPanel = panelFcn(mainContentGrid, 'Title', '');
    maximaPanel.Layout.Row = 5;
    handles.maximaPanels(k) = maximaPanel;
    maximaOuterGrid = uigridlayout(maximaPanel, [2, 1]);
    maximaOuterGrid.RowHeight = {'fit', '1x'};
    maximaOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(maximaOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Local Maxima Detection', 'FontWeight', 'bold');
    handles.maximaEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.maximaEnabled{k});
    % (contents of maximaPanel...)
    maximaGrid = uigridlayout(maximaOuterGrid, [7, 4]);
    maximaGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    maximaGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};

    % --- Row 1: Mode Controls ---
    labelFcn(maximaGrid, 'Text', 'Mode:');
    handles.maximaModeDrops(k) = dropdownFcn(maximaGrid, 'Items', {'3D', '2D (Slice-by-slice)', 'On Z-Projection'}, 'Value', lastUsed.maximaMode{k}, 'Tooltip', 'Choose where to find local maxima.');
    handles.maximaScaleChecks(k) = checkboxFcn(maximaGrid, 'Text', 'Scale', 'Value', lastUsed.maximaScale{k}, 'Tooltip', 'Use physical units (µm) for neighborhood size.');
    handles.maximaProjectionDrops(k) = dropdownFcn(maximaGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.maximaProjection{k}, 'Tooltip', 'Type of Z-projection to perform if mode is On Z-Projection.');
    handles.maximaProjectionDrops(k).Layout.Column = 4; % Last column

    % --- Row 2: Method Controls ---
    labelFcn(maximaGrid, 'Text', 'Method:');
    handles.maximaMethodDrops(k) = dropdownFcn(maximaGrid, 'Items', {'Simple Regional', 'Extended Maxima', 'Laplacian of Gaussian'}, 'Value', lastUsed.maximaMethod{k}, 'Tooltip', 'Algorithm for detecting local maxima.');
    handles.maximaMethodDrops(k).Layout.Column = [2 4];
    handles.maximaMethodDrops(k).Layout.Row = 2;
    
    % --- Row 3: Neighborhood Controls ---
    labelFcn(maximaGrid, 'Text', 'Neighborhood:');
    handles.maximaNeighborhoodInputs(k) = editFcn(maximaGrid, 'numeric', 'Value', lastUsed.maximaNeighborhoodSize{k}, 'Tooltip', 'Neighborhood size in pixels/voxels (or microns if scaled).');
    handles.maximaNeighborhoodInputs(k).Layout.Column = 2;
    handles.maximaNeighborhoodInputs(k).Layout.Row = 3;
    labelFcn(maximaGrid, 'Text', '(units)');

    handles.hMaxPanel(k) = panelFcn(maximaGrid, 'BorderType', 'none');
    handles.hMaxPanel(k).Layout.Row = 4;
    handles.hMaxPanel(k).Layout.Column = [1 4];
    hMaxGrid = uigridlayout(handles.hMaxPanel(k), [1,3]);
    hMaxGrid.ColumnWidth = {'fit', '1x', 'fit'};
    labelFcn(hMaxGrid, 'Text', 'H-Max Value:', 'Tag', 'hmax_control');
    handles.hMaxInputs(k) = editFcn(hMaxGrid, 'numeric', 'Value', lastUsed.hMaxValue{k}, 'Tooltip', 'Minimum intensity drop between a peak and its lower neighbor.', 'Tag', 'hmax_control');
    labelFcn(hMaxGrid, 'Text', '(intensity)', 'Tag', 'hmax_control');
    handles.logPanel(k) = panelFcn(maximaGrid, 'BorderType', 'none', 'Visible', 'on');
    handles.logPanel(k).Layout.Row = 5;
    handles.logPanel(k).Layout.Column = [1 4];
    logGrid = uigridlayout(handles.logPanel(k), [1, 6]);
    labelFcn(logGrid, 'Text', 'Sigma:', 'Tag', 'log_control');
    handles.logSigmaInputs(k) = editFcn(logGrid, 'numeric', 'Value', lastUsed.sigmaValue{k}, 'Tooltip', 'Sigma for LoG filter. Units are µm if "Scale" is checked (3D), otherwise pixels.', 'Tag', 'log_control');
    labelFcn(logGrid, 'Text', '(units)', 'Tag', 'log_control');
    labelFcn(logGrid, 'Text', 'Threshold:', 'Tag', 'log_control');
    handles.logThresholdInputs(k) = editFcn(logGrid, 'numeric', 'Value', lastUsed.peakThresholdValue{k}, 'Tooltip', 'Minimum peak intensity after LoG filtering to be considered a maximum.', 'Tag', 'log_control');
    labelFcn(logGrid, 'Text', 'intensity', 'Tag', 'log_control');
    handles.showMaximaChecks(k) = checkboxFcn(maximaGrid, 'Text', 'Show Maxima on Previews', 'Value', lastUsed.showMaxima{k}, 'Tooltip', 'If checked, detected maxima will be overlaid on the preview images.');
    handles.showMaximaChecks(k).Layout.Row = 6;
    handles.showMaximaChecks(k).Layout.Column = [1 2];
    labelFcn(maximaGrid, 'Text', 'Maxima Color:');
    handles.maximaColorDrops(k) = dropdownFcn(maximaGrid, 'Items', {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan'}, 'Value', lastUsed.maximaColor{k}, 'Tooltip', 'Color used to display maxima overlays for this channel.');
    handles.maximaColorDrops(k).Layout.Row = 6;
    handles.maximaColorDrops(k).Layout.Column = [3 4];
    
    % --- Row 7: Display on All Previews ---
    handles.displayOnAllPreviewsChecks(k) = checkboxFcn(maximaGrid, 'Text', 'Display on All Previews', 'Value', lastUsed.displayOnAllPreviews{k}, 'Tooltip', 'If checked, this channel''s maxima will be displayed on all preview panels, not just its own.');
    handles.displayOnAllPreviewsChecks(k).Layout.Row = 7;
    handles.displayOnAllPreviewsChecks(k).Layout.Column = [1 4];

    % 6. Local Maxima Fitting
    gaussFitPanel = panelFcn(mainContentGrid, 'Title', '');
    gaussFitPanel.Layout.Row = 6;
    handles.gaussFitPanels(k) = gaussFitPanel;
    gaussFitOuterGrid = uigridlayout(gaussFitPanel, [2, 1]);
    gaussFitOuterGrid.RowHeight = {'fit', '1x'};
    gaussFitOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(gaussFitOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Local Maxima Fitting', 'FontWeight', 'bold');
    handles.gaussFitEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', lastUsed.gaussFitEnabled{k});

    gaussFitGrid = uigridlayout(gaussFitOuterGrid);
    gaussFitGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
    gaussFitGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};

    % Voxel Window Size
    labelFcn(gaussFitGrid, 'Text', 'Voxel Window Size:', 'Tooltip', 'Width of the cubic window for Gaussian fitting. Must be an odd number.');
    handles.gaussFitVoxelWindowSlider(k) = uislider(gaussFitGrid, 'Limits', [3, 21], 'MajorTicks', 3:2:21, 'Value', lastUsed.gaussFitVoxelWindowSize{k}, 'MajorTickLabels', arrayfun(@num2str, 3:2:21, 'UniformOutput', false));
    handles.gaussFitVoxelWindowSlider(k).Layout.Column = [2 4];

    % Local Background Correction Method
    labelFcn(gaussFitGrid, 'Text', 'BG Correction:', 'Tooltip', 'Method for local background estimation around each maxima.');
    handles.gaussFitBgCorrMethodDrop(k) = dropdownFcn(gaussFitGrid, 'Items', {'Mean Surrounding Subtraction', 'Local Plane Fitting', 'Local Polynomial Fitting'}, 'Value', lastUsed.gaussFitBgCorrMethod{k});
    handles.gaussFitBgCorrMethodDrop(k).Layout.Column = [2 4];

    % Background Correction Parameters
    labelFcn(gaussFitGrid, 'Text', 'BG Width:', 'Tooltip', 'Width of the surrounding region for background calculation (in voxels).');
    handles.gaussFitBgCorrWidthEdit(k) = uieditfield(gaussFitGrid, 'numeric', 'Value', lastUsed.gaussFitBgCorrWidth{k});
    labelFcn(gaussFitGrid, 'Text', 'Poly Degree:', 'Tag', ['poly_degree_label_' num2str(k)], 'Tooltip', 'Degree of the polynomial for fitting (up to 3).');
    handles.gaussFitPolyDegreeEdit(k) = uieditfield(gaussFitGrid, 'numeric', 'Value', lastUsed.gaussFitPolyDegree{k}, 'Tag', ['poly_degree_edit_' num2str(k)]);

    % Local Maxima Fitting Method
    labelFcn(gaussFitGrid, 'Text', 'Fit Method:', 'Tooltip', 'Mathematical model to fit to detected local maxima.');
    fitMethods = {'1D (X,Y,Z) Gaussian', '2D (XY) + 1D (Z) Gaussian', '3D Gaussian', 'Distorted 3D Gaussian', 'Radial Symmetry'};
    handles.gaussFitMethodDrop(k) = dropdownFcn(gaussFitGrid, 'Items', fitMethods, 'Value', lastUsed.gaussFitMethod{k});
    handles.gaussFitMethodDrop(k).Layout.Column = [2 4];

    % Solver Parameters (Row 5)
    labelFcn(gaussFitGrid, 'Text', 'Max Iterations:', 'Tooltip', 'Maximum number of iterations for the fitting algorithm.');
    handles.gaussFitMaxIterationsEdit(k) = uieditfield(gaussFitGrid, 'numeric', 'Value', lastUsed.gaussFitMaxIterations{k});
    handles.gaussFitMaxIterationsEdit(k).Layout.Row = 5;
    handles.gaussFitMaxIterationsEdit(k).Layout.Column = 2;
    labelFcn(gaussFitGrid, 'Text', 'Tolerance:', 'Tooltip', 'Termination tolerance for the fitting algorithm.');
    handles.gaussFitToleranceEdit(k) = uieditfield(gaussFitGrid, 'numeric', 'Value', lastUsed.gaussFitTolerance{k});
    handles.gaussFitToleranceEdit(k).Layout.Row = 5;
    handles.gaussFitToleranceEdit(k).Layout.Column = 4;
    
    % Radial Symmetry Radius (Row 6)
    labelFcn(gaussFitGrid, 'Text', 'Radial Radius:', 'Tooltip', 'Search radius for radial symmetry method (in pixels).', 'Tag', 'radial_radius_label');
    handles.gaussFitRadialRadiusEdit(k) = uieditfield(gaussFitGrid, 'numeric', 'Value', lastUsed.gaussFitRadialRadius{k});
    handles.gaussFitRadialRadiusEdit(k).Layout.Row = 6;
    handles.gaussFitRadialRadiusEdit(k).Layout.Column = 2;
    handles.gaussFitRadialRadiusEdit(k).Tag = ['radial_radius_edit_' num2str(k)];

    % 7. Fit Filtering
    fitFilterPanel = panelFcn(mainContentGrid, 'Title', '');
    fitFilterPanel.Layout.Row = 7;
    handles.fitFilterPanels(k) = fitFilterPanel;
    fitFilterOuterGrid = uigridlayout(fitFilterPanel, [2, 1]);
    fitFilterOuterGrid.RowHeight = {'fit', '1x'};
    fitFilterOuterGrid.Padding = [2 2 2 2];
    % Title-like row
    titleGrid = uigridlayout(fitFilterOuterGrid, [1, 2]);
    titleGrid.ColumnWidth = {'1x', 'fit'};
    titleGrid.Padding = [0 0 0 0];
    labelFcn(titleGrid, 'Text', 'Fit Filtering', 'FontWeight', 'bold');
    % Safely get fitFilterEnabled value, defaulting to false if not present
    fitFilterEnabledValue = false;
    if isfield(lastUsed, 'fitFilterEnabled') && k <= length(lastUsed.fitFilterEnabled)
        fitFilterEnabledValue = lastUsed.fitFilterEnabled{k};
    end
    handles.fitFilterEnabledChecks(k) = checkboxFcn(titleGrid, 'Text', '', 'Value', fitFilterEnabledValue);
    % Content grid - 6 columns for checkbox, min input, max input, label, spacer, spacer
    fitFilterGrid = uigridlayout(fitFilterOuterGrid, [4, 6]);
    fitFilterGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
    fitFilterGrid.ColumnWidth = {'fit', '1x', '1x', 'fit', 'fit', 'fit'};
    fitFilterGrid.Padding = [5 5 5 5];
    fitFilterGrid.RowSpacing = 5;
    fitFilterGrid.ColumnSpacing = 5;

    % Row 1: R² value filtering
    % Safely get fitFilterRSquaredEnabled value, defaulting to false if not present
    fitFilterRSquaredEnabledValue = false;
    if isfield(lastUsed, 'fitFilterRSquaredEnabled') && k <= length(lastUsed.fitFilterRSquaredEnabled)
        fitFilterRSquaredEnabledValue = lastUsed.fitFilterRSquaredEnabled{k};
    end
    handles.fitFilterRSquaredEnabledChecks(k) = checkboxFcn(fitFilterGrid, 'Text', 'R²/Quality Score:', 'Value', fitFilterRSquaredEnabledValue, 'Tooltip', 'Filter fits by R² (Gaussian) or Quality Score (Radial Symmetry)');
    handles.fitFilterRSquaredEnabledChecks(k).Layout.Row = 1;
    handles.fitFilterRSquaredEnabledChecks(k).Layout.Column = 1;
    % Safely get fitFilterRSquaredMin value, defaulting to 0.8 if not present
    fitFilterRSquaredMinValue = 0.8;
    if isfield(lastUsed, 'fitFilterRSquaredMin') && k <= length(lastUsed.fitFilterRSquaredMin)
        fitFilterRSquaredMinValue = lastUsed.fitFilterRSquaredMin{k};
    end
    handles.fitFilterRSquaredMinInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterRSquaredMinValue, 'Tooltip', 'Minimum R-squared value (0.0-1.0, higher = better fit)');
    handles.fitFilterRSquaredMinInputs(k).Layout.Row = 1;
    handles.fitFilterRSquaredMinInputs(k).Layout.Column = 2;
    % Safely get fitFilterRSquaredMax value, defaulting to 1.0 if not present
    fitFilterRSquaredMaxValue = 1.0;
    if isfield(lastUsed, 'fitFilterRSquaredMax') && k <= length(lastUsed.fitFilterRSquaredMax)
        fitFilterRSquaredMaxValue = lastUsed.fitFilterRSquaredMax{k};
    end
    handles.fitFilterRSquaredMaxInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterRSquaredMaxValue, 'Tooltip', 'Maximum R-squared value (0.0-1.0, higher = better fit)');
    handles.fitFilterRSquaredMaxInputs(k).Layout.Row = 1;
    handles.fitFilterRSquaredMaxInputs(k).Layout.Column = 3;
    tempLabel1 = labelFcn(fitFilterGrid, 'Text', '(0.0-1.0)', 'Tooltip', 'R-squared range: 0.0 (poor fit) to 1.0 (perfect fit)');
    tempLabel1.Layout.Row = 1;
    tempLabel1.Layout.Column = 4;
    tempLabel2 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel2.Layout.Row = 1;
    tempLabel2.Layout.Column = 5;
    tempLabel3 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel3.Layout.Row = 1;
    tempLabel3.Layout.Column = 6;

    % Row 2: Sigma sum filtering
    % Safely get fitFilterSigmaSumEnabled value, defaulting to false if not present
    fitFilterSigmaSumEnabledValue = false;
    if isfield(lastUsed, 'fitFilterSigmaSumEnabled') && k <= length(lastUsed.fitFilterSigmaSumEnabled)
        fitFilterSigmaSumEnabledValue = lastUsed.fitFilterSigmaSumEnabled{k};
    end
    handles.fitFilterSigmaSumEnabledChecks(k) = checkboxFcn(fitFilterGrid, 'Text', 'Sigma Sum:', 'Value', fitFilterSigmaSumEnabledValue, 'Tooltip', 'Filter fits by sum of sigma values range (size constraint)');
    handles.fitFilterSigmaSumEnabledChecks(k).Layout.Row = 2;
    handles.fitFilterSigmaSumEnabledChecks(k).Layout.Column = 1;
    % Safely get fitFilterSigmaSumMin value, defaulting to 0.0 if not present
    fitFilterSigmaSumMinValue = 0.0;
    if isfield(lastUsed, 'fitFilterSigmaSumMin') && k <= length(lastUsed.fitFilterSigmaSumMin)
        fitFilterSigmaSumMinValue = lastUsed.fitFilterSigmaSumMin{k};
    end
    handles.fitFilterSigmaSumMinInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterSigmaSumMinValue, 'Tooltip', 'Minimum sum of sigma values (pixels or µm)');
    handles.fitFilterSigmaSumMinInputs(k).Layout.Row = 2;
    handles.fitFilterSigmaSumMinInputs(k).Layout.Column = 2;
    % Safely get fitFilterSigmaSumMax value, defaulting to 10.0 if not present
    fitFilterSigmaSumMaxValue = 10.0;
    if isfield(lastUsed, 'fitFilterSigmaSumMax') && k <= length(lastUsed.fitFilterSigmaSumMax)
        fitFilterSigmaSumMaxValue = lastUsed.fitFilterSigmaSumMax{k};
    end
    handles.fitFilterSigmaSumMaxInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterSigmaSumMaxValue, 'Tooltip', 'Maximum sum of sigma values (pixels or µm)');
    handles.fitFilterSigmaSumMaxInputs(k).Layout.Row = 2;
    handles.fitFilterSigmaSumMaxInputs(k).Layout.Column = 3;
    tempLabel4 = labelFcn(fitFilterGrid, 'Text', '(units)', 'Tooltip', 'Units depend on image spacing settings');
    tempLabel4.Layout.Row = 2;
    tempLabel4.Layout.Column = 4;
    tempLabel5 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel5.Layout.Row = 2;
    tempLabel5.Layout.Column = 5;
    tempLabel6 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel6.Layout.Row = 2;
    tempLabel6.Layout.Column = 6;

    % Row 3: Amplitude filtering
    % Safely get fitFilterAmplitudeEnabled value, defaulting to false if not present
    fitFilterAmplitudeEnabledValue = false;
    if isfield(lastUsed, 'fitFilterAmplitudeEnabled') && k <= length(lastUsed.fitFilterAmplitudeEnabled)
        fitFilterAmplitudeEnabledValue = lastUsed.fitFilterAmplitudeEnabled{k};
    end
    handles.fitFilterAmplitudeEnabledChecks(k) = checkboxFcn(fitFilterGrid, 'Text', 'Amplitude:', 'Value', fitFilterAmplitudeEnabledValue, 'Tooltip', 'Filter fits by amplitude range (brightness threshold)');
    handles.fitFilterAmplitudeEnabledChecks(k).Layout.Row = 3;
    handles.fitFilterAmplitudeEnabledChecks(k).Layout.Column = 1;
    % Safely get fitFilterAmplitudeMin value, defaulting to 100.0 if not present
    fitFilterAmplitudeMinValue = 100.0;
    if isfield(lastUsed, 'fitFilterAmplitudeMin') && k <= length(lastUsed.fitFilterAmplitudeMin)
        fitFilterAmplitudeMinValue = lastUsed.fitFilterAmplitudeMin{k};
    end
    handles.fitFilterAmplitudeMinInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterAmplitudeMinValue, 'Tooltip', 'Minimum amplitude value (intensity units)');
    handles.fitFilterAmplitudeMinInputs(k).Layout.Row = 3;
    handles.fitFilterAmplitudeMinInputs(k).Layout.Column = 2;
    % Safely get fitFilterAmplitudeMax value, defaulting to 10000.0 if not present
    fitFilterAmplitudeMaxValue = 10000.0;
    if isfield(lastUsed, 'fitFilterAmplitudeMax') && k <= length(lastUsed.fitFilterAmplitudeMax)
        fitFilterAmplitudeMaxValue = lastUsed.fitFilterAmplitudeMax{k};
    end
    handles.fitFilterAmplitudeMaxInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterAmplitudeMaxValue, 'Tooltip', 'Maximum amplitude value (intensity units)');
    handles.fitFilterAmplitudeMaxInputs(k).Layout.Row = 3;
    handles.fitFilterAmplitudeMaxInputs(k).Layout.Column = 3;
    tempLabel7 = labelFcn(fitFilterGrid, 'Text', '(intensity)', 'Tooltip', 'Intensity units');
    tempLabel7.Layout.Row = 3;
    tempLabel7.Layout.Column = 4;
    tempLabel8 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel8.Layout.Row = 3;
    tempLabel8.Layout.Column = 5;
    tempLabel9 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel9.Layout.Row = 3;
    tempLabel9.Layout.Column = 6;

    % Row 4: Intensity filtering
    % Safely get fitFilterIntensityEnabled value, defaulting to false if not present
    fitFilterIntensityEnabledValue = false;
    if isfield(lastUsed, 'fitFilterIntensityEnabled') && k <= length(lastUsed.fitFilterIntensityEnabled)
        fitFilterIntensityEnabledValue = lastUsed.fitFilterIntensityEnabled{k};
    end
    handles.fitFilterIntensityEnabledChecks(k) = checkboxFcn(fitFilterGrid, 'Text', 'Intensity:', 'Value', fitFilterIntensityEnabledValue, 'Tooltip', 'Filter fits by integrated intensity range (total signal)');
    handles.fitFilterIntensityEnabledChecks(k).Layout.Row = 4;
    handles.fitFilterIntensityEnabledChecks(k).Layout.Column = 1;
    % Safely get fitFilterIntensityMin value, defaulting to 1000.0 if not present
    fitFilterIntensityMinValue = 1000.0;
    if isfield(lastUsed, 'fitFilterIntensityMin') && k <= length(lastUsed.fitFilterIntensityMin)
        fitFilterIntensityMinValue = lastUsed.fitFilterIntensityMin{k};
    end
    handles.fitFilterIntensityMinInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterIntensityMinValue, 'Tooltip', 'Minimum integrated intensity value');
    handles.fitFilterIntensityMinInputs(k).Layout.Row = 4;
    handles.fitFilterIntensityMinInputs(k).Layout.Column = 2;
    % Safely get fitFilterIntensityMax value, defaulting to 100000.0 if not present
    fitFilterIntensityMaxValue = 100000.0;
    if isfield(lastUsed, 'fitFilterIntensityMax') && k <= length(lastUsed.fitFilterIntensityMax)
        fitFilterIntensityMaxValue = lastUsed.fitFilterIntensityMax{k};
    end
    handles.fitFilterIntensityMaxInputs(k) = editFcn(fitFilterGrid, 'numeric', 'Value', fitFilterIntensityMaxValue, 'Tooltip', 'Maximum integrated intensity value');
    handles.fitFilterIntensityMaxInputs(k).Layout.Row = 4;
    handles.fitFilterIntensityMaxInputs(k).Layout.Column = 3;
    tempLabel10 = labelFcn(fitFilterGrid, 'Text', '(intensity)', 'Tooltip', 'Integrated intensity units');
    tempLabel10.Layout.Row = 4;
    tempLabel10.Layout.Column = 4;
    tempLabel11 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel11.Layout.Row = 4;
    tempLabel11.Layout.Column = 5;
    tempLabel12 = labelFcn(fitFilterGrid, 'Text', '');
    tempLabel12.Layout.Row = 4;
    tempLabel12.Layout.Column = 6;

    % 8. Classification Panel (SVM)
    classifyPanel = panelFcn(mainContentGrid, 'Title', '');
    classifyPanel.Layout.Row = 8;
    handles.classifyPanels(k) = classifyPanel;
    classifyOuterGrid = uigridlayout(classifyPanel, [2, 1]);
    classifyOuterGrid.RowHeight = {'fit', '1x'};
    classifyOuterGrid.Padding = [2 2 2 2];
    % Title row
    classifyTitleGrid = uigridlayout(classifyOuterGrid, [1, 2]);
    classifyTitleGrid.ColumnWidth = {'1x', 'fit'};
    classifyTitleGrid.Padding = [0 0 0 0];
    labelFcn(classifyTitleGrid, 'Text', 'Signal Classification (SVM)', 'FontWeight', 'bold');
    handles.classifyEnabledChecks(k) = checkboxFcn(classifyTitleGrid, 'Text', '', 'Value', false);
    % Content grid
    classifyGrid = uigridlayout(classifyOuterGrid, [6, 3]);
    classifyGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    classifyGrid.ColumnWidth = {'fit', '1x', 'fit'};
    classifyGrid.Padding = [5 5 5 5];
    classifyGrid.RowSpacing = 5;
    % Row 1: Status
    statusLbl = labelFcn(classifyGrid, 'Text', 'Status:');
    statusLbl.Layout.Row = 1; statusLbl.Layout.Column = 1;
    handles.classifyStatusLabels(k) = labelFcn(classifyGrid, 'Text', 'No classifier loaded', 'FontColor', [0.5 0.5 0.5]);
    handles.classifyStatusLabels(k).Layout.Row = 1; handles.classifyStatusLabels(k).Layout.Column = [2 3];
    % Row 2: Features
    featLbl = labelFcn(classifyGrid, 'Text', 'Features:');
    featLbl.Layout.Row = 2; featLbl.Layout.Column = 1;
    handles.classifyFeatureCountLabels(k) = labelFcn(classifyGrid, 'Text', '0 features', 'FontColor', [0.5 0.5 0.5]);
    handles.classifyFeatureCountLabels(k).Layout.Row = 2; handles.classifyFeatureCountLabels(k).Layout.Column = 2;
    handles.classifySelectFeaturesButtons(k) = buttonFcn(classifyGrid, 'Text', 'Select...', 'Enable', 'off');
    handles.classifySelectFeaturesButtons(k).Layout.Row = 2; handles.classifySelectFeaturesButtons(k).Layout.Column = 3;
    % Row 3: Load classifier
    loadLbl = labelFcn(classifyGrid, 'Text', 'Classifier:');
    loadLbl.Layout.Row = 3; loadLbl.Layout.Column = 1;
    handles.classifyLoadButtons(k) = buttonFcn(classifyGrid, 'Text', 'Load Classifier...', 'Tooltip', 'Load a trained classifier from SNAP_classify');
    handles.classifyLoadButtons(k).Layout.Row = 3; handles.classifyLoadButtons(k).Layout.Column = [2 3];
    % Row 4: Open trainer
    trainLbl = labelFcn(classifyGrid, 'Text', 'Training:');
    trainLbl.Layout.Row = 4; trainLbl.Layout.Column = 1;
    handles.classifyTrainButtons(k) = buttonFcn(classifyGrid, 'Text', 'Open SNAP_classify...', 'Tooltip', 'Open the standalone classifier training application');
    handles.classifyTrainButtons(k).Layout.Row = 4; handles.classifyTrainButtons(k).Layout.Column = [2 3];
    % Row 5: Apply button
    handles.classifyApplyButtons(k) = buttonFcn(classifyGrid, 'Text', 'Apply to Current Results', 'Enable', 'off');
    handles.classifyApplyButtons(k).Layout.Row = 5; handles.classifyApplyButtons(k).Layout.Column = [1 3];
    % Row 6: Filter noise option
    handles.classifyFilterNoiseChecks(k) = checkboxFcn(classifyGrid, 'Text', 'Auto-filter noise predictions', 'Value', true);
    handles.classifyFilterNoiseChecks(k).Layout.Row = 6; handles.classifyFilterNoiseChecks(k).Layout.Column = [1 3];

    % 9. Spacer Panel (to fill remaining vertical space)
    spacerPanel = panelFcn(mainContentGrid, 'BorderType', 'none');
    spacerPanel.Layout.Row = 9;

    % Set initial "scroll" state
    initialState = lastUsed.navPanelIndex{k};
    % Set initial navigation state for channels
    % Grid rows: 1=Spacing, 2=Deconv, 3=Preproc, 4=BG, 5=Maxima, 6=Gaussian Fit, 7=Fit Filter, 8=Classification, 9=Spacer
    if initialState == 0 % Stage 0 - show panels 1-5 (Spacing, Deconv, Preproc, BG, Maxima)
        mainContentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 0, 0, 0, '1x'};
    elseif initialState == 1 % Stage 1 - show panels 2-6 (Deconv, Preproc, BG, Maxima, Gaussian Fit)
        mainContentGrid.RowHeight = {0, 'fit', 'fit', 'fit', 'fit', 'fit', 0, 0, '1x'};
    elseif initialState == 2 % Stage 2 - show panels 3-7 (Preproc, BG, Maxima, Gaussian Fit, Fit Filter)
        mainContentGrid.RowHeight = {0, 0, 'fit', 'fit', 'fit', 'fit', 'fit', 0, '1x'};
    elseif initialState == 3 % Stage 3 - show panels 4-8 (BG, Maxima, Gaussian Fit, Fit Filter, Classification)
        mainContentGrid.RowHeight = {0, 0, 0, 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
    elseif initialState == 4 % Stage 4 - show panels 6-8 (Gaussian Fit, Fit Filter, Classification)
        mainContentGrid.RowHeight = {0, 0, 0, 0, 0, 'fit', 'fit', 'fit', '1x'};
    elseif initialState == 5 % Stage 5 - show panels 7-8 (Fit Filter, Classification) - maximum Classification view
        mainContentGrid.RowHeight = {0, 0, 0, 0, 0, 0, 'fit', 'fit', '1x'};
    else % Default fallback - Stage 0
        mainContentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 0, 0, 0, '1x'};
    end
end

% --- Right Panel for Previews ---
rightPanel = panelFcn(mainGrid, 'Title', 'Previews');
rightPanel.Layout.Row = 1;
rightPanel.Layout.Column = 2;

rightGrid = uigridlayout(rightPanel, [2, 3]);
rightGrid.Padding = [10 10 10 10];
rightGrid.RowSpacing = 10;
rightGrid.ColumnSpacing = 10;

% --- Pre-allocate Handle Arrays for 5 Previews ---
numPreviews = 5;
handles.previewPanels = gobjects(1, numPreviews);
handles.previewAxes = gobjects(1, numPreviews);
handles.previewContentDrops = gobjects(1, numPreviews);
handles.previewModeDrops = gobjects(1, numPreviews);
handles.previewProjectionDrops = gobjects(1, numPreviews);
handles.zSliders = gobjects(1, numPreviews);
handles.zLabels = gobjects(1, numPreviews);
handles.brightnessSliders = gobjects(1, numPreviews);
handles.brightnessResets = gobjects(1, numPreviews);

for i = 1:numPreviews
    p = panelFcn(rightGrid, 'Title', ['Preview ' num2str(i)]);
    handles.previewPanels(i) = p;
    
    previewGrid = uigridlayout(p, [4, 1]);
    previewGrid.RowHeight = {'fit', 'fit', '1x', 'fit'};
    
    % Top controls
    topCtrlGrid = uigridlayout(previewGrid, [1, 4]);
    topCtrlGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
    labelFcn(topCtrlGrid, 'Text', 'Content:');
    handles.previewContentDrops(i) = dropdownFcn(topCtrlGrid, 'Items', {'None'}, 'Value', 'None');
    labelFcn(topCtrlGrid, 'Text', 'Mode:');
    handles.previewModeDrops(i) = dropdownFcn(topCtrlGrid, 'Items', {'Z-Projection', 'Z-Stack'}, 'Value', lastUsed.previewModes{i});

    % Z-Controls
    zCtrlGrid = uigridlayout(previewGrid, [1,3]);
    zCtrlGrid.ColumnWidth = {'fit', '1x', 'fit'};
    handles.previewProjectionDrops(i) = dropdownFcn(zCtrlGrid, 'Items', {'Max', 'Min', 'Mean', 'Median'}, 'Value', lastUsed.previewProjections{i});
    handles.zSliders(i) = sliderFcn(zCtrlGrid, 'Limits', [1,100], 'Value', 1);
    handles.zLabels(i) = labelFcn(zCtrlGrid, 'Text', 'Z: 1/100');

    % Axes
    handles.previewAxes(i) = axesFcn(previewGrid);
    handles.previewAxes(i).XTick = []; handles.previewAxes(i).YTick = [];

    % Brightness controls
    brightCtrlGrid = uigridlayout(previewGrid, [1,3]);
    brightCtrlGrid.ColumnWidth = {'fit', '1x', 'fit'};
    labelFcn(brightCtrlGrid, 'Text', 'Brightness:');
    handles.brightnessSliders(i) = sliderFcn(brightCtrlGrid, 'Limits', [0.1, 5], 'Value', 1);
    handles.brightnessResets(i) = buttonFcn(brightCtrlGrid, 'Text', 'Reset');
end

% --- Actions & Summary Panel ---
actionsPanel = panelFcn(rightGrid, 'Title', 'Actions & Summary');
actionsGrid = uigridlayout(actionsPanel, [5, 1]);
actionsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};

% Global Z Slider - Compact Layout
zPanel = panelFcn(actionsGrid, 'Title', 'Global Z Control');
zGrid = uigridlayout(zPanel, [3,1]);
zGrid.RowHeight = {'fit', 'fit', 'fit'};
zGrid.Padding = [5 5 5 5];
zGrid.RowSpacing = 3;

handles.globalZLabel = labelFcn(zGrid, 'Text', 'Global Z: 1/1');
handles.globalZSlider = sliderFcn(zGrid, 'Limits', [1,100], 'Value', 1, 'Enable', 'off');

% Play/Pause Controls - Compact
playPauseGrid = uigridlayout(zGrid, [1,2]);
playPauseGrid.ColumnWidth = {'1x', '1x'};
playPauseGrid.RowHeight = {'fit'};
playPauseGrid.Padding = [0 0 0 0];
playPauseGrid.ColumnSpacing = 3;

handles.playButton = buttonFcn(playPauseGrid, 'Text', '▶ Play', 'FontWeight', 'bold');
handles.pauseButton = buttonFcn(playPauseGrid, 'Text', '⏸ Pause', 'FontWeight', 'bold');
handles.playButton.Enable = 'off'; % Initially disabled until z-stack data is loaded
handles.pauseButton.Enable = 'off'; % Initially disabled until play starts

% Combined Status and Maxima Counts Panel
statusMaximaPanel = panelFcn(actionsGrid, 'BorderType', 'none');
statusMaximaGrid = uigridlayout(statusMaximaPanel, [3, 1]);
statusMaximaGrid.RowHeight = {'fit', 'fit', 'fit'};
statusMaximaGrid.Padding = [5 5 5 5];
statusMaximaGrid.RowSpacing = 2;

handles.statusLabel = labelFcn(statusMaximaGrid, 'Text', 'Status: Idle', 'FontWeight', 'bold');
maximaTitleLabel = labelFcn(statusMaximaGrid, 'Text', 'Local Maxima Counts:', 'FontWeight', 'bold');
handles.maximaCountLabel = labelFcn(statusMaximaGrid, 'Text', '', 'VerticalAlignment', 'top');

buttonsPanel = panelFcn(actionsGrid, 'BorderType', 'none');
buttonsGrid = uigridlayout(buttonsPanel, [2, 1]);
buttonsGrid.RowHeight = {'fit', 'fit'};
buttonsGrid.RowSpacing = 5;
handles.updateLivePreviewButton = buttonFcn(buttonsGrid, 'Text', 'Update Previews', 'FontWeight', 'bold');
handles.updateLivePreviewButton.Layout.Row = 1;
handles.abortButton = buttonFcn(buttonsGrid, 'Text', 'Abort Processing', 'FontWeight', 'bold', ...
    'BackgroundColor', [0.8 0.2 0.2], 'FontColor', [1 1 1], 'Enable', 'off', ...
    'Tooltip', 'Stop ongoing processing immediately');
handles.abortButton.Layout.Row = 2;

% Export Options Section - Adaptive Checklist
exportPanel = panelFcn(actionsGrid, 'Title', 'Export Options');
exportMainGrid = uigridlayout(exportPanel, [1, 2]);
exportMainGrid.ColumnWidth = {'1x', 'fit'};
exportMainGrid.Padding = [5 5 5 5];
exportMainGrid.RowSpacing = 5;

% Left Column: Scrollable adaptive checklist (will be populated after data is loaded)
exportListPanel = panelFcn(exportMainGrid, 'BorderType', 'none');
exportListPanel.Layout.Row = 1;
exportListPanel.Layout.Column = 1;
exportListGrid = uigridlayout(exportListPanel, [10, 1]); % Max 10 items
exportListGrid.RowHeight = repmat({'fit'}, 1, 10);
exportListGrid.Padding = [0 0 0 0];
exportListGrid.RowSpacing = 2;

% Pre-allocate export checkboxes (will be shown/hidden dynamically)
handles.exportItemChecks = gobjects(1, 20); % Support up to 20 export items
handles.exportItemLabels = cell(1, 20);
handles.exportItemTypes = cell(1, 20);
handles.numExportItems = 0; % Initialize to 0 - will be set by updateExportChecklist
for i = 1:20
    handles.exportItemChecks(i) = checkboxFcn(exportListGrid, 'Text', '', 'Value', false, 'Visible', 'off');
    if i <= 10
        handles.exportItemChecks(i).Layout.Row = i;
    end
end

% Right Column: Export controls and buttons
exportControlsPanel = panelFcn(exportMainGrid, 'BorderType', 'none');
exportControlsPanel.Layout.Row = 1;
exportControlsPanel.Layout.Column = 2;
exportControlsGrid = uigridlayout(exportControlsPanel, [4, 1]);
exportControlsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
exportControlsGrid.Padding = [0 0 0 0];
exportControlsGrid.RowSpacing = 5;

% Image format selection
formatGrid = uigridlayout(exportControlsGrid, [2, 1]);
formatGrid.RowHeight = {'fit', 'fit'};
formatGrid.Padding = [0 0 0 0];
labelFcn(formatGrid, 'Text', 'Image Format:');
handles.exportImageFormatDrop = dropdownFcn(formatGrid, 'Items', {'TIFF', 'PNG', 'JPEG'}, 'Value', lastUsed.exportImageFormat, 'Tooltip', 'Image format for processed image export');

% Set directory button
handles.exportDirButton = buttonFcn(exportControlsGrid, 'Text', 'Set Directory', 'Tooltip', 'Select export directory');

% Select all button
handles.exportSelectAllButton = buttonFcn(exportControlsGrid, 'Text', 'Select All', 'Tooltip', 'Select all available export items');

% Main export button (at bottom)
handles.exportAllSelectedButton = buttonFcn(exportControlsGrid, 'Text', 'Export', 'FontWeight', 'bold', 'Tooltip', 'Export all checked items');


% --- Analysis Panel ---
analysisPanel = panelFcn(mainGrid, 'Title', 'Nuclei Signal Analysis');
analysisPanel.Layout.Row = 1;
analysisPanel.Layout.Column = 3;

analysisMainGrid = uigridlayout(analysisPanel, [2, 1]);
analysisMainGrid.RowHeight = {'fit', '1x'};
analysisMainGrid.Padding = [5 5 5 5];

% Analysis controls
analysisControlsPanel = panelFcn(analysisMainGrid, 'BorderType', 'none');
analysisControlsGrid = uigridlayout(analysisControlsPanel, [2, 1]);
analysisControlsGrid.RowHeight = {'fit', 'fit'};

% Collapse/expand button
handles.analysisCollapseButton = buttonFcn(analysisControlsGrid, 'Text', '◄ Collapse', 'FontSize', 10);

% Info label
handles.analysisInfoLabel = labelFcn(analysisControlsGrid, 'Text', 'Shows signal composition for each nucleus when segmentation is enabled.', ...
    'FontSize', 10, 'WordWrap', 'on');

% Analysis table container
handles.analysisTablePanel = panelFcn(analysisMainGrid, 'BorderType', 'none');
analysisTableGrid = uigridlayout(handles.analysisTablePanel, [1, 1]);

% Create the analysis table
handles.analysisTable = uitable(analysisTableGrid);
handles.analysisTable.ColumnName = {'Nucleus', 'Area/Vol', 'Ch1', 'Ch2', 'Ch3', 'Ratio'};
handles.analysisTable.ColumnWidth = {'auto', 'auto', 'auto', 'auto', 'auto', 'auto'};
handles.analysisTable.Data = {};
handles.analysisTable.RowName = {};

% Store analysis panel state
handles.analysisCollapsed = false;
handles.analysisData = [];

end