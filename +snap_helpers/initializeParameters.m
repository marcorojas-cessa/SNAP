function defaultParams = initializeParameters(Nmax)
    % Initializes and returns a structure with default parameters.

    defaultParams.dicFilePath = '';
    defaultParams.nucFilePath = '';
    defaultParams.channelFilePaths = repmat({''}, 1, Nmax);
    defaultParams.numChannels = 1;
    
    % Image Spacing
    defaultParams.xySpacing = num2cell(0.128866*ones(1, Nmax));
    defaultParams.zSpacing = num2cell(0.3*ones(1, Nmax));
    
    % Deconvolution
    defaultParams.deconvEnabled = num2cell(false(1, Nmax));
    defaultParams.deconvMethod = repmat({'Lucy-Richardson'}, 1, Nmax);
    % Lucy-Richardson parameters
    defaultParams.deconvLRIterations = num2cell(10*ones(1, Nmax));
    defaultParams.deconvLRDamping = num2cell(0*ones(1, Nmax)); % 0 = no damping
    % Wiener parameters
    defaultParams.deconvWienerNSR = num2cell(0.01*ones(1, Nmax)); % Noise-to-Signal Ratio
    % Blind deconvolution parameters
    defaultParams.deconvBlindIterations = num2cell(10*ones(1, Nmax));
    defaultParams.deconvBlindUnderRelax = num2cell(0*ones(1, Nmax)); % Under-relaxation parameter
    % PSF parameters (shared across all methods)
    defaultParams.deconvPSFSource = repmat({'Generate'}, 1, Nmax); % 'Generate' or 'Load File'
    defaultParams.deconvPSFFilePath = repmat({''}, 1, Nmax);
    defaultParams.deconvPSFSigmaXY = num2cell(1.5*ones(1, Nmax)); % PSF sigma in XY (pixels)
    defaultParams.deconvPSFSigmaZ = num2cell(2.0*ones(1, Nmax));  % PSF sigma in Z (slices)
    defaultParams.deconvPSFSizeXY = num2cell(5*ones(1, Nmax));    % PSF half-size in XY (pixels)
    defaultParams.deconvPSFSizeZ = num2cell(3*ones(1, Nmax));     % PSF half-size in Z (slices)
    
    % Pre-processing
    defaultParams.preProcMethod = repmat({'None'}, 1, Nmax); % New dropdown parameter
    defaultParams.preProcMode = repmat({'2D (Slice-by-slice)'}, 1, Nmax);
    defaultParams.preProcScale = num2cell(false(1, Nmax));
    defaultParams.preProcProjection = repmat({'Max'}, 1, Nmax);
    defaultParams.smoothGaussianValues = num2cell(ones(1, Nmax));
    defaultParams.smoothMedianValues = num2cell(ones(1, Nmax));
    defaultParams.nlmFilterStrength = num2cell(0.1*ones(1, Nmax)); % Non-local means filter strength
    defaultParams.nlmSearchWindow = num2cell(21*ones(1, Nmax)); % Non-local means search window size
    defaultParams.nlmComparisonWindow = num2cell(7*ones(1, Nmax)); % Non-local means comparison window size
    defaultParams.waveletName = repmat({'haar'}, 1, Nmax);
    defaultParams.waveletLevel = num2cell(ones(1, Nmax));
    defaultParams.waveletThresholdRule = repmat({'sqtwolog'}, 1, Nmax);
    defaultParams.waveletThresholdMethod = repmat({'soft'}, 1, Nmax);
    defaultParams.preprocClipAtZero = num2cell(true(1, Nmax)); % Section-wide clip
    defaultParams.preprocEnabled = num2cell(false(1, Nmax));
    
    % Background Correction
    defaultParams.bgCorrMode = repmat({'2D (Slice-by-slice)'}, 1, Nmax);
    defaultParams.bgCorrScale = num2cell(false(1, Nmax));
    defaultParams.bgCorrProjection = repmat({'Max'}, 1, Nmax);
    defaultParams.bgMethod = repmat({'None'}, 1, Nmax);
    defaultParams.bgParam = num2cell(10*ones(1, Nmax));
    defaultParams.bgCorrClipAtZero = num2cell(true(1, Nmax)); % Section-wide clip
    defaultParams.bgCorrEnabled = num2cell(false(1, Nmax));

    % Local Maxima Detection
    defaultParams.maximaMode = repmat({'2D (Slice-by-slice)'}, 1, Nmax);
    defaultParams.maximaScale = num2cell(false(1, Nmax));
    defaultParams.maximaProjection = repmat({'Max'}, 1, Nmax);
    defaultParams.maximaMethod = repmat({'Simple Regional'}, 1, Nmax);
    % Local maxima neighborhood size (in pixels/voxels or microns if scaled)
    % 1 = 3x3 window (good for closely packed peaks), 2 = 5x5 window (more spacing)
    defaultParams.maximaNeighborhoodSize = num2cell(1*ones(1, Nmax));
    defaultParams.hMaxValue = num2cell(0*ones(1, Nmax));
    defaultParams.sigmaValue = num2cell(2*ones(1, Nmax));
    defaultParams.peakThresholdValue = num2cell(0.01*ones(1, Nmax));
    defaultParams.showMaxima = num2cell(false(1, Nmax));
    defaultParams.maximaColor = repmat({'Red'}, 1, Nmax);
    defaultParams.maximaEnabled = num2cell(false(1, Nmax));
    defaultParams.displayOnAllPreviews = num2cell(false(1, Nmax));

    % Nuclei Spacing
    defaultParams.nucXYSpacing = 0.128866;
    defaultParams.nucZSpacing = 0.3;

    % Nuclei Deconvolution
    defaultParams.nucDeconvEnabled = false;
    defaultParams.nucDeconvMethod = 'Lucy-Richardson';
    defaultParams.nucDeconvLRIterations = 10;
    defaultParams.nucDeconvLRDamping = 0;
    defaultParams.nucDeconvWienerNSR = 0.01;
    defaultParams.nucDeconvBlindIterations = 10;
    defaultParams.nucDeconvBlindUnderRelax = 0;
    defaultParams.nucDeconvPSFSource = 'Generate';
    defaultParams.nucDeconvPSFFilePath = '';
    defaultParams.nucDeconvPSFSigmaXY = 1.5;
    defaultParams.nucDeconvPSFSigmaZ = 2.0;
    defaultParams.nucDeconvPSFSizeXY = 5;
    defaultParams.nucDeconvPSFSizeZ = 6;

    % Nuclei Pre-processing
    defaultParams.nucPreProcMode = '2D (Slice-by-slice)';
    defaultParams.nucPreProcScale = false;
    defaultParams.nucPreProcProjection = 'Max';
    defaultParams.nucPreProcMethod = 'None';
    defaultParams.nucPreprocClipAtZero = true;
    defaultParams.nucSmoothGaussianValue = 1;
    defaultParams.nucSmoothMedianValue = 1;
    defaultParams.nucNlmFilterStrength = 0.1; % Non-local means filter strength for nuclei
    defaultParams.nucNlmSearchWindow = 21; % Non-local means search window size for nuclei
    defaultParams.nucNlmComparisonWindow = 7; % Non-local means comparison window size for nuclei
    defaultParams.nucWaveletName = 'haar';
    defaultParams.nucWaveletLevel = 1;
    defaultParams.nucWaveletThresholdRule = 'sqtwolog';
    defaultParams.nucWaveletThresholdMethod = 'soft';
    defaultParams.nucPreprocEnabled = false;

    % Nuclei Background Correction
    defaultParams.nucBgCorrMode = '2D (Slice-by-slice)';
    defaultParams.nucBgCorrScale = false;
    defaultParams.nucBgCorrProjection = 'Max';
    defaultParams.nucBgMethod = 'None';
    defaultParams.nucBgParam = 10;
    defaultParams.nucBgCorrClipAtZero = true;
    defaultParams.nucBgCorrEnabled = false;

    % Nuclei Segmentation
    defaultParams.nucSegMode = '2D (Slice-by-slice)';
    defaultParams.nucSegProjection = 'Max';
    
    % New segmentation method parameters
    defaultParams.nucSegMainMethod = 'Absolute';
    defaultParams.nucSegSubMethod = 'Std Multiplier';
    defaultParams.nucSegAbsoluteThreshold = 25;
    defaultParams.nucSegStdMultiplier = 2.0;
    defaultParams.nucSegOffset = 500;
    defaultParams.nucSegLocalAlgorithm = 'Otsu';
    defaultParams.nucSegLocalRadius = 25; % Radius in pixels (equivalent to window size 51)
    
    % Algorithm-specific parameters (dynamic parameter system)
    defaultParams.nucSegAlgParam1 = 25; % Default radius value
    defaultParams.nucSegAlgParam2 = 0;
    
    % Individual algorithm parameters (for backward compatibility)
    defaultParams.nucSegBernsenContrast = 15;
    defaultParams.nucSegMeanC = 0;
    defaultParams.nucSegMedianC = 0;
    defaultParams.nucSegMidGreyC = 0;
    defaultParams.nucSegNiblackK = 0.2;
    defaultParams.nucSegNiblackC = 0;
    defaultParams.nucSegPhansalkarK = 0.25;
    defaultParams.nucSegPhansalkarR = 0.5;
    defaultParams.nucSegSauvolaK = 0.5;
    defaultParams.nucSegSauvolaR = 128;
    
    % Default checkbox states (dynamic parameter system)
    defaultParams.nucSegAlgParam1Default = true;
    defaultParams.nucSegAlgParam2Default = true;
    
    % Individual default states (for backward compatibility)
    defaultParams.nucSegBernsenContrastDefault = true;
    defaultParams.nucSegMeanCDefault = true;
    defaultParams.nucSegMedianCDefault = true;
    defaultParams.nucSegMidGreyCDefault = true;
    defaultParams.nucSegNiblackKDefault = true;
    defaultParams.nucSegNiblackCDefault = true;
    defaultParams.nucSegPhansalkarKDefault = true;
    defaultParams.nucSegPhansalkarRDefault = true;
    defaultParams.nucSegSauvolaKDefault = true;
    defaultParams.nucSegSauvolaRDefault = true;
    defaultParams.nucShowSeg = false;
    defaultParams.nucSegEnabled = false;
    
    % Nuclei Inclusion/Exclusion
    defaultParams.nucInclusionExclusionEnabled = false;
    defaultParams.nucInclusionExclusionMode = 'Include Inside Nuclei';
    defaultParams.nucInclusionExclusionApplyTo = 'All Channels';
    defaultParams.nucExcludeEdges = true;
    defaultParams.nucNavPanelIndex = 0;
    
    % Nuclei Filtering
    defaultParams.nucFilterEnabled = false;
    defaultParams.nucFilterSizeEnabled = false;
    defaultParams.nucFilterMinSize = 100;
    defaultParams.nucFilterSizeUnit = 'pixels';
    defaultParams.nucFilterCircularityEnabled = false;
    defaultParams.nucFilterMinCircularity = 0.3;
    defaultParams.nucFilterSolidityEnabled = false;
    defaultParams.nucFilterMinSolidity = 0.7;

    % Export settings
    defaultParams.exportNucleiEnabled = false;
    defaultParams.exportChannelDataEnabled = false;
    defaultParams.exportClusteredDataEnabled = false;
    defaultParams.exportParametersEnabled = false;
    defaultParams.exportProcessedImagesEnabled = false;
    defaultParams.exportImageFormat = 'TIFF';
    defaultParams.exportDirectory = pwd; % Default to current working directory

    % Preview settings
    defaultParams.previewContents = repmat({'None'}, 1, 5);
    defaultParams.previewModes = repmat({'Z-Stack'}, 1, 5);
    defaultParams.previewProjections = repmat({'Max'}, 1, 5);

    % Index for the currently visible navigation panel in each channel tab
    % 0 = Stage 0 (default), 1 = Stage 1, 2 = Stage 2
    defaultParams.navPanelIndex = num2cell(zeros(1, Nmax));

    % --- Local Maxima Fitting Parameters ---
    defaultParams.gaussFitVoxelWindowSize = num2cell(repmat(7, 1, Nmax));
    defaultParams.gaussFitBgCorrMethod = repmat({'Mean Surrounding Subtraction'}, 1, Nmax);
    defaultParams.gaussFitBgCorrWidth = num2cell(repmat(1, 1, Nmax));
    defaultParams.gaussFitPolyDegree = num2cell(repmat(2, 1, Nmax));
    defaultParams.gaussFitMethod = repmat({'1D (X,Y,Z) Gaussian'}, 1, Nmax);
    defaultParams.gaussFitPlotCheck = num2cell(zeros(1, Nmax));
    defaultParams.gaussFitMaxIterations = num2cell(repmat(200, 1, Nmax));
    defaultParams.gaussFitTolerance = num2cell(repmat(1e-6, 1, Nmax));
    defaultParams.gaussFitRadialRadius = num2cell(repmat(3, 1, Nmax));  % Radius for radial symmetry method
    defaultParams.gaussFitEnabled = num2cell(false(1, Nmax));


    % --- Fit Filtering Parameters ---
    defaultParams.fitFilterEnabled = num2cell(false(1, Nmax));
    defaultParams.fitFilterRSquaredEnabled = num2cell(false(1, Nmax));
    defaultParams.fitFilterRSquaredMin = num2cell(repmat(0.8, 1, Nmax));
    defaultParams.fitFilterRSquaredMax = num2cell(repmat(1.0, 1, Nmax));
    defaultParams.fitFilterSigmaSumEnabled = num2cell(false(1, Nmax));
    defaultParams.fitFilterSigmaSumMin = num2cell(repmat(0.0, 1, Nmax));
    defaultParams.fitFilterSigmaSumMax = num2cell(repmat(10.0, 1, Nmax));
    defaultParams.fitFilterAmplitudeEnabled = num2cell(false(1, Nmax));
    defaultParams.fitFilterAmplitudeMin = num2cell(repmat(100.0, 1, Nmax));
    defaultParams.fitFilterAmplitudeMax = num2cell(repmat(10000.0, 1, Nmax));
    defaultParams.fitFilterIntensityEnabled = num2cell(false(1, Nmax));
    defaultParams.fitFilterIntensityMin = num2cell(repmat(1000.0, 1, Nmax));
    defaultParams.fitFilterIntensityMax = num2cell(repmat(100000.0, 1, Nmax));
    
    % --- Additional Fit Result Parameters ---
    defaultParams.gaussFitRhoXY = num2cell(nan(1, Nmax));
    defaultParams.gaussFitRhoXZ = num2cell(nan(1, Nmax));
    defaultParams.gaussFitRhoYZ = num2cell(nan(1, Nmax));
    defaultParams.gaussFitAlphaX = num2cell(nan(1, Nmax));
    defaultParams.gaussFitAlphaY = num2cell(nan(1, Nmax));
    defaultParams.gaussFitAlphaZ = num2cell(nan(1, Nmax));

    % --- Preview Parameters ---
    defaultParams.previewModes = repmat({'Z-Projection'}, 1, 5);
    defaultParams.previewProjections = repmat({'Max'}, 1, 5);
    
    % --- Z-Playback Parameters ---
    defaultParams.zPlaybackActive = false;
    defaultParams.zPlaybackSpeed = 0.2; % seconds per frame (5 FPS)
end
