function JOBLIST = fmcJobList_fullField()

% Job information
DefaultJob.ImageType = 'experimental';
DefaultJob.SetType = 'micro_piv';
DefaultJob.CaseName = 'piv_benchmark_vortex_images';
% DefaultJob.DataRepository = '~/Desktop/vortex/raw';
DefaultJob.DataRepository = ['~/Documents/School/VT/Research/' ...
    'Aether/spc/analysis/data/micro_piv/ts3_000014/raw'];

% Job options
DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.DataRepositoryIsAbsolute = 1;
DefaultJob.JobOptions.ImageRotationAngle = 0;
DefaultJob.JobOptions.SimulateBeam = 0;
DefaultJob.JobOptions.SimulateNoise = 0;
DefaultJob.JobOptions.NumberOfPasses = 5;
DefaultJob.JobOptions.SkipExisting = 0;
DefaultJob.JobOptions.LeftHanded = 0;
DefaultJob.JobOptions.ComparisonType = 'Eulerian';
DefaultJob.JobOptions.StartFromExistingField = 1;
DefaultJob.JobOptions.StartPass = 1;
DefaultJob.JobOptions.RunCompiled = 1;

% Image parameters
DefaultJob.Parameters.Images.BaseName = 'frame_';
DefaultJob.Parameters.Images.Extension = '.tif';
DefaultJob.Parameters.Images.NumberOfDigits = 7;
DefaultJob.Parameters.Images.CorrelationStep = 1;
DefaultJob.Parameters.Images.ParticleConcentration = 0.1;

% Set options
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.Skip = 1;

% Stard and end images
DefaultJob.Parameters.Images.Start = 1;
DefaultJob.Parameters.Images.End = 1;
DefaultJob.Parameters.Images.FrameStep = 1;
DefaultJob.Parameters.Images.ColorChannel = 3;

% correlation step list
% c_steps = [11, 12];  

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0;
DefaultJob.Parameters.Noise.Std = 0;
DefaultJob.Parameters.Noise.Background = 0;

% Grid parameters
DefaultJob.Parameters.Processing.Grid.Spacing.X = 16; % Horizontal spacing between grid points (pixels)
DefaultJob.Parameters.Processing.Grid.Spacing.Y = 16; % Vertical spacing between grid points (pixels)
DefaultJob.Parameters.Processing.Grid.Buffer.Y = [32 32]; % Top and bottom grid buffer (pixels)
DefaultJob.Parameters.Processing.Grid.Buffer.X = [32 32]; % Left and right grid buffer (pixels)

% Interrogation Region Parameters
DefaultJob.Parameters.Processing.InterrogationRegion.Height = 128;
DefaultJob.Parameters.Processing.InterrogationRegion.Width = 128;
DefaultJob.Parameters.Processing.InterrogationRegion.SpatialWindowFraction = [0.5 0.5];
DefaultJob.Parameters.Processing.InterrogationRegion.FMIWindowSize = [2, 2, 0];
DefaultJob.Parameters.Processing.InterrogationRegion.FMIWindowType = 'hann1';
DefaultJob.Parameters.Processing.InterrogationRegion.ZeroMeanRegion = 1;

% Image resampling parameters
DefaultJob.Parameters.Processing.Resampling.NumberOfRings = 64;
DefaultJob.Parameters.Processing.Resampling.NumberOfWedges = 256;
DefaultJob.Parameters.Processing.Resampling.MinimumRadius = 2;

% FFT Dimensions.
DefaultJob.Parameters.Processing.FFTSize = [128, 128];

% Correlation parameters
DefaultJob.Parameters.Processing.Correlation.Method = 'rpc'; % Correlation type
DefaultJob.Parameters.Processing.Correlation.SpatialRPCDiameter = 2.8; % RPC filter diameter (pixels)
DefaultJob.Parameters.Processing.Correlation.FMCDiameter = 3.3; % FMC filter diameter (pixels)
DefaultJob.Parameters.Processing.Correlation.FMCFilterType = 'relative'; % Relative or fixed FMC filter diameter

% Peak fit method
DefaultJob.Parameters.Processing.Correlation.PeakFitMethod =...
    'least-squares';

DefaultJob.Parameters.Processing.DwoDifferenceMethod = 'central';
DefaultJob.Parameters.Processing.DoDiscreteWindowOffset = 0;
DefaultJob.Parameters.Processing.MultiPeak = 0;
DefaultJob.Parameters.Processing.FmcDifferenceMethod = 'forward';
DefaultJob.Parameters.Processing.DwoConvergenceCriteria = 0.005;
DefaultJob.Parameters.Processing.DwoConverge = 0;
DefaultJob.Parameters.Processing.DwoMaxConvergenceIterations = 1;
DefaultJob.Parameters.Processing.DoImageDeformation = 0;
DefaultJob.Parameters.Processing.DoImageDisparity = 0;
DefaultJob.Parameters.Processing.Smoothing.DoSmoothing = 0;
DefaultJob.Parameters.Processing.Smoothing.KernelDiameter = 7; 
DefaultJob.Parameters.Processing.Smoothing.KernelGaussianStdDev = 1;

% Universal Outlier Detection Parameters
DefaultJob.Parameters.Processing.Validation.UodStencilRadius = 1;
DefaultJob.Parameters.Processing.Validation.UodThreshold = 3;
DefaultJob.Parameters.Processing.Validation.UodExpectedDifference = [0.1, 0.1];

% Default Processing parameters
defaultProcessing = DefaultJob.Parameters.Processing;

% Job 2, Pass 1
SegmentItem = DefaultJob;
SegmentItem.Parameters.Images.CorrelationStep = 1;
SegmentItem.JobOptions.NumberOfPasses = 1;

SegmentItem.JobOptions.StartFromExistingField = 0;
SegmentItem.JobOptions.StartPass = 1;

% Pass 1
SegmentItem.Parameters.Processing(1) = defaultProcessing;
SegmentItem.Parameters.Processing(1).DwoDifferenceMethod = 'central';
SegmentItem.Parameters.Processing(1).Grid.Spacing.X = 16;
SegmentItem.Parameters.Processing(1).Grid.Spacing.Y = 16;
SegmentItem.Parameters.Processing(1).Grid.Buffer.Y = [0, 0];
SegmentItem.Parameters.Processing(1).Grid.Buffer.X = [0, 0];
SegmentItem.Parameters.Processing(1).InterrogationRegion.Height = 128;
SegmentItem.Parameters.Processing(1).InterrogationRegion.Width = 128;
SegmentItem.Parameters.Processing(1).DoDiscreteWindowOffset = 0;
SegmentItem.Parameters.Processing(1).DwoConverge = 0;
SegmentItem.Parameters.Processing(1).DwoMaxConvergenceIterations = 1;
SegmentItem.Parameters.Processing(1).DoImageDeformation = 0;
SegmentItem.Parameters.Processing(1).DoImageDisparity = 0;
SegmentItem.Parameters.Processing(1).Smoothing.DoSmoothing = 0;
SegmentItem.Parameters.Processing(1).Correlation.Method = 'rpc';
SegmentItem.Parameters.Processing(1). ...
    InterrogationRegion.SpatialWindowFraction = [0.50 0.50];

% Add to job list
JOBLIST(1) = SegmentItem;



end

