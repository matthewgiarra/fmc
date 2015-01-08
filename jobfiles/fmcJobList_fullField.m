function JOBLIST = fmcJobList_fullField()

% Job information
DefaultJob.ImageType = 'experimental';
DefaultJob.SetType = 'vortex';
DefaultJob.CaseName = 'vortexring_2013-11-12_d03_f60_t06';
DefaultJob.DataRepository = 'analysis/data';
% DefaultJob.DataRepository = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/stagnation_flow/grid_01_test_01/raw';

% Job options
DefaultJob.JobOptions.NumberOfProcessors = 8;
DefaultJob.JobOptions.DataRepositoryIsAbsolute = 0;
DefaultJob.JobOptions.ImageRotationAngle = 0;
DefaultJob.JobOptions.SimulateBeam = 0;
DefaultJob.JobOptions.SimulateNoise = 0;
DefaultJob.JobOptions.NumberOfPasses = 3;
DefaultJob.JobOptions.SkipExisting = 0;
DefaultJob.JobOptions.LeftHanded = 0;
DefaultJob.JobOptions.ComparisonType = 'Eulerian';
DefaultJob.JobOptions.StartFromExistingField = 1;
DefaultJob.JobOptions.StartPass = 1;
DefaultJob.JobOptions.RunCompiled = 1;

% Image parameters
DefaultJob.Parameters.Images.BaseName = 'vortexring_d03_f60_t06_';
DefaultJob.Parameters.Images.Extension = '.tif';
DefaultJob.Parameters.Images.NumberOfDigits = 6;
DefaultJob.Parameters.Images.CorrelationStep = 1;
DefaultJob.Parameters.Images.ParticleConcentration = 0.1;

% Set options
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.Skip = 1;

% Stard and end images
DefaultJob.Parameters.Images.Start = 980;
DefaultJob.Parameters.Images.End = 980;
DefaultJob.Parameters.Images.FrameStep = 10;

% correlation step list
c_steps = [11, 12];  

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

% Image resampling parameters
DefaultJob.Parameters.Processing.Resampling.NumberOfRings = 64;
DefaultJob.Parameters.Processing.Resampling.NumberOfWedges = 256;
DefaultJob.Parameters.Processing.Resampling.MinimumRadius = 2;

% FFT Dimensions.
DefaultJob.Parameters.Processing.FFTSize = [128, 128];

% Correlation parameters
DefaultJob.Parameters.Processing.Correlation.Method = 'fmc'; % Correlation type
DefaultJob.Parameters.Processing.Correlation.SpatialRPCDiameter = 2.8; % RPC filter diameter (pixels)
DefaultJob.Parameters.Processing.Correlation.FMCDiameter = 3.3; % FMC filter diameter (pixels)
DefaultJob.Parameters.Processing.Correlation.FMCFilterType = 'relative'; % Relative or fixed FMC filter diameter
DefaultJob.Parameters.Processing.DwoDifferenceMethod = 'central';
DefaultJob.Parameters.Processing.DoDiscreteWindowOffset = 0;
DefaultJob.Parameters.Processing.MultiPeak = 0;
DefaultJob.Parameters.Processing.FmcDifferenceMethod = 'forward';
DefaultJob.Parameters.Processing.DwoConvergenceCriteria = 0.005;
DefaultJob.Parameters.Processing.DwoConverge = 0;
DefaultJob.Parameters.Processing.DwoMaxConvergenceIterations = 1;
DefaultJob.Parameters.Processing.DoImageDeformation = 1;
DefaultJob.Parameters.Processing.DoImageDisparity = 0;
DefaultJob.Parameters.Processing.Smoothing.DoSmoothing = 1;
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
SegmentItem.JobOptions.ImageRotationAngle = 0;
SegmentItem.Parameters.Images.CorrelationStep = 1;
SegmentItem.JobOptions.NumberOfPasses = 1;

SegmentItem.JobOptions.StartFromExistingField = 0;
SegmentItem.JobOptions.StartPass = 1;

% Pass 1
SegmentItem.Parameters.Processing(1) = defaultProcessing;
SegmentItem.Parameters.Processing(1).DwoDifferenceMethod = 'central';
SegmentItem.Parameters.Processing(1).FmcDifferenceMethod = 'central';
SegmentItem.Parameters.Processing(1).InterrogationRegion.FMIWindowSize = [1 1 0];
SegmentItem.Parameters.Processing(1).InterrogationRegion.FMIWindowType = 'hann1';
SegmentItem.Parameters.Processing(1).Resampling.NumberOfRings = 64;
SegmentItem.Parameters.Processing(1).Resampling.NumberOfWedges = 256;
SegmentItem.Parameters.Processing(1).Resampling.MinimumRadius = 1;
SegmentItem.Parameters.Processing(1).Grid.Spacing.X = 16;
SegmentItem.Parameters.Processing(1).Grid.Spacing.Y = 16;
SegmentItem.Parameters.Processing(1).Grid.Buffer.Y = [0, 0];
SegmentItem.Parameters.Processing(1).Grid.Buffer.X = [0, 0];
SegmentItem.Parameters.Processing(1).InterrogationRegion.Height = 128;
SegmentItem.Parameters.Processing(1).InterrogationRegion.Width = 128;
SegmentItem.Parameters.Processing(1).DoDiscreteWindowOffset = 0;
SegmentItem.Parameters.Processing(1).DwoConverge = 0;
SegmentItem.Parameters.Processing(1).DwoMaxConvergenceIterations = 1;
SegmentItem.Parameters.Processing(1).DoImageDeformation = 1;
SegmentItem.Parameters.Processing(1).DoImageDisparity = 0;
SegmentItem.Parameters.Processing(1).Smoothing.DoSmoothing = 1;
SegmentItem.Parameters.Processing(1).Correlation.Method = 'rpc_deform';
SegmentItem.Parameters.Processing(1).InterrogationRegion.SpatialWindowFraction = [0.50 0.50];

num_passes = SegmentItem.JobOptions.NumberOfPasses;

% Initialize second pass
SegmentItem.Parameters.Processing(2) = SegmentItem.Parameters.Processing(1);

% Update region size
SegmentItem.Parameters.Processing(2).InterrogationRegion.Height = 64;
SegmentItem.Parameters.Processing(2).InterrogationRegion.Width = 64;
SegmentItem.Parameters.Processing(2).Grid.Spacing.Y = 8;
SegmentItem.Parameters.Processing(2).Grid.Spacing.X = 8;

% Repeat the second pass
for p = 2 : num_passes
    SegmentItem.Parameters.Processing(p) = SegmentItem.Parameters.Processing(2);
end

% Append segment item.
for k = 1 : length(c_steps);
    for p = 1 : num_passes
        SegmentItem.Parameters.Images.CorrelationStep = c_steps(k);
    end
   
    % Turn off smoothing for the last pass.
    SegmentItem.Parameters.Processing(end).Smoothing.DoSmoothing = 0;
    
   JOBLIST(k) = SegmentItem;
end


% Switch to FMC
for p = 1 : num_passes
    SegmentItem.Parameters.Processing(p).Correlation.Method = 'fmc_deform';
end

% Length of current joblist
nJobs = length(JOBLIST);

% Append segment item.
for k = 1 : length(c_steps);
   SegmentItem.Parameters.Images.CorrelationStep = c_steps(k);
   JOBLIST(nJobs + k) = SegmentItem;
end

end

