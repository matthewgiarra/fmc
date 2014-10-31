function JOBLIST = fmcJobList_mc()

DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.BooleanGenerateParticleImages = false;
DefaultJob.JobOptions.BooleanRunAnalysis = true;
DefaultJob.JobOptions.FlipYTranslation = false;
DefaultJob.JobOptions.SkipExistingSets = false;
DefaultJob.JobOptions.RepositoryPathIsAbsolute = 0;
DefaultJob.JobOptions.DoAffineTransform = 0;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = 'FMCtest_2014-03-18_rotation_translation';
DefaultJob.CorrelationType = 'fmc';
DefaultJob.Parameters.RegionHeight = 64;
DefaultJob.Parameters.RegionWidth = 64;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 1000;
DefaultJob.Parameters.RepositoryPath =  '/home/voodoo/Desktop/images';

DefaultJob.Parameters.Processing.SpatialWindowFraction = [0.5 0.5];
DefaultJob.Parameters.Processing.FMIWindowFraction = [1, 1, 0];
DefaultJob.Parameters.Processing.SpatialWindowType = 'fraction';
DefaultJob.Parameters.Processing.FMIWindowType = 'fraction';
DefaultJob.Parameters.Processing.SpatialRPCDiameter = 2.8;
DefaultJob.Parameters.Processing.FmiRPCDiameter = 3.3; % Diameter for 128 wedges
DefaultJob.Parameters.Processing.FFTSize = [256, 256];
DefaultJob.Parameters.Processing.FMIWindowType = 'hann1';

% These are the Log-Polar resampling parameters.
DefaultJob.Parameters.Processing.NumberOfWedges.Min = 128;
DefaultJob.Parameters.Processing.NumberOfWedges.Max = 128;
DefaultJob.Parameters.Processing.NumberOfRings.Min = 64;
DefaultJob.Parameters.Processing.NumberOfRings.Max = 64;
DefaultJob.Parameters.Processing.MinimumRadius = 2;

% This is the mean of the additive gaussian white noise
% as a fraction of the maximum image intensity
DefaultJob.Parameters.Processing.Noise.Mean = 0.00;

% This is the 99.5% confidence interval of the noise
% as a fraction of the maximum image intensity.
DefaultJob.Parameters.Processing.Noise.Std = 0.00;

% % JOB 1
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_scaling_translation';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% 
% % JOB 2
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_rotation_translation';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% 
% % Job 3
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_rotation_scaling';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;
% 
% % JOB 4
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_scaling_translation';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;
% 
% % JOB 5
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_rotation_translation';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;
% % 
% % Job 6
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-03-31_rotation_scaling';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;

% % JOB 7
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-04-29_rotation_translation';
% SegmentItem.CorrelationType = 'lpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% DefaultJob.Parameters.Processing.MinimumRadius = 2;
% JOBLIST(1) = SegmentItem;

% 
% JOB 8
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_2014-04-29_rotation_translation';
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.RegionHeight = 128;
SegmentItem.Parameters.RegionWidth = 128;
SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
SegmentItem.Parameters.Processing.MinimumRadius = 5;
JOBLIST(1) = SegmentItem;


% % % 
% JOB 9
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-04-29_scaling_translation';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;

% % JOB 10
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-04-29_scaling_translation';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;
% 
% JOB 10
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-04-29_rotation_scaling';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% % 
% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-04-29_rotation_scaling';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;


% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-01_shearing_rotation';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 128;
% SegmentItem.Parameters.RegionWidth = 128;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% 
% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-01_shearing_rotation';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;


% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-06_shearing_only';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% 
% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-06_shearing_only';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;


% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-08_rotation_only';
% SegmentItem.CorrelationType = 'lpc';
% SegmentItem.SetType = 'mc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% SegmentItem.Parameters.Processing.MinimumRadius = 2;
% SegmentItem.Parameters.Processing.FMIWindowFraction = [1, 1, 0];
% SegmentItem.Parameters.Processing.FFTSize = [64, 64];
% JOBLIST(1) = SegmentItem;
% 

% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-13_shearing_only';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.SetType = 'mc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% SegmentItem.Parameters.Processing.MinimumRadius = 2;
% SegmentItem.Parameters.Processing.FMIWindowFraction = [1, 1, 0];
% SegmentItem.Parameters.Processing.FFTSize = [64, 64];
% JOBLIST(1) = SegmentItem;


% % JOB 11
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-13_shearing_only';
% SegmentItem.CorrelationType = 'fmc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(1) = SegmentItem;
% 
% SegmentItem = DefaultJob;
% SegmentItem.CaseName = 'FMCtest_2014-05-13_shearing_only';
% SegmentItem.CorrelationType = 'rpc';
% SegmentItem.Parameters.RegionHeight = 64;
% SegmentItem.Parameters.RegionWidth = 64;
% SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
% SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.50 * [1 1];
% SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
% SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
% SegmentItem.Parameters.Processing.NumberOfWedges.Min = 128;
% SegmentItem.Parameters.Processing.NumberOfWedges.Max = 128;
% JOBLIST(end + 1) = SegmentItem;
% 


end







