function JOBLIST = fmcJobList_mc()

DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.BooleanGenerateParticleImages = false;
DefaultJob.JobOptions.BooleanRunAnalysis = true;
DefaultJob.JobOptions.FlipYTranslation = false;
DefaultJob.JobOptions.SkipExistingSets = false;
DefaultJob.JobOptions.RepositoryPathIsAbsolute = 0;
DefaultJob.JobOptions.DoAffineTransform = 0;
DefaultJob.JobOptions.RunCompiled = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = 'full_field_registration';
DefaultJob.CorrelationType = 'fmc';
DefaultJob.Parameters.RegionHeight = 64;
DefaultJob.Parameters.RegionWidth = 64;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 10;
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

% JOB 1
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'full_field_registration';
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.Sets.ImagesPerSet = 10;
SegmentItem.Parameters.RegionHeight = 512;
SegmentItem.Parameters.RegionWidth = 512;
SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
SegmentItem.Parameters.Processing.SpatialWindowFraction = 0.10 * [1 1];
SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 256;
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 256;
SegmentItem.Parameters.Processing.FFTSize = [128, 128];
JOBLIST(1) = SegmentItem;

% JOB 2
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_translation_rotation';
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.Sets.ImagesPerSet = 1000;
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.Parameters.Processing.FmiRPCDiameter = 3.3;
SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 256;
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 256;
SegmentItem.Parameters.Processing.FFTSize = [128, 128];
SegmentItem.Parameters.Processing.SpatialWindowFraction = [0.5 0.5];
JOBLIST(end + 1) = SegmentItem;


end







