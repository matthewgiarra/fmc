function JOBLIST = fmcJobList()


DefaultJob.JobOptions.NumberOfProcessors = 12;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.BooleanGenerateParticleImages = 0;
DefaultJob.JobOptions.BooleanRunAnalysis = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = 'FMCtest_2013-10-21_rotation_scaling_translation';
DefaultJob.CorrelationType = 'fmc';
DefaultJob.Parameters.RegionHeight = 64;
DefaultJob.Parameters.RegionWidth = 64;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1000;
DefaultJob.Parameters.Sets.ImagesPerSet = 1000;

DefaultJob.Parameters.Processing.SpatialWindowFraction = [0.5 0.5];
DefaultJob.Parameters.Processing.FMIWindowFraction = [1 0.8];
DefaultJob.Parameters.Processing.SpatialWindowType = 'fraction';
DefaultJob.Parameters.Processing.FMIWindowType = 'fraction';
DefaultJob.Parameters.Processing.SpatialRPCDiameter = 2.8;
DefaultJob.Parameters.Processing.FmiRPCDiameter = 3.3; % Diameter for 128 wedges

% This is the number of wedges  a multiple of the image diameter.
DefaultJob.Parameters.Processing.NumberOfWedges.Min = 255;
DefaultJob.Parameters.Processing.NumberOfWedges.Max = 255;
DefaultJob.Parameters.Processing.NumberOfRings.Min = 64;
DefaultJob.Parameters.Processing.NumberOfRings.Max = 64;

% Subsequent jobs
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_2013-09-03_smallRotations_largeRange';
SegmentItem.Parameters.RegionHeight = 128;
SegmentItem.Parameters.RegionWidth = 128;
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 16; % 225 wedges for 128x128 or 64x64 regions
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 625; % 225 wedges
SegmentItem.Parameters.Processing.NumberOfRings.Min = 64;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 64;
JOBLIST(1) = SegmentItem;

% Subsequent jobs
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_2013-10-17_smallRotations_largeRange';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 16; % 225 wedges for 128x128 or 64x64 regions
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 625; % 225 wedges
SegmentItem.Parameters.Processing.NumberOfRings.Min = 32;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 32;
JOBLIST(end + 1) = SegmentItem;

% Subsequent jobs
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_2013-10-21_rotation_scaling_translation';
SegmentItem.CorrelationType = 'fmc';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.Parameters.Processing.NumberOfWedges.Min = 112; % 112 wedges
SegmentItem.Parameters.Processing.NumberOfWedges.Max = 112; % 112 wedges
SegmentItem.Parameters.Processing.NumberOfRings.Min = 32;
SegmentItem.Parameters.Processing.NumberOfRings.Max = 32;
JOBLIST(end + 1) = SegmentItem;

% Subsequent jobs
SegmentItem = DefaultJob;
SegmentItem.CaseName = 'FMCtest_2013-10-21_rotation_scaling_translation';
SegmentItem.CorrelationType = 'rpc';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
JOBLIST(end + 1) = SegmentItem;

end

