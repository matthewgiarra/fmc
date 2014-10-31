
% Image gen jobfile
% Job Distribution Parameters

% Start and end sets
startSet = 1  +  10 * (s-1);
endSet = min(startSet + 9, 100);

% Analysis paths
addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fmcDisparity;
addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fmcDisparity/jobfiles;
addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fmcDisparity/testCodes;

% % Analysis paths
% addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_withPrana/jobfiles;
% addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_withPrana;

% Interpolation paths
%addpath ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fastLogPolar/ba_interpolation

cd('~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fmcDisparity/jobfiles');
analysisJobList_fmc = fmcJobList_fullField;

% cd('~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_withPrana/jobfiles');
% analysisJobList = fmcJobList_fullField_deform;

for n = 1 : length(analysisJobList_fmc)
	analysisJobList_fmc(n).Parameters.Sets.Start = startSet;
	analysisJobList_fmc(n).Parameters.Sets.End = endSet;
end


% cd ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fastLogPolar;
% runFmcFullFieldJobList(analysisJobList);

% % Close open pools
% if matlabpool('size') > 0
%     matlabpool close
% end
%  

cd ~/mnt/imac/School/VT/Research/Aether/FMISPOMF/analysis/src/FMC/branch_fmcDisparity;
runFmcFullFieldJobList(analysisJobList_fmc);



