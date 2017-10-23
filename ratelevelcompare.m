%---------------------------------------------------------------------
%---------------------------------------------------------------------
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171004';

% files to compare
file1 = '1155_20171004_01_01_2563_BBN_LEVEL.dat';
file2 = '1155_20171004_01_01_2563_BBN_LEVEL_optoON.dat';

%---------------------------------------------------------------------
% settings for processing data
%---------------------------------------------------------------------
% filter
HPFreq = 350;
LPFreq = 6500;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
% Channel Number (use 8 for single channel data)
channelNumber = 8;
% binSize for PSTH (milliseconds)
binSize = 5;
% SAVE PLOTS?
saveFIG = 0;
savePNG = 0;
savePDF = 1;

%---------------------------------------------------------------------
%% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% read in data files, plot data
%---------------------------------------------------------------------
% build datapath
datapath = fullfile(data_root_path, animalID, dateID);

[~, Dinf1, S1] = optoproc('file', fullfile(datapath, file1), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'plotPSTH', 'savepdf');
							
[~, Dinf2, S2] = optoproc('file', fullfile(datapath, file2), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'plotPSTH', 'savepdf');
							
%---------------------------------------------------------------------
%% construct rate-level function
%---------------------------------------------------------------------
% time window (post-stimulus onset) for analysis
WindowLen = 25;
AnalysisWindow = Dinf1.audio.Delay + [0 WindowLen];

rlf1 = computeRLF(S1.spiketimes, AnalysisWindow);
rlf2 = computeRLF(S2.spiketimes, AnalysisWindow);

figure;
errorbar(Dinf1.audio.Level, rlf1.mean, rlf1.std);
hold on
	errorbar(Dinf2.audio.Level, rlf2.mean, rlf2.std);
hold off

xlabel('dB SPL');
ylabel('mean spikes/trial')
title({file1, file2}, 'Interpreter', 'none')
legend({'ctrl', 'optoOn'})




