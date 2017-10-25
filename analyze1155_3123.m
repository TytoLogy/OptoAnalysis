%
% need to compute mean fr for ch 11 rebound response across illumination
% levels
%

% '1155_20171006_04_03_3123_BBN_optoOFF_ch5_ch11_3.dat'
BBNratelevel =  '1155_20171006_04_03_3123_BBN_optoOFF_ch5ch11.dat';

% no stim (opto or sound)
control = '1155_20171006_04_03_3123_BBN_optoOFF_ch5ch11_2.dat';
% 500 opto only
o500 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11.dat';
% 1000 opto only
o1000 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_2.dat';
% 1500 opto only
o1500 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_3.dat';
% 2000 opto only
o2000 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_4.dat';
% 2500 opto only
o2500 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_5.dat';
% 3000 opto only
o3000 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_6.dat';
% 3500 opto only
o3500 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_7.dat';
% 3750 opto only
o3750 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_8.dat';
% 3750 rate level bbn
oBBNratelevel = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_9.dat';

% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171006';

% Channel Number (use 8 for single channel data)
channelNumber = 11;

%---------------------------------------------------------------------
%% settings for processing data
%---------------------------------------------------------------------

% time window (post-stimulus onset) for analysis : this is added to stim
% delay
WindowLen = [0 75];

% filter
HPFreq = 350;
LPFreq = 6500;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
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

file = o3750;

[~, fname] = fileparts(file);
plotname = sprintf('%s_ch%d', fname, channelNumber);

[~, Dinf, S] = optoproc('file', fullfile(datapath, file), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'plotPSTH', 'savepdf', ...
								'plotFile', plotname);

%---------------------------------------------------------------------
% construct rate-level function
%---------------------------------------------------------------------
AnalysisWindow = Dinf.audio.Delay + WindowLen;

rlf = computeRLF(S.spiketimes, AnalysisWindow);

figure;
if length(Dinf.audio.Level) > 1
	errorbar(Dinf.audio.Level, rlf.mean, rlf.std);
else
	plot(Dinf.audio.Level, rlf.mean, 'o');
end
xlabel('dB SPL');
ylabel('mean spikes/trial +/- s.d.')
title({file}, 'Interpreter', 'none')
