%------------------------------------------------------------------------
% analyze data from animal 1122
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 27 Oct, 2017 (SJS) 
%	- adapted from analyze1122.m
% 
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% setup
%------------------------------------------------------------------------
%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%------------------------------------------------------------------------
%% Data file information
%------------------------------------------------------------------------
% animal id 
animalID = '1122';
% date code (helps locate files in data directory)
dateID = '20171024';
% unit #
unit = '01';
% penetration
penetration = '03';
% recording depth
depth = '3751';
% channel
channelNumber = 8;

% path to data
datapath = fullfile(data_root_path, animalID, dateID);

% datafiles
file_list = {	...
	% 30 reps, 200 ms delay, 1000 acq
	'1122_20171024_01_03_3751_OptoInhibOFF.dat'; ...
	% opto: 200 del, 400 dur, 3000 mV LED, ISI 500
	'1122_20171024_01_03_3751_OptoInhibON.dat'; ...	
	% 30 reps, 200 ms delay, 1000 acq
	'1122_20171024_01_03_3751_OptoInhibOFF_2.dat'; ...
	% opto: 100 del, 400 dur, 1000 mV LED, ISI 500	
	'1122_20171024_01_03_3751_OptoInhibON_2.dat'; ...
	% 30 reps, 250 ISI, 1000 acq
	'1122_20171024_02_03_3918_OptoInhibOFF.dat'; ...
	% 30 reps, opto: 100 del, 400 dur, 3000 mV, aborted
	'1122_20171024_02_03_3918_OptoInhibON.dat'; ...
};

%------------------------------------------------------------------------
%% settings for processing data
%------------------------------------------------------------------------
% bandpass filter frequencies
HPFreq = 250;
LPFreq = 6000;
% threshold
spikeThreshold = 'default';
binSize = 5;

%------------------------------------------------------------------------
%% plot data files 1, 2
%------------------------------------------------------------------------
datafile = file_list{1};
[D1, Dinf1, S1] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', [0 20]);
							
datafile = file_list{2};
[D2, Dinf2, S2] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', [0 20]);

%------------------------------------------------------------------------
%% plot data files 3, 4
%------------------------------------------------------------------------
datafile = file_list{3};
[D3, Dinf3, S3] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', [0 25]);
							
datafile = file_list{4};
[D4, Dinf4, S4] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', [0 25]);
							
							
							
							