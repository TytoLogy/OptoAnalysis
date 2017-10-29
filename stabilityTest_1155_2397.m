% stability check

%-----------------------------------------------------------------------------
%{
idea: for each stimulus/intensity combination, take first n trials, and
last m trials.

Compare firing rate of two samples.



1155_20171025_02_02_2397_BBN
1155_20171025_02_02_2397_BBN_optoON_9 % 3000 mV

%}
%-----------------------------------------------------------------------------

%---------------------------------------------------------------------
% settings and files
%---------------------------------------------------------------------
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171025';
% unit #
unit = '02';
% penetration
penetration = '02';
% recording depth
depth = '2397';
% channel
channelNumber = 8;

% for all files:
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, varied command voltage to ThorLabs LED
optoRLF.files = {	...
	'1155_20171025_02_02_2397_BBN.dat', ...					% rlf, 0:10:60
	'1155_20171025_02_02_2397_BBN_optoON.dat', ...			% rlf, 0:10:60, opto 3500 mV
	'1155_20171025_02_02_2397_BBN_optoON_2.dat', ...		% rlf, 0:10:60, opto 500 mV
	'1155_20171025_02_02_2397_BBN_optoON_3.dat', ...		% rlf, 0:10:60, opto 250 mV
	'1155_20171025_02_02_2397_BBN_optoON_4.dat', ...		% rlf, 0:10:60, opto 750 mV
	'1155_20171025_02_02_2397_BBN_optoON_5.dat', ...		% rlf, 0:10:60, opto 1000 mV
	'1155_20171025_02_02_2397_BBN_optoON_6.dat', ...		% rlf, 0:10:60, opto 1500 mV
	'1155_20171025_02_02_2397_BBN_optoON_7.dat', ...		% rlf, 0:10:60, opto 2000 mV
	'1155_20171025_02_02_2397_BBN_optoON_8.dat', ...		% rlf, 0:10:60, opto 2500 mV
	'1155_20171025_02_02_2397_BBN_optoON_9.dat' ...			% rlf, 0:10:60, opto 3000 mV
};
optoRLF.LEDintensity = [ ...
		0, ...
		3500, ...
		500, ...
		250, ...
		750, ...
		1000, ...
		1500, ...
		2000, ...
		2500, ...
		3000 ...
];
optoRLF.LEDpower = [ ...
	0.0, ...
	6.8, ...
	1.4, ... 
	0.7, ... 
	2.0, ...
	2.5, ...
	3.6, ...
	4.5, ...
	5.3, ...
	6.1 ...
];

%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Process RLF data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% # of data files
nfiles = length(optoRLF.files);
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% get indices to sort by opto intensity
[sortedIntensity, sortedIndices] = sort(optoRLF.LEDintensity);
% data window (absolute time re: sweep onset)
AnalysisWindow = [150 200];
% psth bin size
binSize = 5;
% threshold (RMS)
spikeThreshold = 3;

