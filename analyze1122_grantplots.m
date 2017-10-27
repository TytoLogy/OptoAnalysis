%------------------------------------------------------------------------
% analyze data from animal 1122, plots for grant
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
spikeThreshold = 3.1;
binSize = 5;
ylimits = [0 20];

%------------------------------------------------------------------------
%% plot all PSTH data for files 1, 2
%------------------------------------------------------------------------
datafile = file_list{1};
[D1, Dinf1, S1, T1] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'HPFreq', HPFreq, ...
								'LPFreq', LPFreq, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', ylimits);
							
datafile = file_list{2};
[D2, Dinf2, S2, T2] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'HPFreq', HPFreq, ...
								'LPFreq', LPFreq, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', ylimits);

%------------------------------------------------------------------------
%% plot all PSTH data for files 3, 4
%------------------------------------------------------------------------
datafile = file_list{3};
[D3, Dinf3, S3] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'HPFreq', HPFreq, ...
								'LPFreq', LPFreq, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', ylimits);
							
datafile = file_list{4};
[D4, Dinf4, S4] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'HPFreq', HPFreq, ...
								'LPFreq', LPFreq, ...
								'binsize', binSize, ...
								'threshold', spikeThreshold, ...
								'plotPSTH', ...
								'plotRowCols', [3 2], ...
								'yLimits', ylimits);

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot BBN data
%------------------------------------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%---------------------------------------------------------------------
%% threshold and binsize for grant plots
%---------------------------------------------------------------------
spikeThreshold = 5;
binSize = 20;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot BBN data for file 1
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf1;
T = T1;
[~, fname] = fileparts(file_list{1});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find BBN index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('BBN', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% determine global max
%---------------------------------------------------------------------
% # of reps
% find rms, max vals for each stim
netrmsvals = rms(T{v});
maxvals = max(abs(T{v}));
% compute overall mean rms for threshold
fprintf('Calculating mean and max RMS for data...\n');
mean_rms = mean(netrmsvals);
fprintf('\tMean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(maxvals);
fprintf('\tGlobal max abs value: %.4f\n', global_max);
%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------
% spiketimes = spikeschmitt2(T{v}', spikeThreshold*mean_rms, 1, Fs);
spiketimes = spikeschmitt2(T{v}', -0.35*global_max, 1, Fs);
% convert spike times in seconds to milliseconds
for r = 1:length(spiketimes)
	spiketimes{r} = (1000/Fs)*spiketimes{r};
end
%---------------------------------------------------------------------
% psth
%---------------------------------------------------------------------
[histdata.H, histdata.bins] = psth(spiketimes, binSize, ...
														[0 Dinf.test.AcqDuration]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% new figure
hF = figure(10);
% name figure
set(hF, 'Name', [fname '_' varlist{v}]);
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
grantplot(hF, t, T{v}, global_max, spiketimes, histdata);
% captions
subplot(1, 3, 1)
title({fname, varlist{v}}, 'Interpreter', 'none');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot BBN data for file 2
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf2;
T = T2;
[~, fname] = fileparts(file_list{2});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find BBN index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('BBN', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% determine global max
%---------------------------------------------------------------------
% find rms, max vals for each stim
netrmsvals = rms(T{v});
maxvals = max(abs(T{v}));
% compute overall mean rms for threshold
fprintf('Calculating mean and max RMS for data...\n');
mean_rms = mean(netrmsvals);
fprintf('\tMean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(maxvals);
fprintf('\tGlobal max abs value: %.4f\n', global_max);
%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------
% spiketimes = spikeschmitt2(T{v}', spikeThreshold*mean_rms, 1, Fs);
spiketimes = spikeschmitt2(T{v}', -0.35*global_max, 1, Fs);
% convert spike times in seconds to milliseconds
for r = 1:length(spiketimes)
	spiketimes{r} = (1000/Fs)*spiketimes{r};
end
%---------------------------------------------------------------------
% psth
%---------------------------------------------------------------------
[histdata.H, histdata.bins] = psth(spiketimes, binSize, ...
														[0 Dinf.test.AcqDuration]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% new figure
hF = figure(11);
% name figure
set(hF, 'Name', [fname '_' varlist{v}]);
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
grantplot(hF, t, T{v}, global_max, spiketimes, histdata);
% captions
subplot(1, 3, 1)
title({fname, varlist{v}}, 'Interpreter', 'none');





