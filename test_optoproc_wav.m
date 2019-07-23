%------------------------------------------------------------------------
% analyze data from animal 1108
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 June, 2016 (SJS) 
%	- adapted from viewOptoData.m
% 
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% setup
%------------------------------------------------------------------------
% add opto path
addpath('/Users/sshanbhag/Work/Code/Matlab/dev/TytoLogy/Experiments/Opto')

%------------------------------------------------------------------------
%% settings for processing data
%------------------------------------------------------------------------
% bandpass filter frequencies
HPFreq = 500;
LPFreq = 6000;

% PSTH bin width (ms)
psthBin = 10;

% Channel for data
channelNumber = 8;

%------------------------------------------------------------------------
%% Specify Data File Location/Name
%------------------------------------------------------------------------

% set paths, data file name
if ispc
% 	datapath = 'E:\Data\SJS\1058';
% 	datafile = '1058_20160623_0_02_1500_FREQ.dat';
% 	datapath = 'E:\Data\SJS\1012\20160727';
% 	datafile = '1012_20160727_5_3_1_OPTO.dat';
	datapath = 'E:\Data\SJS';
else 
	% assume this is a mac
% 	datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1012/20160727';
% 	datafile = '1012_20160727_5_3_1_OPTO.dat';
% 	datafile = '1012_20160727_4_3_1_LEVEL.dat';	
% 	datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1110/20170515';
% 	datafile = '1110_20170515_01_01_4038_WAV_OPTO.dat';
% 	datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1108/20170529';
	datapath = '/Users/sshanbhag/Desktop/20190722';
end
% datafile = '1225_20190115_03_02_894_BBN.dat';
datafile = '1323_20190722_01_02_508_WAV.dat';
%% Read Data (don't plot it though)
[D, Dinf, spikes, traces] = optoproc('file', fullfile(datapath, datafile), ...
													'channel', channelNumber);

% save('fratestdata.mat', 'D', 'Dinf', 'spikes', 'traces', 'datapath', 'datafile');

%% Plot Level data
% plotRawLevelData(<trace cell array>, <Dinf struct>) plot raw data traces!
% what a surprise!
if strcmpi(Dinf.test.Type, 'LEVEL')
	plotRawLevelData(traces, Dinf);		
end

%% Plot freq + level data as FRA
frawin = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
freqs = spikes.varlist{1};
levels = spikes.varlist{2};

FRA = computeFRA(spikes.spiketimes, freqs, levels, frawin);
FRA.fname = datafile;
plotFRA(FRA, 'dB', 'checker');

%% Plot raster and psth plot across levels

% set spike holdoff window (ms)
spwindow = 0.75;
tVal = -8; %#ok<NASGU>
% get threshold value
tVal = input('Enter Threshold value: ');
if isempty(tVal) || ~isnumeric(tVal)
	fprintf('Cancelled\n');
	return
end

% time window for counting spikes - use Delay to Delay+Duration interval
stimulusWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
spikeTimes = cell(nlevels, 1); 
for l = 1:nlevels
	spikeTimes{l} = spikeschmitt2(levelData{l}', tVal, spwindow, Dinf.indev.Fs);
	for n = 1:length(spikeTimes{l})
		spikeTimes{l}{n} = spikeTimes{l}{n}*(1000/Dinf.indev.Fs);
	end
	[H, bins] = psth(spikeTimes{l}, psthBin, ...
								[0 Dinf.test.AcqDuration], Dinf.indev.Fs);
	figure
	% rasters
	subplot(211)
	rasterplot(spikeTimes{l});
	title({datafile, sprintf('Channel %d, Level %d', ...
										channelNumber, ...
										Dinf.test.stimcache.vrange(l))}, ...
										'Interpreter', 'none');
	yminmax = ylim;
	line(stimulusWindow(1) * [1 1], yminmax);
	line(stimulusWindow(2) * [1 1], yminmax);

	% PSTH
	subplot(212)
	bar(bins, H, 1);
	xlim([0 Dinf.test.AcqDuration])
	xlabel('Time (ms)')
	ylabel('spikes/bin');
	yminmax = ylim;
	line(stimulusWindow(1) * [1 1], yminmax);
	line(stimulusWindow(2) * [1 1], yminmax);

end


%% rate-level function

% do some posthoc changes on Dinf
if isfield(Dinf.test, 'Name')
	Dinf.test.Name = char(Dinf.test.Name);
end

% time window for counting spikes - use Delay to Delay+Duration interval
analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
% check for level and spiketime data
if strcmpi(Dinf.test.Type, 'LEVEL') && exist('spikeTimes', 'var')
	RLF = computeRLF(spikeTimes, Dinf.audio.Level, analysisWindow);
	figure;
	plotRLF(RLF, 'median');
	% build title string
	tStr = {sprintf('RLF [%d %d] dB SPL', min(RLF.levels), max(RLF.levels)), ...
				[datafile ', ' sprintf('Channel %d', channelNumber)]};
	title(tStr, 'Interpreter', 'none');
end


%%
if strcmpi(Dinf.test.Type, 'WavFile')
	% time vector for plotting - assume all traces are equal length
	t = (1000/Fs)*((1:length(D{1}.datatrace(:, 1))) - 1);
	tracesByStim = cell(nwavs, 1);
	for w = 1:nwavs
		% create temporary array to hold data
		tracesByStim{w} = zeros(length(D{1}.datatrace(:, 1)), Dinf.test.Reps);
		for n = 1:Dinf.test.Reps
			dIndx = stimindex{w}(n);
			tracesByStim{w}(:, n) = filtfilt(filtB, filtA, ...
													D{dIndx}.datatrace(:, channelIndex));
		end
		stackplot(t, tracesByStim{w}, 'colormode', 'black');
		title({	datafile, sprintf('Stimulus: %s', wavlist{w})}, ...
					'Interpreter', 'none');
		xlabel('ms')
		ylabel('Trial')
	end
=======
%% load data
%------------------------------------------------------------------------
dname = '/Users/sshanbhag/Work/Data/Mouse/IC/1302/20190507';
fname = '1302_20190507_03_03_711.9_WAV.dat';

%		optoproc('file', fullfile(dname, fname), 'plotPSTH');

if ~exist('D', 'var')
	if ~exist('wavproc.mat', 'file')
		% optoproc has been modified temporarily to save mat file of 
		% data/info in 'wavproc.mat' when plotPSTH is specified
		optoproc('file', fullfile(dname, fname), 'plotPSTH');
	end
 	load('wavproc.mat')
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% -----------rest is optoproc_plotPSTH_WAV code 
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%  H = optoproc_plotPSTH_WAV(spiketimes, Dinf, binSize, nvars, varlist, timeLimits, ...
%												yLimits, titleString)
% note nvars, varlist might be redundant!
% * change spiketimes to spikesByStim?
% title string not necessary?
%------------------------------------------------------------------------
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%% first, get information about stimuli
%------------------------------------------------------------------------

% get unique stimuli in order they appear in stimList using 'stable' option
% iA will be indices of first occurrence of each unique stim
% iC will identify which of the unique stim is in each row
[uniqueStim, iA, iC]  = unique(Dinf.test.wavlist, 'stable');
nStim = numel(uniqueStim);
% create list of indices into stimList, spiketimes for each stimulus 
% (will use later for many things)
stimIndices = cell(nStim, 1);
for s = 1:nStim
	stimIndices{s} = find(iC == s);
end

% determine if there is a NULL stimulus
if any(strcmpi('NULL', uniqueStim))
	hasNULL = true;
	nullIndex = find(strcmpi('NULL', uniqueStim));
else
	hasNULL = false;
end
%------------------------------------------------------------------------
%% figure out attenuation levels
%------------------------------------------------------------------------
% could do this from:
% a) Dinf.test.Level
% 	BUT, this may not apply to all stimuli (e.g., 'NULL')
% b) look at individual stimList entries:Dinf.stimList(n).audio.Level
% 	should be one entry per level/stimulus combination
% 	laborious?

%------------------------------------------------------------------------
% first way
%------------------------------------------------------------------------
% dbLevels = Dinf.test.Level;

%------------------------------------------------------------------------
% second way
%------------------------------------------------------------------------
% create cell array to store levels for each stim (might be different
% number of levels for each stim, so can't use array)
dbLevelsByStim = cell(nStim, 1);
nLevels = zeros(nStim, 1);
% loop through stim
fprintf('Determining Stimulus Levels...\n');
for s = 1:nStim
	% ASSUME that repeats of same stimulus are at different levels...
	nLevels(s) = length(stimIndices{s});
	dbLevelsByStim{s} = zeros(nLevels(s), 1);
	fprintf('\tStimulus: %s\n', uniqueStim{s});
	fprintf('\t\tnLevels: %d\n', nLevels(s));
	fprintf('\t\tLevels: ');
	for n = 1:nLevels(s)
		dbLevelsByStim{s}(n) = Dinf.stimList(stimIndices{s}(n)).audio.Level;
		fprintf('%d  ', Dinf.stimList(stimIndices{s}(n)).audio.Level);
	end
	fprintf('\n');
end

%------------------------------------------------------------------------
%% setup plots / options
%------------------------------------------------------------------------

% create array to hold figure handles
hPR = cell(nUnique, 1);

% global options for raster and psth matrix
plotopts.timelimits = timeLimits;
% assign y axis limits if provided
if ~isempty(yLimits)
	plotopts.psth_ylimits = yLimits;
end

% options for raster and psth matrix
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 16;
plotopts.raster_color = [0 0 0];
plotopts.psth_binwidth = binSize;
plotopts.plotgap = 0.001;
plotopts.xlabel = 'msec';
plotopts.stimulus_times_plot = 3;
plotopts.stimulus_on_color{1} = [0 0 1];
plotopts.stimulus_off_color{1} = [0 0 1];
plotopts.stimulus_onoff_pct(1) = 60;
% add on off bars for opto stim
if Dinf.opto.Enable
	% add colors for second stimulus
	plotopts.stimulus_on_color{2} = [1 0 0];
	plotopts.stimulus_off_color{2} = [1 0 0];
	plotopts.stimulus_onoff_pct(2) = 80;
end


% plot name
[fpath, fname, fext] = fileparts(strrep(Dinf.filename, '\', '/'));
fname = [fname fext];

% titles for stimuli
varlist = Dinf.test.wavlist;
nvars = length(varlist);
titleString = cell(nvars, 1);
for v = 1:nvars
	if v == 1 
		titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
	else
		titleString{v} = sprintf('wav name: %s', varlist{v});
	end
end



for s = 1:nStim	
	% check on levels and stimIndices
	if nLevels(s) ~= length(stimIndices{s})
		error('%s: mismatch in nLevels and stimIndices', mfilename);
	end
	
	% create list of stimulus times
	plotopts.stimulus_times = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		% need to have [stim_onset stim_offset], so add delay to 
		% [0 duration] to compute proper times. then, multiply by 0.001 to
		% give times in seconds (Dinf values are in milliseconds)
		plotopts.stimulus_times{s, 1} = 0.001 * (Dinf.audio.Delay + ...
															[0 Dinf.audio.Duration]);
		% if opto is Enabled, add it to the array by concatenation
		if Dinf.opto.Enable
			plotopts.stimulus_times{s, 1} = [stimulus_times{s, 1}; ...
												 0.001 * (Dinf.opto.Delay + ...
															[0 Dinf.opto.Dur]) ];
		end
	end
	% create plot titles
	plotopts.plot_titles = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		if l == 1
			plotopts.plot_titles{l} = {fname, sprintf('%s', uniqueStim{s})};
		else
			plotopts.plot_titles{l} = '';
		end
	end

	% assign figure index
	hPR{s} = figure;
	
	% assign spikes
	Spikes = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		Spikes{l} = spiketimes{stimIndices{s}(l)};
	end
	rasterpsthmatrix(Spikes, plotopts);
end
