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
	datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1225/20190115';
end
% datafile = '1225_20190115_03_02_894_BBN.dat';
datafile = '1225_20190115_02_01_744_FRA.dat';
%% Read Data (don't plot it though)
[D, Dinf, spikes, traces] = optoproc('file', fullfile(datapath, datafile), ...
													'channel', channelNumber);

save('tmpdata.mat', 'D', 'Dinf', 'spikes', 'traces', 'datapath', 'datafile');

%% Plot Level data
if strcmpi(Dinf.test.Type, 'LEVEL')
	plotRawLevelData(traces, Dinf);		
end

%% Plot freq + level data

if strcmpi(Dinf.test.Type, 'FREQ+LEVEL')
	
	% time vector for plotting
	t = (1000/Fs)*((1:length(D{1}.datatrace(:, 1))) - 1);
	% get overall max value, so all traces can be plotted on same scale
	levelData = cell(nlevels, 1);
	% initialize maxVals
	maxVals = zeros(nlevels, 1);
	for l = 1:nlevels
		dlist = stimindex{l};
		ntrials = length(dlist);
		levelData{l} = zeros(length(D{1}.datatrace(:, 1)), ntrials);
		for n = 1:ntrials
			levelData{l}(:, n) = filtfilt(filtB, filtA, ...
											D{dlist(n)}.datatrace(:, channelIndex));
		end
		maxVals(l) = max(max(levelData{l}));
	end

	maxVal = max(maxVals);
	% stacked data plot
	for l = 1:nlevels
		stackplot(t, levelData{l}, 'colormode', 'black', 'Ymax', maxVal);
		title(sprintf('Channel %d, Level %d', channelNumber, ...
									Dinf.test.stimcache.vrange(l)));
	end
	
	% overlay plot
	for l = 1:nlevels
		overlayplot(t, levelData{l}, 'Ymax', maxVal);
		title(sprintf('Channel %d, Level %d', channelNumber, ...
									Dinf.test.stimcache.vrange(l)));
	end
		
end



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
end
