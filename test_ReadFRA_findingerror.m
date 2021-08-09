%------------------------------------------------------------------------
% test_optoproc.m
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 Aug 2021 (SJS)
%	- adapted from test_readOptoData.m
% 
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% setup
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% settings for processing data
%------------------------------------------------------------------------
% bandpass filter frequencies
HPFreq = 500;
LPFreq = 6000;

% PSTH bin width (ms)
psthBin = 10;


%------------------------------------------------------------------------
%% Specify Data File Location/Name
%------------------------------------------------------------------------
% filebase = '1395_20200204_01_02_3600';
% filebase = '1429_20200707_01_01_2942';
filebase = '1395_20200204_01_01_3650';
% filebase = '1395_20200204_01_02_3600';
database = '/media/Data/NeuroData/Mouse/Raw';

oname = OptoFileName(filebase);

% build FRA filename
datafile = [oname.fileWithoutOther '_FRA.dat'];
% build path to FRA file
datapath = fullfile(database, oname.animal, oname.datecode);

infile = fullfile(datapath, datafile);

% Channel for data
channelNumber = 5;

%% Read in raw Data 

[d, dinf] = readOptoData(infile);

%%  read in filtered DAta

[d, dinf, tr] = getFilteredOptoData( infile, ...
											'Filter', [HPFreq LPFreq], ...
											'Channel', channelNumber);

%% call optoproc to read in data (don't plot it though)
if ~isempty(infile)
   [D, Dinf, spikes, traces] = ...
                        optoproc('file', infile, ...
													'channel', channelNumber);
else
   [D, Dinf, spikes, traces] = ...
                        optoproc('channel', channelNumber);
end

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
end
