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
outpath_base = fullfile(data_root_path, 'Analyzed');

%------------------------------------------------------------------------
%% settings for processing data
%------------------------------------------------------------------------
% bandpass filter frequencies
HPFreq = 500;
LPFreq = 6000;

%------------------------------------------------------------------------
%% Read Data
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
	'1122_20171024_01_03_3751_OptoInhibOFF.dat', ...
	'1122_20171024_01_03_3751_OptoInhibOFF_2.dat', ...
	'1122_20171024_01_03_3751_OptoInhibON.dat', ...
	'1122_20171024_01_03_3751_OptoInhibON_2.dat', ...
	'1122_20171024_02_03_3918_OptoInhibOFF.dat', ...
	'1122_20171024_02_03_3918_OptoInhibON.dat' ...
};

datafile = file_list{1};

% read in data
[D, Dinf] = readOptoData(fullfile(datapath, datafile));

%% define filter for data
% sampling rate
Fs = Dinf.indev.Fs;
% build bandpass filter, store coefficients in filtB, filtA
fband = [HPFreq LPFreq] ./ (0.5 * Fs);
[filtB, filtA] = butter(5, fband);

%% Get test info

% try to get information from test Type
if isfield(Dinf.test, 'Type')
	% convert ascii characters from binary file 
	Dinf.test.Type = char(Dinf.test.Type);
	fprintf('Test type: %s\n', Dinf.test.Type);
else
	% otherwise, need to find a different way
	if isfield(Dinf.test, 'optovar_name')
		Dinf.test.optovar_name = char(Dinf.test.optovar_name);
	end
	if isfield(Dinf.test, 'audiovar_name')
		Dinf.test.audiovar_name = char(Dinf.test.audiovar_name);
		if strcmpi(Dinf.test.audiovar_name, 'WAVFILE')
			% test is WAVfile
			Dinf.test.Type = Dinf.test.audiovar_name;
		end
	end
	fprintf('Test type: %s\n', Dinf.test.Type);
end

%%
% Some test-specific things...

% for FREQ test, find indices of stimuli with same frequency
if isnumeric(Dinf.test.Type)
	Dinf.test.Type = char(Dinf.test.Type);
end	
switch upper(Dinf.test.Type)
	case 'FREQ'
		% list of frequencies, and # of freqs tested
		freqlist = cell2mat(Dinf.test.stimcache.FREQ);
		nfreqs = length(Dinf.test.stimcache.vrange);
		% locate where trials for each frequency are located in the 
		% stimulus cache list - this will be used to pull out trials of
		% same frequency
		stimindex = cell(nfreqs, 1);
		for f = 1:nfreqs
			stimindex{f} = find(Dinf.test.stimcache.vrange(f) == freqlist);
		end
		
% for LEVEL test, find indices of stimuli with same level (dB SPL)
	case 'LEVEL'
		% list of legvels, and # of levels tested
		levellist = Dinf.test.stimcache.LEVEL;
		nlevels = length(Dinf.test.stimcache.vrange);
		% locate where trials for each frequency are located in the 
		% stimulus cache list - this will be used to pull out trials of
		% same frequency
		stimindex = cell(nlevels, 1);
		for l = 1:nlevels
			stimindex{l} = find(Dinf.test.stimcache.vrange(l) == levellist);
		end

% for OPTO test...
	case 'OPTO'
	
	% for WavFile, need to find indices with same filename.
	case 'WAVFILE'
		% get list of stimuli (wav file names)
		nwavs = length(Dinf.stimList);
		wavlist = cell(nwavs, 1);
		stimindex = cell(nwavs, 1);
		for w = 1:nwavs
			stype = Dinf.stimList(w).audio.signal.Type;
			if strcmpi(stype, 'null')
				wavlist{w} = 'null';
			elseif strcmpi(stype, 'noise')
				wavlist{w} = 'BBN';
			elseif strcmpi(stype, 'wav')
				[~, wavlist{w}] = fileparts(Dinf.stimList(w).audio.signal.WavFile);
			else
				error('%s: unknown type %s', mfilename, stype);
			end
			stimindex{w} = find(Dinf.test.stimIndices == w);
		end

	otherwise
		error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end


%% Pull out trials, apply filter, store in matrix
if isfield(Dinf.channels, 'nRecordChannels')
	nchan = Dinf.channels.nRecordChannels;
	channelList = Dinf.channels.RecordChannelList;
else
	nchan = Dinf.channels.nInputChannels;
	channelList = Dinf.channels.InputChannels;
end


%% find channel data
channelNumber = 8;
channelIndex = find(channelList == channelNumber);
if isempty(channelIndex)
	error('Channel not recorded')
end

%% Plot data for one channel, process will vary depending on stimulus type

if strcmpi(Dinf.test.Type, 'FREQ')
	% time vector for plotting
	t = (1000/Fs)*((1:length(D{1}.datatrace(:, 1))) - 1);
	for f = 1:nfreqs
		dlist = stimindex{f};
		ntrials = length(dlist);
		tmpM = zeros(length(D{1}.datatrace(:, 1)), ntrials);
		for n = 1:ntrials
			tmpM(:, n) = filtfilt(filtB, filtA, ...
											D{dlist(n)}.datatrace(:, channelIndex));
		end
		stackplot(t, tmpM, 'colormode', 'black');
		title(sprintf('Channel %d, Freq %d', channelNumber, ...
									Dinf.test.stimcache.vrange(f)));
	end
end

%% Plot Level data
if strcmpi(Dinf.test.Type, 'LEVEL')
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


%% rasters, psth

% set spike holdoff window (ms)
spwindow = 0.75;
tVal = -8;
% % get threshold value
% tVal = input('Enter Threshold value: ');
% if isempty(tVal) || ~isnumeric(tVal)
% 	fprintf('Cancelled\n');
% 	return
% end

%%
%{
if strcmpi(Dinf.test.Type, 'LEVEL')
	spikeTimes = cell(nlevels, 1);
	for l = 1:nlevels
		spikeTimes{l} = spikeschmitt2(levelData{l}', tVal, spwindow, Dinf.indev.Fs);
		for n = 1:length(spikeTimes{l})
			spikeTimes{l}{n} = spikeTimes{l}{n}*(1000/Dinf.indev.Fs);
		end
		[H, bins] = psth(spikeTimes{l}, 10, [0 1000], Dinf.indev.Fs);
		figure
		subplot(211)
		rasterplot(spikeTimes{l});
		title({datafile, sprintf('Channel %d, Level %d', channelNumber, ...
										Dinf.test.stimcache.vrange(l))});
		subplot(212)
		bar(bins, H, 1);
		xlim([0 1000])
	end
end
%}



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

if strcmpi(Dinf.test.Type, 'WavFile')
	spikeTimes = cell(nwavs, 1);
	for l = 1:nwavs
		spikeTimes{l} = spikeschmitt2(levelData{l}', tVal, spwindow, Dinf.indev.Fs);
		for n = 1:length(spikeTimes{l})
			spikeTimes{l}{n} = spikeTimes{l}{n}*(1000/Dinf.indev.Fs);
		end
		[H, bins] = psth(spikeTimes{l}, 10, [0 1000], Dinf.indev.Fs);
		figure
		subplot(211)
		rasterplot(spikeTimes{l});
		title({datafile, sprintf('Channel %d, Level %d', channelNumber, ...
										Dinf.test.stimcache.vrange(l))});
		subplot(212)
		bar(bins, H, 1);
		xlim([0 1000])
	end
end
