%% settings for processing data
HPFreq = 350;
LPFreq = 6500;
 
%% Read Data
[D, Dinf, tracesByStim] = getFilteredOptoData('', [HPFreq LPFreq]);
if isempty(D)
	return
end

%% determine global RMS and max
nstim = length(tracesByStim);
ntrials = zeros(nstim, 1);
for s = 1:nstim
	ntrials(s) = length(tracesByStim{s});
end
if min(ntrials) ~- max(ntrials)
	error('unequal trials');
else
	ntrials = max(ntrials);
end

rmsvals = zeros(nstim, ntrials);
for s = 1:nstim
	for t = 1:ntrials
		rmsvals(s, t) = rms(tracesByStim{s}{t});
	end
end


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
if strcmpi(Dinf.test.Type, 'LEVEL')
    % time vector for plotting
    t = (1000/Fs)*((1:length(D{1}.datatrace(:, 1))) - 1);
%   for l = 1:nlevels
    for l = 1:nlevels
        dlist = stimindex{l};
        ntrials = length(dlist);
        tmpM = zeros(length(D{1}.datatrace(:, 1)), ntrials);
        for n = 1:ntrials
            tmpM(:, n) = filtfilt(filtB, filtA, ...
                                            D{dlist(n)}.datatrace(:, channelIndex));
        end
        stackplot(t, tmpM, 'colormode', 'black');
        title(sprintf('Channel %d, Level %d', channelNumber, ...
                                    Dinf.test.stimcache.vrange(l)));
    end
end
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
%           tracesByStim{w}(:, n) = D{dIndx}.datatrace(:, channelIndex);
        end
        stackplot(t, tracesByStim{w}, 'colormode', 'black');
        title({ datafile, sprintf('Stimulus: %s', wavlist{w})}, ...
                    'Interpreter', 'none');
        xlabel('ms')
        ylabel('Trial')
    end
end