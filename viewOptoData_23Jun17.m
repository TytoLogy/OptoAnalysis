function [data, datainfo] = viewOptoData(varargin)
%------------------------------------------------------------------------
% [data, datainfo] = viewOptoData(varargin)
%------------------------------------------------------------------------
% % TytoLogy:Experiments:opto Application
%--------------------------------------------------------------------------
% Reads binary data file created by the opto program, plots traces
%
% If a datafile name is provided in varargin (e.g.,
% viewOptoData('c:\mydir\mynicedata.dat'), the program will attempt to
% read from that file. 
%
% Otherwise it will open a dialog window for the user
% to select the data (.dat) file.
%
%------------------------------------------------------------------------
% Output Arguments:
%
% data          contains the read data in a cell structure array.
% datainfo      has the file header information.
%
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 June, 2016 (SJS)
%           - adapted from readHPData.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%   *Documentation!
%--------------------------------------------------------------------------
%% settings for processing data
HPFreq = 350;
LPFreq = 6500;
nargin
 
%% Read Data
% set paths
if ispc
%   datapath = 'E:\Data\SJS\1058';
%   datafile = '1058_20160623_0_02_1500_FREQ.dat';
%   datapath = 'E:\Data\SJS\1012\20160727';
%   datafile = '1012_20160727_5_3_1_OPTO.dat';
  datapath = 'E:\Data\SJS';
%     datapath = 'C:\Users\cgrimsley\Desktop\1116\20170618';
    datafile = '';
else
%   datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1012/20160727';
%   datafile = '1012_20160727_5_3_1_OPTO.dat';
%   datafile = '1012_20160727_4_3_1_LEVEL.dat';
%   datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1110/20170515';
%   datafile = '1110_20170515_01_01_4038_WAV_OPTO.dat';
    datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto';
    datafile = '';
end
if nargin
    [datapath, datafile, dataext] = fileparts(varargin{1});
    datafile = [datafile dataext];
end
if isempty(datafile)
    % get data file from user
    [datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
                                                    fullfile(datapath, datafile));
    if datafile == 0
        return
    end
end
% read in raw data
[D, Dinf] = readOptoData(fullfile(datapath, datafile));
% read in test data (if it exists)
[~, fbase, ~] = fileparts(datafile);
testfile = [fbase '_testdata.mat'];
if exist(fullfile(datapath, testfile), 'file')
    load(fullfile(datapath, testfile), 'testdata');
else
    testdata = [];
end
%%
% read in wav info (if it exists)
wavinfofile = [fbase '_wavinfo.mat'];
if exist(fullfile(datapath, wavinfofile), 'file')
%   load(fullfile(datapath, wavinfofile), 'stimList');
    load(fullfile(datapath, wavinfofile));
else
    stimList = [];
end
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
    nchan = Dinf.channels.nRecordChannels; %#ok<NASGU>
    channelList = Dinf.channels.RecordChannelList;
else
    nchan = Dinf.channels.nInputChannels; %#ok<NASGU>
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
%% work from testdata
% if no testdata, get outta here
if isempty(testdata)
    return
else
    spikesByStim = testdata.SpikeTimes;
    PSTHByStim = testdata.PSTH;
    if length(spikesByStim) ~= nwavs
        error('%s: mismatch in # of stimuli in testdata and # of wavs', mfilename)
    end
end
% Plot data
if strcmpi(Dinf.test.Type, 'WavFile')
    winfoIndx = 0;
    for w = 1:nwavs
        figure
        subplot(211)
        title({ datafile, sprintf('Stimulus: %s', wavlist{w})}, ...
                    'Interpreter', 'none');
        rasterplot(spikesByStim{w}, [0 1000]);
        subplot(212)
        bar(PSTHByStim.bins, PSTHByStim.hvals{w}, 1);
        xlabel('ms')
       xlim([0 1000])
        % draw line for stimulus onset/offset
        onset = [];
        offset = [];
        % stimulus lines depend on stimulus type
        if strcmpi(stimList(w).audio.signal.Type, 'noise')
            onset = stimList(w).audio.Delay;
            offset = onset + stimList(w).audio.Duration;
        elseif strcmpi(stimList(w).audio.signal.Type, 'wav')
            winfoIndx = winfoIndx + 1;
            onset = stimList(w).audio.Delay;
            offset = onset + 1000*wavInfo(winfoIndx).Duration;
        elseif strcmpi(stimList(w).audio.signal.Type, 'null')
        else
            % WTF!!!???
            error('WTF?')
        end
        if ~isempty(onset)
            text(onset, 0, ':', ...
                                'Color', 'r', ...
                                'FontSize', 12, ...
                                'FontWeight', 'bold');
            text(offset, 0, '|', ...
                                'Color', 'r', ...
                                'FontSize', 12, ...
                                'FontWeight', 'bold');
        end
        % draw line for opto stim
        if stimList(w).opto.Enable
            onset = stimList(w).opto.Delay;
            offset = onset + stimList(w).opto.Dur;
            text(onset, 0, ':', ...
                                'Color', 'g', ...
                                'FontSize', 12, ...
                                'FontWeight', 'bold');
            text(offset, 0, '|',  ...
                                'Color', 'g', ...
                                'FontSize', 12, ...
                                'FontWeight', 'bold');
        end
    end
end
%% assign outputs
data = D;
datainfo = Dinf;
 