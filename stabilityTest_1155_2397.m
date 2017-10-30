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
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% get file information 
% defines animalID, dateID, unit
%	datafiles in arrays, structs:
%		optoRLF
% 		bbnRLF
% 		bbn
% 		optoSlidingWin2K
% 		optoSlidingWin3K
% 		freqtuning
% 		optoFreqTuning
% 		optoWav		
%---------------------------------------------------------------------
files_for_1155_2397;

%% select data files

% datafiles{1} = '1155_20171025_02_02_2397_BBN.dat'
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: OFF
datafiles{1} = '1155_20171025_02_02_2397_BBN.dat';

% datafiles{2} = '1155_20171025_02_02_2397_BBN_optoON.dat'
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, 3500 or 6.8 mW
datafiles{2} = '1155_20171025_02_02_2397_BBN_optoON.dat';

% datafiles{2} = '1155_20171025_02_02_2397_BBN_optoON_9.dat'
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, 3000mV or 6.1 mW
datafiles{3} = '1155_20171025_02_02_2397_BBN_optoON_9.dat';


% datafiles{4} = '1155_20171025_02_02_2397_OptoInhibOFF.dat'
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
datafiles{4} = '1155_20171025_02_02_2397_OptoInhibOFF.dat';

% datafiles{4} = '1155_20171025_02_02_2397_OptoInhibON.dat'
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		250 ms duration, 150 ms delay, 3000mV or 6.1 mW
datafiles{5} = '1155_20171025_02_02_2397_OptoInhibON.dat';


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Settings for Processing data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% get indices to sort by opto intensity
[sortedIntensity, sortedIndices] = sort(optoRLF.LEDintensity);
% data window (absolute time re: sweep onset)
AnalysisWindow_BBN = [150 200];
AnalysisWindow_optoWav = [200 250];
% psth bin size
binSize = 5;
% threshold (RMS)
spikeThreshold = 3;

%---------------------------------------------------------------------
%% allocate data in sorted intensity order
%---------------------------------------------------------------------
% 	Dinf = data information struct
% 	S = spikes struct
% 		spiketimes: {ntrials×1 cell}
% 		mean_rms: overall mean rms value
% 		global_max: overal max value
% 		Threshold: threshold used
%	F = RLF function data struct
% 		nlevels	# levels
% 		level		db SPL levels [nlevels 1]
% 		mean		mean value for each level [nlevels 1]
% 		std		std deviation value for each level [nlevels 1]
data = struct(	'datafile', datafiles, ...
						'D', [], ...
						'Dinf', [], ...
						'S', [], ...
						'T', [], ...
						'F', []	);
					
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Load & get spike count data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
for f = 1:length(datafiles)
	
	[data(f).D, data(f).Dinf, data(f).S, data(f).T] = ...
				optoproc('file', fullfile(datapath, data(f).datafile), ...
							'channel', channelNumber, ...
							'binsize', binSize, ...
							'Threshold', spikeThreshold);
						
	% compute RLF (or mean response rate for non RLF data)
	if strcmpi(data(f).Dinf.test.Type, 'LEVEL')
		data(f).F = computeRLF(data(f).S.spiketimes, AnalysisWindow_BBN);
		data(f).F.stimuli = char(data(1).Dinf.test.Name);
	else
		data(f).F = computeRLF(data(f).S.spiketimes, AnalysisWindow_optoWav);
		data(f).F.stimuli = data(f).Dinf.test.wavlist;
	end
	data(f).F.nlevels = length(data(f).Dinf.audio.Level);
	data(f).F.level = data(f).Dinf.audio.Level;
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Process data
%---------------------------------------------------------------------
%---------------------------------------------------------------------

% locations of level or stimuli in F arrays (60 db for BBN rlfs, 65 for
% opto files)
dLocations = [7 7 7 2 2];
for f = 1:3
	d = dLocations(f);
	fprintf('Datafile: %s\n', data(f).datafile);
	tstart = data(f).Dinf.time_start;
	tend = data(f).Dinf.time_end;
	telapsed = tend - tstart;
	fprintf('\tStart, end, elapsed time (mm:ss.xxx): %s\t%s\t%s\n', ...
															datestr(tstart, 'HH:MM:SS.FFF'), ...
															datestr(tend, 'HH:MM:SS.FFF'), ... 
															datestr(telapsed, 'MM:SS.FFF'));
	F = data(f).F;
	fprintf('\tMean ± sd: %.2f ± %.2f\n', F.mean(d), F.std(d));
	fprintf('\tMedian [ci]: %.2f [%.2f %.2f]\n', F.median(d), ...
																F.median_ci{d}(1), ...
																F.median_ci{d}(2));
	fprintf('\n');
	data(f).F.times = struct(	'tstart', tstart, ...
										'tend', tend, ...
										'telapsed', telapsed);
end
% 
% opto1.mean = data(2).F.mean(7);
% opto1.sd = data(2).F.std(7);
% opto1.spikeCount = data(2).F.spikeCount{7};
% 
% ctrl2.mean = data(3).F.mean(2);
% ctrl2.sd = data(3).F.std(2);
% ctrl2.spikeCount = data(3).F.spikeCount{2};
% 
% opto2.mean = data(4).F.mean(2);
% opto2.sd = data(4).F.std(2);
% opto2.spikeCount = data(4).F.spikeCount{2};
% 
% 
