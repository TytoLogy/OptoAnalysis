% stability check

%-----------------------------------------------------------------------------
%{
idea: for each stimulus/intensity combination, take first n trials, and
last m trials.

Compare firing rate of two samples.
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
files_for_1155_2616;

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
% out for analysis
analysispath = fullfile(outpath_base, animalID, dateID, ...
									'/OptoAnalysisPlots/Stability');
if ~exist(analysispath, 'dir')
	mkdir(analysispath);
end

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
							'Threshold', spikeThreshold, ...
							'plotPath', analysispath, ...
							'plotPSTH', ...
							'savePDF');
	
% 	[data(f).D, data(f).Dinf, data(f).S, data(f).T] = ...
% 				optoproc('file', fullfile(datapath, data(f).datafile), ...
% 							'channel', channelNumber, ...
% 							'binsize', binSize, ...
% 							'Threshold', spikeThreshold, ...
% 							'plotPath', analysispath, ...
% 							'plotPSTH', ...
% 							'savePDF');
						
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

%% look at rates of three BBN data with/without opto inact
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

%% plot data
title_id = [animalID '-' dateID '-' unit '-' penetration '-' depth];
figure(1)
xdata = [1 2 3];
xnames = {	['0 mW, ' datestr(data(1).F.times.tstart, 'HH:MM')]; ...
				['6.8 mW, ' datestr(data(2).F.times.tstart, 'HH:MM')]; ... 
				['6.1 mW, ' datestr(data(3).F.times.tstart, 'HH:MM')]	};
ydata = zeros(1, 3);
yerr = ydata;
for f = 1:3
	ydata(f) = data(f).F.mean(dLocations(f));
	yerr(f) =  data(f).F.std(dLocations(f));
end
bar(xdata, ydata)
hold on
	errorbar(xdata, ydata, yerr, '.');
hold off
set(gca, 'XTickLabel', xnames)
title({title_id, 'Avg # Spikes'});
xlabel('Light level, Exp. Time');
ylabel('Avg Spikes ± s.d.');

%% bootstrap estimates
nboot = 1000;
boot_values = cell(3, 1);
boot_samples = cell(3, 1);
mean_ks = cell(3, 1);
std_ks = cell(3, 1);
mstats = zeros(3, 2);
figure(2)
for f = 1:3
	bdata = data(f).F.spikeCount{dLocations(f)};
	% bootstrapped estimates of mean and std deviation
	[boot_values{f}, boot_samples{f}] = bootstrp(	nboot, ...
																	@(x)[mean(x) std(x)], ...
																	bdata);
	% get probability density estimates
	% mean
	mean_ks{f} = struct('prob', [], 'values', []);
	[mean_ks{f}.prob, mean_ks{f}.values] = ksdensity(boot_values{f}(:, 1));
	% find peak of probability, use maxindx to get the max value
	[maxprob, maxindx] = max(mean_ks{f}.prob);
	maxval = mean_ks{f}.values(maxindx);
	% plot ks kurve
	subplot(2, 3, f)
	plot(mean_ks{f}.values, mean_ks{f}.prob);
	if f == 1
		xlabel('avg spikes');
		ylabel('Mean: P(x)');
	elseif f == 2
		title({title_id, 'Bootstrapped Estimates; pdfs'});
	end
	% place point at max on curve
	textoffset = 0.1 * diff(xlim);
	hold on
		plot(maxval, maxprob, '.r')
		text(maxval+textoffset, maxprob, sprintf('%.4f', maxval));
	hold off
	mstats(f, 1) = maxval;
	
	% std
	std_ks{f} = struct('prob', [], 'values', []);
	[std_ks{f}.prob, std_ks{f}.values] = ksdensity(boot_values{f}(:, 2));
		% find peak of probability, use maxindx to get the max value
	[maxprob, maxindx] = max(std_ks{f}.prob);
	maxval = std_ks{f}.values(maxindx);
	% plot ks kurve
	subplot(2, 3, f+3)
	plot(std_ks{f}.values, std_ks{f}.prob);
	if f == 1
		xlabel('stdev spikes');
		ylabel('StdDev: P(x)');
	end
	% place point at max on curve
	textoffset = 0.1 * diff(xlim);
	hold on
		plot(maxval, maxprob, '.r')
		text(maxval+textoffset, maxprob, sprintf('%.4f', maxval));
	hold off
	mstats(f, 2) = maxval;
end

%%
figure(3)
y2data = zeros(1, 3);
y2err = y2data;
for f = 1:3
	y2data(f) = mstats(f, 1);
	y2err(f) =  mstats(f, 2);
end
bar(xdata, y2data)
hold on
	errorbar(xdata, y2data, y2err, '.');
hold off
set(gca, 'XTickLabel', xnames)
title({title_id, 'Bootstrapped Avg # Spikes'});
xlabel('Light level, Exp. Time');
ylabel('Avg Spikes ± s.d.');




