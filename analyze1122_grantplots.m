%------------------------------------------------------------------------
% analyze data from animal 1122, plots for grant
%	1122 has optical cannuli implanted over left hemisphere MG
%	data were recorded from BLA in left hemisphere using tungsten electrode
%  animal under urethane anesthesia, head fixed
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
% path for plots
plotpath = fullfile(plotpath_base, animalID, dateID, 'GrantFigs');
if ~exist(plotpath, 'dir')
	mkdir(plotpath);
end

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
% settings for grant plots
%---------------------------------------------------------------------
spikeThreshold = -0.35;
binSize = 20;
psthlim = [0 15];
% initialize stimulus onset offset line structs
sound_onoff = struct('onset', 0, 'offset', 1, 'color', 'b');
opto_onoff = struct('onset', 0, 'offset', 1, 'color', 'g');
% subplots?
subplotONOFF = 0;
% save figure?
saveFIG = 1;

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
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end

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
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% *******************LFH**********************
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%---------------------------------------------------------------------
% settings for grant plots
%---------------------------------------------------------------------
spikeThreshold = -0.35;
binSize = 20;
psthlim = [0 20];
% initialize stimulus onset offset line structs
sound_onoff = struct('onset', 0, 'offset', 1, 'color', 'b');
opto_onoff = struct('onset', 0, 'offset', 1, 'color', 'g');
% subplots?
subplotONOFF = 0;
% save figure?
saveFIG = 1;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot LFH data for file 1
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf1;
T = T1;
[~, fname] = fileparts(file_list{1});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find LFH index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('P100_9_LFH', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot LFH data for file 2
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf2;
T = T2;
[~, fname] = fileparts(file_list{2});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find LFH index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('P100_9_LFH', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% *******************NULL**********************
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%---------------------------------------------------------------------
% settings for grant plots
%---------------------------------------------------------------------
spikeThreshold = -0.35;
binSize = 20;
psthlim = [0 20];
% initialize stimulus onset offset line structs
sound_onoff = struct('onset', 0, 'offset', 1, 'color', 'b');
opto_onoff = struct('onset', 0, 'offset', 1, 'color', 'g');
% subplots?
subplotONOFF = 0;
% save figure?
saveFIG = 1;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot LFH data for file 1
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf1;
T = T1;
[~, fname] = fileparts(file_list{1});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find NULL index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('NULL', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot NULL data for file 2
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Dinf = Dinf2;
T = T2;
[~, fname] = fileparts(file_list{2});
%---------------------------------------------------------------------
% get list of stimuli (wav file names) and find NULL index (store in v)
%---------------------------------------------------------------------
varlist = Dinf.test.wavlist;
v = find(strcmpi('NULL', varlist));
% sample rate
Fs = Dinf.indev.Fs;
%---------------------------------------------------------------------
% stimulus on/off times
%---------------------------------------------------------------------
sound_onoff.onset = Dinf.audio.Delay;
sound_onoff.offset = Dinf.audio.Delay + Dinf.audio.Duration;
if Dinf.opto.Enable == 1
	opto_onoff.onset = Dinf.opto.Delay;
	opto_onoff.offset = Dinf.opto.Delay + Dinf.opto.Dur;
	stim_onoff = [sound_onoff opto_onoff];
else
	stim_onoff = sound_onoff;
end
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
spiketimes = spikeschmitt2(T{v}', spikeThreshold*global_max, 1, Fs);
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
% compute the mean rate
%---------------------------------------------------------------------
rdata = computeRLF({spiketimes}, [sound_onoff.onset sound_onoff.offset]);
%---------------------------------------------------------------------
% plot data
%---------------------------------------------------------------------
% new figure if subplotONOFF
if subplotONOFF
	figure; %#ok<UNRCH>
end
% time vector for plotting
t = (1000/Fs)*((1:length(T{v}(:, 1))) - 1);
[aH, fH] = grantplot(t, T{v}, global_max, spiketimes, histdata, stim_onoff, ...
								psthlim, subplotONOFF);
% captions
title(aH(1), {fname, varlist{v}}, 'Interpreter', 'none');
% name figure
if subplotONOFF == 1
	set(fH, 'Name', [fname '_' varlist{v} '_plots']);
else
	set(fH(1), 'Name', [fname '_' varlist{v} '_traces']);
	set(fH(2), 'Name', [fname '_' varlist{v} '_raster']);
	set(fH(3), 'Name', [fname '_' varlist{v} '_psth']);
end
% save fig?
if saveFIG
	if subplotONOFF ==1
		plotfile = get(fH, 'Name');
		savefig(fH, fullfile(plotpath, plotfile));
	else
		for f = 1:3
			plotfile = get(fH(f), 'Name');
			savefig(fH(f), fullfile(plotpath, plotfile));
			print(fH(f), fullfile(plotpath, plotfile), '-dpdf');
			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
				print(fH(f), fullfile(plotpath, plotfile), '-dmeta');
			end
		end
	end
	matfile = fullfile(plotpath, [fname '_' varlist{v} '_figdata.mat']);
	save(matfile, 'rdata', '-MAT');	
end
