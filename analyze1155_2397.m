%---------------------------------------------------------------------
% settings and files
%---------------------------------------------------------------------
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171025';
% unit #
unit = '02';
% penetration
penetration = '02';
% recording depth
depth = '2397';
% channel
channelNumber = 8;

% for all files:
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, varied command voltage to ThorLabs LED
optoRLF.files = {	...
	'1155_20171025_02_02_2397_BBN.dat', ...					% rlf, 0:10:70
	'1155_20171025_02_02_2397_BBN_optoON.dat', ...			% rlf, 0:10:70, opto 3500 mV
	'1155_20171025_02_02_2397_BBN_optoON_2.dat', ...		% rlf, 0:10:70, opto 500 mV
	'1155_20171025_02_02_2397_BBN_optoON_3.dat', ...		% rlf, 0:10:70, opto 250 mV
	'1155_20171025_02_02_2397_BBN_optoON_4.dat', ...		% rlf, 0:10:70, opto 750 mV
	'1155_20171025_02_02_2397_BBN_optoON_5.dat', ...		% rlf, 0:10:70, opto 1000 mV
	'1155_20171025_02_02_2397_BBN_optoON_6.dat', ...		% rlf, 0:10:70, opto 1500 mV
	'1155_20171025_02_02_2397_BBN_optoON_7.dat', ...		% rlf, 0:10:70, opto 2000 mV
	'1155_20171025_02_02_2397_BBN_optoON_8.dat', ...		% rlf, 0:10:70, opto 2500 mV
	'1155_20171025_02_02_2397_BBN_optoON_9.dat' ...			% rlf, 0:10:70, opto 3000 mV
};
optoRLF.LEDintensity = [ ...
		0, ...
		3500, ...
		500, ...
		250, ...
		750, ...
		1000, ...
		1500, ...
		2000, ...
		2500, ...
		3000 ...
];
optoRLF.LEDpower = [ ...
	0.0, ...
	6.8, ...
	1.4, ... 
	0.7, ... 
	2.0, ...
	2.5, ...
	3.6, ...
	4.5, ...
	5.3, ...
	6.1 ...
];

bbnRLF.files = '1155_20171025_02_02_2397_BBN_3.dat';	% rlf, 0:10:60

bbn.files = '1155_20171025_02_02_2397_BBN_4.dat';	% BBN, 60 dB, 15 reps


% sliding "pulse" of opto stimulation
% for all files:
% 	acq duration = 600 ms, ISI = 250 ms
% 	15 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], [0 30 60] dB, 150ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		50 ms duration, delay varied, 2000 mV command voltage to ThorLabs LED
optoSlidingWin2K.files = {	...
	'1155_20171025_02_02_2397_BBN_5.dat', ...				% no opto
	'1155_20171025_02_02_2397_BBN_optoON_10.dat', ...	% opto: 50 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_11.dat', ...	% opto: 100 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_12.dat', ...	% opto: 150 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_13.dat', ...	% opto: 175 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_14.dat', ...	% opto: 200 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_15.dat', ...	% opto: 250 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_16.dat' ...	% opto: 300 ms delay, 2000 mV, 50ms dur
};

optoSlidingWin3K.files = {	...
	'1155_20171025_02_02_2397_BBN_optoON_17.dat', ...	% opto: 100 ms delay, 3000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_18.dat', ...	% opto: 150 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_19.dat', ...	% opto: 200 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_20.dat' ...	% opto: 250 ms delay, 2000 mV, 50ms dur
};


% 	acq duration = 400 ms, ISI = 100 ms
% 	15 reps
% 	sound:
% 		tone, freq:5k:5k:80k, 35dB, 100ms delay, 100ms duration, 5ms ramp
freqtuning.files = 	'1155_20171025_02_02_2397_FREQ_TUNING.dat';

% for all files:
% 	acq duration = 600 ms, ISI = 300 ms
% 	15 reps
% 	sound:
% 		tone, freq:5k:1250:25k, 35dB, 150ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, varied command voltage to ThorLabs LED
optoFreqTuning.files = { ...
	'1155_20171025_02_02_2397_FREQ.dat', ...				% freq,5k:1250:25k opto 0 mV
	'1155_20171025_02_02_2397_FREQ_optoON.dat', ...		% freq, opto 3000 mV
	'1155_20171025_02_02_2397_FREQ_optoON_2.dat', ...	% freq, opto 500 mV
	'1155_20171025_02_02_2397_FREQ_optoON_3.dat', ...	% freq, opto 1000 mV
	'1155_20171025_02_02_2397_FREQ_optoON_4.dat' ...	% freq, opto 1500 mV
};

% for all files:
% 	20 reps
% 	opto: 
% 		250 ms duration, 150 ms delay, 3000 mV to ThorLabs LED
optoWav.files = {	...
	'1155_20171025_02_02_2397_OptoInhibOFF.dat', ...
	'1155_20171025_02_02_2397_OptoInhibON.dat' ...
};

%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Process RLF data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% # of data files
nfiles = length(optoRLF.files);
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% get indices to sort by opto intensity
[sortedIntensity, sortedIndices] = sort(optoRLF.LEDintensity);
% data window (absolute time re: sweep onset)
AnalysisWindow = [150 200];
% psth bin size
binSize = 5;
% threshold (RMS)
spikeThreshold = 3;

%---------------------------------------------------------------------
%% allocate RLFdata in sorted intensity order
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
RLFdata = struct(	'datafile', optoRLF.files(sortedIndices), ...
						'LEDintensity', num2cell(sortedIntensity), ...
						'LEDpower', num2cell(optoRLF.LEDpower(sortedIndices)), ...
						'D', [], ...
						'Dinf', [], ...
						'S', [], ...
						'T', [], ...
						'F', []	);

%---------------------------------------------------------------------
%% compute RLFs
%---------------------------------------------------------------------
% loop through files in sorted order
for f = 1:nfiles
	fprintf('Analyzing %s\n', RLFdata(f).datafile);
% 	[~, RLFdata(f).Dinf, RLFdata(f).S] = ...
% 					optoproc('file', fullfile(datapath, RLFdata(f).datafile), ...
% 								'channel', channelNumber, ...
% 								'binsize', binSize, ...
% 								'Threshold', spikeThreshold, 'plotpsth', 'savepdf');
	[RLFdata(f).D, RLFdata(f).Dinf, RLFdata(f).S, RLFdata(f).T] = ...
					optoproc('file', fullfile(datapath, RLFdata(f).datafile), ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'Threshold', spikeThreshold);
	% compute RLF
	RLFdata(f).F = computeRLF(RLFdata(f).S.spiketimes, AnalysisWindow);
	RLFdata(f).F.nlevels = length(RLFdata(f).Dinf.audio.Level);
	RLFdata(f).F.level = RLFdata(f).Dinf.audio.Level;
end

%---------------------------------------------------------------------
%% save RLF data
%---------------------------------------------------------------------
rlfdatapath = fullfile(outpath_base, animalID, dateID); 
if ~exist(rlfdatapath, 'dir')
	mkdir(rlfdatapath);
end
rlfdatafile = [animalID '_' dateID '_' unit '_' penetration '_' depth '_' ...
						'RLFdata.mat'];
save(fullfile(rlfdatapath, rlfdatafile), 'RLFdata', 'optoRLF', ...
			'optoFreqTuning', '-MAT');

%---------------------------------------------------------------------
%% plot RLF
%---------------------------------------------------------------------
ptitle_base = [animalID '-' dateID '-' unit '-' penetration '-' depth ];...
legend_str = cell(nfiles, 1);

figure;
for f = 1:nfiles

	rlf = RLFdata(f).F;
	opto_mW = RLFdata(f).LEDpower;

% 	ebneg = zeros(rlf.nlevels, 1);
% 	ebpos = zeros(rlf.nlevels, 1);
% 	for l = 1:rlf.nlevels
% 		ebneg(l) = rlf.mean_ci{l}(1);
% 		ebpos(l) = rlf.mean_ci{l}(2);
% 	end
% 	errorbar(rlf.level, rlf.mean, ebneg, ebpos);
% 	for l = 1:rlf.nlevels
% 		ebneg(l) = rlf.median_ci{l}(1);
% 		ebpos(l) = rlf.median_ci{l}(2);
% 	end
% 	errorbar(rlf.level, rlf.median, ebneg, ebpos);
% 	if f == 1
% 		errorbar(rlf.level, rlf.mean, rlf.std);
% 		title( {	ptitle, sprintf('window: [%d %d] ms', AnalysisWindow) }, ...
% 					'Interpreter', 'none');
% 		xlabel('Sound Level (dB SPL)');
% 		ylabel('mean spikes/trial +/- s.d.');
% 	else
% 		hold on
% 			errorbar(rlf.level, rlf.mean, rlf.std);
% 		hold off
% 	end
	if f == 1
		plot(rlf.level, rlf.mean, '.-');
		title( {	ptitle_base, sprintf('window: [%d %d] ms', AnalysisWindow) }, ...
					'Interpreter', 'none');
		xlabel('Sound Level (dB SPL)');
		ylabel('mean spikes/trial');
		set(gcf, 'Name', [ptitle_base '-OptoLevels']);
		set(gcf, 'FileName', [ptitle_base '-OptoLevels.fig']);
	else
		hold on
			plot(rlf.level, rlf.mean, '.-');
		hold off
	end
	legend_str{f} = sprintf('%.1f mW', RLFdata(f).LEDpower);
end
legend(legend_str, 'Location', 'best');
% set markers
lH = get(gca, 'Children');
for l = 1:length(lH)
	set(lH(l), 'MarkerSize', 12);
end
%---------------------------------------------------------------------
%% plot RLF of [0,250,500,750,1000,1500,2000,2500,3000]
%---------------------------------------------------------------------
figure;
findx = [1 2 3 4 5 6 7 8 9];
for f = findx
	rlf = RLFdata(f).F;
	opto_mV = RLFdata(f).LEDintensity;

	if f == 1
		plot(rlf.level, rlf.mean, '.-');
		title( {	ptitle_base, sprintf('window: [%d %d] ms', AnalysisWindow) }, ...
					'Interpreter', 'none');
		xlabel('Sound Level (dB SPL)');
		ylabel('mean spikes per trial');
		set(gcf, 'Name', [ptitle_base '-OptoLevels']);
		set(gcf, 'FileName', [ptitle_base '-OptoLevels.fig']);
	else
		hold on
			plot(rlf.level, rlf.mean, '.-');
		hold off
	end
	legend_str{f} = sprintf('%.1f mW', RLFdata(f).LEDpower);
end
legend(legend_str(findx), 'Location', 'best');
% set markers
lH = get(gca, 'Children');
for l = 1:length(lH)
	set(lH(l), 'MarkerSize', 12);
end
%---------------------------------------------------------------------
%% plot RLF of [0, 0.7, 2.0, 5.3]]
%---------------------------------------------------------------------
figure;
findx = [1 2 4 8];
for f = findx
	rlf = RLFdata(f).F;
	% either create new plot or add lines to plot
	if f == 1
		plot(rlf.level, rlf.mean, '.-');
		title( {	ptitle_base, sprintf('window: [%d %d] ms', AnalysisWindow) }, ...
					'Interpreter', 'none');
		xlabel('Sound Level (dB SPL)');
		ylabel('mean spikes per trial');
		set(gcf, 'Name', [ptitle_base '-OptoLevels']);
		set(gcf, 'FileName', [ptitle_base '-OptoLevels.fig']);
	else
		hold on
			plot(rlf.level, rlf.mean, '.-');
		hold off
	end
	legend_str{f} = sprintf('%.1f mW', RLFdata(f).LEDpower);
end
legend(legend_str(findx), 'Location', 'northwest');
% set markers, line width
lH = get(gca, 'Children');
for l = 1:length(lH)
	set(lH(l), 'MarkerSize', 12);
	set(lH(l), 'LineWidth', 1.25);
end
set(gca, 'Tickdir', 'out');
set(gca, 'Box', 'off');
ylim([0 10])
set(gca, 'YTick', 0:5:10);
rlfplotfile = [animalID '_' dateID '_' unit '_' penetration '_' depth '_RLF'];
set(gcf, 'Name', rlfplotfile);
savefig(gcf, fullfile(rlfdatapath , rlfplotfile));

%---------------------------------------------------------------------
%% get select PSTH data
%---------------------------------------------------------------------
findx = [1 2 4 8];
tgtlevel = 50;
subplotONOFF = 0;
% save figure?
saveFIG = 1;
if tgtlevel == 30
	psthminmax = [0 25];
elseif tgtlevel == 50
	psthminmax = [0 40];
else
	psthminmax = [0 20];
end
fname = [animalID '_' dateID '_' unit '_' penetration '_' depth];
% initialize stimulus onset offset line structs
sound_onoff = struct('onset', 0, 'offset', 1, 'color', 'b');
opto_onoff = struct('onset', 0, 'offset', 1, 'color', 'g');

for f = findx
	rlf = RLFdata(f);
	tgtindx = (tgtlevel == rlf.F.level);
	spiketimes = rlf.S.spiketimes{tgtindx};
	[histdata.H, histdata.bins] = psth(spiketimes, binSize, ...
														[0 rlf.Dinf.test.AcqDuration]);
	% sample rate
	Fs = rlf.Dinf.indev.Fs;
	traces = rlf.T{tgtindx};
	% time vector for plotting
	t = (1000/Fs)*((1:length(traces(:, 1))) - 1);
	% stimulus on/off times
	sound_onoff.onset = rlf.Dinf.audio.Delay;
	sound_onoff.offset = rlf.Dinf.audio.Delay + rlf.Dinf.audio.Duration;
	if rlf.Dinf.opto.Enable == 1
		opto_onoff.onset = rlf.Dinf.opto.Delay;
		opto_onoff.offset = rlf.Dinf.opto.Delay + rlf.Dinf.opto.Dur;
		stim_onoff = [sound_onoff opto_onoff];
	else
		stim_onoff = sound_onoff;
	end
	if subplotONOFF
		figure
	end
	[aH, fH] = grantplot(t, traces, max(max(abs(traces))), ...
									spiketimes, histdata, stim_onoff, ...
									psthminmax, subplotONOFF);
	% captions
	titlestr = {fname, [num2str(tgtlevel) ' dB SPL'], ...
						[num2str(rlf.LEDpower) ' mW'] };
	title(aH(1), titlestr, 'Interpreter', 'none');
	% name figure
	basename = [fname '_' num2str(tgtlevel) 'dB' num2str(10*rlf.LEDpower)];
	if subplotONOFF == 1
		set(fH, 'Name',  [basename '_plots']);
	else
		set(fH(1), 'Name', [basename '_T']);
		set(fH(2), 'Name', [basename '_R']);
		set(fH(3), 'Name', [basename '_P']);
	end
	% save fig?
	if saveFIG
		if subplotONOFF == 1
			plotfile = get(fH, 'Name');
			savefig(fH, fullfile(rlfdatapath, plotfile));
		else
			for p = 1:3
				plotfile = get(fH(p), 'Name');
				savefig(fH(p), fullfile(rlfdatapath, plotfile));
				print(fH(p), fullfile(rlfdatapath, plotfile), '-dpdf');
	% 			print(fH(f), fullfile(plotpath, plotfile), '-dpng', '-r300');
	% 			if strcmpi(computer, 'PCWIN') || strcmpi(computer, 'PCWIN64')
	% 				print(fH(f), fullfile(rlfdatapath, plotfile), '-dmeta');
	% 			end
			end
		end
	% 	matfile = fullfile(rlfdatapath, [fname '_' varlist{v} '_figdata.mat']);
	% 	save(matfile, 'rdata', '-MAT');	
	end
	
end

%---------------------------------------------------------------------




%---------------------------------------------------------------------
%% plot sliding window, 2000 mV data
%---------------------------------------------------------------------
% loop through files in sorted order
for f = 1:length(optoSlidingWin2K.files)
	fprintf('plotting %s\n', optoSlidingWin2K.files{f});
	optoproc('file', fullfile(datapath, optoSlidingWin2K.files{f}), ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'Threshold', spikeThreshold, 'plotpsth', 'savepdf');
end

