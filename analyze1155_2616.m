%---------------------------------------------------------------------
% settings and files
%---------------------------------------------------------------------
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171025';
% unit #
unit = '05';
% penetration
penetration = '03';
% recording depth
depth = '2616';
% channel
channelNumber = 8;

% files
'1155_20171025_05_03_2616_BBN.dat';					% rate level function BBN 0:10:70
'1155_20171025_05_03_2616_BBN_optoON.dat';		% rlf, opto on, 3500 mV

% for all files:
% 	acq duration = 500 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:70dB, 100ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 50 ms delay, varied command voltage to ThorLabs LED
optoRLF.files = {	...
	'1155_20171025_05_03_2616_BBN_2.dat', ...				% rlf, 0:10:70
	'1155_20171025_05_03_2616_BBN_optoON_2.dat', ...	% rlf, 0:10:70, opto 3500 mV
	'1155_20171025_05_03_2616_BBN_optoON_3.dat', ...	% rlf, 0:10:70, opto 100 mV
	'1155_20171025_05_03_2616_BBN_optoON_4.dat', ...	% rlf, 0:10:70, opto 500 mV
	'1155_20171025_05_03_2616_BBN_optoON_5.dat', ...	% rlf, 0:10:70, opto 1000 mV
	'1155_20171025_05_03_2616_BBN_optoON_6.dat', ...	% rlf, 0:10:70, opto 750 mV
	'1155_20171025_05_03_2616_BBN_optoON_7.dat'  ...	% rlf, 0:10:70, opto 1500 mV
};
optoRLF.LEDintensity = [ ...
	0, ...
	3500, ...
	100, ...
	500, ...
	1000, ...
	750, ...
	1500 ...
];

% sliding "pulse" of opto stimulation
% for all files:
% 	acq duration = 500 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], [0 40 60] dB, 100ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		25 ms duration, delay varied, 3000 mV command voltage to ThorLabs LED
optoSlidingWin_files = {	...
	'1155_20171025_05_03_2616_BBN_3.dat', ...				% sliding opto, control, BBN [0,40,60]dB; 
	'1155_20171025_05_03_2616_BBN_optoON_9.dat', ...	% opto: 50 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_10.dat', ...	% opto: 100 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_11.dat', ...	% opto: 150 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_12.dat'  ...	% opto: 200 ms delay, 3000 mV, 25ms dur
};

misc_files = { ...
	'1155_20171025_05_03_2616_BBN_optoON_13.dat', ...	% opto: 50 ms delay, 3000 mV, 50ms dur
	'1155_20171025_05_03_2616_BBN_optoON_14.dat', ...	% opto: 100 ms delay, 3000 mV, 50ms dur
	'1155_20171025_05_03_2616_BBN_optoON_15.dat' ...	% opto: 150 ms delay, 3000 mV, 50ms dur
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
nfiles = length(optoRLF_files);
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% get indices to sort by intensity
[sortedIntensity, sortedIndices] = sort(optoRLF.LEDintensity);
% data window (absolute time re: sweep onset)
AnalysisWindow = [100 200];
% psth bin size
binSize = 5;
% threshold (RMS)
spikeThreshold = 3.25;

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
						'Dinf', [], ...
						'S', [], ...
						'F', []	);

%---------------------------------------------------------------------
%% compute RLFs
%---------------------------------------------------------------------
% loop through files in sorted order
for f = 1
	fprintf('Analyzing %s\n', RLFdata(f).datafile);
	[~, RLFdata(f).Dinf, RLFdata(f).S] = ...
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
%% plot RLF
%---------------------------------------------------------------------
ptitle_base = [animalID '-' dateID '-' unit '-' penetration '-' depth '-' ...
						'opto_mV'];...
						
for f = 1
	figure;
	rlf = RLFdata(f).F;
	opto_mV = RLFdata(f).LEDintensity;
	ptitle = [ptitle_base num2str(opto_mV)];
	ebneg = zeros(rlf.nlevels, 1);
	ebpos = zeros(rlf.nlevels, 1);
% 	for l = 1:rlf.nlevels
% 		ebneg(l) = rlf.mean_ci{l}(1);
% 		ebpos(l) = rlf.mean_ci{l}(2);
% 	end
% 	errorbar(rlf.level, rlf.mean, ebneg, ebpos);
	for l = 1:rlf.nlevels
		ebneg(l) = rlf.median_ci{l}(1);
		ebpos(l) = rlf.median_ci{l}(2);
	end
	errorbar(rlf.level, rlf.median, ebneg, ebpos);
	title( {	ptitle, sprintf('window: [%d %d] ms', AnalysisWindow) }, ...
				'Interpreter', 'none');
	xlabel('Sound Level (dB SPL)');
	ylabel('mean spikes/trial +/- s.d.')
end
