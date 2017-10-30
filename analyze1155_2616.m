%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% settings and files
%---------------------------------------------------------------------
files_for_1155_2616;


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Process RLF data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% # of data files
nfiles = length(optoRLF.files);
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% get indices to sort by intensity
[sortedIntensity, sortedIndices] = sort(optoRLF.LEDintensity);
% data window (absolute time re: sweep onset)
AnalysisWindow = [100 200];
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
						'Dinf', [], ...
						'S', [], ...
						'F', []	);

%---------------------------------------------------------------------
%% compute RLFs
%---------------------------------------------------------------------
% loop through files in sorted order
for f = 1:nfiles
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
ptitle_base = [animalID '-' dateID '-' unit '-' penetration '-' depth ];...
legend_str = cell(nfiles, 1);

figure;
for f = 1:nfiles

	rlf = RLFdata(f).F;
	opto_mV = RLFdata(f).LEDintensity;

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
	legend_str{f} = sprintf('%d mV', RLFdata(f).LEDintensity);
end

legend(legend_str)

%---------------------------------------------------------------------
%% plot RLF of 0, 500, 1000, 1500, 3500
%---------------------------------------------------------------------
figure;

for f = [1 3 5 6 7]

	rlf = RLFdata(f).F;
	opto_mV = RLFdata(f).LEDintensity;

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
	legend_str{f} = sprintf('%d mV', RLFdata(f).LEDintensity);
end

legend(legend_str([1 3 5 6 7]))