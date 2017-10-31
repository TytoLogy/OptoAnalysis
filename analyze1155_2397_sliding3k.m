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

%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');
% out for analysis
analysispath = fullfile(outpath_base, animalID, dateID, ...
									'/OptoAnalysisPlots/SlidingWin');
if ~exist(analysispath, 'dir')
	mkdir(analysispath);
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Process data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% # of data files
% nfiles = length(optoRLF.files);
nfiles = length(optoSlidingWin3K.files);
% path to data
datapath = fullfile(data_root_path, animalID, dateID);
% data window (absolute time re: sweep onset)
AnalysisWindow = [150 200];
% psth bin size
binSize = 5;
% threshold (RMS)
spikeThreshold = 3.1;

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
data = struct(	'datafile', optoSlidingWin3K.files, ...
						'LEDintensity', optoSlidingWin3K.LEDintensity, ...
						'LEDpower', optoSlidingWin3K.LEDpower, ...
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
	fprintf('Analyzing %s\n', data(f).datafile);

	[~, data(f).Dinf, data(f).S] = ...
					optoproc('file', fullfile(datapath, data(f).datafile), ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'Threshold', spikeThreshold, 'plotpsth');

% 	[~, data(f).Dinf, data(f).S] = ...
% 					optoproc('file', fullfile(datapath, data(f).datafile), ...
% 								'channel', channelNumber, ...
% 								'binsize', binSize, ...
% 								'Threshold', spikeThreshold, 'plotpsth', 'savepdf', ...
% 								'plotPath', analysispath);
% 	[data(f).D, data(f).Dinf, data(f).S, data(f).T] = ...
% 					optoproc('file', fullfile(datapath, data(f).datafile), ...
% 								'channel', channelNumber, ...
% 								'binsize', binSize, ...
% 								'Threshold', spikeThreshold);
	% compute RLF
	% need to figure out analysis window
	data(f).F = computeRLF(data(f).S.spiketimes, AnalysisWindow);
	data(f).F.nlevels = length(data(f).Dinf.audio.Level);
	data(f).F.level = data(f).Dinf.audio.Level;
end

%---------------------------------------------------------------------
%% save RLF data
%---------------------------------------------------------------------
% rlfdatapath = fullfile(outpath_base, animalID, dateID); 
% if ~exist(rlfdatapath, 'dir')
% 	mkdir(rlfdatapath);
% end
rlfdatafile = [animalID '_' dateID '_' unit '_' penetration '_' depth '_' ...
						'slidingdata2k.mat'];
save(fullfile(analysispath, rlfdatafile), 'data', 'optoSlidingWin3K', '-MAT');



%---------------------------------------------------------------------
%% plot data across files for each audio intensity level
%---------------------------------------------------------------------

% make sure nlevels match
tmp = zeros(nfiles, 1);
for f = 1:nfiles
	tmp(f) = data(f).F.nlevels;
end
if any(diff(tmp))
	error('nlevels mismatch!');
else
	nlevels = tmp(1);
end
clear tmp;
clear plotopts

% plot options
plotopts.timelimits = [0 data(f).Dinf.test.AcqDuration];
plotopts.timelimits = [50 350];
plotopts.psth_binwidth = binSize;
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 24;
plotopts.raster_color = 'k';
plotopts.psth_color = 'k';
plotopts.stimulus_times_plot = 3;
plotopts.stimulus_on_color = { 'b', 'g' };
plotopts.stimulus_off_color = { 'b', 'g' };
plotopts.stimulus_onoff_pct = [0 50];
plotopts.stimulus_times = cell(nfiles, 1);
plotopts.psth_ylimits = [0 30];
plotopts.plotgap = 0.001;
plotopts.vertgap = 0.035;
plotopts.heightscale = 1.05;
plotopts.widthscale = 0.5;
% loop through levels
for l = 1:nlevels
	% create/switch to figure
	figure(l)
	% allocate spikes
	S = cell(nfiles, 1);
	% loop through files
	for f = 1:nfiles
		S{f} = data(f).S.spiketimes{l};
		if data(f).Dinf.opto.Enable
			plotopts.stimulus_times{f}(1, :) = data(f).Dinf.audio.Delay + ...
															[0 data(f).Dinf.audio.Duration];
			plotopts.stimulus_times{f}(2, :) = data(f).Dinf.opto.Delay + ...
															[0 data(f).Dinf.opto.Dur];
		else
			plotopts.stimulus_times{f} = data(f).Dinf.audio.Delay + ...
															[0 data(f).Dinf.audio.Duration];
		end
		% convert onset/offset times to seconds (odd quirk)
		plotopts.stimulus_times{f} = 0.001 * plotopts.stimulus_times{f};
	end
	
	[H(l), pO(l)] = rasterpsthmatrix(S, plotopts); %#ok<SAGROW>
	set(figure(l), 'Position',  [680     6   560   972]);
	set(figure(l), 'Name', sprintf('Sliding_3K_%ddB', data(f).F.level(l)));
	set(figure(l), 'FileName', sprintf('Sliding_3K_%ddB.fig', data(f).F.level(l)));
end


