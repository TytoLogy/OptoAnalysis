function varargout = optoproc(varargin)
%------------------------------------------------------------------------
% optoproc
%------------------------------------------------------------------------
% % TytoLogy:Experiments:opto Application
%--------------------------------------------------------------------------
% Processes data collected by the opto program
%
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ???
%
% Revisions:
%	see git!
%------------------------------------------------------------------------
% TO DO:
%	- Document
%	- Functionalize
%--------------------------------------------------------------------------

%---------------------------------------------------------------------
% settings for processing data
%---------------------------------------------------------------------
datafile = '';
plotpath_base = ''; %#ok<NASGU>
% filter
HPFreq = 350;
LPFreq = 6500;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
% Channel Number (use 8 for single channel data)
channelNumber = 8;
% binSize for PSTH (milliseconds)
binSize = 5;
% Plot Traces?
plotTraces = 0;
% plotPSTH?
plotPSTH = 0;
% SAVE PLOTS?
saveFIG = 0;
savePNG = 0;
savePDF = 0;
timeLimits = [];
yLimits = [];

%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
if nargin
	argIndx = 1;
	while argIndx <= nargin
		fprintf('%s\n', upper(varargin{argIndx}))
		switch upper(varargin{argIndx})
			case {'DATAFILE', 'FILE'}
				datafile = varargin{argIndx + 1};
				[datapath, dfile, dext] = fileparts(datafile);
				datafile = [dfile dext];
				argIndx = argIndx + 2;
			case 'PLOTPATH'
				plotpath_base = varargin{argIndx + 1}; %#ok<NASGU>
				argIndx = argIndx + 2;
			case 'HPFREQ'
				HPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'LPFREQ'
				LPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'THRESHOLD'
				Threshold = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'CHANNEL'
				channelNumber = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'BINSIZE'
				binSize = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'PLOT_TRACES'
				plotTraces = 1;
				argIndx = argIndx + 1;
			case 'PLOT_PSTH'
				plotPSTH = 1;
				argIndx = argIndx + 1;
			case 'SAVEFIG'
				saveFIG = 1;
				argIndx = argIndx + 1;
			case 'SAVEPNG'
				savePNG = 1;
				argIndx = argIndx + 1;
			case 'SAVEPDF'
				savePDF = 1;
				argIndx = argIndx + 1;
			case 'TIMELIMITS'
				timeLimits = 1;
				argIndx = argIndx + 1;
			case 'YLIMITS'
				yLimits = 1;
				argIndx = argIndx + 1;
			otherwise
				error('%s: unknown input arg %s', mfilename, varargin{argIndx});
		end
	end
end

%---------------------------------------------------------------------
% need to get information about system
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths; %#ok<ASGLU>
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% data file things
%---------------------------------------------------------------------
if isempty(datafile)
	% get data file from user
	[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
														data_root_path);
	% abort if cancelled
	if datafile == 0
		fprintf('Cancelled\n');
		return
	end	
end

%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(datapath, datafile), ...
											'Filter', [HPFreq LPFreq], ...
											'Channel', channelNumber);
if isempty(D)
	return
end

%---------------------------------------------------------------------
% get info from filename
%---------------------------------------------------------------------
[~, fname] = fileparts(datafile);
usc = find(fname == '_');
endusc = usc - 1;
startusc = usc + 1;
animal = fname(1:endusc(1));
datecode = fname(startusc(1):endusc(2));
penetration = fname(startusc(2):endusc(3)); %#ok<NASGU>
unit = fname(startusc(3):endusc(4)); %#ok<NASGU>
other = fname(startusc(end):end); %#ok<NASGU>

%---------------------------------------------------------------------
% create plot output dir
%---------------------------------------------------------------------
if any([saveFIG savePDF savePNG])
	plotpath = fullfile(plotpath_base, animal, datecode); 
	fprintf('Files will be written to:\n\t%s\n', plotpath);
	if ~exist(plotpath, 'dir')
		mkdir(plotpath);
	end
end

%---------------------------------------------------------------------
% determine global RMS and max
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps
if strcmpi(Dinf.test.Type, 'WavFile')
	nstim = Dinf.test.nCombinations;
	nreps = Dinf.test.Reps;
else
	nstim = Dinf.test.stimcache.ntrials;
	nreps = Dinf.test.stimcache.nreps;
end
% allocate matrices
netrmsvals = zeros(nstim, nreps);
maxvals = zeros(nstim, nreps);
% find rms, max vals for each stim
for s = 1:nstim
	netrmsvals(s, :) = rms(tracesByStim{s});
	maxvals(s, :) = max(abs(tracesByStim{s}));
end
% compute overall mean rms fir threshold
fprintf('Calculating mean and max RMS for data...\n');
mean_rms = mean(reshape(netrmsvals, numel(netrmsvals), 1));
fprintf('\tMean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(max(maxvals));
fprintf('\tGlobal max abs value: %.4f\n', global_max);

%---------------------------------------------------------------------
% Some test-specific things...
%---------------------------------------------------------------------
switch upper(Dinf.test.Type)
	case 'FREQ'
		% list of frequencies, and # of freqs tested
		varlist = Dinf.test.stimcache.vrange;
		nvars = length(varlist);
		titleString = cell(nvars, 1);
		for v = 1:nvars
			if v == 1
				titleString{v} = {fname, ...
										sprintf('Frequency = %.0f kHz', 0.001*varlist(v))};
			else
				titleString{v} = sprintf('Frequency = %.0f kHz', 0.001*varlist(v));
			end
		end
	case 'LEVEL'
		% list of legvels, and # of levels tested
		varlist = Dinf.test.stimcache.vrange;
		nvars = length(varlist);
		titleString = cell(nvars, 1);
		for v = 1:nvars
			if v == 1
				titleString{v} = {fname, sprintf('Level = %d dB SPL', varlist(v))};
			else
				titleString{v} = sprintf('Level = %d dB SPL', varlist(v));
			end
		end
	case 'OPTO'
		% not yet implemented
	case 'WAVFILE'
		% get list of stimuli (wav file names)
		varlist = Dinf.test.wavlist;
		nvars = length(varlist);
		titleString = cell(nvars, 1);
		for v = 1:nvars
			if v == 1 
				titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
			else
				titleString{v} = sprintf('wav name: %s', varlist{v});
			end
		end
	otherwise
		error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end

%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------
Fs = Dinf.indev.Fs;
spiketimes = cell(nvars, 1);
for v = 1:nvars
	spiketimes{v} = spikeschmitt2(tracesByStim{v}', Threshold*mean_rms, 1, Fs);
	for r = 1:length(spiketimes{v})
		spiketimes{v}{r} = (1000/Fs)*spiketimes{v}{r};
	end
end

%---------------------------------------------------------------------
% Plot raw data
%---------------------------------------------------------------------
if plotTraces
	% new figure
	hF = figure;
	% name figure
	set(hF, 'Name', [fname '_sweeps']);
	% determine # of columns of plots
	if nvars <= 6
		prows = nvars;
		pcols = 1;
	elseif iseven(nvars)
		prows = nvars/2;
		pcols = 2;
	else
		prows = ceil(nvars/2);
		pcols = 2;
	end
	% loop through variable
	for v = 1:nvars
		% time vector for plotting
		t = (1000/Fs)*((1:length(tracesByStim{v}(:, 1))) - 1);
		subplot(prows, pcols, v);
		% flip tracesByStim in order to have sweeps match the raster plots
		stackplot(t, fliplr(tracesByStim{v}), 'colormode', 'black', ...
												'ymax', global_max, ...
												'figure', hF, 'axes', gca);
		title(titleString{v}, 'Interpreter', 'none');
		xlabel('ms')
		ylabel('Trial')
	end
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [fname '_traces']); 
		if saveFIG
			savefig(hF, pname, 'compact');
		end
		if savePDF
			print(hF, pname, '-dpdf');
		end
		if savePNG
			print(hF, pname, '-dpng', '-r300');
		end
	end
end

%---------------------------------------------------------------------
% raster, psths
%---------------------------------------------------------------------
if plotPSTH
	hPR = figure;
	if isempty(timeLimits)
		plotopts.timelimits = [0 ceil(max(t))];
	else
		plotopts.timelimits = timeLimits;
	end
	if ~isempty(yLimits)
		plotopts.psth_ylimits = yLimits;
	end
	plotopts.raster_tickmarker = '.';
	plotopts.raster_ticksize = 16;
	plotopts.raster_color = [0 0 0];
	plotopts.psth_binwidth = binSize;
	plotopts.plotgap = 0.001;
	plotopts.xlabel = 'msec';
	plotopts.stimulus_times_plot = 3;
	plotopts.stimulus_on_color{1} = [0 0 1];
	plotopts.stimulus_off_color{1} = [0 0 1];
	plotopts.stimulus_onoff_pct(1) = 60;
	if Dinf.opto.Enable
		% add colors for second stimulus
		plotopts.stimulus_on_color{2} = [1 0 0];
		plotopts.stimulus_off_color{2} = [1 0 0];
		plotopts.stimulus_onoff_pct(2) = 80;
	end

	% create times to indicate stimuli
	stimulus_times = cell(nvars, 1);
	for v = 1:nvars
		% need to have [stim_onset stim_offset], so add delay to 
		% [0 duration] to compute proper times. then, multiply by 0.001 to
		% give times in seconds (Dinf values are in milliseconds)
		stimulus_times{v, 1} = 0.001 * (Dinf.audio.Delay + ...
															[0 Dinf.audio.Duration]);
		% if opto is Enabled, add it to the array by concatenation
		if Dinf.opto.Enable
			stimulus_times{v, 1} = [stimulus_times{v, 1}; ...
												 0.001 * (Dinf.opto.Delay + ...
															[0 Dinf.opto.Dur]) ];
		end
	end

	% adjust depending on # of columns of plots
	if nvars <= 5
		plotopts.plot_titles = titleString;
		plotopts.stimulus_times = stimulus_times;
		rasterpsthmatrix(spiketimes, plotopts);
	elseif iseven(nvars)
		plotopts.plot_titles = reshape(titleString, [prows pcols]);
		plotopts.stimulus_times = reshape(stimulus_times, [prows pcols]);
		rasterpsthmatrix(reshape(spiketimes, [prows pcols]), plotopts);
	else
		% need to add 'dummy' element to arrays
	% 	titleString = [titleString; {''}];
	% 	spiketimes = [spiketimes; {{}}];
		plotopts.plot_titles = reshape([titleString; {''}], [prows pcols]);
		plotopts.stimulus_times = reshape([stimulus_times; stimulus_times{end}], [prows pcols]);
		rasterpsthmatrix(reshape([spiketimes; {{}}], [prows pcols]), plotopts);
	end
	% set plot name
	set(hPR, 'Name', fname)

	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [fname '_rp']);
		if saveFIG
			savefig(hPR, pname, 'compact');
		end
		if savePDF
			print(hPR, pname, '-dpdf');
		end
		if savePNG
			print(hPR, pname, '-dpng', '-r300');
		end
	end
end

%---------------------------------------------------------------------
% outputs
%---------------------------------------------------------------------
if nargout
	varargout{1} = D;
	varargout{2} = Dinf;
	varargout{3} = struct('spiketimes', {spiketimes}, 'mean_rms', mean_rms, ...
									'global_max', global_max, 'Threshold', Threshold);
	if plotPSTH
		varargout{4} = plotopts;
	end
end
