function varargout = optoproc(varargin)
%------------------------------------------------------------------------
% optoproc
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% Processes data collected by the opto program
%
%------------------------------------------------------------------------
%  Input Options:
%   DATAFILE or FILE		path and name of .dat file
%   PLOTPATH				path (directory) for output of plot files
%   HPFREQ					high pass filter cutoff frequency for neural data (Hz)
%   LPFREQ					low pass filter cutoff frequency for neural data (Hz)
%   THRESHOLD				RMS spike threshold (# RMS)
%   CHANNEL					Input data channel (1-16)
%   BINSIZE					binsize for PSTH (ms)
%   PLOT_TRACES or		draw plots of all individual traces (1 = yes, 2 = no)
%		PLOTTRACES
%   PLOT_PSTH or			draw plots of PSTHs (1 = yes, 2 = no)
%		PLOTPSTH	
%   PLOT_RLF				plot Rate-Level Function
%	 PLOT_FTC				plot Frequency Tuning Curve
%   PLOT_FRA				plot FrequencyResponseArea data
%   SAVEFIG					save .fig plots (1 = yes, 2 = no)
%   SAVEPNG					save plots as .png files (1 = yes, 2 = no)
%   SAVEPDF					save plots as .pdf files (1 = yes, 2 = no)
%   TIMELIMITS				limit time range to specified interval
%									(e.g., [0 1000], in ms);
%   YLIMITS					limit y axis to specified interval (e.g., [-1 1])
%   PLOT_FILENAME			user-specified plot file base name (e.g., '1142_unit1')
%   PLOTROWCOLS			# of rows and columns for plots	
%									([2 3] is 2 rows, 3 cols)
%	 SHOW_DEFAULTS			show default values for options
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
%	22 Jan 19 (SJS): added input arg documentation
%	7 Feb 19 (SJS): working to streamline processing of data - pulling
%						 operations into separate functions as much as possible
%	26 Mar 19 (SJS): function additions
% 		- fixed issue in autoRowsCols for nvars
% 		- added PLOT_FRA
%	27 Mar 19 (SJS):
%		- added PLOT_RLF, PLOT_FTC
%--------------------------------------------------------------------------

%---------------------------------------------------------------------
% settings for processing data
%---------------------------------------------------------------------
datafile = '';
plotpath_base = ''; 
userPLOTPATH = 0;
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
% Plot rate level function?
plotRateLevelFun = 0;
% Plot frequency tuning curve?
plotFreqTuningCrv = 0;
% Plot FRA?
plotFreqRespArea = 0;
% SAVE PLOTS?
saveFIG = 0;
savePNG = 0;
savePDF = 0;
timeLimits = [];
yLimits = [];
% plot file name
plotFileName = '';
% plot rows/cols
autoRowCols = 1;
prows = [];
pcols = [];

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
				plotpath_base = varargin{argIndx + 1};
				userPLOTPATH = 1;
				argIndx = argIndx + 2;
			case 'HPFREQ'
				HPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'LPFREQ'
				LPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'THRESHOLD'
				tmp = varargin{argIndx + 1};
				if ischar(tmp)
					if strcmpi(tmp, 'DEFAULT')
						fprintf('%s: using default threshold: %d\n', ...
															mfilename, Threshold);
					else
						error('%s: unknown threshold command: %s', mfilename, tmp);
					end
				elseif isnumeric(tmp)
					if tmp > 0
						Threshold = tmp;
					else
						error('%s: invalid threshold value: %.4f', mfilename, tmp)
					end
				else
					error('%s: invalid threshold value: %s', mfilename, tmp)
				end
				argIndx = argIndx + 2;
			case 'CHANNEL'
				channelNumber = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'BINSIZE'
				binSize = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case {'PLOT_TRACES', 'PLOTTRACES'}
				plotTraces = 1;
				argIndx = argIndx + 1;
			case {'PLOT_PSTH', 'PLOTPSTH'}
				plotPSTH = 1;
				argIndx = argIndx + 1;
			case {'PLOT_RLF', 'PLOTRLF', 'RLF'}
				plotRateLevelFun = 1;
				argIndx = argIndx + 1;
			case {'PLOT_FTC', 'PLOTFTC', 'FTC'}
				plotFreqTuningCrv = 1;
				argIndx = argIndx + 1;
			case {'PLOT_FRA', 'PLOTFRA', 'FRA'}
				plotFreqRespArea = 1;
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
				yLimits = varargin{argIndx + 1};
				if numel(yLimits) ~= 2
					error('%s: yLimits must be 2 element numerical vector', ...
												mfilename);
				end
				argIndx = argIndx + 2;
			case {'PLOT_FILENAME', 'PLOTFILENAME', 'PLOTFILE', 'PLOT_FILE'}
				plotFileName = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case {'PLOTROWCOLS'}
				autoRowCols = 0;
				tmp = varargin{argIndx + 1};
				if isnumeric(tmp)
					if length(tmp) == 2
						prows = tmp(1);
						pcols = tmp(2);
					else
						error('%s: need to provide [plot_rows plot_cols]', mfilename);
					end
				else
					error('%s: invalid argument to plotRowsCols %s', tmp);
				end
				argIndx = argIndx + 2;
			case 'SHOW_DEFAULTS'
				fprintf('%s: Default values:\n', mfilename)
				fprintf('\tDATAFILE: %s\n', datafile);
				fprintf('\tPLOTPATH: %s\n', plotpath_base);
				fprintf('\tHPFREQ: %d\n', HPFreq);
				fprintf('\tLPFREQ: %d\n', LPFreq);
				fprintf('\tTHRESHOLD: %d\n', Threshold);
				fprintf('\tCHANNEL: %d\n', channelNumber);
				fprintf('\tBINSIZE: %d\n', binSize);
				fprintf('\tPLOT_TRACES: %d\n', plotTraces);
				fprintf('\tPLOT_PSTH: %d\n', plotPSTH);
				fprintf('\tPLOT_RLF: %d\n', plotRateLevelFun);
				fprintf('\tPLOT_FTC: %d\n', plotFreqTunCrv);
				fprintf('\tPLOT_FRA: %d\n', plotPSTH);
				fprintf('\tSAVEFIG: %d\n', saveFIG);
				fprintf('\tSAVEPNG: %d\n', savePNG);
				fprintf('\tSAVEPDF: %d\n', savePDF);
				fprintf('\tTIMELIMITS: '); 
					fprintf('%.2f ', timeLimits);	fprintf('\n');
				fprintf('\tYLIMITS: ');
					fprintf('%.2f ', yLimits); fprintf('\n');
				fprintf('\tPLOTROWCOLS: ');
				if autoRowCols
					fprintf('autoRowCols = 1\n');
				else
					fprintf('%d rows %d cols\n', prows, pcols);
				end
				return
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
if ~userPLOTPATH
	plotpath_base = fullfile(data_root_path, 'Analyzed');
end

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
	warning('%s: D is empty???!!!??!!', mfilename);
	return
end

%---------------------------------------------------------------------
% get info from filename - this makes some assumptions about file
% name structure!
% <animal id #>_<date>_<penetration #>_<unit #>_<other info>.dat
%---------------------------------------------------------------------
% break up file name into <fname>.<ext> (~ means don't save ext info)
[~, fname] = fileparts(datafile);
% locate underscores in fname
usc = find(fname == '_');
% location of start and end underscore indices
%    abcde_edcba
%        ^ ^
%        | |
%        | ---- endusc index
%        ---startusc index
endusc = usc - 1;
startusc = usc + 1;
animal = fname(1:endusc(1));
datecode = fname(startusc(1):endusc(2));
penetration = fname(startusc(2):endusc(3)); %#ok<NASGU>
unit = fname(startusc(3):endusc(4)); %#ok<NASGU>
other = fname(startusc(end):end); %#ok<NASGU>

if isempty(plotFileName)
	plotFileName = fname;
end

%---------------------------------------------------------------------
% create plot output dir if plots will be saved to files
%---------------------------------------------------------------------
if any([saveFIG savePDF savePNG]) && any([plotTraces plotPSTH])				
	if userPLOTPATH
		plotpath = plotpath_base;
	else
		plotpath = fullfile(plotpath_base, animal, datecode); 
		fprintf('Files will be written to:\n\t%s\n', plotpath);
		if ~exist(plotpath, 'dir')
			mkdir(plotpath);
		end
	end
end

%---------------------------------------------------------------------
% determine global RMS and max - used for thresholding
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
% compute overall mean rms for threshold
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
		% list of levels, and # of levels tested
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
	case 'FREQ+LEVEL'
		% list of freq, levels
		varlist = cell(2, 1);
		% # of freqs in nvars(1), # of levels in nvars(2)
		nvars = zeros(2, 1);
		for v = 1:2
			varlist{v} = unique(Dinf.test.stimcache.vrange(v, :), 'sorted');
			nvars(v) = length(varlist{v});
		end
		titleString = fname;

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
% local copy of sample rate
Fs = Dinf.indev.Fs;

% different approaches to storage depending on test type
switch upper(Dinf.test.Type)
	case 'FREQ+LEVEL'
		% for FRA data, nvars has values [nfreqs nlevels];
		spiketimes = cell(nvars(2), nvars(1));
		for v1 = 1:nvars(1)
			for v2 = 1:nvars(2)
				% use rms threshold to find spikes
				spiketimes{v2, v1} = ...
						spikeschmitt2(tracesByStim{v2, v1}', Threshold*mean_rms, ...
																			1, Fs, 'ms');
			end
		end
	otherwise
		% if test is not FREQ+LEVEL (FRA), nvars will be a single number
		spiketimes = cell(nvars, 1);
		for v = 1:nvars
			% use rms threshold to find spikes
			spiketimes{v} = spikeschmitt2(tracesByStim{v}', Threshold*mean_rms, ...
																				1, Fs, 'ms');
		end
end

%---------------------------------------------------------------------
% determine # of columns of plots
%---------------------------------------------------------------------
if autoRowCols
	if numel(nvars) == 1
		% for data that are not "2D" (e.g., FRA), adjust # of columns based
		% on the number of variable levels or types (nvars)
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
	else
		% otherwise, plot levels in rows, freqs in columns
		% rows = levels, cols = freqs
		prows = nvars(2);
		pcols = nvars(1);
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
		pname = fullfile(plotpath, [plotFileName '_traces']); 
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
		% time vector for plotting
		t = (1000/Fs)*((1:length(tracesByStim{1}(:, 1))) - 1);
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
	if numel(nvars) == 1
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
			plotopts.stimulus_times = ...
								reshape(	[stimulus_times; stimulus_times{end}], ...
											[prows pcols]);
			rasterpsthmatrix(reshape([spiketimes; {{}}], [prows pcols]), plotopts);
		end
	else
		% for 2D data... might be overkill....
		stimulus_times = cell(nvars(2), nvars(1));
 		plotopts.plot_titles = cell(nvars(2), nvars(1));
		for f = 1:nvars(1)
			for l = 1:nvars(2)
				% need to have [stim_onset stim_offset], so add delay to 
				% [0 duration] to compute proper times. then, multiply by 0.001 to
				% give times in seconds (Dinf values are in milliseconds)
				stimulus_times{l, f} = 0.001 * (Dinf.audio.Delay + ...
																	[0 Dinf.audio.Duration]);
				% if opto is Enabled, add it to the array by concatenation
				if Dinf.opto.Enable
					stimulus_times{l, f} = [stimulus_times{l, f}; ...
														 0.001 * (Dinf.opto.Delay + ...
																	[0 Dinf.opto.Dur]) ];
				end
				plotopts.plot_titles{l, f} = sprintf('%d kHz', varlist{1}(f));
			end
		end

		plotopts.stimulus_times = stimulus_times;
		rasterpsthmatrix(spiketimes, plotopts);
		
	end
	% set plot name
	set(hPR, 'Name', plotFileName)

	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_rp']);
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
% Rate-Level function plot
%---------------------------------------------------------------------
if plotRateLevelFun && strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not FRA data
	warning('optoproc: Test type (%s) is FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotRateLevelFun && any(strcmpi(Dinf.test.Type, {'LEVEL', 'BBN'}))
	% time window for counting spikes - use Delay to Delay+Duration interval
	analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	RLF = computeRLF(spiketimes, Dinf.audio.Level, analysisWindow);
	hRLF = plotRLF(RLF, 'median');
	% build title string
	tStr = {sprintf('RLF [%d %d] dB SPL', min(RLF.levels), max(RLF.levels)), ...
					[datafile ', ' sprintf('Channel %d', channelNumber)]};
	title(tStr, 'Interpreter', 'none');
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_RLF']);
		if saveFIG
			savefig(hRLF, pname, 'compact');
		end
		if savePDF
			print(hRLF, pname, '-dpdf');
		end
		if savePNG
			print(hRLF, pname, '-dpng', '-r300');
		end
	end	
end

%---------------------------------------------------------------------
% Frequency tuning curve plot
%---------------------------------------------------------------------
if plotFreqTunCrv && strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not useful for FRA data (at moment - this could be used to plot a
	% "slice" of the FRA!!!!
	warning('optoproc: Test type (%s) is FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotFreqTunCrv && ~strcmpi(Dinf.test.Type, 'FREQ')
	% other non freq test
	warning('optoproc: Test type (%s) is not compatible with FTC plot!', ...
					Dinf.test.Type);
elseif plotFreqTunCrv && strcmpi(Dinf.test.Type, 'FREQ')
	% time window for counting spikes - use Delay to Delay+Duration interval
	analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	FTC = computeFTC(spiketimes, Dinf.audio.Freqs, analysisWindow);
	hFTC = plotFTC(FTC, 'median');
	% build title string
	tStr = {sprintf('FTC [%d %d] kHz', min(FTC.freqs), max(FTC.freqs)), ...
					[datafile ', ' sprintf('Channel %d', channelNumber)]};
	title(tStr, 'Interpreter', 'none');
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_FTC']);
		if saveFIG
			savefig(hFTC, pname, 'compact');
		end
		if savePDF
			print(hFTC, pname, '-dpdf');
		end
		if savePNG
			print(hFTC, pname, '-dpng', '-r300');
		end
	end	
end


%---------------------------------------------------------------------
% FRA plot
%---------------------------------------------------------------------
if plotFreqRespArea && ~strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not FRA data
	warning('optoproc: Test type (%s) is not FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotFreqRespArea
	% they are FRA data
	% window for spike count
	frawin = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	% calculate FRA stored in struct FRA
	FRA = computeFRA(spiketimes, varlist{1}, varlist{2}, frawin);
	% set fname to data file name
	FRA.fname = datafile;
	hFRA = plotFRA(FRA, 'dB');
	
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_FRA']);
		if saveFIG
			savefig(hFRA, pname, 'compact');
		end
		if savePDF
			print(hFRA, pname, '-dpdf');
		end
		if savePNG
			print(hFRA, pname, '-dpng', '-r300');
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
									'global_max', global_max, 'Threshold', Threshold, ...
									'nvars', nvars, 'varlist', {varlist});
	varargout{4} = tracesByStim;
	if plotPSTH
		varargout{5} = plotopts;
	end
end
