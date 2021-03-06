function varargout = optoproc_plotPSTH_WAVbyLevel(spikesByStim, ...
										Dinf, binSize, ...
										nrowcols, timeLimits, yLimits, varargin)
%------------------------------------------------------------------------
%  H = optoproc_plotPSTH_WAVbyLevel(spikesByStim, Dinf, binSize, 
% 												nrowcols, timeLimits, yLimits, 
%												stimulus_times)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
%  Input Args:
%	 spikesByStim		cell array of spiketimes, organized by stimulus type
%							and level
%	 Dinf					data information struct
%	 binSize				size of histogram bins, milliseconds
%	 nrowcols			[nrows ncols] for plotting
%	 timeLimits			xaxis limits for plots [min max]
%	 yLimits				yaxis limits for plots [min max]
%	Optional Inputs:
%	 stimulusTimes		{nstimuli, 1} cell array with 
%							[stim onset, stim offset] times in milliseconds
%	 offset]
%  Output Args:
%	 H		handle to figure
%------------------------------------------------------------------------
% See Also: computeFRA, optoproc
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 24 May 2019 (SJS), adapting from plotPSTHMATRIX & plotPSTH_WAV
%
% Revisions:
%  title string not necessary?
% 7 Jul 2020 (SJS): adding input stimulus times option
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% definitions, input args
%------------------------------------------------------------------------
if ~isempty(varargin)
	stimulusTimes = varargin{1};
else
	stimulusTimes = [];
end

%------------------------------------------------------------------------
%% first, get information about stimuli
%------------------------------------------------------------------------

% get unique stimuli in order they appear in stimList using 'stable' option
% iA will be indices of first occurrence of each unique stim
% iC will identify which of the unique stim is in each row
[uniqueStim, ~, iC]  = unique(Dinf.test.wavlist, 'stable');
nStim = numel(uniqueStim);
% create list of indices into stimList, spiketimes for each stimulus 
% (will use later for many things)
stimIndices = cell(nStim, 1);
for s = 1:nStim
	stimIndices{s} = find(iC == s);
end

% determine if there is a NULL stimulus
if any(strcmpi('NULL', uniqueStim))
	hasNULL = true;
	nullIndex = find(strcmpi('NULL', uniqueStim));
else
	hasNULL = false;
end

%------------------------------------------------------------------------
% first way
%------------------------------------------------------------------------
% dbLevels = Dinf.test.Level;

%------------------------------------------------------------------------
% second way
%------------------------------------------------------------------------
% create cell array to store levels for each stim (might be different
% number of levels for each stim, so can't use array)
dbLevelsByStim = cell(nStim, 1);
nLevels = zeros(nStim, 1);
% loop through stim
fprintf('Determining Stimulus Levels...\n');
for s = 1:nStim
	% ASSUME that repeats of same stimulus are at different levels...
	nLevels(s) = length(stimIndices{s});
	dbLevelsByStim{s} = zeros(nLevels(s), 1);
	fprintf('\tStimulus: %s\n', uniqueStim{s});
	fprintf('\t\tnLevels: %d\n', nLevels(s));
	fprintf('\t\tLevels: ');
	for l = 1:nLevels(s)
		dbLevelsByStim{s}(l) = Dinf.stimList(stimIndices{s}(l)).audio.Level;
		fprintf('%d  ', Dinf.stimList(stimIndices{s}(l)).audio.Level);
	end
	fprintf('\n');
end

% need to check nlevels
if hasNULL && length(unique(nLevels(nullIndex))) > 1
	warning('odd mismatch in expected levels for NULL stimulus');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% setup plots / options
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% create array to hold figure handles - assume max(nLevels) is correct
hPR = cell(max(nLevels), 1);

% determine rows, columns if necessary
if isempty(nrowcols)
	nrowcols = optoproc_figure_plotrowcol(nStim);
end

%------------------------------------------------------------------------
% global options for raster and psth matrix
%------------------------------------------------------------------------
plotopts.timelimits = timeLimits;
% assign y axis limits if provided
if ~isempty(yLimits)
	plotopts.psth_ylimits = yLimits;
end

% options for raster and psth matrix
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 16;
plotopts.raster_color = [0 0 0];
plotopts.psth_binwidth = binSize;
plotopts.plotgap = 0.001;
plotopts.vertgap = 0.035;
plotopts.xlabel = 'msec';
plotopts.stimulus_times_plot = 3;
plotopts.stimulus_on_color{1} = [1 0.4 0.4];
plotopts.stimulus_off_color{1} = [1 0.4 0.4];
plotopts.stimulus_onoff_pct(1) = 0;
% add on off bars for opto stim
if Dinf.opto.Enable
	% add colors for second stimulus
	plotopts.stimulus_on_color{2} = [1 0 0];
	plotopts.stimulus_off_color{2} = [1 0 0];
	plotopts.stimulus_onoff_pct(2) = 80;
end

% plot name
% need to replace back slash by front slash for UNIXy OSes
[~, fbase, fext] = fileparts(strrep(Dinf.filename, '\', '/'));
fname = [fbase fext];

% titles for stimuli
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

% create list of stimulus times
if isempty(stimulusTimes)
	% build list of stimulus times if one is not provided
	plotopts.stimulus_times = cell(nrowcols(1), nrowcols(2));
	for r = 1:nrowcols(1)
		for c = 1:nrowcols(2)
			% need to have [stim_onset stim_offset], so add delay to 
			% [0 duration] to compute proper times. then, multiply by 0.001 to
			% give times in seconds (Dinf values are in milliseconds)
			plotopts.stimulus_times{r, c} = 0.001 * (Dinf.audio.Delay + ...
															[0 Dinf.audio.Duration]);
			% if opto is Enabled, add it to the array by concatenation
			if Dinf.opto.Enable
				plotopts.stimulus_times{r, c} = [plotopts.stimulus_times{r, c}; ...
																0.001 * (Dinf.opto.Delay + ...
																[0 Dinf.opto.Dur]) ];
			end
		end
	end
else
	% assign stimulus times from provided list
	plotopts.stimulus_times = cell(nrowcols(1), nrowcols(2));
	s = 1;
	for r = 1:nrowcols(1)
		for c = 1:nrowcols(2)
			if s <= nStim
				% need to have [stim_onset stim_offset], so add delay to 
				% [0 duration] to compute proper times. then, multiply by 0.001 to
				% give times in seconds (Dinf values are in milliseconds)
				plotopts.stimulus_times{r, c} = 0.001 * stimulusTimes{s};
				% if opto is Enabled, add it to the array by concatenation
				if Dinf.opto.Enable
					plotopts.stimulus_times{r, c} = ...
										[plotopts.stimulus_times{r, c}; ...
											0.001 * (Dinf.opto.Delay + ...
													[0 Dinf.opto.Dur]) ];
				end
			end
			s = s + 1;
		end
	end	
end


% create plot titles
plotopts.plot_titles = cell(nrowcols(1), nrowcols(2));
s = 1;
for r = 1:nrowcols(1)
	for c = 1:nrowcols(2)
		if s <= nStim
			plotopts.plot_titles{r, c} = {sprintf('%s', uniqueStim{s})};
		end
		s = s + 1;
	end
end

% loop through levels
for l = 1:max(nLevels)
	% reset stimulus index
	s = 1;
	
	% add level to 1st plot titles (will be odd for NULL stim, but no easy
	% alternative... also, this assumes that levels are uniform
	plotopts.plot_titles{1, 1} = ...
					{	sprintf('%s  %d dB SPL', fname, Dinf.test.Level(l)) ...
						sprintf('%s', uniqueStim{s})};

	% allocate Spikes cell array to hold spiketimes for stimuli
	Spikes = cell(nrowcols(1), nrowcols(2));

	% loop through stimuli (in rows and cols)
	for r = 1:nrowcols(1)
		for c = 1:nrowcols(2)
			if s <= nStim
				% assign spikes - need to account for NULL stimulus if present,
				% since there will be only 1 level for the NULL stim...
				if hasNULL && (s == nullIndex)
					% only 1 level for null stimulus
					Spikes{r, c} = spikesByStim{stimIndices{s}(1)};
				else
					% check on levels and stimIndices
					if nLevels(s) ~= length(stimIndices{s})
						error('%s: mismatch in nLevels and stimIndices', mfilename);
					end
					Spikes{r, c} = spikesByStim{stimIndices{s}(l)};
				end
			end
			% increment stimulus index
			s = s + 1;
		end
	end
	
	% assign figure handle to hPR
	hPR{l} = figure;
	% plot!
	rasterpsthmatrix(Spikes, plotopts);
	% replace '.' in fbase with 'p' for similar reasons
	set(hPR{l}, 'Name', sprintf('%s_%ddB', strrep(fbase, '.', 'p'), ...
												Dinf.test.Level(l)));
	set(hPR{l}, 'FileName', sprintf('%s_%ddB', strrep(fbase, '.', 'p'), ...
												Dinf.test.Level(l)));
end

if nargout
	varargout{1} = hPR;
end
