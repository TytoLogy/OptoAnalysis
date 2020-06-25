function varargout = optoproc_plotPSTH_byWAV(spikesByStim, Dinf, binSize, ...
												timeLimits, yLimits)
%------------------------------------------------------------------------
%  H = optoproc_plotPSTH_WAV(spikesByStim, Dinf, binSize, timeLimits, ...
%												yLimits)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
%  Input Args:
%	 
%
%  Output Args:
%	 H		handle to figure
%------------------------------------------------------------------------
% See Also: computeFRA, optoproc
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 22 May 2019 (SJS), adapting from plotPSTHMATRIX
%
% Revisions:
% 23 May 2019 (SJS): working... ???
%------------------------------------------------------------------------


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

%------------------------------------------------------------------------
%% figure out attenuation levels
%------------------------------------------------------------------------
% could do this from:
% a) Dinf.test.Level
% 	BUT, this may not apply to all stimuli (e.g., 'NULL')
% b) look at individual stimList entries:Dinf.stimList(n).audio.Level
% 	should be one entry per level/stimulus combination
% 	laborious?

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
	for n = 1:nLevels(s)
		dbLevelsByStim{s}(n) = Dinf.stimList(stimIndices{s}(n)).audio.Level;
		fprintf('%d  ', Dinf.stimList(stimIndices{s}(n)).audio.Level);
	end
	fprintf('\n');
end

%------------------------------------------------------------------------
%% setup plots / options
%------------------------------------------------------------------------

% create array to hold figure handles
hPR = cell(nStim, 1);

% global options for raster and psth matrix
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
plotopts.xlabel = 'msec';
plotopts.stimulus_times_plot = 3;
plotopts.stimulus_on_color{1} = [0 0 1];
plotopts.stimulus_off_color{1} = [0 0 1];
plotopts.stimulus_onoff_pct(1) = 60;
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

% Loop through stim
for s = 1:nStim	
	% check on levels and stimIndices
	if nLevels(s) ~= length(stimIndices{s})
		error('%s: mismatch in nLevels and stimIndices', mfilename);
	end
	
	% create list of stimulus times
	plotopts.stimulus_times = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		% need to have [stim_onset stim_offset], so add delay to 
		% [0 duration] to compute proper times. then, multiply by 0.001 to
		% give times in seconds (Dinf values are in milliseconds)
		plotopts.stimulus_times{s, 1} = 0.001 * (Dinf.audio.Delay + ...
															[0 Dinf.audio.Duration]);
		% if opto is Enabled, add it to the array by concatenation
		if Dinf.opto.Enable
			plotopts.stimulus_times{s, 1} = [plotopts.stimulus_times{s, 1}; ...
												 0.001 * (Dinf.opto.Delay + ...
															[0 Dinf.opto.Dur]) ];
		end
	end
	% create plot titles
	plotopts.plot_titles = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		if l == 1
			plotopts.plot_titles{l} = {fname, sprintf('%s', uniqueStim{s})};
		else
			plotopts.plot_titles{l} = '';
		end
	end

	% assign figure index
	hPR{s} = figure;

	% assign spikes
	Spikes = cell(nLevels(s), 1);
	for l = 1:nLevels(s)
		Spikes{l} = spikesByStim{stimIndices{s}(l)};
	end

	rasterpsthmatrix(Spikes, plotopts);
	% replace '.' in fbase with 'p' for similar reasons
	set(hPR{s}, 'Name', sprintf('%s_%s', strrep(fbase, '.', 'p'), ...
												uniqueStim{s}));
	set(hPR{s}, 'FileName', sprintf('%s_%s', strrep(fbase, '.', 'p'), ...
												uniqueStim{s}));
end


if nargout
	varargout{1} = hPR;
end
