dname = '/Users/sshanbhag/Work/Data/Mouse/IC/1302/20190507';
fname = '1302_20190507_03_03_711.9_WAV.dat';

if ~exist('D', 'var')
	if ~exist('wavproc.mat', 'file')
		optoproc('file', fullfile(dname, fname), 'plotPSTH');
	end
	load('wavproc.mat')
end

%% -----------rest is optoproc_plotPSTH_WAV code 
%  H = optoproc_plotPSTH_WAV(spiketimes, Dinf, binSize, nvars, varlist, timeLimits, ...
%												yLimits, titleString)

% get unique stimuli in order they appear in stimList
unique_stim = unique(Dinf.test.wavlist, 'stable');
nUnique = numel(unique_stim);

% create array to hold figure handles
hPR = cell(nUnique, 1);


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

for p = 1:nUnique
	% local copy of stimuli
	
	% get # of levels for this stimulus
	
	
	
	
	
	
end
