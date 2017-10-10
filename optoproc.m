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
clear all %#ok<CLSCR>

%---------------------------------------------------------------------
%% settings for processing data
%---------------------------------------------------------------------
% filter
HPFreq = 350;
LPFreq = 6500;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
% Channel Number (use 8 for single channel data)
channelNumber = 5;

%---------------------------------------------------------------------
%% need to get information about system
%---------------------------------------------------------------------
if ~exist('username', 'file')
	warning('Cannot find <username.m> function... assuming mac for os');
	uname = 'sshanbhag';
	os_type = 'MACI64';
	hname = 'parvati';
else
	[uname, os_type, hname] = username;
end

switch os_type
	case {'PCWIN', 'PCWIN64'}
		% assume we are using the opto computer (optocom)
		data_root_path = 'E:\Data\SJS';
		tytology_root_path = 'C:\TytoLogy';
	
	case {'MAC', 'MACI', 'GLNXA64', 'MACI64'}
		data_root_path = '/Users/sshanbhag/Work/Data/Mouse/Opto';
		tytology_root_path = ...
								'/Users/sshanbhag/Work/Code/Matlab/dev/TytoLogy';
end
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
%% Read Data
%---------------------------------------------------------------------
% add animal and datestring if desired
animal = '1155';
datestring = '20171006';
datafile = '1155_20171006_04_03_3123_FREQoptoON_ch5ch11_3.dat';

% build datapath
datapath = fullfile(data_root_path, animal, datestring);

% datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1157/20170707/';
% datafile = '1157_20170707_01_01_639_BBN_LEVEL_dur100.dat';
% datafile = '1157_20170707_01_01_639_BBN_LEVEL.dat';

% get data file from user
[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
													fullfile(datapath, datafile));
% abort if cancelled
if datafile == 0
	fprintf('Cancelled\n');
	return
end

% add path to opto - needed for 
if ~exist('readOptoData.m', 'file')
	addpath(fullfile(tytology_root_path, 'Experiments', 'Opto'));
end
% 
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(datapath, datafile), ...
											'Filter', [HPFreq LPFreq], ...
											'Channel', channelNumber);
if isempty(D)
	return
end

%---------------------------------------------------------------------
%% get info from filename
%---------------------------------------------------------------------
[~, fname] = fileparts(datafile);
usc = find(fname == '_');
endusc = usc - 1;
startusc = usc + 1;
animal = fname(1:endusc(1));
datecode = fname(startusc(1):endusc(2));
penetration = fname(startusc(2):endusc(3));
unit = fname(startusc(3):endusc(4));
other = fname(startusc(end):end);

%---------------------------------------------------------------------
% create plot output dir
%---------------------------------------------------------------------
plotpath = fullfile(plotpath_base, animal, datecode);
fprintf('Files will be written to:\n\t%s\n', plotpath);
if ~exist(plotpath, 'dir')
	mkdir(plotpath);
end

%---------------------------------------------------------------------
%% determine global RMS and max
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps
nstim = Dinf.test.stimcache.ntrials;
nreps = Dinf.test.stimcache.nreps;
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
%% Some test-specific things...
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
%% find spikes!
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
%% Plot raw data
%---------------------------------------------------------------------
hF = figure(1);
set(hF, 'Name', [fname '_sweeps']);

% determine # of columns of plots
if nvars <= 5
	prows = nvars;
	pcols = 1;	
elseif iseven(nvars)
	prows = nvars/2;
	pcols = 2;
else
	prows = ceil(nvars/2);
	pcols = 2;
end

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
pname = fullfile(plotpath, [fname '_traces']);
savefig(hF, pname, 'compact');
print(hF, pname, '-dpng', '-r300');

%---------------------------------------------------------------------
%% raster, psths
%---------------------------------------------------------------------
hPR = figure(2);
plotopts.timelimits = [0 ceil(max(t))];
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 16;
plotopts.psth_binwidth = 10;
plotopts.plotgap = 0.001;
plotopts.xlabel = 'msec';

% adjust depending on # of columns of plots
if nvars <= 5
	plotopts.plot_titles = titleString;
elseif iseven(nvars)
	plotopts.plot_titles = reshape(titleString, [prows pcols]);
	spiketimes = reshape(spiketimes, [prows pcols]);
else
	titleString = [titleString; {''}];
	spiketimes = [spiketimes; {{}}];
	plotopts.plot_titles = reshape(titleString, [prows pcols]);
	spiketimes = reshape(spiketimes, [prows pcols]);
end

% plot!
rasterpsthmatrix(spiketimes, plotopts)
	
set(hPR, 'Name', fname)

% save plot
pname = fullfile(plotpath, [fname '_rp']);
savefig(hPR, pname, 'compact');
print(hPR, pname, '-dpng', '-r300');
