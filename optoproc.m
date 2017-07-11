%% settings for processing data
HPFreq = 350;
LPFreq = 6500;
 
%% Read Data
datapath = '/Users/sshanbhag/Work/Data/Mouse/Opto/1157/20170707/';
datafile = '1157_20170707_01_01_639_BBN_LEVEL_dur100.dat';
[D, Dinf, tracesByStim] = getFilteredOptoData(fullfile(datapath, datafile), ...
																	[HPFreq LPFreq]);
if isempty(D)
	return
end

[~, fname] = fileparts(datafile);

%% determine global RMS and max
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
% compute overall mean rms
mean_rms = mean(reshape(rmsvals, numel(rmsvals), 1));
fprintf('Mean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(max(maxvals));
fprintf('Global max abs value: %.4f\n', global_max);


%% Some test-specific things...
switch upper(Dinf.test.Type)
    case 'FREQ'
        % list of frequencies, and # of freqs tested
        varlist = Dinf.test.stimcache.vrange;
        nvars = length(varlist);
		  titleString = cell(nvars, 1);
		  for v = 1:nvars
			  titleString{v} = {fname, ...
										sprintf('Frequency = %d kHz', 0.001*varlist(v))};
		  end
    case 'LEVEL'
        % list of legvels, and # of levels tested
        varlist = Dinf.test.stimcache.vrange;
        nvars = length(varlist);
		  titleString = cell(nvars, 1);
		  for v = 1:nvars
			  titleString{v} = {fname, sprintf('Level = %d dB SPL', varlist(v))};
		  end
    case 'OPTO'
		 % not yet implemented
    case 'WAVFILE'
        % get list of stimuli (wav file names)
        varlist = Dinf.test.wavlist;
        nvars = length(varlist);
		  titleString = cell(nvars, 1);
		  for v = 1:nvars
			  titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
		  end
    otherwise
        error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end
 
%% Plot raw data
Fs = Dinf.indev.Fs;
hF = figure(1);
set(hF, 'Name', [fname '_sweeps']);
for v = 1:nvars
	% time vector for plotting
	t = (1000/Fs)*((1:length(tracesByStim{v}(:, 1))) - 1);
	subplot(nvars, 1, v);
	% flip tracesByStim in order to have sweeps match the raster plots
	stackplot(t, fliplr(tracesByStim{v}), 'colormode', 'black', ...
											'ymax', global_max, ...
											'figure', hF, 'axes', gca);
	title(titleString{v}, 'Interpreter', 'none');
   xlabel('ms')
   ylabel('Trial')
end

%% find spikes!
spiketimes = cell(nvars, 1);
for v = 1:nvars
	spiketimes{v} = spikeschmitt2(tracesByStim{v}', 4.5*mean_rms, 1, Fs);
	for r = 1:length(spiketimes{v})
		spiketimes{v}{r} = (1000/Fs)*spiketimes{v}{r};
	end
end


%% raster, psths
figure(2)
plotopts.timelimits = [0 ceil(max(t))];
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 16;
plotopts.psth_binwidth = 10;
plotopts.plotgap = 0.001;
plotopts.xlabel = 'msec';
plotopts.plot_titles = titleString;
plotopts.plot_titles{1} = {fname, plotopts.plot_titles{1}};
rasterpsthmatrix(spiketimes, plotopts)
set(gcf, 'Name', fname)

