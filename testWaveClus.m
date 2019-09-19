% script to test wave_clus spike sorting

%% make sure it is on path
if ~exist('wave_clus', 'file')
	% if not, add paths
	fprintf('Adding wave_clus paths...\n');
	CODEBASEDIR = '~/Work/Code/Matlab/dev/Analysis/';
	WAVECLUS_BASEDIR = 'wave_clus';
	WAVECLUS_SUBDIR = {	'Batch_files', ...
								'Gui_functions', ...
								'Raw_data_readers', ...
								'SPC', ...
								'figs', ...
								'tools', ...
								'wave_clus_font' ...
							};
	addpath(fullfile(CODEBASEDIR, WAVECLUS_BASEDIR));
	for sd = WAVECLUS_SUBDIR
		addpath(	fullfile(CODEBASEDIR, WAVECLUS_BASEDIR, sd{1}));
	end
	fprintf('...done\n');
end

%% test method to convert 1 channel of data
% this should be included in exportChannelForSorting...

% input dir and file
datadir = '/Users/sshanbhag/Work/Data/Mouse/IC/1344';
datafile = '1344_20190916_04_02_1532_FREQ_TUNING.dat';

% output dir
outputdir = pwd;

% channel(s) of data to obtain
% channel = 8;
channel = [1 2 3];

% can probably just use readOptoData
% % the data will be further processed during sorting. for now, use a fairly
% % broad filter
% filtband = [5000 10000];
% % get the data - only need the raw data (stored in D cell array of each
% % stimulus presentation) and the information in Dinf about the test.
% [D, Dinf] = getFilteredOptoData(fullfile(datadir, datafile), ...
% 											'filter', filtband, ...
% 											'channel', channel);
										
[D, Dinf] = readOptoData(fullfile(datadir, datafile));
% Fix test info
Dinf = correctTestType(Dinf);
fprintf('Test type: %s\n', Dinf.test.Type);

%% quick plot of data for one trial

% need colorspace to distinguish channels
lineStyles = linspecer(16);
figure(10)
axes('NextPlot','replacechildren', 'ColorOrder', lineStyles);
plot(D(1).datatrace);
for c = 1:16
	tl{c} = num2str(c);
end
legend(tl);



%% Get data for channel(s) and write to file(s)

% get basename for output file
[~, outputfile_base] = fileparts(datafile);

% number of traces/stimulus presentations
nstims = length(D);
% use length of first response to determine length of traces
nsamples = length(D{1}.datatrace);
% number of channels to loop
nchan = length(channel);
% sample rate
sr = Dinf.indev.Fs;

% loop through channels
for cIndx = 1:nchan

	% allocate temporary data storage for current channel
	tmpC = zeros(nsamples, nstims);
	for s = 1:nstims
		tmpC(:, s) = D{s}.datatrace(:, channel(cIndx));
	end

	% convert to single-row 
	data = reshape(tmpC, 1, []);

	% write to file
	outputfile = sprintf('%s_C%d.mat', outputfile_base, channel(cIndx));
	fprintf('Writing %d points to %s\n', length(data), outputfile)
	save(fullfile(outputdir, outputfile), '-MAT', 'data', 'sr');
end


