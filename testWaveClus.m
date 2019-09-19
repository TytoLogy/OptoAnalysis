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
datadir = '/Users/sshanbhag/Work/Data/Mouse/IC/1344';
datafile = '1344_20190916_04_02_1532_FREQ_TUNING.dat';

% channel of data to obtain
channel = 8;

% can probably just use readOptoData
% % the data will be further processed during sorting. for now, use a fairly
% % broad filter
% filtband = [5000 10000];
% % get the data - only need the raw data (stored in D cell array of each
% % stimulus presentation) and the information in Dinf about the test.
% [D, Dinf] = getFilteredOptoData(fullfile(datadir, datafile), ...
% 											'filter', filtband, ...
% 											'channel', channel);
										
readOptoData(fullfile(datadir, datafile));
%%
