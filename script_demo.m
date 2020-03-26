%------------------------------------------------------------------------
% example optoproc script
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% define files and directory
%------------------------------------------------------------------------

% path to data
DataPath = '~/Work/Data/TestData/MT';

% path to save plots
OutputPath = '~/tmp';

% list of data files
DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'; ...
 					'1382_20191212_02_02_3200_BBN.dat'; ...
 					'1382_20191212_02_02_3200_WAV.dat';	};
				
% activity on channels = [4, 5, 11, 14];
Channels = [4, 5, 11, 14];

% define threshold
RMSThreshold = 5;

%% Loop through channels

% for c = 1:length(Channels)
for c = 1:2
	
	plotpath = fullfile(OutputPath, sprintf('Channel%d', Channels(c)));
	mkdir(plotpath);
	
	%------------------------------------------------------------------------
	% plot freq tuning data!
	%------------------------------------------------------------------------
	optoproc('plotFTC', 'Channel', Channels(c), 'Threshold', RMSThreshold, ...
				'savePNG', ...
				'file', fullfile(DataPath, DataFile{1}), 'Plotpath', plotpath)

	%------------------------------------------------------------------------
	% plot BBN rate level data!
	%------------------------------------------------------------------------
	optoproc('plotRLF', 'Channel', Channels(c), 'Threshold', RMSThreshold, ...
				'savePNG', ...
				'file', fullfile(DataPath, DataFile{2}), 'Plotpath', plotpath)

	%------------------------------------------------------------------------
	% plot WAV by level data!
	%------------------------------------------------------------------------
	optoproc('plot_psth_by_level', 'Channel', Channels(c), 'Threshold', RMSThreshold, ...
				'savePNG', ...
				'file', fullfile(DataPath, DataFile{3}), 'Plotpath', plotpath)

end