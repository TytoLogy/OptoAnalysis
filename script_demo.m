%------------------------------------------------------------------------
% script_demo.m
%------------------------------------------------------------------------
% TytoLogy:OptoAnalysis
%--------------------------------------------------------------------------
% example script for batch-processing data files from the opto experiment
% program.
%
%------------------------------------------------------------------------
% See Also: optoproc, opto (TytoLogy:opto program)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ? (SJS)
%
% Revisions:
%	12 Jun 2020 (SJS): added FRA per MT changes, updated comments
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% define files and directory
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% sample
%------------------------------------------------------------------------
%{
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
%}
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Working from MT
%------------------------------------------------------------------------
% path to data
DataPath = 'D:\Toolkit\Backup\Documents\IC+MG-data\1425\20200604';
% path to save plots
OutputPath = 'D:\Toolkit\Backup\Documents\IC+MG-data\1425\20200604';
% list of data files
DataFile = {	'1425_20200604_01_01_550_FREQ_TUNING.dat'; ...
 					'1425_20200604_01_01_550_BBN.dat'; ...
 					'1425_20200604_01_01_550_WAV.dat'; ...
                '1425_20200604_01_01_550_FRA.dat';	};
% activity on channels = [1:16];
Channels = [1, 2];
% define threshold
RMSThreshold = 5;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Loop through channels and create plots for each file
%------------------------------------------------------------------------
%------------------------------------------------------------------------
for c = 1:length(Channels)
	%------------------------------------------------------------------------
	% specify and create a directory for each channel's data
	%------------------------------------------------------------------------
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
