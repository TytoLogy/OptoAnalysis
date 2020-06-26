%------------------------------------------------------------------------
% wavproc.m
%------------------------------------------------------------------------
% opto:scripts
%--------------------------------------------------------------------------
% example script for plotting opto WAV data (from MTwav_standalone.m)
%
%------------------------------------------------------------------------
% See Also: optoproc, opto (TytoLogy:opto program)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 26 Jun 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% define files and settings
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% path to data
DataPath = '/Users/sshanbhag/Work/Data/TestData/RigTests/000/20200626';
% path to save plots
OutputPath = DataPath;
% list of data files
DataFile = {	'000_20200626_0_0_0_WAV.dat';	...
	'000_20200626_0_0_0_FRA.dat';};
% activity on channels = [1:16];
Channels = [1, 2];
% define threshold
RMSThreshold = 1;



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
	
	%{
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
	%}
	%{
	%------------------------------------------------------------------------
	% plot WAV by level data!
	%------------------------------------------------------------------------
	optoproc('plot_psth_by_level', 'Channel', Channels(c), ...
				'Threshold', RMSThreshold, ...
				'savePNG', ...
				'file', fullfile(DataPath, DataFile{1}), 'Plotpath', plotpath)
	%}		
	%------------------------------------------------------------------------
	% plot FRA !
	%------------------------------------------------------------------------
	% cread plot filename: add channel number to filename
	% first, break up filename
	F = parse_opto_filename(DataFile{2});
	plotfile = [F.animal '_' F.datecode '_' F.unit '_' F.penetration '_' ...
					F.depth '_Ch' num2str(Channels(c)) '_' F.other '.png'];	
	% assemble file
	optoproc('plot_fra', 'Channel', Channels(c), ...
				'Threshold', RMSThreshold, ...
				'savePNG', ...
				'file', fullfile(DataPath, DataFile{2}), ...
				'Plotpath', plotpath, 'Plotfile', plotfile);
			
end