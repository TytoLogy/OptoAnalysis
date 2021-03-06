function varargout = optoplot(varargin)
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
clearvars

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
channelNumber = 8;
% binSize for PSTH (milliseconds)
binSize = 5;
% SAVE PLOTS?
saveFIG = 0;
savePNG = 0;
savePDF = 1;

%---------------------------------------------------------------------
% need to get information about system
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths; %#ok<ASGLU>
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');
	
%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
% add animal and datestring if desired
animal = '1155';
datestring = '20171009';
% build datapath
datapath = fullfile(data_root_path, animal, datestring);
% get data file from user
[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', datapath);
% abort if cancelled
if datafile == 0
	fprintf('Cancelled\n');
	return
end

[~, Dinf] = optoproc('file', fullfile(datapath, datafile), ...
								'plotpath', plotpath_base, ...
								'channel', 5, 'binsize', 5, 'savepdf');

