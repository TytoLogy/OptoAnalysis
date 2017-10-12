function varargout = optoproc(varargin)
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


%---------------------------------------------------------------------
% data file things
%---------------------------------------------------------------------
% need to get information about system
if ~exist('username', 'file')
	warning('Cannot find <username.m> function... assuming mac for os');
	os_type = 'MACI64';
else
	[~, os_type, ~] = username;
end
switch os_type
	case {'PCWIN', 'PCWIN64'}
		% assume we are using the opto computer (optocom)
		data_root_path = 'E:\Data\SJS';
	case {'MAC', 'MACI', 'GLNXA64', 'MACI64'}
		data_root_path = '/Users/sshanbhag/Work/Data/Mouse/Opto';
end
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

