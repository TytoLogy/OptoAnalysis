

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

% add path to opto - needed for 
if ~exist('readOptoData.m', 'file')
	addpath(fullfile(tytology_root_path, 'Experiments', 'Opto'));
end


%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
% add animal and datestring if desired
animal = '1155';
datestring = '20171006';
datafile = '';
% datafile = '1155_20171006_04_03_3123_FREQoptoON_ch5ch11_3.dat';

% build datapath
datapath = fullfile(data_root_path, animal, datestring);

% get data file from user
[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
													fullfile(datapath, datafile));
% abort if cancelled
if datafile == 0
	fprintf('Cancelled\n');
	return
end

P = get_test_properties(fullfile(datapath, datafile));

