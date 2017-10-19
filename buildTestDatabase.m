%------------------------------------------------------------------------
% buildTestDatabase
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% builds database of test properties
%
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 October, 2017
%
% Revisions:
%	see git!
%------------------------------------------------------------------------
% TO DO:
%	- Document
%--------------------------------------------------------------------------


%---------------------------------------------------------------------
% need to get information about system
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% select data directory
%---------------------------------------------------------------------
% add animal and datestring if desired
animal = '1155';
datestring = '';


% build datapath
datapath = fullfile(data_root_path, animal, datestring);

datapath = uigetdir(datapath, 'Select data directory');

if isempty(datafile)
	% get data file from user
	[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
														fullfile(datapath, datafile));
	% abort if cancelled
	if datafile == 0
		fprintf('Cancelled\n');
		return
	end
end

%---------------------------------------------------------------------
%% Read Data
%---------------------------------------------------------------------
[tabledatum, Dinf] = get_test_properties(fullfile(datapath, datafile));

%---------------------------------------------------------------------
% update comment
%---------------------------------------------------------------------
newTxt = uiaskvalue(	'QuestionText', datafile, ...
							'FigureName', 'Additional comments', ...
							'ValueType', 'char', ...
							'Value', '', ...
							'ValueText', '' );






