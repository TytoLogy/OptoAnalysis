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

% query user for desired directory
datapath = uigetdir(datapath, 'Select data directory');

% abort if cancelled
if datapath == 0
	fprintf('Cancelled\n');
	return
end

%---------------------------------------------------------------------
% list files
%---------------------------------------------------------------------
dList = dir(fullfile(datapath, '*.dat'));

nfiles = length(dList);

if nfiles == 0
	error('%s: no .dat files found in directory %s', mfilename, datapath);
end
%---------------------------------------------------------------------
%% Read Data
%---------------------------------------------------------------------

P = [];

for f = 1:nfiles
	datafile = dList(f).name;
	structdatum = get_test_properties(fullfile(datapath, datafile));

% 	%---------------------------------------------------------------------
% 	% update comment
% 	%---------------------------------------------------------------------
% 	newTxt = uiaskvalue(	'QuestionText', datafile, ...
% 								'FigureName', 'Additional comments', ...
% 								'ValueType', 'char', ...
% 								'Value', '', ...
% 								'ValueText', '' );
% 
% 	if ~isempty(newTxt)
% 		structdatum.comment = newTxt;
% 	end

	if f == 1
		P = structdatum;
	else
		P(f) = structdatum;
	end

end

