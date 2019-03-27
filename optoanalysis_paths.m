function [data_root_path, tytology_root_path] = optoanalysis_paths()
%------------------------------------------------------------------------
% optoanalysis_paths
%------------------------------------------------------------------------
% Opto Analysis
%--------------------------------------------------------------------------
% sets up paths for opto analysis
%
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 19 October, 2017
%
% Revisions:
%	see git!
%------------------------------------------------------------------------

%---------------------------------------------------------------------
% need to get information about system
%---------------------------------------------------------------------
if ~exist('username', 'file')
	warning('Cannot find <username.m> function... assuming mac for os');
	uname = 'sshanbhag'; %#ok<NASGU>
	os_type = 'MACI64';
	hname = 'parvati'; %#ok<NASGU>
else
	[uname, os_type, hname] = username; %#ok<ASGLU>
end

switch os_type
	case {'PCWIN', 'PCWIN64'}
		% assume we are using the opto computer (optocom)
		data_root_path = 'E:\Data';
		tytology_root_path = 'C:\TytoLogy';

	case {'MAC', 'MACI', 'GLNXA64', 'MACI64'}
		data_root_path = '/Users/sshanbhag/Work/Data/Mouse/Opto';
		tytology_root_path = ...
								'/Users/sshanbhag/Work/Code/Matlab/dev/TytoLogy';
end

%---------------------------------------------------------------------
% add path to opto - needed for 
%---------------------------------------------------------------------
if ~exist('readOptoData.m', 'file')
	path_to_add = fullfile(tytology_root_path, 'Experiments', 'Opto');
	fprintf('%s:\n\tAdding\n\t\t%s\n\tto caller path\n', ...
						mfilename, path_to_add);
	expr = [	'addpath(''' path_to_add ''')' ];
	evalin('base', expr);
end

% %---------------------------------------------------------------------
% % data file things
% %---------------------------------------------------------------------
% % need to get information about system
% if ~exist('username', 'file')
% 	warning('Cannot find <username.m> function... assuming mac for os');
% 	os_type = 'MACI64';
% else
% 	[~, os_type, ~] = username;
% end
% switch os_type
% 	case {'PCWIN', 'PCWIN64'}
% 		% assume we are using the opto computer (optocom)
% 		data_root_path = 'E:\Data\SJS';
% 	case {'MAC', 'MACI', 'GLNXA64', 'MACI64'}
% 		data_root_path = '/Users/sshanbhag/Work/Data/Mouse/Opto';
% end
