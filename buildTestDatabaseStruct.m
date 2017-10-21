function P = buildTestDatabaseStruct(varargin)
%------------------------------------------------------------------------
% P = buildTestDatabaseStruct(varargin)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% builds database of test properties as a struct, stores in 
% ~/Work/Data/Mouse/Opto/Analyzed/db (default location)
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
% settings
%---------------------------------------------------------------------
% add comment (query user for it)
ADD_COMMENT = 0;
% save db file
SAVE_DB = 0;
FORCE_SAVE = 0;
datapath = [];
%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
if nargin
	argIndx = 1;
	while argIndx <= nargin
		switch upper(varargin{argIndx})
			case {'PATH', 'DIR'}
				datapath = varargin{argIndx+1};
				argIndx = argIndx + 2;
			case {'ADDCOMMENT', 'COMMENT'}
				if ischar(varargin{argIndx+1})
					if upper(varargin{argIndx+1}(1)) == 'Y'
						ADD_COMMENT = 1;
					elseif upper(varargin{argIndx+1}(1)) == 'N'
						ADD_COMMENT = 0;
					else
						error('%s: invalid AddComment value %s', mfilename, ...
																varargin{argIndx+1});
					end
				else
					if varargin{argIndx+1} == 1
						ADD_COMMENT = 1;
					else
						ADD_COMMENT = 0;
					end
				end
				argIndx = argIndx + 2;
			case {'SAVEDB', 'SAVE_DB'}
				if ischar(varargin{argIndx+1})
					if upper(varargin{argIndx+1}(1)) == 'Y'
						SAVE_DB = 1;
					elseif upper(varargin{argIndx+1}(1)) == 'N'
						SAVE_DB = 0;
					else
						error('%s: invalid AddComment value %s', mfilename, ...
																varargin{argIndx+1});
					end
				else
					if varargin{argIndx+1} == 1
						SAVE_DB = 1;
					else
						SAVE_DB = 0;
					end
				end
			case {'FORCE_SAVE', 'FORCE', 'SAVE'}
				if ischar(varargin{argIndx+1})
					if upper(varargin{argIndx+1}(1)) == 'Y'
						FORCE_SAVE = 1;
					elseif upper(varargin{argIndx+1}(1)) == 'N'
						FORCE_SAVE = 0;
					else
						error('%s: invalid AddComment value %s', mfilename, ...
																varargin{argIndx+1});
					end
				else
					if varargin{argIndx+1} == 1
						FORCE_SAVE = 1;
					else
						FORCE_SAVE = 0;
					end
				end
				if FORCE_SAVE
					SAVE_DB = 1;
				end
				argIndx = argIndx + 2;
			otherwise
				error('%s: invalid option %s', mfilename, varargin{argIndx});
		end
	end
end
%---------------------------------------------------------------------
% need to get information about system
%---------------------------------------------------------------------
data_root_path = optoanalysis_paths;

%---------------------------------------------------------------------
% select data directory if not provided
%---------------------------------------------------------------------
if isempty(datapath)
	% build datapath
	datapath = data_root_path;
	% query user for desired directory
	datapath = uigetdir(datapath, 'Select data directory');
	% abort if cancelled
	if datapath == 0
		fprintf('Cancelled\n');
		P = [];
		return
	end
end
%---------------------------------------------------------------------
% get list of files
%---------------------------------------------------------------------
dList = dir(fullfile(datapath, '*.dat'));
nfiles = length(dList);
if nfiles == 0
	error('%s: no .dat files found in directory %s', mfilename, datapath);
end
%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
P = [];
for f = 1:nfiles
	datafile = dList(f).name;
	structdatum = get_test_properties(fullfile(datapath, datafile));
	
	%---------------------------------------------------------------------
	% update comment
	%---------------------------------------------------------------------
	if ADD_COMMENT
		newTxt = uiaskvalue(	'QuestionText', datafile, ...
									'FigureName', 'Additional comments', ...
									'ValueType', 'char', ...
									'Value', '', ...
									'ValueText', '' );

		if ~isempty(newTxt)
			structdatum.comment = newTxt;
		end
	else
		structdatum.comment = '';
	end
	%---------------------------------------------------------------------
	% create or append to struct array
	%---------------------------------------------------------------------
	if f == 1
		P = structdatum;
	else
		P(f) = structdatum; %#ok<AGROW>
	end
end

%---------------------------------------------------------------------
% Save Data
%---------------------------------------------------------------------
if SAVE_DB
	% find location(s) of file separators in datapath
	sloc = find(datapath == filesep);
	% chars after last filesep should be date
	datestr = datapath( (sloc(end) + 1):end );
	% chars between next to last filesep and last filesep are animal
	animalstr = datapath( (sloc(end-1) + 1):(sloc(end)-1) );
	% build proposed file name
	dfile = fullfile(data_root_path, 'Analyzed', 'db', ...
													[animalstr '_' datestr '_db.mat']);
	% force save?
	if ~FORCE_SAVE
		[dfile, fpath] = uiputfile('*.mat', 'Save test properties to', dfile);
	end
	
	if dfile == 0
		return
	else
		save(fullfile(fpath, dfile), 'P', '-MAT');
	end
end


