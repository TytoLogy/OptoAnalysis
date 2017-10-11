function P = get_test_properties(varargin)
%------------------------------------------------------------------------
% get_test_properties
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% extracts information about data collected by the opto program
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
%	- Functionalize
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% define some things
%--------------------------------------------------------------------------
datafile = [];
P = [];

%--------------------------------------------------------------------------
% check inputs
%--------------------------------------------------------------------------
if nargin
	if exist(varargin{1}, 'file') == 2
		datafile = varargin{1};
	else
		error([mfilename ': datafile ' varargin{1} ' not found.']);
	end
end

%--------------------------------------------------------------------------
% have user select data file if one was not provided
%--------------------------------------------------------------------------
if isempty(datafile)
	[datafile, datapath] = uigetfile('*.dat','Select data file');
	if datafile == 0
		disp('user cancelled datafile load')
		return
	end
	datafile = fullfile(datapath, datafile);
end

%--------------------------------------------------------------------------
% read in data
%--------------------------------------------------------------------------
fprintf('Reading test data from %s\n', datafile);

% open file
fp = fopen(datafile, 'r');

% read the header
try
	P.Dinf = readOptoDataFileHeader(fp);
catch errMsg
	errMsg.message
	fclose(fp);
	return;
end
fclose(fp);

% check if stimcache exists
if isfield(P.Dinf.test, 'stimcache')
	% if so, get # of reps and trials from size of trialRandomSequence
	[numreps, numtrials] = size(P.Dinf.test.stimcache.trialRandomSequence);
else
	% otherwise, compute from Reps and nCombinations
	numreps = P.Dinf.test.Reps;
	numtrials = P.Dinf.test.nCombinations;
end
fprintf('%s: had %d reps, %d trials\n', datafile, numreps, numtrials);

%---------------------------------------------------------------------
% get info from filename
%---------------------------------------------------------------------
P.filename = datafile;
P.fname_props = parse_opto_filename(datafile);

%---------------------------------------------------------------------
% determine global RMS and max
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps
P.nstim = P.Dinf.test.stimcache.ntrials;
P.nreps = P.Dinf.test.stimcache.nreps;




