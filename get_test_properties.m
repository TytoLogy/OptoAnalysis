function [P, varargout] = get_test_properties(varargin)
%--------------------------------------------------------------------------
% get_test_properties
%--------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% extracts information about data collected by the opto program, ret
%
%--------------------------------------------------------------------------
% See Also:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 11 October, 2017
%
% Revisions:
%	see git!
%--------------------------------------------------------------------------
% TO DO:
%	- Document
%--------------------------------------------------------------------------

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% define some things
%---------------------------------------------------------------------
%---------------------------------------------------------------------
datafile = [];

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% check inputs
%---------------------------------------------------------------------
%---------------------------------------------------------------------
if nargin
	if exist(varargin{1}, 'file') == 2
		datafile = varargin{1};
	else
		error([mfilename ': datafile ' varargin{1} ' not found.']);
	end
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% have user select data file if one was not provided
%---------------------------------------------------------------------
%---------------------------------------------------------------------
if isempty(datafile)
	[datafile, datapath] = uigetfile('*.dat','Select data file');
	if datafile == 0
		disp('user cancelled datafile load')
		return
	end
	datafile = fullfile(datapath, datafile);
end

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% read in data
%---------------------------------------------------------------------
%---------------------------------------------------------------------
fprintf('Reading test data from %s\n', datafile);
%---------------------------------------------------------------------
% open file
%---------------------------------------------------------------------
fp = fopen(datafile, 'r');
%---------------------------------------------------------------------
% read the header
%---------------------------------------------------------------------
try
	Dinf = readOptoDataFileHeader(fp);
catch errMsg
	errMsg.message
	fclose(fp);
	return;
end
fclose(fp);

%--------------------------------------------------------------------------
% get stimulus information
%--------------------------------------------------------------------------
% build wavinfo_matfile name
[datapath, basename] = fileparts(datafile);
wavinfo_matfile = fullfile(datapath, [basename '_wavinfo.mat']);
% check if wavinfo exists
if exist(wavinfo_matfile, 'file')
	fprintf('Loading stimList from %s\n', wavinfo_matfile);
	load(wavinfo_matfile, 'stimList');
	Dinf.stimList = stimList;
else
	Dinf.stimList = [];
end

%---------------------------------------------------------------------
% check if stimcache exists
%---------------------------------------------------------------------
if isfield(Dinf.test, 'stimcache')
	%-----------------------------------------------------------------
	% if so, get # of reps and trials from size of trialRandomSequence
	%-----------------------------------------------------------------
	[numreps, numtrials] = size(Dinf.test.stimcache.trialRandomSequence);
else
	%-----------------------------------------------------------------
	% otherwise, compute from Reps and nCombinations
	%-----------------------------------------------------------------
	numreps = Dinf.test.Reps;
	numtrials = Dinf.test.nCombinations;
end
fprintf('%s: had %d reps, %d trials\n', datafile, numreps, numtrials);

% correct Dinf test type
Dinf = correctTestType(Dinf);

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% assign data to output struct
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% get info from filename
%---------------------------------------------------------------------
[fpath, fname, fext] = fileparts(datafile);
P.filepath = fpath; 
P.filename = [fname fext];
fname_props = parse_opto_filename(datafile);
fnames = fieldnames(fname_props);
for n = 1:length(fnames)
	P.(fnames{n}) = fname_props.(fnames{n});
end

%---------------------------------------------------------------------
% test things
%---------------------------------------------------------------------
P.nstim = numtrials;
P.nreps = numreps;
% not all tests have Type as field (e.g., standalone OptoInhibxxx),
% so check that it is present
if isfield(Dinf.test, 'Type')
	P.TestType = char(Dinf.test.Type);
else
	% if not, just use the test name
	P.TestType = char(Dinf.test.Name);
end
P.TestName = char(Dinf.test.Name);
testfields = {'Reps', 'Randomize', 'Block', 'saveStim', ...
					'AcqDuration', 'SweepPeriod'};
for n = 1:length(testfields)
	if isfield(Dinf.test, testfields{n})
		P.(testfields{n}) = Dinf.test.(testfields{n});
	else
		P.(testfields{n}) = [];
	end
end
if isfield(Dinf.test, 'stimseq')
	P.stimseq = reshape(Dinf.test.stimseq', [1, numel(Dinf.test.stimseq)]);
elseif isfield(Dinf.test, 'stimIndices')
	if isfield(Dinf, 'stimList')
		P.stimseq = {Dinf.stimList, force_row(Dinf.test.stimIndices)};
	else
		P.stimseq = force_row(Dinf.test.stimIndices);
	end
else
	disp('cannot find stimulus sequence');
	P.stimseq = [];
end
%---------------------------------------------------------------------
% audio
%---------------------------------------------------------------------
P.Audio_Signal = Dinf.audio.signal;
P.Audio_Signal.Type = char(P.Audio_Signal.Type);
P.Audio_Delay = Dinf.audio.Delay;
P.Audio_Dur = Dinf.audio.Duration;
P.Audio_Level = Dinf.audio.Level;
P.Audio_Ramp = Dinf.audio.Ramp;
P.Audio_Frozen = Dinf.audio.Frozen;
P.Audio_ISI = Dinf.audio.ISI;
%---------------------------------------------------------------------
% opto
%---------------------------------------------------------------------
P.Opto_Enable = Dinf.opto.Enable;
P.Opto_Delay = Dinf.opto.Delay;
P.Opto_Dur = Dinf.opto.Dur;
P.Opto_Amp = Dinf.opto.Amp;
%---------------------------------------------------------------------
% record channels
%---------------------------------------------------------------------
P.nChannels = Dinf.channels.nRecordChannels;
% need to make sure all elements of P are in a single row or else
% conversion from struct to table will barf
P.ChannelList = force_row(Dinf.channels.RecordChannelList);
%---------------------------------------------------------------------
% I/O 
%---------------------------------------------------------------------
P.iFs = Dinf.indev.Fs;
P.oFs = Dinf.outdev.Fs;
%---------------------------------------------------------------------
% comment is a place holder
%---------------------------------------------------------------------
P.comment = 'NaN';
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% assign to outputs
%---------------------------------------------------------------------
%---------------------------------------------------------------------

if nargout > 1
	if nargout >= 2
		varargout{1} = Dinf;
	end
	if nargout == 3
		varargout{2} = struct2table(P);
	end
	if nargout > 3
		error('incorrect outputs');
	end
end


