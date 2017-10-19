function [T, varargout] = get_test_properties(varargin)
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
T = [];

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

%---------------------------------------------------------------------
%---------------------------------------------------------------------
% assign data to output struct
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% get info from filename
%---------------------------------------------------------------------
P.filename = datafile;
fname_props = parse_opto_filename(datafile);
fnames = fieldnames(fname_props);
for n = 1:length(fnames)
	P.(fnames{n}) = fname_props.(fnames{n});
end
%---------------------------------------------------------------------
% get  # of stimuli (called ntrials by opto) as well as # of reps
%---------------------------------------------------------------------
P.nstim = Dinf.test.stimcache.ntrials;
P.nreps = Dinf.test.stimcache.nreps;
%---------------------------------------------------------------------
% test things
%---------------------------------------------------------------------
% P.Test = rmfield(Dinf.test, 'stimcache');
% P.Test.Type = char(Dinf.test.Type);
% P.Test.Name = char(Dinf.test.Name);


P.TestType = char(Dinf.test.Type);
P.TestName = char(Dinf.test.Name);
testfields = {'Reps', 'Randomize', 'Block', 'saveStim', ...
					'AcqDuration', 'SweepPeriod'};
for n = 1:length(testfields)
	P.(testfields{n}) = Dinf.test.(testfields{n});
end
P.stimseq = reshape(Dinf.test.stimseq', [1, numel(Dinf.test.stimseq)]);
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
P.comment = ' ';
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% assign to outputs
%---------------------------------------------------------------------
%---------------------------------------------------------------------
T = struct2table(P);

if nargout > 1
	if nargout >= 2
		varargout{1} = Dinf;
	end
	if nargout == 3
		varargout{2} = P;
	end
	if nargout > 3
		error('incorrect outputs');
	end
end


