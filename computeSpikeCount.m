function varargout = computeSpikeCount(spikesByStim, Dinf, Winf)
%------------------------------------------------------------------------
% [C, dBLevelsByStim] = computeSpikeCount(spikesByStim, Dinf, Winf, countWindow)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
%  Input Args:
%	 spikesByStim		cell array of spiketimes, organized by stimulus type
%							and level
%	 Dinf					data information struct
%	 Winf					WAV stimulus info (from the test's wavinfo.mat file)
%	
%  Output Args:
%------------------------------------------------------------------------
% See Also: computeFRA, optoproc, optoproc_plotPSTH_WAVbyLevel
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 August 2019 (SJS), adapting from optoproc_plotPSTH_WAVbyLevel
%
% Revisions:
% 
%------------------------------------------------------------------------

% hard wav stimulus onset, milliseconds
WAV_ONSET = 45;
% fixed latency for adding to stimulus offset for spike count, milliseconds
FIXED_LATENCY = 0;

%------------------------------------------------------------------------
% first, get information about stimuli
%------------------------------------------------------------------------
% get unique stimuli in order they appear in stimList using 'stable' option
% iA will be indices of first occurrence of each unique stim
% iC will identify which of the unique stim is in each row
[uniqueStim, ~, iC]  = unique(Dinf.test.wavlist, 'stable');
nStim = numel(uniqueStim);
% create list of indices into stimList, spiketimes for each stimulus 
% (will use later for many things)
stimIndices = cell(nStim, 1);
for s = 1:nStim
	stimIndices{s} = find(iC == s);
end

% determine if there is a NULL stimulus - this will be important when
% dealing with the different stimulus levels (only 1 level for NULL)
if any(strcmpi('NULL', uniqueStim))
	hasNULL = true;
	nullIndex = find(strcmpi('NULL', uniqueStim));
else
	hasNULL = false;
end

%------------------------------------------------------------------------
% need to determine stimulus level information
%------------------------------------------------------------------------
% dBLevelsByStim will be a cell array of size {nStim, 1}
% each element will be a vector of db SPL levels for each stimulus
% nLevels will be the # of levels for each stimulus
[dBLevelsByStim, nLevels] = opto_find_stimlevels(uniqueStim, stimIndices, Dinf);
% need to check nlevels
if hasNULL && length(unique(nLevels(nullIndex))) > 1
	warning('odd mismatch in expected levels for NULL stimulus');
end

%------------------------------------------------------------------------
% setup 
%------------------------------------------------------------------------

% create array to hold raw spike counts is correct
C = cell(nStim, max(nLevels));

% build list of stimulus onset/offset times
stimOnsetOffset = zeros(nStim, 2);
for s = 1:nStim
	switch uniqueStim{s}
		case 'null'
			% for null stimulus use whole window
			onset = 0;
			offset = Dinf.test.AcqDuration;
		case 'BBN'
			% for BBN, use [delay delay+duration]
			onset = Winf.noise.Delay;
			offset = Winf.noise.Delay + Winf.noise.Duration + FIXED_LATENCY;
		otherwise
			% for WAV, need to be specific so find index to wavInfo
			wx = strcmp([uniqueStim{s} '.wav'], {Winf.wavInfo.Filename});
			if isempty(wx)
				error('%s: could not find %s in wavInfo!', ...
										mfilename, uniqueStim{s})
			end
			W = Winf.wavInfo(wx);
			onset = Dinf.audio.Delay + WAV_ONSET;
			offset = onset + bin2ms(W.OffsetBin, W.SampleRate) + FIXED_LATENCY;
	end
	% store onset/offset
	stimOnsetOffset(s, :) = [onset ceil(offset)];
end

%% process data

% loop through levels
for l = 1:max(nLevels)
	for s = 1:nStim
		fprintf('%s\t\t', uniqueStim{s});

		% get spikes - need to account for NULL stimulus if present,
		% since there will be only 1 level for the NULL stim...
		if hasNULL && (s == nullIndex)
			% only 1 level for null stimulus
			spiket = spikesByStim{stimIndices{s}(1)};
			fprintf('%d (%d) dB SPL\n', dBLevelsByStim{s}(1), Dinf.test.Level(l));
		else
			% check on levels and stimIndices
			if nLevels(s) ~= length(stimIndices{s})
				error('%s: mismatch in nLevels and stimIndices', mfilename);
			end
			spiket = spikesByStim{stimIndices{s}(l)};
			fprintf('%d (%d) dB SPL\n', dBLevelsByStim{s}(l), Dinf.test.Level(l));
		end

		% now count spikes for each rep
		nReps = length(spiket);
		% allocate count storate
		spikeCount = zeros(nReps, 1);
		for r = 1:nReps
			spikeCount(r) = sum(between(spiket{r}, stimOnsetOffset(s, 1), ...
																stimOnsetOffset(s, 2) ) );
		end
		C{s, l} = spikeCount;
	end
end

if nargout
	varargout{1} = C;
	varargout{2} = dBLevelsByStim;
	varargout{3} = nLevels;
end

