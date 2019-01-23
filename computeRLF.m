function varargout = computeRLF(spikeTimes, levels, analysisWindow)
%------------------------------------------------------------------------
% process Spiketimes into rate level function
%------------------------------------------------------------------------
% assumes spikeTimes is in format:
% 		spikeTimes{nLevels, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 23 Oct 2017 (SJS) 
%	- adapted from viewOptoData.m
% 
% Revisions:
%	22 Jan 2019 (SJS): added levels as input
%------------------------------------------------------------------------


nLevels = length(levels);
spikeCount = cell(nLevels, 1);
rlf.levels = levels;
rlf.window = analysisWindow;
rlf.mean = zeros(nLevels, 1);
rlf.std = zeros(nLevels, 1);
rlf.mean_ci = cell(nLevels, 1);
rlf.median = zeros(nLevels, 1);
rlf.median_ci = cell(nLevels, 1);

for n = 1:nLevels
	nReps = length(spikeTimes{n});
	spikeCount{n} = zeros(nReps, 1);
	for r = 1:nReps
		spikeCount{n}(r) = sum(between(spikeTimes{n}{r}, analysisWindow(1), ...
																analysisWindow(2) ) );
	end
	rlf.mean(n) = mean(spikeCount{n}); 
	rlf.std(n) = std(spikeCount{n});
	rlf.mean_ci{n} = bootci(2000, @mean, spikeCount{n});
	rlf.median(n) = median(spikeCount{n});
	rlf.median_ci{n} = bootci(2000, @median, spikeCount{n});
end
rlf.spikeCount = spikeCount;

varargout{1} = rlf;

if nargout > 1
	varargout{2} = spikeCount;
end