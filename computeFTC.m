function varargout = computeFTC(spikeTimes, freqs, analysisWindow)
%------------------------------------------------------------------------
% process Spiketimes into frequency tuning curve
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% assumes spikeTimes is in format:
% 		spikeTimes{nFreqs, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%------------------------------------------------------------------------
% See Also: plotFTC, optoproc, opto
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 27 Mar 2019 (SJS) 
%	- adapted from computeRLF.m
% 
% Revisions:
%------------------------------------------------------------------------


nFreqs = length(freqs);
spikeCount = cell(nFreqs, 1);
ftc.freqs = freqs;
ftc.window = analysisWindow;
ftc.mean = zeros(nFreqs, 1);
ftc.std = zeros(nFreqs, 1);
ftc.mean_ci = cell(nFreqs, 1);
ftc.median = zeros(nFreqs, 1);
ftc.median_ci = cell(nFreqs, 1);

for n = 1:nFreqs
	nReps = length(spikeTimes{n});
	spikeCount{n} = zeros(nReps, 1);
	for r = 1:nReps
		spikeCount{n}(r) = sum(between(spikeTimes{n}{r}, analysisWindow(1), ...
																analysisWindow(2) ) );
	end
	ftc.mean(n) = mean(spikeCount{n}); 
	ftc.std(n) = std(spikeCount{n});
	ftc.mean_ci{n} = bootci(2000, @mean, spikeCount{n});
	ftc.median(n) = median(spikeCount{n});
	ftc.median_ci{n} = bootci(2000, @median, spikeCount{n});
end
ftc.spikeCount = spikeCount;

varargout{1} = ftc;

if nargout > 1
	varargout{2} = spikeCount;
end