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
% Output:
%  curveStruct			struct with fields:
%		curveStruct.fname			filename
% 		curveStruct.spikeCount	cell array of spike counts/trial at each x
% 		curveStruct.xdata			stimulus freqs
%		curveStruct.xlabel		label for x axis
%		curveStruct.window		time window [tstart tend] in ms used for analysis
% 		curveStruct.mean			mean values at each freq
% 		curveStruct.std			std. dev. at each freq
% 		curveStruct.mean_ci		cell array of 95% conf. intervals for mean
% 		curveStruct.median		median spike counts at each freq
% 		curveStruct.median_ci	cell array of 95% conf interval for median
%------------------------------------------------------------------------
% See Also: plotCurveAndCI, optoproc, opto
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 27 Mar 2019 (SJS) 
%	- adapted from computeRLF.m
% 
% Revisions:
%	28 Mar 2019 (SJS): updated to output curveStruct for use with
%								plotCurveAndCI()
%------------------------------------------------------------------------


nFreqs = length(freqs);
spikeCount = cell(nFreqs, 1);
% convert frequency from Hz to kHz
curveStruct.xdata = 0.001*freqs;
curveStruct.xlabel = 'Frequency (kHz)';
curveStruct.window = analysisWindow;
curveStruct.mean = zeros(nFreqs, 1);
curveStruct.std = zeros(nFreqs, 1);
curveStruct.mean_ci = cell(nFreqs, 1);
curveStruct.median = zeros(nFreqs, 1);
curveStruct.median_ci = cell(nFreqs, 1);

for n = 1:nFreqs
	nReps = length(spikeTimes{n});
	spikeCount{n} = zeros(nReps, 1);
	for r = 1:nReps
		spikeCount{n}(r) = sum(between(spikeTimes{n}{r}, analysisWindow(1), ...
																analysisWindow(2) ) );
	end
	curveStruct.mean(n) = mean(spikeCount{n}); 
	curveStruct.std(n) = std(spikeCount{n});
	curveStruct.mean_ci{n} = bootci(2000, @mean, spikeCount{n});
	curveStruct.median(n) = median(spikeCount{n});
	curveStruct.median_ci{n} = bootci(2000, @median, spikeCount{n});
end
curveStruct.spikeCount = spikeCount;

varargout{1} = curveStruct;

if nargout > 1
	varargout{2} = spikeCount;
end