function [outdata, outstamps] = exportChannelAsVector(D, Dinf, channel)
%------------------------------------------------------------------------
% [outdata, outstamps] = exportChannelAsVector(D, Dinf, channel)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% function to combine indvidual trials of neural data into a single stream
%--------------------------------------------------------------------------
% 	outdata will be a cell array of size {# channels, 1} with each row 
% 		containing a vector of all sweeps for that channel
% 	outstamps will be cell array of size {# channels, 2}
% 			column 1: traceIndices(# sweeps, 2) holding start index (sample #)
% 							and end index for each sweep within the long combined
% 							channel data vector. This information will be used
% 							for 
% 			column 2: # of samples for each sweep (more for diagnostic use)
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: September, 2019 (SJS)
%
% Revisions:
%	26 Sep 2019 (SJS): renamed, added comments
%--------------------------------------------------------------------------

% number of traces/stimulus presentations
nstims = length(D);
% number of channels to loop through
nchan = length(channel);
% allocate outdata and outstamps
outdata = cell(nchan, 1);
outstamps = cell(nchan, 2);

% loop through channels
for cIndx = 1:nchan

	% allocate temporary data storage for current channel
	tmpD = cell(nstims, 1);
	% allocate array to hold sweep start and end
	traceIndices = zeros(nstims, 2);
	% allocate array to hold # of samples in trace
	nsamples = zeros(nstims, 1);
	for s = 1:nstims
		% store start index for trace
		if s == 1
			% if this is firs
			traceIndices(s, 1) = 1;
		else
			traceIndices(s, 1) = traceIndices(s-1, 2) + 1;
		end
		% store trace in cell array
		tmpD{s} = D{s}.datatrace(:, channel(cIndx));
		% store # of samples for this trace
		nsamples(s) = length(D{s}.datatrace(:, channel(cIndx)));
		% store end index for trace - this will be the start index + 
		% length of trace (nsamples) and need to subtract 1 
		traceIndices(s, 2) = traceIndices(s, 1) + nsamples(s) - 1;
	end

	% convert to single-row vector and store in outdata
	outdata{cIndx} = cell2mat(reshape(tmpD, [], 1));
	% store stop and start indices
	outstamps{cIndx, 1} = traceIndices;
	% store # of samples per trace
	outstamps{cIndx, 2} = nsamples;
end