function FRA = computeFRA(spiketimes, freqs, levels, frawin)
%------------------------------------------------------------------------
% FRA = computeFRA(spikes, freqs, levels, frawin)
%------------------------------------------------------------------------
% TytoLogy:OptoAnalysis:computeFRA
%------------------------------------------------------------------------
% 
% Using spiketimes, list of stimulus frequencies and levels and
% an analysis window, calculate mean firing rates for each
% stimulus freq and level combination and create FRA response structure
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	spiketimes  {nlevels, nfreqs} cell array of spike times 
%                 each element is an {ntrials, 1} cell array of
%                 spiketime vectors
%  freqs       [nfreqs, 1] vector of stimulus frequency values (Hz)
%  levels      [levels, 1] vector of stimulus level values (dB SPL)
%  frawin      [start, end] start and end time of analysis window, msec
% 
% Output Arguments:
% 	FRA   FRA data struct
%   SpikeTimes  spiketimes that are inside the analysis 
%                 window defined in frawin
%   window      frawin
%   SpikeCount  {nlevels, nfreqs} cell array with vectors of # of spikes for 
%               each stimulus presentation
%   MeanCount   [nlevels, nfreqs] matrix of mean spike count
%   StdDevCount [nlevels, nfreqs] matrix of spike count std. deviation
%   Freqs       list of frequencies
%   Levels      list of levels
%   nfreqs      # frequencies
%   nlevels     # of levels
%   fname       source data file name
%
%------------------------------------------------------------------------
% See also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ? (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------
% get dimensions of data
%------------------------------------------------
nlevels = length(levels);
nfreqs = length(freqs);
[nl, nf] = size(spiketimes);
if (nl ~= nlevels) || (nf ~= nfreqs)
	error('mismatch in levels and/or freqs')
end

%------------------------------------------------
% get the spikes that are inside the analysis window, 
% store in SpikeTimes
%------------------------------------------------
FRA.SpikeTimes = cell(nlevels, nfreqs);
FRA.window = frawin;
for f = 1:nfreqs
	for l = 1:nlevels
		nreps = length(spiketimes{l, f});
		FRA.SpikeTimes{l, f} = cell(nreps, 1);
		for r = 1:nreps
			FRA.SpikeTimes{l, f}{r} = find_valid_timestamps( ...
																spiketimes{l, f}{r}, ...
																frawin(1), ...
																frawin(2), 0);
		end
	end
end

%------------------------------------------------
% compute mean and std dev spike count
%------------------------------------------------
FRA.SpikeCount = cell(nlevels, nfreqs);
FRA.MeanCount = zeros(nlevels, nfreqs);
FRA.StdDevCount = zeros(nlevels, nfreqs);

% loop through frequencies (sorted)
for f = 1:nfreqs
	% loop through levels from low to high
	for l = 1:nlevels
		% add up number of spikes for this freq and atten combination
		% and store each rep in obj.SpikeCount cell matrix
		nreps = length(FRA.SpikeTimes{l, f});
		FRA.SpikeCount{l, f} = zeros(nreps, 1);
		for r = 1:nreps
			FRA.SpikeCount{l, f}(r) = length(FRA.SpikeTimes{l, f}{r});
		end	% END r LOOP
		% compute mean and std dev spike count
		FRA.MeanCount(l, f) = mean(FRA.SpikeCount{l, f});
		FRA.StdDevCount(l, f) = std(FRA.SpikeCount{l, f});
	end	% END l LOOP
end	% END f LOOP

FRA.Freqs = freqs;
FRA.Levels = levels;
FRA.nfreqs = nfreqs;
FRA.nlevels = nlevels;
FRA.fname = '';

