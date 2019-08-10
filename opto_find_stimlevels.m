function [dbLevelsByStim, nLevels] = opto_find_stimlevels(uniqueStim, stimIndices, Dinf)

% create cell array to store levels for each stim (might be different
% number of levels for each stim, so can't use array)
nStim = numel(uniqueStim);
dbLevelsByStim = cell(nStim, 1);
nLevels = zeros(nStim, 1);
% loop through stim
fprintf('Determining Stimulus Levels...\n');
for s = 1:nStim
	% ASSUME that repeats of same stimulus are at different levels...
	nLevels(s) = length(stimIndices{s});
	dbLevelsByStim{s} = zeros(nLevels(s), 1);
	fprintf('\tStimulus: %s\n', uniqueStim{s});
	fprintf('\t\tnLevels: %d\n', nLevels(s));
	fprintf('\t\tLevels: ');
	for l = 1:nLevels(s)
		dbLevelsByStim{s}(l) = Dinf.stimList(stimIndices{s}(l)).audio.Level;
		fprintf('%d  ', Dinf.stimList(stimIndices{s}(l)).audio.Level);
	end
	fprintf('\n');
end