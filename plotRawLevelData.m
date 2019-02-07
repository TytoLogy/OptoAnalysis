function varargout = plotRawLevelData(traces, Dinf)


if ~strcmpi(Dinf.test.Type, 'LEVEL')
	warning('%s: test type (%s) is not LEVEL!', mfilename, Dinf.test.Type);
end

% get overall max value, so all traces can be plotted on same scale
levels = unique(Dinf.test.stimcache.LEVEL, 'sorted');
nlevels = length(levels);
if length(traces) ~= nlevels
	error('%s: length of traces (%d) ~= nlevels (%d)', length(traces), nlevels);
end

%% find overall max value
% initialize maxVals
maxVals = zeros(nlevels, 1);
for l = 1:nlevels
	maxVals(l) = max(max(traces{l}));
end
maxVal = max(maxVals);

% time vector for plotting - use length of 1st column of data
t = (1000/Dinf.indev.Fs)*((1:length(traces{1}(:, 1))) - 1);
if ~ispc
	[~, fname] = fileparts(replace(Dinf.filename, '/', filesep));
else
	[~, fname] = fileparts(Dinf.filename);
end

%% stacked data plot
sH = zeros(nlevels, 1);
for l = 1:nlevels
	[~, sH(l)] = stackplot(t, traces{l}, 'colormode', 'black', 'Ymax', maxVal);
	title({fname, sprintf('Level %d', levels(l))}, 'Interpreter', 'none');
end

%% overlay plot
oH = zeros(nlevels, 1);
for l = 1:nlevels
	[~, oH(l)] = overlayplot(t, traces{l}, 'Ymax', maxVal);
	title({fname, sprintf('Level %d', levels(l))}, 'Interpreter', 'none');
end

%% outputs
if nargout > 0
	if nargout >= 1
		varargout{1} = sH;
	end
	if nargout == 2
		varargout{2} = oH;
	end
end
