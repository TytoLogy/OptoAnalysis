function lH = draw_stim_onoffset(stimulus_times)

default_color = 'b';

if ~isstruct(stimulus_times)
	error('%s: stimulus_times must be a struct!', mfilename);
elseif ~all(isfield(stimulus_times, {'onset', 'offset'}))
	error('%s: stimulus_times must have onset and offset fields', mfilename);
else
	nstim = length(stimulus_times);
	lH = cell(nstim, 1);
	% set color if needed
	if ~isfield(stimulus_times, 'color')
		for n = 1:nstim
			stimulus_times(n).color = default_color;
		end
	end
end

% plot stimulus onset/offset lines
for t = 1:nstim
	% get ylimits
	L = ylim;
	% onset line
	lH{t}(1) = line(stimulus_times(t).onset .* [1 1], L, ...
							'Color', stimulus_times(t).color);
	% offset line
	lH{t}(2) = 	line(stimulus_times(t).offset .* [1 1], L, ...
							'Color', stimulus_times(t).color);
end
