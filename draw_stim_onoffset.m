function lH = draw_stim_onoffset(stimulus_times)
%------------------------------------------------------------------------
% lH = draw_stim_onoffset(stimulus_times)
%------------------------------------------------------------------------
% OptoAnalysis
%------------------------------------------------------------------------
% given stimulus_times input struct, adds lines to plot that
% indicate stimulus onset and offset
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Created: ??? (SJS)
% 
% Revisions:
%  28 Nov 2023 (SJS): added some documentation
%------------------------------------------------------------------------

% default line color
default_color = 'b';

% check input data
if ~isstruct(stimulus_times)
	error('%s: stimulus_times must be a struct!', mfilename);
elseif ~all(isfield(stimulus_times, {'onset', 'offset'}))
	error('%s: stimulus_times must have onset and offset fields', mfilename);
else
   % get some information about inputs
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
