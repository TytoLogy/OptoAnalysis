function H = draw_stim_onoffset(stimulus_times, varargin)

if numel(varargin)
	if ishghandle(varargin{1})
		axH = varargin{1};
	else
		axH = gca;
	end
else
	axH = gca;
end

if ~
% plot stimulus onset/offset lines if stimulus_times provided
if isfield(plotopts, 'stimulus_times')
	if any(plotopts.stimulus_times_plot == [2 3])
		if iscell(plotopts.stimulus_times)
			[nstim, ~] = size(plotopts.stimulus_times{row, col});
			if (nstim > 1) && ...
						( length(plotopts.stimulus_onoff_pct) == 1 )
				plotopts.stimulus_onoff_pct = ...
							plotopts.stimulus_onoff_pct * ones(nstim, 1);
			end
			for t = 1:nstim
				% get ylimits
				L = ylim;
				L(1) = 0.01*plotopts.stimulus_onoff_pct(t)*L(2);
				onset = 1000*plotopts.stimulus_times{row, col}(t, 1);
				offset = 1000*plotopts.stimulus_times{row, col}(t, 2);
				% onset line
				line(onset.*[1 1], L, 'Color', ...
											plotopts.stimulus_on_color{t});
				% offset line
				line(offset.*[1 1], L, 'Color', ...
											plotopts.stimulus_off_color{t});	
			end
		else
			% get ylimits
			L = ylim;
			L(1) = 0.01*plotopts.stimulus_onoff_pct*L(2);
			if any(length(plotopts.stimulus_times) == [1 2])
				%draw onset line
				onset = 1000*plotopts.stimulus_times(1);
				line(onset.*[1 1], L, 'Color', ...
										plotopts.stimulus_on_color);
			end
			if length(plotopts.stimulus_times) == 2
				% draw offset line
				offset = 1000*plotopts.stimulus_times(2);
				line(offset.*[1 1], L, 'Color', ...
										plotopts.stimulus_off_color);
			end
		end
	end
end	% END OF if isfield(plotopts, 'stimulus_times')
