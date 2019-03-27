function varargout = plotPSTHMATRIX(spiketimes, Dinf, binSize, nvars, varlist, ...
												nrowcols, timeLimits, yLimits, titleString)
%------------------------------------------------------------------------
%  H = plotPSTHMATRIX(spiketimes, Dinf, nvars, varlist, timeLimits, ...
%												nrowcols, yLimits, titleString)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
%  Input Args:
%	 
%
%  Output Args:
%	 H		handle to figure
%------------------------------------------------------------------------
% See Also: computeFRA, optoproc
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 29 March 2019 (SJS), pulled code from optoproc
%
% Revisions:
%------------------------------------------------------------------------

% create figure
hPR = figure;

prows = nrowcols(1);
pcols = nrowcols(2);

plotopts.timelimits = timeLimits;
% assign y axis limits if provided
if ~isempty(yLimits)
	plotopts.psth_ylimits = yLimits;
end
% options for raster and psth matrix
plotopts.raster_tickmarker = '.';
plotopts.raster_ticksize = 16;
plotopts.raster_color = [0 0 0];
plotopts.psth_binwidth = binSize;
plotopts.plotgap = 0.001;
plotopts.xlabel = 'msec';
plotopts.stimulus_times_plot = 3;
plotopts.stimulus_on_color{1} = [0 0 1];
plotopts.stimulus_off_color{1} = [0 0 1];
plotopts.stimulus_onoff_pct(1) = 60;
% add on off bars for opto stim
if Dinf.opto.Enable
	% add colors for second stimulus
	plotopts.stimulus_on_color{2} = [1 0 0];
	plotopts.stimulus_off_color{2} = [1 0 0];
	plotopts.stimulus_onoff_pct(2) = 80;
end

% create times to indicate stimuli
if numel(nvars) == 1
	stimulus_times = cell(nvars, 1);
	for v = 1:nvars
		% need to have [stim_onset stim_offset], so add delay to 
		% [0 duration] to compute proper times. then, multiply by 0.001 to
		% give times in seconds (Dinf values are in milliseconds)
		stimulus_times{v, 1} = 0.001 * (Dinf.audio.Delay + ...
															[0 Dinf.audio.Duration]);
		% if opto is Enabled, add it to the array by concatenation
		if Dinf.opto.Enable
			stimulus_times{v, 1} = [stimulus_times{v, 1}; ...
												 0.001 * (Dinf.opto.Delay + ...
															[0 Dinf.opto.Dur]) ];
		end
	end

	% adjust depending on # of columns of plots
	if nvars <= 5
		plotopts.plot_titles = titleString;
		plotopts.stimulus_times = stimulus_times;
		rasterpsthmatrix(spiketimes, plotopts);
	elseif iseven(nvars)
		plotopts.plot_titles = reshape(titleString, [prows pcols]);
		plotopts.stimulus_times = reshape(stimulus_times, [prows pcols]);
		rasterpsthmatrix(reshape(spiketimes, [prows pcols]), plotopts);
	else
		% need to add 'dummy' element to arrays
	% 	titleString = [titleString; {''}];
	% 	spiketimes = [spiketimes; {{}}];
		plotopts.plot_titles = reshape([titleString; {''}], [prows pcols]);
		plotopts.stimulus_times = ...
							reshape(	[stimulus_times; stimulus_times{end}], ...
										[prows pcols]);
		rasterpsthmatrix(reshape([spiketimes; {{}}], [prows pcols]), plotopts);
	end
else
	% for 2D data (FRA)... might be overkill....
	stimulus_times = cell(nvars(2), nvars(1));
	plotopts.plot_titles = cell(nvars(2), nvars(1));
	for f = 1:nvars(1)
		for l = 1:nvars(2)
			% need to have [stim_onset stim_offset], so add delay to 
			% [0 duration] to compute proper times. then, multiply by 0.001 to
			% give times in seconds (Dinf values are in milliseconds)
			stimulus_times{l, f} = 0.001 * (Dinf.audio.Delay + ...
																[0 Dinf.audio.Duration]);
			% if opto is Enabled, add it to the array by concatenation
			if Dinf.opto.Enable
				stimulus_times{l, f} = [stimulus_times{l, f}; ...
													 0.001 * (Dinf.opto.Delay + ...
																[0 Dinf.opto.Dur]) ];
			end
			plotopts.plot_titles{l, f} = sprintf('%d kHz', varlist{1}(f));
		end
	end

	plotopts.stimulus_times = stimulus_times;
	rasterpsthmatrix(spiketimes, plotopts);

end

if nargout
	varargout{1} = hPR;
end
