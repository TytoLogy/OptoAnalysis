function [axH, figH] = grantplot(t, traces, trace_ymax, spiketimes, ...
										histdata, stimulus_times, ...
										psth_yLim, subOnOff)
%---------------------------------------------------------------------
% Plot things for grant
%---------------------------------------------------------------------

tMinMax = [min(t) max(t)];

%---------------------------------------------------------------------
% Plot raw data
%---------------------------------------------------------------------
if subOnOff
	subplot(1, 3, 1);
	figH = gcf;
else
	figH(1) = figure;
 	axes;
end
% flip tracesByStim in order to have sweeps match the raster plots
[~, ~, axH(1)] = stackplot(t, fliplr(traces), 'colormode', 'black', ...
										'ymax', trace_ymax, ...
										'figure', gcf, 'axes', gca);
xlabel('ms')
ylabel('Trial')
% adjust y tick labels
ytl = get(gca, 'YTickLabels');
set(gca, 'YTickLabels', flipud(ytl));
% draw stimulus on/off
draw_stim_onoffset(stimulus_times);
config_plots(axH(1));

%---------------------------------------------------------------------
% Plot raster data
%---------------------------------------------------------------------
if subOnOff
	subplot(1, 3, 2);
	figH = gcf;
else
	figH(2) = figure;
	axes;
end
axH(2) = rasterplot(spiketimes, tMinMax, '.', 20, 'k');
xlabel('Time (ms)');
ylabel('Trial');
% adjust y tick labels
ytl = get(gca, 'YTickLabels');
set(gca, 'YTickLabels', flipud(ytl));
% draw stimulus on/off
draw_stim_onoffset(stimulus_times);
config_plots(axH(2));

%---------------------------------------------------------------------
% Plot psth data
%---------------------------------------------------------------------
if subOnOff
	subplot(1, 3, 3);
	figH = gcf;
else
	figH(3) = figure;
	axes;
end
bar(histdata.bins, histdata.H, 1, 'EdgeColor', 'k', 'FaceColor', 'k');
axH(3) = gca;
xlabel('Time (ms)');
ylabel('Spike Count');
% psth limits
if ~isempty(psth_yLim)
	ylim(psth_yLim);
end
% draw stimulus on/off
draw_stim_onoffset(stimulus_times);
config_plots(axH(3));

end

%---------------------------------------------------------------------
% config plots
%---------------------------------------------------------------------
function config_plots(h)
	set(h, 'TickDir', 'out')
	set(h, 'Box', 'off');
end

