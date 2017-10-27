function figH = grantplot(figH, t, traces, trace_ymax, spiketimes, histdata)

if ~isempty(figH)
	figure(figH);
else
	figH = figure;
end

tMinMax = [min(t) max(t)];

%---------------------------------------------------------------------
% Plot raw data
%---------------------------------------------------------------------
subplot(1, 3, 1)
% flip tracesByStim in order to have sweeps match the raster plots
stackplot(t, fliplr(traces), 'colormode', 'black', ...
										'ymax', trace_ymax, ...
										'figure', figH, 'axes', gca);
xlabel('ms')
ylabel('Trial')
% adjust y tick labels
ytl = get(gca, 'YTickLabels');
set(gca, 'YTickLabels', flipud(ytl));









%---------------------------------------------------------------------
% Plot raster data
%---------------------------------------------------------------------
subplot(1, 3, 2);
rasterplot(spiketimes, tMinMax, '.', 18, 'k');
xlabel('Time (ms)');
ylabel('Trial');
% adjust y tick labels
ytl = get(gca, 'YTickLabels');
set(gca, 'YTickLabels', flipud(ytl));

%---------------------------------------------------------------------
% Plot psth data
%---------------------------------------------------------------------
subplot(1, 3, 3)
bar(histdata.bins, histdata.H, 1, 'EdgeColor', 'b', 'FaceColor', 'b');
xlabel('Time (ms)');
ylabel('Spike Count');


