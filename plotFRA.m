function plotFRA(FRA, varargin)
%------------------------------------------------
% plot as color patch
%------------------------------------------------

if nargin == 2
	if strcmpi(varargin{1}, 'ATTEN')
		yStr = 'Attenuation (dB)';
		yAtten = 1;
	elseif strcmpi(varargin{1}, 'dB')
		yStr = 'dB SPL';
		yAtten = 0;
	else
		error('invalid xStr');
	end
else
	yStr = 'Attenuation (dB)';	
	yAtten = 1;
end

% create figure and subplot
figure
subplot(211);
% plot as color patch.  need to flip sorted atten to get
% higher atten (lower amplitude tones) at bottom of plot and
% lower atten (higher amp) at top, per FRA plot convention
% the y axis labels will be in reverse order, but we'll take care
% of that later
xdata = log10(FRA.Freqs);
if yAtten
	ydata = fliplr(FRA.Levels);
else
	ydata = FRA.Levels;
end
pcolor(xdata, ydata, FRA.MeanCount);
% show color legend
colorbar
% deal with labels and title
xlabel('Log Frequency (kHz)');
ylabel(yStr);
title(	{	FRA.fname, ...
				sprintf('Avg Spike Count, [%d-%d] ms window', ...
											FRA.window(1), FRA.window(2)), ...
			}, ...
			'Interpreter', 'none');
% re-do X tick labels to that they're more readable
xt = get(gca, 'XTick');
xtl = cell(length(xt), 1);
get(gca, 'XTickLabel');
for n = 1:length(xt)
	xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
end
set(gca, 'XTickLabel', xtl);
% correct the Y tick labels, as promised
if yAtten
	set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
end
%------------------------------------------------
% create subplot and plot waterfall 
% (in fashion similar to color patch)
%------------------------------------------------
subplot(212);
waterfall(xdata, ydata, FRA.MeanCount);
% labels, again deal with log labels
xlabel('Log Frequency (kHz)');
ylabel(yStr);
zlabel('Spike Count');
xt = get(gca, 'XTick');
xtl = cell(length(xt), 1);
get(gca, 'XTickLabel');
for n = 1:length(xt)
	xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
end
set(gca, 'XTickLabel', xtl);
if yAtten
	set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')))
end
