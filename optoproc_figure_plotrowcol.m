function nrowcols = optoproc_figure_plotrowcol(nStim)

% max number of plots (rows) per column
MAX_PLOTS_PER_COL = 6;

% for data that are not "2D" (e.g., FRA), adjust # of columns based
% on the number of variable levels or types (nvars)
if numel(nStim) ~= 1
	error('%s: not written to handle 2 dim nStim', mfilename);
elseif nStim == 0
	error('%s: nStim is 0!', mfilename);
end

% need to determine how many columns to use. 
% do this by dividing # of stim by max number of plots per col, and
% rounding up
ncols = ceil(nStim / MAX_PLOTS_PER_COL);

% then compute rows
if nStim <= MAX_PLOTS_PER_COL
	nrows = nStim;
else
	nrows = ceil(nStim/ncols);
end

nrowcols = [nrows ncols];