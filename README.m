Data from MTwav:

PSTH + Raster plots for each stimulus (all levels):

	optoproc('PLOT_PSTH')

PSTH + Raster plots for all stimuli at one level per page:
	optoproc('PLOT_PSTH_BY_LEVEL');
	

Other data:

for BBN, Freq tuning:
	
optoproc('PLOT_PSTH_MATRIX')
	
specify bin size for psth:

	for 5 ms bins:
	optoproc('PLOT_PSTH', 'binsize', 5)


tip: you can specify file to save time:

	optoproc('file', <<file with path>>, ....other options....)
	e.g. optoproc('file', 'D:\dirname\anotherdir\data\9999\9999_20210701_234_WAV.dat', ...)
	
