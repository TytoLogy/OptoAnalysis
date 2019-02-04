% fratest.m
% testing things to be incorporated into getFilteredOptoData.m
% 4 Feb 2019

load fratestdata.mat

%----------------------------------------------------------------------
%% Some test-specific things...
%----------------------------------------------------------------------
% for FREQ test, find indices of stimuli with same frequency
if isnumeric(Dinf.test.Type)
	Dinf.test.Type = char(Dinf.test.Type);
end

switch upper(Dinf.test.Type)

	% for FRA (FREQ+LEVEL) test, find indices of stimuli with
	% freq and same level (dB SPL)
	case 'FREQ+LEVEL'
		% convert stimtype and curvetype to strings
		if isnumeric(Dinf.test.stimcache.stimtype)
			Dinf.test.stimcache.stimtype = char(Dinf.test.stimcache.stimtype);
		end
		if isnumeric(Dinf.test.stimcache.curvetype)
			Dinf.test.stimcache.curvetype = char(Dinf.test.stimcache.curvetype);
		end
		% if necessary, convert cells to matrices
		testcell = {'splval', 'rmsval', 'atten', 'FREQ', 'LEVEL'};
		for c = 1:length(testcell)
			if iscell(Dinf.test.stimcache.(testcell{c}))
				Dinf.test.stimcache.(testcell{c}) = ...
										cell2mat(Dinf.test.stimcache.(testcell{c}));
			end
		end
		% list of stimulus freqs, # of freqs tested
		freqlist = unique(Dinf.test.stimcache.FREQ, 'sorted');
		nfreqs = length(freqlist);
		% list of stimulus levels, # of levels tested
		levellist = unique(Dinf.test.stimcache.LEVEL, 'sorted');
		nlevels = length(levellist);
	otherwise
		error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end

%----------------------------------------------------------------------
%% Pull out trials, apply filter, store in matrix
%----------------------------------------------------------------------
%{
 Dinf.test.stimcache
               stimtype: 'tone'
              curvetype: 'FREQ+LEVEL'
             freezeStim: 0
                  nreps: 10
               saveStim: 0
                ntrials: 126
                 nstims: 1260
                 repnum: [1260×1 double]
               trialnum: [1260×1 double]
                 splval: [1260×1 double]
                 rmsval: [1260×1 double]
                  atten: [1260×1 double]
                   FREQ: [1260×1 double]
                  LEVEL: [1260×1 double]
                   opto: {1260×1 cell}
                radvary: 1
    trialRandomSequence: [10×126 double]
                  vname: [70 82 69 81 43 76 69 86 69 76]
                 vrange: [2×126 double]
                stimvar: {1×1260 cell}
%}

if strcmpi(Dinf.test.Type, 'FREQ+LEVEL')
	tracesByStim = cell(nfreqs, nlevels);
	
	for f = 1:nfreqs
		for l = 1:nlevels
			
			
		end
	end
	



end