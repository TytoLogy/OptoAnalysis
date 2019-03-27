function varargout = plotFTC(ftcStruct, varargin)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% plot frequency tuning curve
%------------------------------------------------------------------------
% assumes ftcStruct is in format:
%		ftcStruct.fname		filename
% 		ftcStruct.spikeCount	cell array of spike counts/trial at each freq
% 		ftcStruct.freqs		stimulus freqs
%		ftcStruct.window		time window [tstart tend] in ms used for analysis
% 		ftcStruct.mean			mean values at each freq
% 		ftcStruct.std			std. dev. at each freq
% 		ftcStruct.mean_ci		cell array of 95% conf. intervals for mean
% 		ftcStruct.median		median spike counts at each freq
% 		ftcStruct.median_ci	cell array of 95% conf interval for median
%------------------------------------------------------------------------
% See Also: computeFTC, optoproc, opto
%------------------------------------------------------------------------
 		
%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 27 Mar 2019  (SJS) 
%	- adapted from plotRLF.m
% Revisions:
%------------------------------------------------------------------------

%---------------------------------------------------------------------
% Defaults
%---------------------------------------------------------------------

dataToPlot = 'MEAN';

%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
if nargin > 1
	argIndx = 1;
	while argIndx <= (nargin - 1)
		switch upper(varargin{argIndx})
			case {'MEAN', 'AVERAGE', 'AVG'}
				dataToPlot = 'MEAN';
				argIndx = argIndx + 1;
			case {'MEDIAN', 'MED'}
				dataToPlot = 'MEDIAN';
				argIndx = argIndx + 1;
			otherwise
				error('%s: unknown option %s', mfilename, varargin{argIndx});
		end
	end
end

% compute confidence intervals
cimatrix = zeros(length(ftcStruct.freqs), 2);
for l = 1:length(ftcStruct.freqs)
	if strcmpi(dataToPlot, 'MEAN')
		cimatrix(l, :) = ftcStruct.mean_ci{l}';
	else
		cimatrix(l, :) = ftcStruct.median_ci{l}';
	end
end

% create figure
H = figure;

% plot mean
if strcmpi(dataToPlot, 'MEAN')
	ploterrea(ftcStruct.freqs, ftcStruct.mean, cimatrix);
	ylabel('Mean Spike Count');
% plot median
else
	ploterrea(ftcStruct.freqs, ftcStruct.median, cimatrix);
	ylabel('Median Spike Count');	
end
% x label
xlabel('Frequency')
grid on
	
if nargout
	varargout{1} = H;
end

