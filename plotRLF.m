function varargout = plotRLF(rlfStruct, varargin)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%------------------------------------------------------------------------
% plot rate level function
%------------------------------------------------------------------------
% assumes rlfStruct is in format:
%		rlfStruct.fname		filename
% 		rlfStruct.spikeCount	cell array of spike counts/trial at each level
% 		rlfStruct.levels		stimulus levels
%		rlfStruct.window		time window [tstart tend] in ms used for analysis
% 		rlfStruct.mean			mean values at each level
% 		rlfStruct.std			std. dev. at each level
% 		rlfStruct.mean_ci		cell array of 95% conf. intervals for mean
% 		rlfStruct.median		median spike counts at each level
% 		rlfStruct.median_ci	cell array of 95% conf interval for median
%------------------------------------------------------------------------
% See Also: computeRLF, optoproc, opto
%------------------------------------------------------------------------
 		
%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 23 Oct 2017 (SJS) 
%	- adapted from viewOptoData.m
% Revisions:
%	27 Mar 2019 (SJS): changes for inclusion into optoproc
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
cimatrix = zeros(length(rlfStruct.levels), 2);
for l = 1:length(rlfStruct.levels)
	if strcmpi(dataToPlot, 'MEAN')
		cimatrix(l, :) = rlfStruct.mean_ci{l}';
	else
		cimatrix(l, :) = rlfStruct.median_ci{l}';
	end
end

% create figure
H = figure;

% plot mean
if strcmpi(dataToPlot, 'MEAN')
	ploterrea(rlfStruct.levels, rlfStruct.mean, cimatrix);
	ylabel('Mean Spike Count');
% plot median
else
	ploterrea(rlfStruct.levels, rlfStruct.median, cimatrix);
	ylabel('Median Spike Count');	
end
% x label
xlabel('dB SPL')
grid on
	
if nargout
	varargout{1} = H;
end

