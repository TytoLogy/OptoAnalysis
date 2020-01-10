function F = parse_opto_filename(filename)
%------------------------------------------------------------------------
% parse_opto_filename
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% extracts information from opto filename
%
%------------------------------------------------------------------------
% Input Arguments:
% 	filename		data file from opto program
%					e.g., 1372_20191126_03_01_1500_FREQ_TUNING.dat
% 
% Output Arguments:
% 	F	struct with fields:
% 		animal: '1372'
% 		datecode: '20191126'
% 		unit: '03'
% 		penetration: '01'
% 		depth: '1500'
% 		other: 'FREQ_TUNING'
%------------------------------------------------------------------------
% See Also: opto_createDataFileName
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 October, 2017
%
% Revisions:
%	9 Aug 2019 (SJS): 
%		penetration and unit are reversed: fixed
%		added rec depth
%	10 Jan 2019 (SJS): 
%	 - fixed issue with final startusc and test name
%	 - added base to F struct to hold filename
%	 - added input/output arg info to help
%------------------------------------------------------------------------

if isempty(filename)
	F = [];
	return
end

%---------------------------------------------------------------------
% get info from filename
%---------------------------------------------------------------------
% only need base filename (no path or ext)
[~, fname] = fileparts(filename);
F.base = fname;
% locate underscores
usc = find(fname == '_');
% last underscore index
endusc = usc - 1;
% first underscore index
startusc = usc + 1;
% animal #
F.animal = fname(1:endusc(1));
% date (<year><month><date> => YYYYMMDD, e.g., 20170401 for April 1, 2017)
F.datecode = fname(startusc(1):endusc(2));
% unit number
F.unit = fname(startusc(2):endusc(3));
% penetration number
F.penetration = fname(startusc(3):endusc(4));
% recording depth
F.depth = fname(startusc(4):endusc(5));
% this should be test name
if startusc(4) == startusc(end)
	F.other = fname(startusc(end):end);
else
	F.other = fname(startusc(5):end);
end
