function F = parse_opto_filename(filename)
%------------------------------------------------------------------------
% parse_opto_filename
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% extracts information from opto filename
%
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 October, 2017
%
% Revisions:
%	see git!
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
% penetration number
F.penetration = fname(startusc(2):endusc(3));
% unit number
F.unit = fname(startusc(3):endusc(4));
F.other = fname(startusc(end):end);
