function finfo = parse_dat_filename(datafile)
%---------------------------------------------------------------------
% get info from filename - this makes some assumptions about file
% name structure!
% <animal id #>_<date>_<penetration #>_<unit #>_<other info>.dat
%---------------------------------------------------------------------
% break up file name into <fname>.<ext> (~ means don't save ext info)
[~, fname] = fileparts(datafile);
% locate underscores in fname
usc = find(fname == '_');
% location of start and end underscore indices
%    abcde_edcba
%        ^ ^
%        | |
%        | ---- endusc index
%        ---startusc index
endusc = usc - 1;
startusc = usc + 1;
finfo.animal = fname(1:endusc(1));
finfo.datecode = fname(startusc(1):endusc(2));
finfo.penetration = fname(startusc(2):endusc(3));
finfo.unit = fname(startusc(3):endusc(4)); 
finfo.other = fname(startusc(end):end); 
finfo.fname = fname;