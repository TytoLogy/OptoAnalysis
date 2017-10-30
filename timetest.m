%---------------------------------------------------------------------
% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
outpath_base = fullfile(data_root_path, 'Analyzed');
%---------------------------------------------------------------------
% settings and files
%---------------------------------------------------------------------
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171025';
% unit #
unit = '02';
% penetration
penetration = '02';
% recording depth
depth = '2397';
% channel
channelNumber = 8;
% data file
datafile = '1155_20171025_02_02_2397_OptoInhibOFF.dat';
% data path
datapath = fullfile(data_root_path, animalID, dateID);
%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
[~, Dinf] = getFilteredOptoData( ...
											fullfile(datapath, datafile), ...
											'Channel', channelNumber);
%---------------------------------------------------------------------
% Display time info, show how to compute elapsed time
%---------------------------------------------------------------------
fprintf('Start time: %s\n', datestr(Dinf.time_start));
fprintf('End time: %s\n', datestr(Dinf.time_end));
% elapsed time is calculated by simple subtraction of the time number
% values
fprintf('Elapsed time: %s\n',  ...
							datestr(Dinf.time_end-Dinf.time_start, 'MM:SS.FFF'));
