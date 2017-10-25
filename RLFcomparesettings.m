%---------------------------------------------------------------------
%---------------------------------------------------------------------
% compare settings
%---------------------------------------------------------------------
%---------------------------------------------------------------------


%---------------------------------------------------------------------
% data sources
%---------------------------------------------------------------------

%{
% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171004';
% files to compare
file1 = '1155_20171004_01_01_2563_BBN_LEVEL.dat';
file2 = '1155_20171004_01_01_2563_BBN_LEVEL_optoON.dat';
% Channel Number (use 8 for single channel data)
channelNumber = 8;
%}

% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171006';
% files to compare
% file1 = '1155_20171006_03_03_2988_BBN_optoOFF_ch5.dat';
% file2 = '1155_20171006_03_03_2988_BBN_optoON_ch5.dat';
file1 = '1155_20171006_04_03_3123_BBN_optoOFF_ch5ch11.dat';
file2 = '1155_20171006_04_03_3123_BBN_optoON_ch5ch11_9.dat';

% Channel Number (use 8 for single channel data)
channelNumber = 5;

%---------------------------------------------------------------------
% settings for processing data
%---------------------------------------------------------------------

% time window (post-stimulus onset) for analysis : this is added to stim
% delay
WindowLen = [0 75];

% filter
HPFreq = 350;
LPFreq = 6500;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
% binSize for PSTH (milliseconds)
binSize = 5;
% SAVE PLOTS?
saveFIG = 0;
savePNG = 0;
savePDF = 1;