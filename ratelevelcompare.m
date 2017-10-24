%---------------------------------------------------------------------
%---------------------------------------------------------------------
RLFcomparesettings

%---------------------------------------------------------------------
%% set paths to things:
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths;
% output path for plots
plotpath_base = fullfile(data_root_path, 'Analyzed');

%---------------------------------------------------------------------
% read in data files, plot data
%---------------------------------------------------------------------
% build datapath
datapath = fullfile(data_root_path, animalID, dateID);

[~, Dinf1, S1] = optoproc('file', fullfile(datapath, file1), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'plotPSTH', 'savepdf');
							
[~, Dinf2, S2] = optoproc('file', fullfile(datapath, file2), ...
								'plotpath', plotpath_base, ...
								'channel', channelNumber, ...
								'binsize', binSize, ...
								'plotPSTH', 'savepdf');
							
%---------------------------------------------------------------------
%% construct rate-level function
%---------------------------------------------------------------------
AnalysisWindow = Dinf1.audio.Delay + WindowLen;

rlf1 = computeRLF(S1.spiketimes, AnalysisWindow);
rlf2 = computeRLF(S2.spiketimes, AnalysisWindow);

figure;
errorbar(Dinf1.audio.Level, rlf1.mean, rlf1.std);
hold on
	errorbar(Dinf2.audio.Level, rlf2.mean, rlf2.std);
hold off

xlabel('dB SPL');
ylabel('mean spikes/trial +/- s.d.')
title({file1, file2}, 'Interpreter', 'none')
legend({'ctrl', 'optoOn'})




