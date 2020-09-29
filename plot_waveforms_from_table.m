function hWF = plot_waveforms_from_table(st, Fs, varargin)
%------------------------------------------------------------------------
% h = plot_waveforms_from_table(spiketable)
%------------------------------------------------------------------------
% TytoLogy:OptoAnalysis
%--------------------------------------------------------------------------
% given input table of spikes, plots waveforms
%------------------------------------------------------------------------
% See Also: optosort, optoproc, opto (TytoLogy:opto program)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 23 September 2020 (SJS)
%	 Uses code from optosort:plot_curve_demo.m
% Revisions:
%
%------------------------------------------------------------------------

% (1) vertically concatenate data stored in cell array of sweeps in
% st.spiketable into a single table
tmpT = vertcat(st.spiketable{:});

% (2) extract just the wave field 
tmpwav = tmpT.Wave;

% (3) plot overlaid waveforms

% plot in new figure
figure;

% need time vector for x-axis
[~, nBins] = size(tmpwav);
t_ms = (1000/Fs) * (0:(nBins - 1));

% plot waveforms, with mean and legend
% need tmpwav to be in column format - time in rows, indv. units by column
% so send the function the transpose of the tmpwav matrix.
hWF = plot_spike_waveforms(t_ms, tmpwav', 'MEAN', true, 'LEGEND', true);

% add title to plot
% create title string with 2 rows:
%	filename (either from st struct or S.Info.FileInfo{findx}.F.file
%	channel and unit
tstr = {	st.fileName, ...
			sprintf('Channel %d Unit %d', st.channel, st.unit)};
title(tstr, 'Interpreter', 'none');	

% set figure filename - use the base from the FreqTuningInfo.F object
% set(gcf, 'Name', sprintf('%s_Ch%d_Un%d', ...
%					S.Info.FileInfo{findx}.F.base, st.channel, st.unit));