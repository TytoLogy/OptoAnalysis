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
% for all files:
% 	acq duration = 600 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:60dB, 100ms delay, 150ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, varied command voltage to ThorLabs LED
optoRLF.files = {	...
	'1155_20171025_02_02_2397_BBN.dat', ...					% rlf, 0:10:60
	'1155_20171025_02_02_2397_BBN_optoON.dat', ...			% rlf, 0:10:60, opto 3500 mV
	'1155_20171025_02_02_2397_BBN_optoON_2.dat', ...		% rlf, 0:10:60, opto 500 mV
	'1155_20171025_02_02_2397_BBN_optoON_3.dat', ...		% rlf, 0:10:60, opto 250 mV
	'1155_20171025_02_02_2397_BBN_optoON_4.dat', ...		% rlf, 0:10:60, opto 750 mV
	'1155_20171025_02_02_2397_BBN_optoON_5.dat', ...		% rlf, 0:10:60, opto 1000 mV
	'1155_20171025_02_02_2397_BBN_optoON_6.dat', ...		% rlf, 0:10:60, opto 1500 mV
	'1155_20171025_02_02_2397_BBN_optoON_7.dat', ...		% rlf, 0:10:60, opto 2000 mV
	'1155_20171025_02_02_2397_BBN_optoON_8.dat', ...		% rlf, 0:10:60, opto 2500 mV
	'1155_20171025_02_02_2397_BBN_optoON_9.dat' ...			% rlf, 0:10:60, opto 3000 mV
};
optoRLF.LEDintensity = [ ...
		0, ...
		3500, ...
		500, ...
		250, ...
		750, ...
		1000, ...
		1500, ...
		2000, ...
		2500, ...
		3000 ...
];
optoRLF.LEDpower = [ ...
	0.0, ...
	6.8, ...
	1.4, ... 
	0.7, ... 
	2.0, ...
	2.5, ...
	3.6, ...
	4.5, ...
	5.3, ...
	6.1 ...
];

bbnRLF.files = '1155_20171025_02_02_2397_BBN_3.dat';	% rlf, 0:10:60

bbn.files = '1155_20171025_02_02_2397_BBN_4.dat';	% BBN, 60 dB, 15 reps

% sliding "pulse" of opto stimulation
% for all files:
% 	acq duration = 600 ms, ISI = 250 ms
% 	15 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], [0 30 60] dB, 150ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		50 ms duration, delay varied, 2000 mV command voltage to ThorLabs LED
optoSlidingWin2K.files = {	...
	'1155_20171025_02_02_2397_BBN_5.dat', ...				% no opto
	'1155_20171025_02_02_2397_BBN_optoON_10.dat', ...	% opto: 50 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_11.dat', ...	% opto: 100 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_12.dat', ...	% opto: 150 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_13.dat', ...	% opto: 175 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_14.dat', ...	% opto: 200 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_15.dat', ...	% opto: 250 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_16.dat' ...	% opto: 300 ms delay, 2000 mV, 50ms dur
};

optoSlidingWin3K.files = {	...
	'1155_20171025_02_02_2397_BBN_optoON_17.dat', ...	% opto: 100 ms delay, 3000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_18.dat', ...	% opto: 150 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_19.dat', ...	% opto: 200 ms delay, 2000 mV, 50ms dur
	'1155_20171025_02_02_2397_BBN_optoON_20.dat' ...	% opto: 250 ms delay, 2000 mV, 50ms dur
};

% 	acq duration = 400 ms, ISI = 100 ms
% 	15 reps
% 	sound:
% 		tone, freq:5k:5k:80k, 35dB, 100ms delay, 100ms duration, 5ms ramp
freqtuning.files = 	'1155_20171025_02_02_2397_FREQ_TUNING.dat';

% for all files:
% 	acq duration = 600 ms, ISI = 300 ms
% 	15 reps
% 	sound:
% 		tone, freq:5k:1250:25k, 35dB, 150ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 100 ms delay, varied command voltage to ThorLabs LED
optoFreqTuning.files = { ...
	'1155_20171025_02_02_2397_FREQ.dat', ...				% freq,5k:1250:25k opto 0 mV
	'1155_20171025_02_02_2397_FREQ_optoON.dat', ...		% freq, opto 3000 mV
	'1155_20171025_02_02_2397_FREQ_optoON_2.dat', ...	% freq, opto 500 mV
	'1155_20171025_02_02_2397_FREQ_optoON_3.dat', ...	% freq, opto 1000 mV
	'1155_20171025_02_02_2397_FREQ_optoON_4.dat' ...	% freq, opto 1500 mV
};

% for all files:
% 	20 reps
% sound: null, wav files, bbn (65 dB)
%	del 200 ms
% 	opto: 
% 		250 ms duration, 150 ms delay, 3000 mV to ThorLabs LED
optoWav.files = {	...
	'1155_20171025_02_02_2397_OptoInhibOFF.dat', ...
	'1155_20171025_02_02_2397_OptoInhibON.dat' ...
};
optoWav.LEDintensity = [ ...
	0, ...
	3000 ...
];
optoWav.LEDpower = [ ...
	0, ...
	6.1 ...
];

