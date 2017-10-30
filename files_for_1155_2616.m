% animal id 
animalID = '1155';
% date code (helps locate files in data directory)
dateID = '20171025';
% unit #
unit = '05';
% penetration
penetration = '03';
% recording depth
depth = '2616';
% channel
channelNumber = 8;
% files

% sound: dur 100 del 100 ramp 5 acqdur 300 isi 250 0:10:70 dB
% opto: dur 200, del 50, 3500 mV (6.8 mW)
bbn.files = { ...
	'1155_20171025	_05_03_2616_BBN.dat', ...			% rate level function BBN 0:10:70
	'1155_20171025_05_03_2616_BBN_optoON.dat' ...	% rlf, opto on, 3500 mV
};
	
% for all files:
% 	acq duration = 500 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], 0:10:70dB, 100ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		200 ms duration, 50 ms delay, varied command voltage to ThorLabs LED
optoRLF.files = {	...
	'1155_20171025_05_03_2616_BBN_2.dat', ...				% rlf, 0:10:70
	'1155_20171025_05_03_2616_BBN_optoON_2.dat', ...	% rlf, 0:10:70, opto 3500 mV
	'1155_20171025_05_03_2616_BBN_optoON_3.dat', ...	% rlf, 0:10:70, opto 100 mV
	'1155_20171025_05_03_2616_BBN_optoON_4.dat', ...	% rlf, 0:10:70, opto 500 mV
	'1155_20171025_05_03_2616_BBN_optoON_5.dat', ...	% rlf, 0:10:70, opto 1000 mV
	'1155_20171025_05_03_2616_BBN_optoON_6.dat', ...	% rlf, 0:10:70, opto 750 mV
	'1155_20171025_05_03_2616_BBN_optoON_7.dat'  ...	% rlf, 0:10:70, opto 1500 mV
};
optoRLF.LEDintensity = [ ...
	0, ...
	3500, ...
	100, ...
	500, ...
	1000, ...
	750, ...
	1500 ...
];
optoRLF.LEDpower = [ ...
	0.0, ...
	6.8, ...
	0.04, ...
	1.4, ...
	2.5, ...
	2.0, ...
	3.6 ...
}

% 0		0.0
% 250		0.7
% 500		1.4
% 750		2.0
% 1000		2.5
% 1500		3.6
% 2000		4.5
% 2500		5.3
% 3000		6.1
% 3500		6.8	

% sliding "pulse" of opto stimulation
% for all files:
% 	acq duration = 500 ms, ISI = 250 ms
% 	20 reps
% 	sound:
% 		noise, BW:[4kHz 80kHz], [0 40 60] dB, 100ms delay, 100ms duration, 5ms ramp
% 	opto: 
% 		25 ms duration, delay varied, 3000 mV command voltage to ThorLabs LED
optoSlidingWin_files = {	...
	'1155_20171025_05_03_2616_BBN_3.dat', ...				% sliding opto, control, BBN [0,40,60]dB; 
	'1155_20171025_05_03_2616_BBN_optoON_9.dat', ...	% opto: 50 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_10.dat', ...	% opto: 100 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_11.dat', ...	% opto: 150 ms delay, 3000 mV, 25ms dur
	'1155_20171025_05_03_2616_BBN_optoON_12.dat'  ...	% opto: 200 ms delay, 3000 mV, 25ms dur
};

misc_files = { ...
	'1155_20171025_05_03_2616_BBN_optoON_13.dat', ...	% opto: 50 ms delay, 3000 mV, 50ms dur
	'1155_20171025_05_03_2616_BBN_optoON_14.dat', ...	% opto: 100 ms delay, 3000 mV, 50ms dur
	'1155_20171025_05_03_2616_BBN_optoON_15.dat' ...	% opto: 150 ms delay, 3000 mV, 50ms dur
};
