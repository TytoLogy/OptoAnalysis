function varargout = optexplore(varargin)
% OPTEXPLORE MATLAB code for optexplore.fig
%      OPTEXPLORE, by itself, creates a new OPTEXPLORE or raises the existing
%      singleton*.
%
%      H = OPTEXPLORE returns the handle to a new OPTEXPLORE or the handle to
%      the existing singleton*.
%
%      OPTEXPLORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTEXPLORE.M with the given input arguments.
%
%      OPTEXPLORE('Property','Value',...) creates a new OPTEXPLORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before optexplore_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to optexplore_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help optexplore

% Last Modified by GUIDE v2.5 15-May-2019 15:24:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optexplore_OpeningFcn, ...
                   'gui_OutputFcn',  @optexplore_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%--------------------------------------------------------------------------


%------------------------------------------------------------------------------
% OPENING/INIT FUNCTIONS
%------------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes just before optexplore is made visible.
function optexplore_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optexplore (see VARARGIN)

	% Choose default command line output for optexplore
	handles.output = hObject;

	% Update handles structure
	guidata(hObject, handles);

	% assign data from varargin
	if ~isempty(varargin)
		handles.inArgs = varargin{1};
	else
		handles.inArgs = [];
	end
	guidata(hObject, handles);
	
	optexplore_Init(hObject, handles);
	
	% UIWAIT makes optexplore wait for user response (see UIRESUME)
	% uiwait(handles.figOptExplore);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function optexplore_Init(hObject, handles)

	H = struct( ...
			'data_root_path', '', ...
			'tytology_root_path', '', ...
			'dataLoaded', 0, ...
			'D', [], ...
			'Dinf', [], ...
			'traces', [], ...
			'spikes', [], ...
			'varlist', [], ...
			'nvars', [], ...
			'threshold', [] ...
			);
		
	% need to get information about system
	[H.data_root_path, H.tytology_root_path] = optoanalysis_paths;
	
	% check input args
	if ~isempty(handles.inArgs)
		% check if input args were provided
		if isstruct(handles.inArgs{1})
			H.D = [];
			H.Dinf = handles.inArgs{1}.Dinf;
			H.traces = handles.inArgs{1}.traces;
			H.spikes = handles.inArgs{1}.spikes;
			H.varlist = handles.inArgs{1}.varlist;
			H.nvars = handles.inArgs{1}.nvars;
			H.dataLoaded = 1;
		end
	end
	
	% store H in handles
	handles.H = H;
	guidata(hObject, handles);
	
	% update gui
	% updateGui(hObject, handles);
	
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = optexplore_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;
%--------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% CALLBACKS
%------------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function listStim_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns listStim
%			contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listStim
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function popSweep_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popSweep contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popSweep



%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% Menu FUNCTIONS
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%--------------------------------------------------------------------
%--------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)



%--------------------------------------------------------------------
%--------------------------------------------------------------------
function menuOpenFile_Callback(hObject, eventdata, handles)
	% get data file from user
	[datafile, datapath] = ...
									uigetfile('*.dat', 'Select opto data file', ...
														handles.H.data_root_path);
	% abort if cancelled...
	if datafile == 0
		fprintf('Cancelled\n');
		return
	% or store info
	else
		handles.H.datafile = datafile;
		handles.H.datapath = datapath;
		handles.H.data_root_path = datapath;
		update_ui_str(handles.fileInfo, ...
								fullfile(handles.H.datapath, handles.H.datafile));
		guidata(hObject, handles);
	end	

	% Read Data and store in handles
	[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(datapath, datafile), ...
											'Filter', [HPFreq LPFreq], ...
											'Channel', channelNumber);
	if isempty(D)
		warning('%s: D is empty???!!!??!!', mfilename);
		handles.H.dataLoaded = false;
		return
	end
	handles.H.D = D;
	handles.H.Dinf = Dinf;
	handles.H.traces = tracesByStim;
	handles.H.dataLoaded = true;
	handles.H.finfo = parse_dat_filename(handles.H.datafile);
	guidata(hObject, handles);

	% Some test-specific things...
	switch upper(Dinf.test.Type)
		case 'FREQ'
			% list of frequencies, and # of freqs tested
			varlist = Dinf.test.stimcache.vrange;
			nvars = length(varlist);
			titleString = cell(nvars, 1);
			for v = 1:nvars
				if v == 1
					titleString{v} = {fname, ...
											sprintf('Frequency = %.0f kHz', 0.001*varlist(v))};
				else
					titleString{v} = sprintf('Frequency = %.0f kHz', 0.001*varlist(v));
				end
			end
		case 'LEVEL'
			% list of levels, and # of levels tested
			varlist = Dinf.test.stimcache.vrange;
			nvars = length(varlist);
			titleString = cell(nvars, 1);
			for v = 1:nvars
				if v == 1
					titleString{v} = {fname, sprintf('Level = %d dB SPL', varlist(v))};
				else
					titleString{v} = sprintf('Level = %d dB SPL', varlist(v));
				end
			end
		case 'FREQ+LEVEL'
			% list of freq, levels
			varlist = cell(2, 1);
			% # of freqs in nvars(1), # of levels in nvars(2)
			nvars = zeros(2, 1);
			for v = 1:2
				varlist{v} = unique(Dinf.test.stimcache.vrange(v, :), 'sorted');
				nvars(v) = length(varlist{v});
			end
			titleString = fname;
		case 'OPTO'
			% not yet implemented
		case 'WAVFILE'
			% get list of stimuli (wav file names)
			varlist = Dinf.test.wavlist;
			nvars = length(varlist);
			titleString = cell(nvars, 1);
			for v = 1:nvars
				if v == 1 
					titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
				else
					titleString{v} = sprintf('wav name: %s', varlist{v});
				end
			end
		otherwise
			error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
	end

	handles.H.varlist = varlist;
	handles.H.titleString = titleString;
	handles.H.nvars = nvars;


%--------------------------------------------------------------------


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% CREATE FUNCTIONS
%------------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
%--------------------------------------------------------------------------
function listStim_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
												get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
% --- Executes during object creation, after setting all properties.
function popSweep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
