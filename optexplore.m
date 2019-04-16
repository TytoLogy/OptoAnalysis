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

% Last Modified by GUIDE v2.5 16-Apr-2019 16:37:48

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
		handles.H = varargin{1};
	else
		handles.H = [];
	end
	guidata(hObject, handles);
	
	optexplore_Init(hObject, handles);
	
	% UIWAIT makes optexplore wait for user response (see UIRESUME)
	% uiwait(handles.figure1);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function optexplore_Init(hObject, handles)
	

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
% --- Executes on selection change in listStim.
function listStim_Callback(hObject, eventdata, handles)
% hObject    handle to listStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listStim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listStim
%--------------------------------------------------------------------------

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
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
