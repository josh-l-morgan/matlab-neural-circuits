function varargout = ImportData(varargin)
% IMPORTDATA MATLAB code for ImportData.fig
%      IMPORTDATA, by itself, creates a new IMPORTDATA or raises the existing
%      singleton*.
%
%      H = IMPORTDATA returns the handle to a new IMPORTDATA or the handle to
%      the existing singleton*.
%
%      IMPORTDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORTDATA.M with the given input arguments.
%
%      IMPORTDATA('Property','Value',...) creates a new IMPORTDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImportData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImportData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImportData

% Last Modified by GUIDE v2.5 17-Jun-2020 09:47:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImportData_OpeningFcn, ...
                   'gui_OutputFcn',  @ImportData_OutputFcn, ...
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


% --- Executes just before ImportData is made visible.
function ImportData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImportData (see VARARGIN)

% Choose default command line output for ImportData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImportData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImportData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_Close.
function pushbutton_Close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1)


% --- Executes on button press in pushbutton_runVastExport.
function pushbutton_runVastExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_runVastExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


