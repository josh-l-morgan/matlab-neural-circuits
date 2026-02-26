function varargout = ImportData(varargin)
%IMPORTDATA MATLAB code file for ImportData.fig
%      IMPORTDATA, by itself, creates a new IMPORTDATA or raises the existing
%      singleton*.
%
%      H = IMPORTDATA returns the handle to a new IMPORTDATA or the handle to
%      the existing singleton*.
%
%      IMPORTDATA('Property','Value',...) creates a new IMPORTDATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ImportData_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      IMPORTDATA('CALLBACK') and IMPORTDATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in IMPORTDATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImportData

% Last Modified by GUIDE v2.5 17-Jun-2020 20:49:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImportData_OpeningFcn, ...
                   'gui_OutputFcn',  @ImportData_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

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


% --- Executes on button press in pushbutton_runVastExport.
function pushbutton_runVastExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_runVastExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_MPN.
function pushbutton_MPN_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MPN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_exportDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exportDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_exportDirectory as text
%        str2double(get(hObject,'String')) returns contents of edit_exportDirectory as a double


% --- Executes during object creation, after setting all properties.
function edit_exportDirectory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exportDirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_exportName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_exportName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_exportName as text
%        str2double(get(hObject,'String')) returns contents of edit_exportName as a double


% --- Executes during object creation, after setting all properties.
function edit_exportName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_exportName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textOut_Callback(hObject, eventdata, handles)
% hObject    handle to textOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textOut as text
%        str2double(get(hObject,'String')) returns contents of textOut as a double


% --- Executes during object creation, after setting all properties.
function textOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_makeCellLibrary.
function pushbutton_makeCellLibrary_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_makeCellLibrary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

makeCellLibrary
