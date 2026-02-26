function varargout = waterWatchFig(varargin)
% WATERWATCHFIG MATLAB code for waterWatchFig.fig
%      WATERWATCHFIG, by itself, creates a new WATERWATCHFIG or raises the existing
%      singleton*.
%
%      H = WATERWATCHFIG returns the handle to a new WATERWATCHFIG or the handle to
%      the existing singleton*.
%
%      WATERWATCHFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WATERWATCHFIG.M with the given input arguments.
%
%      WATERWATCHFIG('Property','Value',...) creates a new WATERWATCHFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before waterWatchFig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to waterWatchFig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help waterWatchFig

% Last Modified by GUIDE v2.5 26-Jul-2018 10:02:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @waterWatchFig_OpeningFcn, ...
                   'gui_OutputFcn',  @waterWatchFig_OutputFcn, ...
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

global watProps

% --- Executes just before waterWatchFig is made visible.
function waterWatchFig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to waterWatchFig (see VARARGIN)
global watProps

watProps = startWatProps(handles);

% if ~strcmp(get(watProps.cam,'running'),'on')
%     startCam
% end
set(handles.text_Status,'backgroundcolor',[.6 .6 .6 ]);
set(handles.text_Status,'string','OFF and Waiting')

% Choose default command line output for waterWatchFig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes waterWatchFig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = waterWatchFig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_threshUp.
function pushbutton_threshUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_threshUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global watProps
watProps.thresh1 = watProps.thresh1 + .25;
updateFields(handles)


% --- Executes on button press in pushbutton_threshDown.
function pushbutton_threshDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_threshDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global watProps
watProps.thresh1 = watProps.thresh1 - .25;
updateFields(handles)

% --- Executes on button press in pushbutton_defineWindow.
function pushbutton_defineWindow_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_defineWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global watProps
[inX inY] = ginput

function edit_integrationTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_integrationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_integrationTime as text
%        str2double(get(hObject,'String')) returns contents of edit_integrationTime as a double
global watProps


val = str2num(eventdata.Source.String);
if val
    watProps.intTime = val;
else
    eventdata.Source.String = num2str(watProps.intTime);
end


% --- Executes during object creation, after setting all properties.
function edit_integrationTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_integrationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ON.
function pushbutton_ON_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global watProps
watProps.watching = 1;
set(handles.text_Status,'backgroundcolor',[1 .6 .2 ]);
set(handles.text_Status,'string','Starting camera')
if ~strcmp(get(watProps.cam,'running'),'on')
    startCam
end
if isempty(watProps.cam)
    set(handles.text_Status,'backgroundcolor',[1 .1 .1 ]);
    set(handles.text_Status,'string','Camera failed. Restart camera and try again')
    pause(3)
end
updateFields(handles)
runWatWatch(handles)

% --- Executes on button press in pushbutton_OFF.
function pushbutton_OFF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global watProps
watProps.watching = 0;
watProps.status = .1;
set(handles.text_Status,'backgroundcolor',[.6 .6 .6 ]);
set(handles.text_Status,'string','OFF and Waiting')
try, stop(watProps.cam),end
%stop(watProps.cam)

% --- Executes on button press in pushbutton_getSource.
function pushbutton_getSource_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_getSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global watProps

sourceVal = uigetdir;
if ~isempty(sourceVal)
    watProps.sourceDir = sourceVal;
    set(handles.text_sourceDirectory,'string',watProps.sourceDir);

end


% --- Executes on slider movement.
function slider_upDown_Callback(hObject, eventdata, handles)
% hObject    handle to slider_upDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global watProps
watProps.y1 = round((1-eventdata.Source.Value) * watProps.iHeight);
watProps.y2 = watProps.y1 + watProps.boxHeight;


% --- Executes during object creation, after setting all properties.
function slider_upDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_upDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_leftRight_Callback(hObject, eventdata, handles)
% hObject    handle to slider_leftRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global watProps
watProps.x1 = round(eventdata.Source.Value * watProps.iWidth);
watProps.x2 = watProps.x1 + watProps.boxWidth;

% --- Executes during object creation, after setting all properties.
function slider_leftRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_leftRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_win1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_win1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global watProps
watProps.win1 = eventdata.Source.Value

% --- Executes during object creation, after setting all properties.
function slider_win1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_win1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_win2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_win2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global watProps
watProps.win2 = eventdata.Source.Value


% --- Executes during object creation, after setting all properties.
function slider_win2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_win2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_boxWidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_boxWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_boxWidth as text
%        str2double(get(hObject,'String')) returns contents of edit_boxWidth as a double

global watProps

val = str2num(eventdata.Source.String);
if val
    watProps.boxWidth = val;
    watProps.x2 = watProps.x1 + val;
else
    eventdata.Source.String = num2str(watProps.boxWidth);
end



% --- Executes during object creation, after setting all properties.
function edit_boxWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_boxWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


global watProps
%eventdata.Source.String = num2str(watProps.boxWidth);

function edit_boxHeight_Callback(hObject, eventdata, handles)
% hObject    handle to edit_boxHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_boxHeight as text
%        str2double(get(hObject,'String')) returns contents of edit_boxHeight as a double

global watProps

val = str2num(eventdata.Source.String);
if val
    watProps.boxHeight = val;
    watProps.y2 = watProps.y1 + val;
else
    eventdata.Source.String = num2str(watProps.boxHeight);
end


% --- Executes during object creation, after setting all properties.
function edit_boxHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_boxHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_thresh1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_thresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global watProps
watProps.thresh1 = eventdata.Source.Value * 256
updateFields(handles)


% --- Executes during object creation, after setting all properties.
function slider_thresh1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_thresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_thresh1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresh1 as text
%        str2double(get(hObject,'String')) returns contents of edit_thresh1 as a double
global watProps

val = str2num(eventdata.Source.String);
if val
    watProps.thresh1 = val;
else
    eventdata.Source.String = watProps.thresh1;
end
updateFields(handles)


% --- Executes during object creation, after setting all properties.
function edit_thresh1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton_threshUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_threshUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function pushbutton_threshUp_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_threshUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton_threshDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_threshDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Pump.
function pushbutton_Pump_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Pump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global watProps

watProps.manHist = cat(1,watProps.manHist,datenum(datetime));
set(handles.text_Status,'backgroundcolor',[1 1 .3 ]);
set(handles.text_Status,'string','Manual PUMPING')
triggerPump


% --- Executes on button press in togglebutton_threshPump.
function togglebutton_threshPump_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_threshPump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_threshPump
global watProps

val = eventdata.Source.Value;
watProps.autoOn = val;


function edit_pumpDuration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pumpDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pumpDuration as text
%        str2double(get(hObject,'String')) returns contents of edit_pumpDuration as a double

global watProps
val = str2num(eventdata.Source.String);
if val
    watProps.pumpDuration = val;
else
    eventdata.Source.String = watProps.pumpDuration;
end
updateFields(handles)


% --- Executes during object creation, after setting all properties.
function edit_pumpDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pumpDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pumpInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pumpInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pumpInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_pumpInterval as a double

global watProps
val = str2num(eventdata.Source.String);
if val
    watProps.pumpInterval = val;
else
    eventdata.Source.String = watProps.pumpInterval;
end
updateFields(handles)


% --- Executes during object creation, after setting all properties.
function edit_pumpInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pumpInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_histUnit.
function listbox_histUnit_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_histUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_histUnit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_histUnit


% --- Executes during object creation, after setting all properties.
function listbox_histUnit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_histUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton_sound.
function togglebutton_sound_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_sound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_sound

global watProps

watProps.soundOn = get(handles.togglebutton_sound,'value')


% --- Executes on selection change in listbox_colorType.
function listbox_colorType_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_colorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_colorType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_colorType

global watProps

watProps.colorMode = get(handles.listbox_colorType,'Value');

% --- Executes during object creation, after setting all properties.
function listbox_colorType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_colorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveDat.
function pushbutton_saveDat_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveDat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_record.
function togglebutton_record_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_record
global watProps
watProps.record = get(handles.togglebutton_record,'Value');


% --- Executes on button press in pushbutton_recordDir.
function pushbutton_recordDir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recordDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global watProps

watProps.recordDir = uigetdir;
set(handles.text_sourceDirectory,'String',watProps.recordDir);
