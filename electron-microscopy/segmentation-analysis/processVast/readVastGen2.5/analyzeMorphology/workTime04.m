function varargout = workTime02(varargin)
% WORKTIME02 MATLAB code for workTime02.fig
%      WORKTIME02, by itself, creates a new WORKTIME02 or raises the existing
%      singleton*.
%
%      H = WORKTIME02 returns the handle to a new WORKTIME02 or the handle to
%      the existing singleton*.
%
%      WORKTIME02('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORKTIME02.M with the given input arguments.
%
%      WORKTIME02('Property','Value',...) creates a new WORKTIME02 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before workTime02_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to workTime02_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help workTime02

% Last Modified by GUIDE v2.5 19-Dec-2014 22:15:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @workTime02_OpeningFcn, ...
                   'gui_OutputFcn',  @workTime02_OutputFcn, ...
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


% --- Executes just before workTime02 is made visible.
function workTime02_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to workTime02 (see VARARGIN)

% Choose default command line output for workTime02
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes workTime02 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%loadDay(handles)
%pause(1)
%countTime(handles)


timeFolder = [pwd '\times\'];
if ~exist(timeFolder,'dir'),mkdir(timeFolder);end

figVar.targ.work = 10;
figVar.targ.play = 2;
figVar.targ.life = 3;
figVar.targ.other =1; 
figVar.tolerance = 15; %minutes
guidata(hObject,figVar);


% --- Outputs from this function are returned to the command line.
function varargout = workTime02_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function currentStatus_CreateFcn(hObject, eventdata, handles)

% --- Executes when selected object is changed in currentStatus.
function currentStatus_SelectionChangeFcn(hObject, eventdata, handles)
newChoice = get(eventdata.NewValue,'Tag');
saveDay(handles);
%countTime(handles,newChoice)
if ~exist('newChoice','var')
    newChoice = 'otherButton';
end
switch newChoice
    case 'workButton'
        textTag = handles.workTime;
    case 'playButton'
         textTag = handles.playTime;
    case 'lifeButton'
        textTag = handles.lifeTime;
    case 'otherButton'
       textTag = handles.otherTime;
    otherwise
        disp('not sure')
end

lastTime = datenum(clock);
while 1
    pause(1)
    
    %%Change times (should be fast)
    timeVal = get(textTag,'Value');
    currentTime = datenum(clock);
    elapsedTime = (currentTime - lastTime);
    newTime = timeVal + elapsedTime;
    set(textTag,'Value',newTime);
    lastTime = currentTime;
    
    [Y, M, D, H, MN, S]  = datevec(newTime);
    timeString = sprintf('%d:%02.0f:%02.0f',D*24+H,MN,S);
    set(textTag,'String',timeString)
    if ~S
        try saveDay(handles);
        catch err
            err
        end
    end
    
end


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)

loadDay(handles)

% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)

saveDay(handles);


% --- Executes during object creation, after setting all properties.
function loadButton_CreateFcn(hObject, eventdata, handles)

function saveDay(handles)


times.workString = get(handles.workTime,'String')
times.playString = get(handles.playTime,'String')
times.lifeString = get(handles.lifeTime,'String')
times.otherString = get(handles.otherTime,'String')

times.workValue = get(handles.workTime,'Value')
times.playValue = get(handles.playTime,'Value')
times.lifeValue = get(handles.lifeTime,'Value')
times.otherValue = get(handles.otherTime,'Value')

timeFolder = [pwd '\times\'];
[Y, M, D, H, MN, S]  = datevec(datenum(clock));
timeName = sprintf('times_%04.0f+%02.0f+%02.0f',Y,M,D);
save([timeFolder 'times_' timeName '.mat'],'times');

function loadDay(handles)

timeFolder = [pwd '\times\'];
[Y, M, D, H, MN, S]  = datevec(datenum(clock));
timeName = sprintf('times_%04.0f+%02.0f+%02.0f',Y,M,D);
load([timeFolder 'times_' timeName '.mat']);


set(handles.workTime,'Value',times.workValue )
set(handles.playTime,'Value',times.playValue)
set(handles.lifeTime,'Value',times.lifeValue)
set(handles.otherTime,'Value',times.otherValue)


set(handles.workTime,'String',times.workString)
set(handles.playTime,'String',times.playString)
set(handles.lifeTime,'String',times.lifeString)
set(handles.otherTime,'String',times.otherString)

function countTime(handles,newChoice)



function workTarget_Callback(hObject, eventdata, handles)
% hObject    handle to workTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of workTarget as text
%        str2double(get(hObject,'String')) returns contents of workTarget as a double

% 
% % --- Executes during object creation, after setting all properties.
% function workTarget_CreateFcn(hObject, eventdata, handles)
% % figVar = guidata(hObject);
% % set(handles.workTarget,'String',num2str(figVar.targ.work));
% 
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes during object creation, after setting all properties.
function workTarget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to workTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
