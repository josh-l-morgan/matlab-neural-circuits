function varargout = AnaGui4(varargin)
% ANAGUI4 M-file for AnaGui4.fig
%      ANAGUI4, by itself, creates a new ANAGUI4 or raises the existing
%      singleton*.
%
%      H = ANAGUI4 returns the handle to a new ANAGUI4 or the handle to
%      the existing singleton*.
%
%      ANAGUI4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANAGUI4.M with the given input arguments.
%
%      ANAGUI4('Property','Value',...) creates a new ANAGUI4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnaGui4_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnaGui4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnaGui4

% Last Modified by GUIDE v2.5 04-Jan-2007 16:25:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AnaGui4_OpeningFcn, ...
                   'gui_OutputFcn',  @AnaGui4_OutputFcn, ...
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


% --- Executes just before AnaGui4 is made visible.
function AnaGui4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnaGui4 (see VARARGIN)

% Choose default command line output for AnaGui4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global DPN DFN TPN 
colormap gray(255)
htext = uicontrol('Style','text','String','Select a File','Position',[625,480,100,40]);
hrtext = uicontrol('Style','text','String','Ready to Begin','Position',[625,40,100,20]);
pause(.01)





% UIWAIT makes AnaGui4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = AnaGui4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in FindThresh.
function FindThresh_Callback(hObject, eventdata, handles)
% hObject    handle to FindThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FindThresh


% --- Executes on button press in Skeletonize.
function Skeletonize_Callback(hObject, eventdata, handles)
% hObject    handle to Skeletonize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Skeletonize


% --- Executes on button press in Ratio.
function Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ratio


% --- Executes on button press in Mask.
function Mask_Callback(hObject, eventdata, handles)
% hObject    handle to Mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mask


% --- Executes on button press in Draw.
function Draw_Callback(hObject, eventdata, handles)
% hObject    handle to Draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Draw


% --- Executes on button press in DotDend.
function DotDend_Callback(hObject, eventdata, handles)
% hObject    handle to DotDend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DotDend


% --- Executes on button press in Sholl.
function Sholl_Callback(hObject, eventdata, handles)
% hObject    handle to Sholl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sholl


% --- Executes on button press in FixData.
function FixData_Callback(hObject, eventdata, handles)
% hObject    handle to FixData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixData


% --- Executes on button press in FindFile.
function FindFile_Callback(hObject, eventdata, handles)
% hObject    handle to FindFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DPN DFN TPN
GetFolder
htext = uicontrol('Style','text','String',DPN,'Position',[625,480,100,40]);


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DPN DFN TPN

hrtext = uicontrol('Style','text','String','Running','Position',[625,40,100,20]);
pause(.01)


if get(handles.FindThresh,'Value'); anaFT; end
if get(handles.Skeletonize,'Value'); anaSk; end
if get(handles.DotFind,'Value'); anaDF; end
if get(handles.FixData,'Value'); anaFix;end
if get(handles.Ratio,'Value'); anaRa; end
if get(handles.Mask,'Value'); anaMa; end
if get(handles.Draw,'Value'); anaDr; end
if get(handles.DotDend,'Value'); anaDD; end
if get(handles.Sholl,'Value'); anaSh; end
if exist([TPN 'data\DotStats']), anaNN; end

hrtext = uicontrol('Style','text','String','Done','Position',[625,40,100,20]);
pause(.01)

% --- Executes on button press in DotFind.
function DotFind_Callback(hObject, eventdata, handles)
% hObject    handle to DotFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DotFind


