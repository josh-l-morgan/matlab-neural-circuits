function varargout = cellNavGuide_14(varargin)
% CELLNAVGUIDE_14 MATLAB code for cellNavGuide_14.fig
%      CELLNAVGUIDE_14, by itself, creates a new CELLNAVGUIDE_14 or raises the existing
%      singleton*.
%
%      H = CELLNAVGUIDE_14 returns the handle to a new CELLNAVGUIDE_14 or the handle to
%      the existing singleton*.
%
%      CELLNAVGUIDE_14('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLNAVGUIDE_14.M with the given input arguments.
%
%      CELLNAVGUIDE_14('Property','Value',...) creates a new CELLNAVGUIDE_14 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cellNavGuide_14_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cellNavGuide_14_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cellNavGuide_14

% Last Modified by GUIDE v2.5 05-Jun-2020 21:12:34




% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cellNavGuide_14_OpeningFcn, ...
    'gui_OutputFcn',  @cellNavGuide_14_OutputFcn, ...
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

% --- Executes just before cellNavGuide_14 is made visible.
function cellNavGuide_14_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellNavGuide_14 (see VARARGIN)

% Choose default command line output for cellNavGuide_14
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cellNavGuide_14 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global glob tis
glob.handles = handles;



function makeGlob

global glob tis

glob.start = 1;
glob.fvDir = ['..\fvLibrary\'];
glob.swc = ['..\swc\'];
glob.p = [];
clear gca
load([glob.fvDir 'obI.mat'])
% load([MPN 'dsObj.mat'])
load([glob.fvDir 'tis.mat'])
%
% glob.MPN = MPN;
% glob.WPN = WPN;

glob.pickIdx = 1;
glob.pickCID = tis.cids(glob.pickIdx);

glob.param.markerSize = 100;
glob.fvRes = 0.1;
glob.em = obI.em;

glob.cellNum = length(tis.cells.cids)
glob.listCellidx = 1:glob.cellNum;
glob.typeID = 0;
glob.g.idx = {};

if 1;%~exist(glob.fvDir,'dir')
    
    glob.fvDir = uigetdir('..','select fvLibrary folder');
    glob.fvDir = [glob.fvDir '\'];
    
    
end






% --- Outputs from this function are returned to the command line.
function varargout = cellNavGuide_14_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in toggleBorders.
function toggleBorders_Callback(hObject, eventdata, handles)
% hObject    handle to toggleBorders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleBorders

global glob

val = get(hObject,'Value');
try
    isHand = isvalid(glob.p.pRefINL);
catch err
    isHand = 0;
end

if val
    if isHand
        glob.p.pRefGCL.Visible = 'on';
        glob.p.pRefINL.Visible = 'on';
    else
        load([glob.fvDir 'ref_gcl nucEdge.mat'])
        glob.p.pRefGCL = renderFVnav(fv,[1 1 0],.5)
        
        load([glob.fvDir 'ref_inl nucEdge.mat'])
        glob.p.pRefINL = renderFVnav(fv,[1 0 1],.5)
    end
else
    if isHand
        glob.p.pRefGCL.Visible = 'off'
        glob.p.pRefINL.Visible = 'off'
    end
    
end




% --- Executes on button press in toggleScaleBar.
function toggleScaleBar_Callback(hObject, eventdata, handles)
% hObject    handle to toggleScaleBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleScaleBar
global glob

val = get(hObject,'Value')
try
    isHand = isvalid(glob.p.pRefMicronBar);
catch err
    isHand = 0;
end

if val
    if isHand
        glob.p.pRefMicronBar.Visible = 'on';
        glob.p.pRefMicronPt1.Visible = 'on';
        glob.p.pRefMicronPt2.Visible = 'on';
    else
        load([glob.fvDir 'ref_10 micron bar.mat'])
        glob.p.pRefMicronBar = renderFVnav(fv,[1 1 1],1)
        
        load([glob.fvDir 'ref_10 micron point 1.mat'])
        glob.p.pRefMicronPt1 = renderFVnav(fv,[1 1 1],.5)
        
        load([glob.fvDir 'ref_10 micron point 2.mat'])
        glob.p.pRefMicronPt2 = renderFVnav(fv,[1 1 1],.5)
    end
else
    if isHand
        glob.p.pRefMicronBar.Visible = 'off';
        glob.p.pRefMicronPt1.Visible = 'off';
        glob.p.pRefMicronPt2.Visible = 'off';
    end
end



% --- Executes on button press in toggleBipolars.
function toggleBipolars_Callback(hObject, eventdata, handles)
% hObject    handle to toggleBipolars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleBipolars

global glob tis

val = get(hObject,'Value')

try
    isHand = isvalid(glob.p.pBipolars);
catch err
    isHand = 0;
end

if val
    if isHand
        glob.p.pBipolars.Visible = 'on';
    else
        bipCids = tis.cells.type.typeLists.bpc;
        vert = [];
        fac = [];
        for i = 1:length(bipCids)
            fileName = sprintf('%s%d.mat',glob.fvDir,bipCids(i));
            load(fileName);
            fac = cat(1,fac,fv.faces+size(vert,1));
            vert = cat(1,vert,fv.vertices);
        end
        
        fv.vertices = vert;
        fv.faces = fac;
        glob.p.pBipolars = renderFVnav(fv,[0 .8 0],.1)
        
    end
else
    if isHand
        glob.p.pBipolars.Visible = 'off';
    end
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global glob

str = get(hObject,'string')
val = str2num(str);

glob.pickCID = val;
showPickedCell








% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupCellType.
function popupCellType_Callback(hObject, eventdata, handles)
% hObject    handle to popupCellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupCellType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupCellType

global glob tis


val = get(hObject,'Value');

if val == 1
    glob.listCellidx = 1:glob.cellNum;
else
    glob.typeID =glob.typeIDs(val);
    glob.listCellidx = find(tis.cells.type.typeID == glob.typeID);
    
end

set(handles.listCIDs,'String',glob.cellStr(glob.listCellidx))
glob.pickIdx = 1;
glob.pickCID = tis.cids(glob.pickIdx);



if glob.typeID < 1
    subTypeNames = {'all'};
else
    subTypeNames = tis.cells.type.subTypeNames{glob.typeID};
    subTypeNames = cat(2,{'all'},subTypeNames);
end

set(handles.subType_popupMenue,'String',subTypeNames)
set(handles.subType_popupMenue,'Value',1)





% --- Executes during object creation, after setting all properties.
function popupCellType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupCellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global glob tis

typeID = tis.cells.type.typeID;
hasType = unique(typeID);
typeStrings = {'all'};
if sum(hasType==0)
    typeStrings{2} = 'unassigned';
end

isType = hasType(hasType>0);
glob.typeIDs = [0 0 isType];
typeStrings = cat(2,typeStrings,...
    {tis.cells.type.typeNames{isType}});

glob.typeStrings = typeStrings;
set(hObject,'String',glob.typeStrings);







% --- Executes on selection change in listCIDs.
function listCIDs_Callback(hObject, eventdata, handles)
% hObject    handle to listCIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listCIDs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listCIDs

global glob tis

val = get(hObject,'Value');
glob.pickIdx = glob.listCellidx(val);
glob.pickCID = tis.cids(glob.pickIdx);
showPickedCell





% --- Executes during object creation, after setting all properties.
function listCIDs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listCIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


global glob tis

% tis.cells.cids
% tis.cells.label

for i = 1:glob.cellNum
    cellStr{i} = sprintf('%d    %s',tis.cells.cids(i),tis.cells.label{i}) ;
end

glob.cellStr = cellStr;

glob.listCellidx = 1:glob.cellNum;


set(hObject,'String',glob.cellStr(glob.listCellidx) )





% --- Executes during object creation, after setting all properties.
function mainAx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainAx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainAx

global glob

view([0 0])
axis off

set(gca,'color',[0 0 0])
%set(gcf,'color',[0 0 0])
set(gca,'visible','off')



% --- Executes on mouse press over axes background.
function mainAx_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainAx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in viewXY.
function viewXY_Callback(hObject, eventdata, handles)
% hObject    handle to viewXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

view([0 270])


% --- Executes on button press in viewXZ.
function viewXZ_Callback(hObject, eventdata, handles)
% hObject    handle to viewXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


view([0 5 ])

% --- Executes during object creation, after setting all properties.
function viewXZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

a = 3;

% --- Executes on button press in viewYZ.
function viewYZ_Callback(hObject, eventdata, handles)
% hObject    handle to viewYZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

view([92 1])


% --- Executes on button press in prevCell.
function prevCell_Callback(hObject, eventdata, handles)
% hObject    handle to prevCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global glob tis
glob.pickIdx = max(glob.pickIdx - 1, 1);
glob.pickCID = tis.cids(glob.pickIdx);
set(handles.listCIDs,'Value',glob.pickIdx)
showPickedCell

% --- Executes on button press in nextCell.
function nextCell_Callback(hObject, eventdata, handles)
% hObject    handle to nextCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob tis

glob.pickIdx = min(glob.pickIdx + 1, length(glob.listCellidx));
glob.pickCID = tis.cids(glob.pickIdx);
set(handles.listCIDs,'Value',glob.pickIdx)

showPickedCell


function showPickedCell

global glob tis

set(glob.handles.edit1,'String',num2str(glob.pickCID))

pushPostSyn_Callback(glob.handles.pushPostSyn);
pushPreSyn_Callback(glob.handles.pushPreSyn);
pushPostCells_Callback(glob.handles.pushPostCells);
pushPreCell_Callback(glob.handles.pushPreCell);

fileName = sprintf('%s%d.mat',glob.fvDir,glob.pickCID);
load(fileName);
%     set(handles.mainAx)
%     patch(handles.mainAx,fv);

try
    glob.p.cell.delete
end
glob.p.cell = renderFVnav(fv,[1 1 1],1,glob.pickCID);
axis off
pause(.01)

toggleSkeleton_Callback(glob.handles.toggleSkeleton);






% --- Executes on button press in raiseLight.
function raiseLight_Callback(hObject, eventdata, handles)
% hObject    handle to raiseLight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob

a = rand * 360 - 180;
e = rand * 180 - 90;
lightangle(glob.light,a,e);
%
% [a e] = view(glob.handles.mainAx)
% a = mod(a+180 + 180,360)-180
% lightangle(glob.light,a,e);
% [a e] = lightangle(glob.light)

% lightangle(glob.light,-179, -89);
% for i = 1:1000
% [a e] = lightangle(glob.light)
% lightangle(glob.light,mod(a+197,179)-179, mod(e + 92,89)-89);
% pause(.04)
% end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

global glob tis

str = get(hObject,'string');

if isempty(str)
    val = [];
else
    val = str2num(str{1});
end


glob.pickCIDref = val;



try
    for r = 1:length(glob.p.refCell)
        glob.p.refCell(r).delete
    end
end

for i = 1:length(val)
    
    fileName = sprintf('%s%d.mat',glob.fvDir,val(i));
    load(fileName);
    
    glob.p.refCell(i) = renderFVnav(fv,[1 0 0],1,val(i));
    
end




% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

makeGlob


% --- Executes on button press in pushPreCell.
function pushPreCell_Callback(hObject, eventdata, handles)
% hObject    handle to pushPreCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global glob tis

val = get(hObject,'Value');



try
    glob.p.preCells.delete
end


if val
    
    targ = glob.pickCID;
    
    isSyn = find(tis.syn.post == targ);
    isPre = unique(tis.syn.pre(isSyn));
    isPre = setdiff(isPre,0);
    showCells = isPre;
    vert = [];
    fac = [];
    for i = 1:length(showCells)
        fileName = sprintf('%s%d.mat',glob.fvDir,showCells(i));
        if exist(fileName,'file')
            load(fileName);
            fac = cat(1,fac,fv.faces+size(vert,1));
            vert = cat(1,vert,fv.vertices);
        end
    end
    
    fv.vertices = vert;
    fv.faces = fac;
    glob.p.preCells = renderFVnav(fv,[.2 .4 1],.6)
    
    preText = sprintf('pre to %d = %s',targ,num2str(isPre'));
    showText(preText);
    
    
end

function showText(textStr)
global glob
set(glob.handles.textOut,'String',textStr );


% --- Executes on button press in pushPostCells.
function pushPostCells_Callback(hObject, eventdata, handles)
% hObject    handle to pushPostCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global glob tis

val = get(hObject,'Value');

try
    glob.p.postCells.delete
end

if val
    
    targ = glob.pickCID;
    
    isSyn = find(tis.syn.pre == targ);
    isPost = unique(tis.syn.post(isSyn));
    isPost = setdiff(isPost,0);
    showCells = isPost;
    vert = [];
    fac = [];
    for i = 1:length(showCells)
        fileName = sprintf('%s%d.mat',glob.fvDir,showCells(i));
        if exist(fileName,'file')
            load(fileName);
            fac = cat(1,fac,fv.faces+size(vert,1));
            vert = cat(1,vert,fv.vertices);
        end
    end
    
    fv.vertices = vert;
    fv.faces = fac;
    glob.p.postCells = renderFVnav(fv,[1 .4 .2],.6)
    
    postText = sprintf('post to %d = %s',targ,num2str(isPost'));
    showText(postText);
    
    
end


% --- Executes on button press in pushPreSyn.
function pushPreSyn_Callback(hObject, eventdata, handles)
% hObject    handle to pushPreSyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global glob tis

val = get(hObject,'Value');

try glob.p.preSyn1.delete
end
try glob.p.preSyn2.delete
end
if val
    
    targ = glob.pickCID;
    isSyn = find(tis.syn.post == targ);
    pos = tis.syn.synPosDS(isSyn,[2 1 3]);
    preClass = tis.syn.preClass(isSyn);
    col = repmat([1 0 0],[length(isSyn) 1]);
    isBip = find(preClass==7);
    hold on
    glob.p.preSyn1 = scatter3(pos(:,1),pos(:,2),pos(:,3),glob.param.markerSize,'markerfacecolor','r',...
        'marker','o','markeredgecolor','w');
    glob.p.preSyn2 = scatter3(pos(isBip,1),pos(isBip,2),pos(isBip,3),glob.param.markerSize,'markerfacecolor','g',...
        'marker','o','markeredgealpha',0);
    set(glob.p.preSyn1,'clipping','off')
    set(glob.p.preSyn2,'clipping','off')
    
    hold off
    
    preSyn = tis.syn.pre(isSyn);
    %postText = sprintf('preSyn to %d = %s',targ,num2str(preSyn'));
    postText = sprintf('all Pre = %d.  Bip post = %d',length(isSyn),length(isBip));
    
    showText(postText);
end


% --- Executes on button press in pushPostSyn.
function pushPostSyn_Callback(hObject, eventdata, handles)
% hObject    handle to pushPostSyn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob tis

val = get(hObject,'Value');

try glob.p.postSyn1.delete
end
try glob.p.postSyn2.delete
end
if val
    
    targ = glob.pickCID;
    isSyn = find(tis.syn.pre == targ);
    pos = tis.syn.synPosDS(isSyn,[2 1 3]);
    postClass = tis.syn.postClass(isSyn);
    col = repmat([1 0 0],[length(isSyn) 1]);
    isRGC = find(postClass==1);
    hold on
    glob.p.postSyn1 = scatter3(pos(:,1),pos(:,2),pos(:,3),glob.param.markerSize,'markerfacecolor','m',...
        'marker','o','markeredgecolor','w');
    glob.p.postSyn2 = scatter3(pos(isRGC,1),pos(isRGC,2),pos(isRGC,3),glob.param.markerSize,'markerfacecolor','b',...
        'marker','o','markeredgealpha',0);
    set(glob.p.postSyn1,'clipping','off')
    set(glob.p.postSyn2,'clipping','off')
    
    hold off
    
    
    preSyn = tis.syn.post(isSyn);
    %postText = sprintf('postSyn to %d = %s',targ,num2str(preSyn'));
    
    postText = sprintf('allPost = %d.  RGC post = %d',length(isSyn),length(isRGC));
    
    showText(postText);
    
end


% --- Executes on button press in pushShrinkMarker.
function pushShrinkMarker_Callback(hObject, eventdata, handles)
% hObject    handle to pushShrinkMarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob

glob.param.markerSize = glob.param.markerSize ^ .9;
glob.param.markerSize = max(1,glob.param.markerSize);
pushPostSyn_Callback(glob.handles.pushPostSyn)
pushPreSyn_Callback(glob.handles.pushPreSyn)


% --- Executes on button press in pushGrowMarker.
function pushGrowMarker_Callback(hObject, eventdata, handles)
% hObject    handle to pushGrowMarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global glob

glob.param.markerSize = glob.param.markerSize ^ 1.1;
glob.param.markerSize = max(1,glob.param.markerSize);
pushPostSyn_Callback(glob.handles.pushPostSyn)
pushPreSyn_Callback(glob.handles.pushPreSyn)


% --------------------------------------------------------------------
function uitoolbar1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitoolbar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uiDataTip_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uiDataTip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushDataTip.
function pushDataTip_Callback(hObject, eventdata, handles)
% hObject    handle to pushDataTip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushDataTip
global glob

val = get(hObject,'Value');

% try
%     glob.dcm.Enable = 'off'
%     %glob.dcm.delete
% end

if val
    
    glob.dcm = datacursormode(handles.figure1);
    glob.dcm.UpdateFcn = @displayCoordinates;
    
    set(glob.dcm,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
else
    rotate3d(handles.figure1)
    delete(findall(handles.figure1,'Type','hggroup'))
end


function txt = displayCoordinates(~,info)
global glob
X = info.Position(2);
Y = info.Position(1);
Z = info.Position(3);
%c_info = getCursorInfo(glob.dcm)

anc2sub = (glob.em.res/ 1000)./ glob.em.dsRes;
allPos2 = [X Y Z] ./ anc2sub;
VASTpos = sprintf('%.0f  %.0f  %.0f',allPos2(end,2),allPos2(end,1),allPos2(end,3));
disp(VASTpos)

anc2um =  1./ glob.em.dsRes;
umPos = [X Y Z] ./ anc2um;
umPosStr = sprintf('%.0f \r%.0f \r%.0f um',umPos(end,2),umPos(end,1),umPos(end,3));

set(glob.handles.textXYZ,'String',VASTpos);
tag = get(info.Target,'Tag');
set(glob.handles.patchTagText,'String',['tag = ' tag]);


txt = umPosStr;

clipboard('copy',VASTpos)






% --- Executes during object creation, after setting all properties.
function raiseLight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to raiseLight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global glob
'bark'
glob.light = lightangle(180,45)


% --- Executes on button press in toggleSkeleton.
function toggleSkeleton_Callback(hObject, eventdata, handles)
% hObject    handle to toggleSkeleton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleSkeleton

global glob

try
    glob.p.skel.delete    ;
    glob.p.plotSkel.delete;
    glob.p.scatterSkel.delete;
end
glob.p.cell.FaceAlpha = 1;


val = get(hObject,'Value');

if val
    fileName = sprintf('%sskelFV_%d.mat',glob.fvDir,glob.pickCID)
    if exist(fileName,'file')
        load(fileName)
        fv.vertices = fv.vertices(:,[3 1 2]) * 10;
        glob.p.skel = renderFVnav(fv,[1 1 0],.6);
    end
    
    glob.p.cell.FaceAlpha = .25;
    
    if 0 % show skel with edges and nodes
        fileName = sprintf('%snep_%d.mat',glob.fvDir,glob.pickCID)
        if exist(fileName,'file')
            load(fileName)
            
            hold on
            pos = nep.pos*10;
            edges = nep.edges;
            glob.p.scatterSkel = scatter3(pos(:,2),pos(:,1),pos(:,3),15,'o',...
                'w','filled');
            glob.p.plotSkel = plot3([pos(edges(:,1),2) pos(edges(:,2),2)]',...
                [pos(edges(:,1),1) pos(edges(:,2),1)]',...
                [pos(edges(:,1),3) pos(edges(:,2),3)]','color',[1 0 0],'linewidth',1);
            
            
            hold off
        end
    end
    
else
    try        glob.p.cell.FaceAlpha = 1;
    end
    
end


% --- Executes on button press in toggleRadius.
function toggleRadius_Callback(hObject, eventdata, handles)
% hObject    handle to toggleRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleRadius

global glob

if 0
    
    try
        glob.p.radius.delete    ;
    end
    glob.p.cell.FaceAlpha = 1;
    
    
    val = get(hObject,'Value');
    
    if val
        fileName = sprintf('%sskelFV_%d.mat',glob.fvDir,glob.pickCID)
        if exist(fileName,'file')
            load(fileName)
            fv.vertices = fv.vertices(:,[3 1 2]) * 10;
            glob.p.radius = renderFVnav(fv,[1 1 0],.6);
        end
        
        glob.p.cell.FaceAlpha = .25;
        
    else
        try        glob.p.cell.FaceAlpha = 1;
        end
        
    end
    
else
    disp('disabled')
end


% --- Executes on selection change in subType_popupMenue.
function subType_popupMenue_Callback(hObject, eventdata, handles)
% hObject    handle to subType_popupMenue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subType_popupMenue contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subType_popupMenue
global glob tis


val = get(hObject,'Value');


glob.subTypeID = val-1;

if val == 1;
    glob.listCellidx = find(tis.cells.type.typeID == glob.typeID);
else
    glob.listCellidx = find((tis.cells.type.typeID == glob.typeID) & ...
        (tis.cells.type.subTypeID == glob.subTypeID) );
    
end

set(handles.listCIDs,'String',glob.cellStr(glob.listCellidx) )
pickCids = tis.cids(glob.listCellidx);
cellText = sprintf('selected cells  = %s ',num2str(pickCids,'%d '))
showText(cellText);





% --- Executes during object creation, after setting all properties.
function subType_popupMenue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subType_popupMenue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global glob tis

if glob.typeID < 1
    subTypeNames = {'all'};
else
    subTypeNames = tis.cell.type.subTypeNames{glob.typeID};
end

set(hObject,'String',subTypeNames)



% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showSub_toggle.
function showSub_toggle_Callback(hObject, eventdata, handles)
% hObject    handle to showSub_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showSub_toggle


global glob tis

if isfield(glob,'g')
    L = length(glob.g)+1;
else
    L = 1;
end

glob.g(L).idx = glob.listCellidx;
glob.g(L).col = repmat([rand rand rand],[length(glob.g(L).idx) 1]);
glob.g(L).alph = .7;

if glob.typeID>0
    typeName = tis.cells.type.typeNames{glob.typeID};
    
    if glob.subTypeID>0
        subTypeName = tis.cells.type.subTypeNames{glob.typeID}{glob.subTypeID};
    else
        subTypeName = '-';
    end
    
else
    typeName = '-';
    subTypeName = '-';
end


glob.g(L).name = sprintf('%d %s %s',L,typeName,subTypeName);
showCellGroup(L)

groupNames ={glob.g.name};
set(handles.renderGroups_popUp,'String',groupNames)
set(handles.renderGroups_popUp,'Value',length(groupNames))








% --- Executes on selection change in subColor_popUp.
function subColor_popUp_Callback(hObject, eventdata, handles)
% hObject    handle to subColor_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subColor_popUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subColor_popUp

global glob tis
G = get(handles.renderGroups_popUp,'Value');

L = length(glob.g(G).idx);


colOption = get(hObject,'Value');

if colOption == 2
    
    col = hsv(L);
    col = col(randperm(L),:);
    
elseif colOption == 1
    
    col = [rand rand rand];
    col = repmat(col,[L 1]);
    
else
    colIdx = colOption - 2;
    colVals = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; .3 .3 .3];
    colPick = colVals(colIdx,:);
    col = repmat(colPick,[L 1]);
    
end

glob.g(G).col = col;
 showCellGroup(G)




% --- Executes during object creation, after setting all properties.
function subColor_popUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subColor_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

colorOptions = {'rand','rainbow','red','green','blue','yellow','magenta','cyan',...
    'white','grey'};
set(hObject,'String',colorOptions);



% --- Executes on selection change in renderGroups_popUp.
function renderGroups_popUp_Callback(hObject, eventdata, handles)
% hObject    handle to renderGroups_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns renderGroups_popUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from renderGroups_popUp

global glob tis






% --- Executes during object creation, after setting all properties.
function renderGroups_popUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to renderGroups_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addRef_toggleButton.
function addRef_toggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to addRef_toggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addRef_toggleButton



global glob tis

if isfield(glob,'g')
    L = length(glob.g)+1;
else
    L = 1;
end


if isfield(glob,'pickCIDref')
    
    refCids = glob.pickCIDref;
    
    
    idx = refCids * 0;
    for i = 1 : length(idx)
        idx(i) = find(tis.cids==refCids(i),1);
    end
    
    
    
    glob.g(L).idx = idx;
    glob.g(L).col = repmat([rand rand rand],[length(glob.g(L).idx) 1]);
    glob.g(L).alph = .7;
    
    if glob.typeID>0
        typeName = tis.cells.type.typeNames{glob.typeID};
        
        if glob.subTypeID>0
            subTypeName = tis.cells.type.subTypeNames{glob.typeID}{glob.subTypeID};
        else
            subTypeName = '-';
        end
        
    else
        typeName = '-';
        subTypeName = '-';
    end
    
    
    glob.g(L).name = sprintf('%d ref',L);
    showCellGroup(L)
    
    groupNames ={glob.g.name};
    set(handles.renderGroups_popUp,'String',groupNames)
    set(handles.renderGroups_popUp,'Value',length(groupNames))
    set(handles.edit2,'String','')
    edit2_Callback(handles.edit2, eventdata, handles)    
    
   
    
end





% --- Executes on button press in deleteGroup.
function deleteGroup_Callback(hObject, eventdata, handles)
% hObject    handle to deleteGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deleteGroup

global glob tis

G = get(handles.renderGroups_popUp,'Value');

if G>0
    
    for p = 1:length(glob.g(G).patch)
        glob.g(G).patch(p).delete
    end
    
    
    keep = setdiff([1:length(glob.g)],G);
    glob.g = glob.g(keep);
    
    groupNames = {glob.g.name};
    if isempty(groupNames)
        groupNames = {' -no groups- '};
    end
    
    newVal = max(min(G,length(groupNames)),1);
    set(handles.renderGroups_popUp,'Value',newVal);
    set(handles.renderGroups_popUp,'String',groupNames)
    
end



function showCellGroup(V)

if ~exist('V','var')
    V = 1:length(glob.g);
end


global glob tis

for v = V
    
    val = tis.cells.cids(glob.g(v).idx);
    
        
    try
        for r = 1:length(glob.g(v).patch)
            glob.g(v).patch(r).delete
        end
    end
    
    
    for i = 1:length(val)
        fileName = sprintf('%s%d.mat',glob.fvDir,val(i));
        load(fileName);
        glob.g(v).patch(i) = renderFVnav(fv,glob.g(v).col(i,:),glob.g(v).alph,val(i));
        
    end
    
end


% --- Executes on selection change in subAlpha_popUp.
function subAlpha_popUp_Callback(hObject, eventdata, handles)
% hObject    handle to subAlpha_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns subAlpha_popUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subAlpha_popUp

global glob tis
G = get(handles.renderGroups_popUp,'Value');

L = length(glob.g(G).idx);

alphOption = get(hObject,'Value');

alph = (alphOption -1) /11;

glob.g(G).alph = alph;
showCellGroup(G)




% --- Executes during object creation, after setting all properties.
function subAlpha_popUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subAlpha_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

alphaOptions = {'0' '.1' '.2' '.3' '.4' '.5' '.5' '.6' '.7' '.8' '.9' '1'};
set(hObject,'String',alphaOptions);
