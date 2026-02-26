function varargout = cellNavGuide_23(varargin)
% CELLNAVGUIDE_23 MATLAB code for cellNavGuide_23.fig
%      CELLNAVGUIDE_23, by itself, creates a new CELLNAVGUIDE_23 or raises the existing
%      singleton*.
%
%      H = CELLNAVGUIDE_23 returns the handle to a new CELLNAVGUIDE_23 or the handle to
%      the existing singleton*.
%
%      CELLNAVGUIDE_23('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLNAVGUIDE_23.M with the given input arguments.
%
%      CELLNAVGUIDE_23('Property','Value',...) creates a new CELLNAVGUIDE_23 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cellNavGuide_23_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cellNavGuide_23_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cellNavGuide_23

% Last Modified by GUIDE v2.5 13-Jun-2020 18:34:53




% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cellNavGuide_23_OpeningFcn, ...
    'gui_OutputFcn',  @cellNavGuide_23_OutputFcn, ...
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

% --- Executes just before cellNavGuide_23 is made visible.
function cellNavGuide_23_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellNavGuide_23 (see VARARGIN)

% Choose default command line output for cellNavGuide_23
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cellNavGuide_23 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global glob tis
glob.handles = handles;
glob.data.path.Parent = glob.handles.figure1;



% --- Outputs from this function are returned to the command line.
function varargout = cellNavGuide_23_OutputFcn(hObject, eventdata, handles)
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
        bipCids = glob.cids(tis.cells.type.typeID == 7);
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

str = get(hObject,'string');
if iscell(str)
    str = str{1};
end
cid = str2num(str);
idx = find(glob.cids==cid,1);

if isempty(idx)
    set(handles.textOut,'String','cell(cid) not found');
else


glob.pickIdx = idx;
glob.pickCID = cid;
showPickedCell
end




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
    glob.subTypeID = 0;
else
    glob.typeID =glob.typeIDs(val);
    glob.listCellidx = find(tis.cells.type.typeID == glob.typeID);
    glob.subTypeID = 0;
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
disp('Showing selected cell')
set(glob.handles.edit1,'String',num2str(glob.pickCID))

pushPostSyn_Callback(glob.handles.pushPostSyn,[],glob.handles);
pushPreSyn_Callback(glob.handles.pushPreSyn,[],glob.handles);
pushPostCells_Callback(glob.handles.pushPostCells,[],glob.handles);
pushPreCell_Callback(glob.handles.pushPreCell,[],glob.handles);

fileName = sprintf('%s%d.mat',glob.fvDir,glob.pickCID);
load(fileName);
%     set(handles.mainAx)
%     patch(handles.mainAx,fv);

try
    glob.p.cell.delete
end
glob.p.cell = renderFVnav(fv,[1 1 1],1,glob.pickCID);
%axis off
pause(.01)

toggleSkeleton_Callback(glob.handles.toggleSkeleton);

toggle_CellData_Callback(glob.handles.toggle_CellData,[],glob.handles) %activate show data for cell





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




% --- Executes on button press in pushPreCell.
function pushPreCell_Callback(hObject, eventdata, handles)
% hObject    handle to pushPreCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global glob tis

val = get(hObject,'Value');

try
    for i = 1:length(glob.p.preCells)
        glob.p.preCells.delete;
    end
end


if val
    
    targ = glob.pickCID;
    
    isSyn = find(tis.syn.post == targ);
    isPre = unique(tis.syn.pre(isSyn));
    isPre = setdiff(isPre,0);
    showCells = isPre;
    
    colOption = get(handles.listbox_preCol,'Value');
    col = getCol(colOption,length(showCells));
    for i = 1:length(showCells)
        fileName = sprintf('%s%d.mat',glob.fvDir,showCells(i));
        if exist(fileName,'file')
            load(fileName);
            glob.p.preCells(i) = renderFVnav(fv,col(i,:),.6)
            glob.p.preCells(i).Tag = num2str(showCells(i));
        end
    end
    
    
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
    for i = 1:length(glob.p.postCells)
        glob.p.postCells(i).delete
    end
end

if val
    
    targ = glob.pickCID;
    
    isSyn = find(tis.syn.pre == targ);
    isPost = unique(tis.syn.post(isSyn));
    isPost = setdiff(isPost,0);
    showCells = isPost;
       
    colOption = get(handles.listbox3,'Value');
    col = getCol(colOption,length(showCells));
    for i = 1:length(showCells)
        fileName = sprintf('%s%d.mat',glob.fvDir,showCells(i));
        if exist(fileName,'file')
            load(fileName);
             glob.p.postCells(i) = renderFVnav(fv,col(i,:),.6)
             glob.p.postCells(i).Tag = num2str(showCells(i));
        end
    end
    
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
    pos = tis.syn.synPosDS(isSyn,[2 1 3]) * glob.em.dsRes(1);
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
    pos = tis.syn.synPosDS(isSyn,[2 1 3])* glob.em.dsRes(1);
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

set(glob.handles.textOut,'String',VASTpos);
tag = get(info.Target,'Tag');

txt = umPosStr;

clipboard('copy',VASTpos)






% --- Executes during object creation, after setting all properties.
function raiseLight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to raiseLight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global glob
glob.light = lightangle(180,45);


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
    subTypeNames = tis.cells.type.subTypeNames{glob.typeID};
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
glob.g(G).col = getCol(colOption,L);
showCellGroup(G)


function[col] = getCol(colOption,L)

if ~exist('L','var')
    L = 1;
end
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

global glob
set(hObject,'String', glob.colorOptions);



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


% --- Executes on button press in toggle_dataOn.
function toggle_dataOn_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_dataOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_dataOn

global glob tisDat

load([glob.fvDir 'tisDat.mat']);
val = get(hObject,'Value');
set(handles.panelData,'Visible',val)


% --- Executes on button press in toggle_CellData.
function toggle_CellData_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_CellData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_CellData

global glob tisDat

val = get(hObject,'Value');

try
    glob.dat.h.cell.delete;
end

if val
    cid = glob.pickCID
    cellIdx = glob.pickIdx;
    cellDat = tisDat.depthHist(cellIdx);
    depth = cellDat.histDat(:,1);
    surface = cellDat.histDat(:,2);
    surface = surface/max(surface);
    
    glob.dat.h.cell = plot(handles.axesDepth,surface,depth,...
        'Color',[0 0 0],'linewidth',2);
    set(handles.axesDepth,'XLim',[0 1],'YLim',[-.5 1.5],'ydir','reverse',...
        'YTick',[-.5:.1:1.5])
end



% --- Executes on button press in toggle_groupData1.
function toggle_groupData1_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_groupData1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_groupData1

global glob tisDat

val = get(hObject,'Value');

if val
    groupID = get(handles.renderGroups_popUp,'Value');
    
    idxs = glob.g(groupID).idx;
    hold on
    for i = 1:length(idxs)
        cellIdx = idxs(i);
        col  = glob.g(groupID).col(i,:);
        cellDat = tisDat.depthHist(cellIdx);
        if cellDat.isDat
            depth = cellDat.histDat(:,1);
            surface = cellDat.histDat(:,2);
            surface = surface/max(surface);
            set(handles.axesDepth, 'NextPlot', 'add')
            glob.dat.h.g1(i) = plot(handles.axesDepth,surface,depth,...
                'Color',col','linewidth',2);
        end
        set(handles.axesDepth,'XLim',[0 1],'YLim',[-.5 1.5],'ydir','reverse',...
            'YTick',[-.5:.1:1.5])
    end
    
else
    try
        for i = 1:length(glob.dat.h.g1)
            glob.dat.h.g1(i).delete
        end
    end
end





% --- Executes on button press in toggle_groupData2.
function toggle_groupData2_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_groupData2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_groupData2



global glob tisDat

val = get(hObject,'Value');

if val
    groupID = get(handles.renderGroups_popUp,'Value');
    
    idxs = glob.g(groupID).idx;
    hold on
    for i = 1:length(idxs)
        cellIdx = idxs(i);
        col  = glob.g(groupID).col(i,:);
        cellDat = tisDat.depthHist(cellIdx);
        if cellDat.isDat
            depth = cellDat.histDat(:,1);
            surface = cellDat.histDat(:,2);
            surface = surface/max(surface);
            set(handles.axesDepth, 'NextPlot', 'add')
            glob.dat.h.g2(i) = plot(handles.axesDepth,surface,depth,...
                'Color',col','linewidth',2);
        end
        set(handles.axesDepth,'XLim',[0 1],'YLim',[-.5 1.5],'ydir','reverse',...
            'YTick',[-.5:.1:1.5])
    end
    
else
    try
        for i = 1:length(glob.dat.h.g1)
            glob.dat.h.g2(i).delete
        end
    end
end




function toggle_dataOn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subAlpha_popUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in toggle_refDat.
function toggle_refDat_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_refDat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_refDat
global glob tisDat

val = get(hObject,'Value');

if val
    
    if isfield(glob,'pickCIDref')
        refCids = glob.pickCIDref;
        
        for i = 1:length(refCids)
            idxs(i) = find(glob.cids==refCids(i),1);
        end
        
        for i = 1:length(idxs)
            cellIdx = idxs(i);
            col  = [1 0 0];
            cellDat = tisDat.depthHist(cellIdx);
            if cellDat.isDat
                depth = cellDat.histDat(:,1);
                surface = cellDat.histDat(:,2);
                surface = surface/max(surface);
                set(handles.axesDepth, 'NextPlot', 'add')
                glob.dat.h.ref(i) = plot(handles.axesDepth,surface,depth,...
                    'Color',col','linewidth',2);
            end
            set(handles.axesDepth,'XLim',[0 1],'YLim',[-.5 1.5],'ydir','reverse',...
                'YTick',[-.5:.1:1.5])
        end
    end
    
else
    try
        for i = 1:length(glob.dat.h.ref)
            glob.dat.h.ref(i).delete
        end
    end
end


% --- Executes on button press in toggle_highlightCell.
function toggle_highlightCell_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_highlightCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_highlightCell

global glob
val = get(hObject,'Value');


if val
    
    glob.dcmH = datacursormode(handles.figure1);
    glob.dcmH.UpdateFcn = @highlightCell;
    
    set(glob.dcmH,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
else
    rotate3d(handles.figure1)
    delete(findall(handles.figure1,'Type','hggroup'))
    try
       set(glob.highlight.oldPatch,...
        'FaceColor', glob.highlight.oldColor,...
        'FaceAlpha',glob.highlight.oldAlpha )
    end
end


function txt = highlightCell(~,info)
global glob
'highlighting'
try
    set(glob.highlight.oldPatch,...
        'FaceColor', glob.highlight.oldColor,...
        'FaceAlpha',glob.highlight.oldAlpha )
end
glob.highlight.oldPatch = info.Target;
glob.highlight.oldAlpha = get(info.Target,'FaceAlpha');
glob.highlight.oldColor =  get(info.Target,'FaceColor');

tag = get(info.Target,'Tag');
set(glob.handles.textOut,'String',['tag = ' tag]);

txt = tag;
set(info.Target,'FaceAlpha',1,'FaceColor',[1 1 .5])


clipboard('copy',tag)


% --- Executes on button press in textCID.
function textCID_Callback(hObject, eventdata, handles)
% hObject    handle to textCID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of textCID

global glob


val = get(hObject,'Value');


% --- Executes on button press in toggle_selectCell.
function toggle_selectCell_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_selectCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_selectCell

global glob
val = get(hObject,'Value');

if val
    
    glob.dcmS = datacursormode(handles.figure1);
    glob.dcmS.UpdateFcn = @cursorSelectCell;
    
    set(glob.dcmS,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
else
    rotate3d(handles.figure1)
    delete(findall(handles.figure1,'Type','hggroup'))
    try
       set(glob.highlight.oldPatch,...
        'FaceColor', glob.highlight.oldColor,...
        'FaceAlpha',glob.highlight.oldAlpha )
    end
end


function txt = cursorSelectCell(~,info)

global glob
'cursor Select'
tag = get(info.Target,'Tag');
set(glob.handles.textOut,'String',['tag = ' tag]);
txt = tag;
clipboard('copy',tag);

cid = str2num(tag);
idx = find(glob.cids==cid,1);
if cid ~= glob.pickCID
    if isempty(idx)
        set(glob.handles.textOut,'String','cell(cid) not found');
    else
        glob.pickIdx = idx;
        glob.pickCID = cid;
        showPickedCell
    end
end
pause(.01)


% --- Executes on selection change in listbox_preCol.
function listbox_preCol_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_preCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_preCol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_preCol

global glob
pushPreCell_Callback(glob.handles.pushPreCell,[],glob.handles);

% --- Executes during object creation, after setting all properties.
function listbox_preCol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_preCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global glob
set(hObject,'String', glob.colorOptions);

% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
global glob
pushPostCells_Callback(glob.handles.pushPostCells,[],glob.handles);

% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global glob
set(hObject,'String', glob.colorOptions);

% --- Executes on button press in pushbutton_preGroup.
function pushbutton_preGroup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_preGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_postGroup.
function pushbutton_postGroup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_postGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton_measurePath.
function togglebutton_measurePath_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_measurePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_measurePath

global glob


val = get(hObject,'Value');

try
    glob.data.path.plotH.delete;
end

if val
    
    glob.data.path.dcm = datacursormode(handles.figure1);
    glob.data.path.dcm.UpdateFcn = @makePath;
    
    set(glob.data.path.dcm,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on');
    plotPath
    
else
    rotate3d(handles.figure1);
    delete(findall(handles.figure1,'Type','hggroup'));
    
end


function txt = makePath(~,info)
global glob
X = info.Position(2);
Y = info.Position(1);
Z = info.Position(3);
%c_info = getCursorInfo(glob.dcm)

newPt = [Y X Z];

if ~isempty(glob.data.path.pts)
    lastPt = glob.data.path.pts(end,:);
else
    lastPt = [0 0 0];
end

if sum(abs(newPt-lastPt))

glob.data.path.pts = cat(1,glob.data.path.pts,[Y X Z]);
glob.data.path.Parent = info.Target.Parent;

anc2sub = (glob.em.res/ 1000)./ glob.em.dsRes;
allPos2 = [X Y Z] ./ anc2sub;
VASTpos = sprintf('%.0f  %.0f  %.0f',allPos2(end,2),allPos2(end,1),allPos2(end,3));
disp(VASTpos)

anc2um =  1./ glob.em.dsRes;
umPos = [X Y Z];
umPosStr = sprintf('%.2f %.2f %.2f um',umPos(end,2),umPos(end,1),umPos(end,3));

tag = get(info.Target,'Tag');


glob.data.path.pos = cat(1,glob.data.path.pos,umPos);

pos = glob.data.path.pos;

L = sqrt((pos(2:end,1)-pos(1:end-1,1)).^2 + (pos(2:end,2)-pos(1:end-1,2)).^2 + ...
    (pos(2:end,3)-pos(1:end-1,3)).^2);
glob.data.path.lengths = L;
txt = sprintf('%0.2f um',sum(L));

clipboard('copy',sum(L))

plotPath
end

function[] =  plotPath()

global glob

pos = glob.data.path.pos;
pts = glob.data.path.pts;
textOutStr = num2str(round(pos));

set(glob.handles.textOut,'String',textOutStr);

try
    glob.data.path.plotH.delete;
end

if ~isempty(pos)
set(glob.data.path.Parent, 'NextPlot', 'add')

glob.data.path.plotH = plot3(glob.data.path.Parent,...
    pts(:,1),pts(:,2),pts(:,3),...
    'linewidth',10,'color','g');
set(glob.data.path.Parent,'Clipping','Off')
end



if get(glob.handles.togglebutton_pathProperties,'Value')
    updatePathPanel
end


% --- Executes on button press in togglebutton_pathProperties.
function togglebutton_pathProperties_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_pathProperties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_pathProperties

val = get(hObject,'Value');

if val
    set(handles.uipanel_pathProperties,'Visible','on')
    updatePathPanel
else
    set(handles.uipanel_pathProperties,'Visible','off')
end


function[] = updatePathPanel()

global glob

h = glob.handles;


pos = glob.data.path.pos;

posStr = num2str(pos,'%0.1f  ');
%posStr2 = mat2strFormat(glob.data.path.pos,1);
set(h.edit_pathPanelPoints,'String',posStr)

lengthsStr = num2str(glob.data.path.lengths,'%0.1f  ');
set(h.edit_pathPanelLengths,'String',lengthsStr)

lengthStr = compose('%0.1f',sum(glob.data.path.lengths));
set(h.edit5,'String',lengthStr)

area = 0;
convexArea = 0;
convexVolume = 0;

if ~isempty(pos)
    area = polyarea(pos(:,1),pos(:,2));
end
if size(pos,1)>2
    [hull2D convexArea] = convhull(pos(:,1),pos(:,2));
end
if size(pos,1)>3
    [hull3D convexVolume] = convhull(pos(:,1),pos(:,2),pos(:,3));   
end

areaStr = compose('%0.1f',area);
set(h.edit_pathPanelXYArea,'String',areaStr)

convAreaStr = compose('%0.1f',convexArea);
set(h.edit_pathPanel_ConvexHullArea,'String',convAreaStr)

volStr = compose('%0.1f',sum(glob.data.path.lengths));
set(h.edit_pathPanelVolume,'String',volStr)


function[Mstr] =  mat2strFormat(M,n)

maxL = length(num2str(max(round(M(:)))))+n+3;
ys = size(M,1);
xs = size(M,2) * maxL-2;
Mstr(ys,xs) = ' ';
for y = 1:size(M,1)
    for x = 1:size(M,2)
        str1 = sprintf('%%0.%df',n);
        str2 = sprintf(str1,M(y,x));
        Mstr(y,maxL*x-length(str2)+1-2:maxL*x-2) = str2;
    end
end




function edit_pathPanelPoints_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathPanelPoints as text
%        str2double(get(hObject,'String')) returns contents of edit_pathPanelPoints as a double


% --- Executes during object creation, after setting all properties.
function edit_pathPanelPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathPanelLengths_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelLengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathPanelLengths as text
%        str2double(get(hObject,'String')) returns contents of edit_pathPanelLengths as a double


% --- Executes during object creation, after setting all properties.
function edit_pathPanelLengths_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelLengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathPanelXYArea_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelXYArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathPanelXYArea as text
%        str2double(get(hObject,'String')) returns contents of edit_pathPanelXYArea as a double


% --- Executes during object creation, after setting all properties.
function edit_pathPanelXYArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelXYArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathPanelVolume_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathPanelVolume as text
%        str2double(get(hObject,'String')) returns contents of edit_pathPanelVolume as a double


% --- Executes during object creation, after setting all properties.
function edit_pathPanelVolume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathPanelVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathPanel_ConvexHullArea_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathPanel_ConvexHullArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathPanel_ConvexHullArea as text
%        str2double(get(hObject,'String')) returns contents of edit_pathPanel_ConvexHullArea as a double


% --- Executes during object creation, after setting all properties.
function edit_pathPanel_ConvexHullArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathPanel_ConvexHullArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_undoPath.
function pushbutton_undoPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_undoPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob

L = length(glob.data.path.pos);

if L>1
    glob.data.path.pos = glob.data.path.pos(1:end-1,:);
    glob.data.path.pts = glob.data.path.pts(1:end-1,:);
    
    pos = glob.data.path.pos;
    
    L = sqrt((pos(2:end,1)-pos(1:end-1,1)).^2 + (pos(2:end,2)-pos(1:end-1,2)).^2 + ...
        (pos(2:end,3)-pos(1:end-1,3)).^2);
    glob.data.path.lengths = L;
    
    pts = glob.data.path.pts;
    textOutStr = num2str(round(pos));
    
    try
        glob.data.path.plotH.delete;
    end
    set(glob.data.path.Parent, 'NextPlot', 'add')
    
    glob.data.path.plotH = plot3(glob.data.path.Parent,...
        pts(:,1),pts(:,2),pts(:,3),...
        'linewidth',10,'color','g');
    
    clipboard('copy',sum(L))
    txt = sprintf('%0.2f um',sum(L));
    
    
end

    updatePathPanel
plotPath
    glob.data.path

% --- Executes on button press in pushbutton_clearPath.
function pushbutton_clearPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob

glob.data.path.pts = [];
glob.data.path.pos = [];
glob.data.path.lengths = [];
delete(findall(handles.figure1,'Type','hggroup'));
updatePathPanel
try
    glob.data.path.plotH.delete;
end


% --- Executes on selection change in popupmenu_funcTypes.
function popupmenu_funcTypes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_funcTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_funcTypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_funcTypes

global glob

val = get(hObject,'Value');
funcIDs = glob.com.typeFunctions{val};
if val == 1
    funcStrs = glob.com.functionFiles;
else
    
    if ~isempty(funcIDs)
        funcStrs = glob.com.functions{funcIDs};
    else
        funcStrs = ' ';
    end
end
set(handles.popupmenu_functions,'String',funcStrs,'Value',1);


% --- Executes during object creation, after setting all properties.
function popupmenu_funcTypes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_funcTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global glob
set(hObject,'String',glob.com.typeStrings)
set(hObject,'Value',1);

% --- Executes on selection change in popupmenu_functions.
function popupmenu_functions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_functions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_functions

% 
% str = get(hObject,'String');
% str = str{get(hObject,'Value')};
% set(handles.edit_eval,'String',str)


% --- Executes during object creation, after setting all properties.
function popupmenu_functions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global glob
set(hObject,'String',glob.com.functionFiles);




function edit_eval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eval as text
%        str2double(get(hObject,'String')) returns contents of edit_eval as a double

'hi'

% --- Executes during object creation, after setting all properties.
function edit_eval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_result_Callback(hObject, eventdata, handles)
% hObject    handle to edit_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_result as text
%        str2double(get(hObject,'String')) returns contents of edit_result as a double


% --- Executes during object creation, after setting all properties.
function edit_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Function.
function pushbutton_Function_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Function (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob
str = get(handles.popupmenu_functions,'String');
if iscell(str)
    str = str{get(handles.popupmenu_functions,'Value')};
end

glob.com.result = '';
try
    output = eval(str);
catch err
    output = err.message 
end
set(handles.edit_result,'String',output)
pause(.5)
if length(glob.com.result)>0
    set(handles.edit_result,'String',glob.com.result)
end

% --- Executes on button press in pushbutton_eval.
function pushbutton_eval_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_eval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global glob
str = get(handles.edit_eval,'String')

glob.com.result = '';
try
    output = eval(str);
catch err
    output = err.message 
end
set(handles.edit_result,'String',output)
pause(.5)
if length(glob.com.result)>0
    set(handles.edit_result,'String',glob.com.result)
end


% --- Executes on button press in togglebutton_showComPanel.
function togglebutton_showComPanel_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_showComPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_showComPanel

val = get(hObject,'Value');
if val
    set(handles.uipanel_Com ,'Visible','On')
else
    set(handles.uipanel_Com , 'Visible','Off')
end

% --- Executes on button press in togglebutton_ShowPathPanel.
function togglebutton_ShowPathPanel_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_ShowPathPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_ShowPathPanel


val = get(hObject,'Value');
if val
    set(handles.uipanel_pathProperties ,'Visible','On')
else
    set(handles.uipanel_pathProperties , 'Visible','Off')
end
