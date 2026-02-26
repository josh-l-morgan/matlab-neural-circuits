function[watProps] = startWatProps(handles);

global watProps


curDir = pwd;
sourceDir = [curDir '\source\'];
%recordDir = [curDir '\record\'];
recordDir  = 'C:\ATUMvid\record\'

if ~exist(sourceDir,'dir'),mkdir(sourceDir),end
if ~exist(recordDir,'dir'),mkdir(recordDir),end

%% set variables
watProps.intTime = 100;
watProps.fadeRate = 0.01;
watProps.x1 = 100;
watProps.x2 = 200;
watProps.y1 = 100;
watProps.y2 = 400;
watProps.win1 = .2;
watProps.win2 = .4;
watProps.thresh1 = 1;
watProps.contrast = 1;
watProps.brightness = 0;
watProps.sourceDir = sourceDir;
watProps.recordDir = recordDir;
watProps.pumpInterval = 30;
watProps.pumpDuration = 1;
watProps.lastPump = tic;
watProps.manHist = [];
watProps.threshHist = [];
watProps.watching = 0;
watProps.startTime = datenum(datetime);
watProps.soundOn = 1;
watProps.status = 1;
watProps.autoOn = 1;
watProps.iHeight = 512;
watProps.iWidth = 512;
watProps.iCol = 3;
watProps.colorMode = 3;
watProps.record = 0;

iDir = dir(watProps.sourceDir);

watProps.useCam = 1;
watProps.cam = [];


%% derived
watProps.boxHeight = watProps.y2-watProps.y1;
watProps.boxWidth = watProps.x2-watProps.x1;

set(handles.axes_Vid,'xtick',[],'ytick',[]);
set(handles.axes_Plot,'xtick',[],'ytick',[]);
set(handles.axes_History,'ytick',[]);

updateFields(handles)








% 
% %% SEt figure
% 
% set(handles.text_sourceDirectory,'string',watProps.sourceDir);
% set(handles.edit_boxHeight,'string',watProps.boxHeight)
% set(handles.edit_boxWidth,'string',watProps.boxWidth)
% set(handles.edit_integrationTime,'string',watProps.intTime)
% set(handles.slider_win1,'Value',watProps.win1)
% set(handles.slider_win2,'Value',watProps.win2)
% set(handles.slider_upDown,'Value',1-watProps.y1/watProps.iHeight)
% set(handles.slider_leftRight,'Value',watProps.x1/watProps.iWidth)
% set(handles.slider_thresh1,'Value',watProps.thresh1/256);
% set(handles.edit_thresh1,'string',watProps.thresh1)
% set(handles.edit_pumpInterval,'string',watProps.pumpInterval)
% set(handles.edit_pumpDuration,'string',watProps.pumpDuration)



