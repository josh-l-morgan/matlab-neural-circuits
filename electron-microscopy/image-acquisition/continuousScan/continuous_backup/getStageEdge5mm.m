%%Take fast image of entire stage area
%clear all
TPN = GetMyDir; %% Get directory to place all images
stripWidth =5; % width of each strip in mm

%% Activate ActiveX
if ~exist('sm')
    sm = startFibics;
end

%%  Set imaging parameters




%%Fibics variables
requestFOV = stripWidth * 1000;
W = 14000;
H = 14000;
xSamp = 20; % Down Sample Xdimension
dwell = .1;
rawDir = [TPN 'raw\'];
if ~exist(rawDir),mkdir(rawDir); end
isoDir = [TPN 'isoDir\'];
if ~exist(isoDir),mkdir(isoDir),end

sm.Fibics_WriteFOV(requestFOV);
FOV = sm.Fibics_ReadFOV;


%%Book end variables
BEdwell = 1;
BEW = length([1:xSamp:W]);
BEH = BEW;

%%Stage variables
Xstart =      110/1000;
Xstop =     0 / 1000;
Ystop =      0/1000;
Ystart =      110/1000;

% Xstart = 40/1000;
% Xstop = 30/1000;
%%Set for bidirectional scanning
fixX = [ 0 -50]/1000000
fixX = [0 0];

Xstep = FOV / -1000000;
Xpos = Xstart: Xstep  :Xstop;
Ypos = Ystart: Xstep : Ystop;

startTime = clock;

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);
smwait(sm,'DP_STAGE_IS');
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);

sm.Set_PassedTypeSingle('DP_STAGE_BACKLASH',1);

clear previousWrite
picNum = 0;
for x = 1:length(Xpos)
    picNum = picNum +1;
    sprintf('imaging col %d  of %d cols', x  ,length(Xpos))
    WriteTo = [rawDir num2str(picNum) '.tif'];
    
    %%Go to start position
    flip = 1 - mod(x,2);
    
    %%Get Bottom
    if flip
    WriteTo = [rawDir 'bottom' num2str(picNum) '.tif'];
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    input('move to bottom of stage then press return')
    tempY = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',tempY);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
   
    
    %%Start Acquisition
    stageEdge.bottom(picNum) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,WriteTo);
    while(sm.Fibics_IsBusy),  pause(.1),  end
    end
     
    
    %%Get Top
    WriteTo = [rawDir 'top' num2str(picNum) '.tif'];
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    input('move to top of stage then press return')
    tempY = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',tempY);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    
    smwait(sm,'DP_STAGE_IS');
    %%Start Acquisition
    stageEdge.top(picNum) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,WriteTo);
    while(sm.Fibics_IsBusy),  pause(.1),  end
    
        %%Get Bottom
    if ~flip
    WriteTo = [rawDir 'bottom' num2str(picNum) '.tif'];
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    input('move to bottom of stage then press return')
    tempY = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',tempY);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
   
    
    %%Start Acquisition
    stageEdge.bottom(picNum) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,WriteTo);
    while(sm.Fibics_IsBusy),  pause(.1),  end
    end
end % Xpos

stopTime= clock;
duration = stopTime - startTime


save([TPN 'stageEdge5.mat'],'stageEdge')








