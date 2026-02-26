%%Take fast image of entire stage area
%clear all
TPN = GetMyDir; %% Get directory to place all images

load('stageEdge.mat')
stripWidth = stageEdge.top(1).FOV;
WD = stageEdge.top(1).WD;

%% Activate ActiveX
if ~exist('sm')
sm = startFibics;
end

%%  Set imaging parameters

%%Fibics variables
requestFOV = stripWidth * 1000;
W = 14000;
H = 14000;
xSamp = 31; % Down Sample Xdimension
dwell = .1;
rawDir = [TPN 'raw\'];
if ~exist(rawDir),mkdir(rawDir); end

isoDir = [TPN 'isoDir\'];
if ~exist(isoDir),mkdir(isoDir),end

BEDir = [TPN 'BEDir\'];
if ~exist(BEDir),mkdir(BEDir),end


sm.Fibics_WriteFOV(requestFOV);
FOV = sm.Fibics_ReadFOV;


%%Book end variables
BEdwell = 1;
BEW = fix(W/xSamp)+1;
BEH = BEW;

%%Stage variables
Xstart =      110/1000;
Xstop =     0 / 1000;
Ystop =      0/1000;
Ystart =      110/1000;
bufY = [-.004 .004];

% Xstart = 40/1000;
% Xstop = 30/1000;
%%Set for bidirectional scanning
fixX = [ 0 -50]/1000000
fixX = [0 0];

startTime = clock;

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);
smwait(sm,'DP_STAGE_IS');
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);
sm.Set_PassedTypeSingle('AP_WD',WD);

sm.Set_PassedTypeSingle('DP_STAGE_BACKLASH',0);

clear previousWrite
stripNum = length(stageEdge.top)-2;
for x = 1:stripNum
    
    flip = 2 - mod(x,2);
    sprintf('imaging strip %d of %d',x,stripNum)
    WriteTo = [rawDir 'strip' num2str(x) '.tif'];
    
    startStop(flip) = stageEdge.top(x);
    startStop(3-flip) = stageEdge.bottom(x);
    
    %%Go to start position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',startStop(1).X);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',startStop(1).Y);
    smwait(sm,'DP_STAGE_IS');
    fibSet.BEstart(x) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,[BEDir 'BEstart' num2str(x) '.tif']);
    while(sm.Fibics_IsBusy),  pause(.1),  end
    
    %%Set startY
    startY = (startStop(1).Y + bufY(flip));
    startY = min(startY, Ystart); startY = max(startY, Ystop);
    
    stopY = (startStop(2).Y + bufY(3-flip));
    stopY = min(stopY, Ystart); stopY = max(stopY, Ystop);
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',startY);
    smwait(sm,'DP_STAGE_IS');

    fibSet.strips(x) = getSettings(sm);
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    
    %%Start Moving
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',stopY);
    
    %% Free time
    previousWrite = [rawDir 'strip' num2str(x-1) '.tif'];
    if exist(previousWrite,'file')
        I = imread(previousWrite,'PixelRegion',{[1,1,H],[1,xSamp,W]});
        imwrite(I,[isoDir num2str(x-1) '.tif'],'Compression','none');
    end
    tic
    %% Wait for stage
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_STAGE_IS'));
        sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
        pause(.1)
    end
    
    %% Wait for image
    while(sm.Fibics_IsBusy),  pause(.1),  end
    imagedX(x) =  sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
    toc;
        
    %% BookEnd stop
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',startStop(2).X);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',startStop(2).Y);
    smwait(sm,'DP_STAGE_IS');
    fibSet.BEstop(x) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,[BEDir 'BEstop' num2str(x) '.tif']);
    while(sm.Fibics_IsBusy),  pause(.1),  end

end

while ~exist(WriteTo,'file'),pause(.1),end
pause(5)
I = imread(WriteTo,'PixelRegion',{[1,1,H],[1,xSamp,W]});
imwrite(I,[isoDir num2str(x) '.tif'],'Compression','none');

imagedX;
realSteps = imagedX(1:end-1)-imagedX(2:end)-FOV/1000000;
sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000;

stopTime= clock;
duration = stopTime - startTime

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',mean([Xstart Xstop]));
smwait(sm,'DP_STAGE_IS');

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',mean(Ystop));

save([TPN 'fibSettings.mat'],'fibSet')
%% Make image
buildSnapshot(TPN)








