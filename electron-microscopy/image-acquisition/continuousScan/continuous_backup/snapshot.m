%%Take fast image of entire stage area
%clear all
TPN = GetMyDir; %% Get directory to place all images
stripWidth = 4.5; % width of each strip in mm

%% Activate ActiveX
if ~exist('sm')
sm = startFibics;
end

%%  Set imaging parameters



%%Fibics variables
requestFOV = stripWidth * 1000;
W = 14000;
H = 14000;
dwell = .1;
rawDir = [TPN 'raw\'];
if ~exist(rawDir),mkdir(rawDir); end
isoDir = [TPN 'isoDir\'];
if ~exist(isoDir),mkdir(isoDir),end

sm.Fibics_WriteFOV(requestFOV);
FOV = sm.Fibics_ReadFOV;

%%Stage variables
Xstart =      110/1000;
Xstop =     11.6591 / 1000;
Ystop =      4.4261/1000;
Ystart =      110/1000;

% Xstart = 40/1000;
% Xstop = 30/1000;
%%Set for bidirectional scanning
Ystop(2) = Ystart(1);
Ystart(2) = Ystop(1);
fixX = [ 0 -50]/1000000
fixX = [0 0];

Xstep = FOV / -1000000;
Xpos = Xstart: Xstep  :Xstop;

startTime = clock;

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);
smwait(sm,'DP_STAGE_IS');
sm.Set_PassedTypeSingle('AP_SCANROTATION',0);

sm.Set_PassedTypeSingle('DP_STAGE_BACKLASH',0);

clear previousWrite
for x = 13:13%length(Xpos)
    
    flip = 2 - mod(x,2);
    sprintf('imaging strip %d of %d',x,length(Xpos))
    WriteTo = [rawDir num2str(x) '.tif'];
    
    %%Go to start position
         
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x)+fixX(flip));
    smwait(sm,'DP_STAGE_IS');
        
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x)+fixX(flip));
    smwait(sm,'DP_STAGE_IS');
    sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
    Ystart(flip);
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystart(flip));
    smwait(sm,'DP_STAGE_IS');
    
    %%Start Acquisition
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    
    %%Start Moving
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystop(flip));
    
    %%Free time
    previousWrite = [rawDir num2str(x-1) '.tif'];
    if exist(previousWrite,'file')
        xSamp = 31;
        I = 255-imread(previousWrite,'PixelRegion',{[1,1,H],[1,xSamp,W]});
        imwrite(I,[isoDir num2str(x-1) '.tif'],'Compression','none');
    end
    tic
    %%Wait for stage
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_STAGE_IS'));
        sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
        pause(.1)
    end
    
    %%Wait for image
    while(sm.Fibics_IsBusy),  pause(.1),  end
    imagedX(x) =  sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
    toc;
end

while ~exist(WriteTo,'file'),pause(.1),end
pause(1)
I = 255-imread(WriteTo,'PixelRegion',{[1,1,H],[1,xSamp,W]});
imwrite(I,[isoDir num2str(x) '.tif'],'Compression','none');

imagedX;
realSteps = imagedX(1:end-1)-imagedX(2:end)-FOV/1000000;
sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000;

stopTime= clock;
duration = stopTime - startTime

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',mean([Xstart Xstop]));
smwait(sm,'DP_STAGE_IS');

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',mean(Ystop));

%% Make image
buildSnapshot(TPN)








