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
xSamp = 31; % Down Sample Xdimension
dwell = .1;
rawDir = [TPN 'raw\'];
if ~exist(rawDir),mkdir(rawDir); end
isoDir = [TPN 'isoDir\'];
if ~exist(isoDir),mkdir(isoDir),end

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

sm.Set_PassedTypeSingle('DP_STAGE_BACKLASH',0);

clear previousWrite
picNum = 0;
for y = 1:length(Ypos)
for x = 1:length(Xpos)
    picNum = picNum +1;
    rc(picNum,:) = [y x];
    sprintf('imaging row %d col %d  of %d rows and %d cols',y, x ,length(Ypos) ,length(Xpos))
    WriteTo = [rawDir num2str(picNum) '.tif'];
    
    %%Go to start position
         
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ypos(y));
    smwait(sm,'DP_STAGE_IS');
    
    %%Start Acquisition
    fibSet(picNum) = getSettings(sm);
    sm.Fibics_AcquireImage(BEW,BEH,BEdwell,WriteTo);
    
    while(sm.Fibics_IsBusy),  pause(.1),  end
    

end % Xpos
end % Ypos

stopTime= clock;
duration = stopTime - startTime


save([TPN 'fibSettings.mat'],'fibSet')








