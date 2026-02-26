%%Take fast image of entire stage area
clear all
stripWidth = 4; % width of each strip in mm
TPN = GetMyDir; %% Get directory to place all images

%% Activate ActiveX

if ~exist('FibicsOn','var') %activate fibics
    sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
    sm.InitialiseRemoting
    sm.Set_PassedTypeSingle('AP_MAG',25);
    sm.Fibics_Initialise();
    sprintf(' Fibics Initializing, pausing 15 seconds...')
    pause(15)
    FibicsOn = 1;
end

%%  Set imaging parameters

%%Fibics variables
FOV = stripWidth * 1000;
W = 14000;
H = 14000;
dwell = .1;
rawDir = [TPN 'aspect\'];
if ~exist(rawDir),mkdir(rawDir); end

%%Stage variables
Xstart =      109.9999/1000;
Xstop =     11.6591 / 1000;
Ystop =      4.4261/1000;
Ystart =      109.9999/1000;

%%Set for bidirectional scanning
Ystop(2) = Ystart(1);
Ystart(2) = Ystop(1);

Xstep = FOV / -1000000;
Xpos = ones(1,4) * mean([Xstart Xstop]);


startTime = clock;

sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);


for x = 1:length(Xpos)
    
    flip = 2 - mod(x,2);
    sprintf('imaging strip %d of %d',x,length(Xpos))
    WriteTo = [rawDir num2str(x) '.tif'];
    
    %%Go to start position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')
    Ystart(flip)
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystart(flip));
    smwait(sm,'DP_STAGE_IS');
    
    %%Start Acquisition
    sm.Fibics_WriteFOV(FOV);
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    
    %%Start Moving
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystop(flip));
    
    %%Wait for stage
    while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_STAGE_IS'));
        sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
        pause(.1)
    end
    
    %%Wait for image
    while(sm.Fibics_IsBusy),  pause(.1),  end
    
end

sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000

stopTime= clock;
duration = stopTime - startTime







