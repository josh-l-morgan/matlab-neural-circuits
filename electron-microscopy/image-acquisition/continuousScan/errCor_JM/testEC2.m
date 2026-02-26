%% Activate ActiveX
sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
sm.InitialiseRemoting
sm.Fibics_Initialise();
%This should give enough time for initialization process
sprintf(' Fibics Initializing, pausing 15 seconds...')
pause(15)

%% Generate fake targets
oldX = sm.Get_ReturnTypeSingle('AP_STAGE_AT_X');
oldY = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
oldWD = sm.Get_ReturnTypeSingle('AP_WD');
oldSX = sm.Get_ReturnTypeSingle('AP_STIG_X');
oldSY = sm.Get_ReturnTypeSingle('AP_STIG_Y');

for t = 1:10
    
    tile(t).imTargX = oldX - t * .00001;    
    tile(t).imTargY = oldY;
    tile(t).focTargX = tile(t).imTargX + [1:10] * .000002; 
    tile(t).focTargY = tile(t).imTargY + [1:10] * .000002;
    tile(t).startWD = oldWD;
    tile(t).startSX = oldSX;
    tile(t).startSY = oldSY;
    
end
retake.tiles = tile;


%% Generate fake imaging parameters
FOV = 4;
W = 1000;
H = 1000;
DwellMicroSec = .5;
ImageFileNameStr = 'C:\MyTestImage_JM.tif';

%% Set auto parameters
focusMag = 25000;



%% Run retakes

tileNum = length(retake.tiles);
for t = 1:tileNum
    
    %%Backlash
    sm.Execute('CMD_STAGE_BACKLASH');
    smwait(sm,'DP_STAGE_IS')
    
    %%Set probable imaging parameters
    param = retake.tiles(t);
    sm.Set_PassedTypeSingle('AP_WD',param.startWD);
    sm.Set_PassedTypeSingle('AP_STIG_X',param.startSX);
    sm.Set_PassedTypeSingle('AP_STIG_Y',param.startSY);
    %sm.Set_PassedTypeSingle('AP_STIG_X',0);  %replace in targeting program
    %sm.Set_PassedTypeSingle('AP_STIG_Y',0);
    
    %%Goto focus position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.focTargX(1));
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',param.focTargY(1));
    smwait(sm,'DP_STAGE_IS');
    
    %%Autofocus
    sm.Set_PassedTypeSingle('AP_MAG',focusMag);
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    smwait(sm,'DP_AUTO_FUNCTION');
    sm.Execute('CMD_AUTO_STIG');
    smwait(sm,'DP_AUTO_FUNCTION');
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    smwait(sm,'DP_AUTO_FUNCTION');
    
    autoWD = sm.Get_ReturnTypeSingle('AP_WD');
    autoSX = sm.Get_ReturnTypeSingle('AP_STIG_X');
    autoSY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
    
    %%Goto final position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.imTargX);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',param.imTargY);
    smwait(sm,'DP_STAGE_IS');
    
    sm.Execute('CMD_STAGE_BACKLASH');
    smwait(sm,'DP_STAGE_IS');
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.imTargX);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',param.imTargY);
    smwait(sm,'DP_STAGE_IS');
    
    
     tic;
    sm.Fibics_WriteFOV(FOV)
    sm.Fibics_AcquireImage(W,H,DwellMicroSec,ImageFileNameStr);
    while(sm.Fibics_IsBusy),  pause(.2),  end

    sprintf('Image took %d seconds',toc)
    
    
    
   sprintf('Finished tile %d of %d',t,tileNum)
end




