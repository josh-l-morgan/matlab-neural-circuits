function[] = retakeControl(sm,wif,retake)
if ~exist('wif','var')
wif = GetMyWafer;
end
if ~exist('retake','var')
load([wif.dir 'retake.mat'])
end



%% Run retakes
tileNum = length(retake.tiles)
for t = 1:tileNum
    param = retake.tiles(t);
    
    FOV = retake.tiles(t).FOV;
    W = retake.tiles(t).TileWidth;
    H = retake.tiles(t).TileHeight;
    DwellMicroSec = retake.tiles(t).DwellTime;
    
    %%Find focus parameters
    %sm.Set_PassedTypeSingle('AP_PIXEL_SIZE',param.PixelSize);
    focusMag = 10000;
    
    %%Image Name
    retakeDir = [wif.dir(1:end-1) '_retake']
    if ~exist(retakeDir,'dir'),mkdir(retakeDir),end
    nam = retake.tileInfo(t).path;
    slashes = find(nam == '\');
    sectionDir = [retakeDir nam(slashes(end-1):slashes(end))];
    if ~exist(sectionDir,'dir'),mkdir(sectionDir),end
    ImageFileNameStr = [sectionDir retake.tileInfo(t).name]; 
    
    %%Backlash
    while(sm.Fibics_IsBusy),  pause(.2),  end
    sm.Execute('CMD_STAGE_BACKLASH');
    smwait(sm,'DP_STAGE_IS')
    
    %%Set probable imaging parameters
    sm.Set_PassedTypeSingle('AP_WD',param.startWD);
    sm.Set_PassedTypeSingle('AP_STIG_X',param.startSX);
    sm.Set_PassedTypeSingle('AP_STIG_Y',param.startSY);
    sm.Set_PassedTypeSingle('AP_CONTRAST', param.Contrast);
    sm.Set_PassedTypeSingle('AP_BRIGHTNESS',param.Brightness);
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.ScanRot);
    
    %%Goto focus position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.focTargX);
    smwait(sm,'DP_STAGE_IS');
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',param.focTargY-20 * 10 ^-6);
    smwait(sm,'DP_STAGE_IS');
    
%% Autofocus
    %sm.Execute('CMD_AUTO_SATURATION');
    %sm.Execute('CMD_AUTO_FOCUS_COARSE');
    
    sm.Set_PassedTypeSingle('AP_MAG',focusMag);
    magwait(sm,focusMag);
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    pause(.1)
    sm.Execute('CMD_ABORT_AUTO');
    sm.Execute('CMD_AUTO_FOCUS_FINE'); 
    smwait(sm,'DP_AUTO_FUNCTION');
    sm.Execute('CMD_AUTO_STIG');
    smwait(sm,'DP_AUTO_FUNCTION');
    sm.Execute('CMD_AUTO_FOCUS_FINE');
    smwait(sm,'DP_AUTO_FUNCTION');
    
    autoWD = sm.Get_ReturnTypeSingle('AP_WD');
    autoSX = sm.Get_ReturnTypeSingle('AP_STIG_X');
    autoSY = sm.Get_ReturnTypeSingle('AP_STIG_Y');
    
    %% Goto final position
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',param.imTargX);
    smwait(sm,'DP_STAGE_IS')
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
    sm.Fibics_AcquireImage(W,H,1,ImageFileNameStr);
    while(sm.Fibics_IsBusy),  pause(.2),  end

    sprintf('Image took %d seconds',toc)
    
    
    
   sprintf('Finished tile %d of %d',t,tileNum)
end



%%

% 
% sm.Set_PassedTypeSingle('AP_AUTO_CONTRAST',10)
% sm.Get_ReturnTypeSingle('AP_AUTO_CONTRAST')
% sm.Get_ReturnTypeSingle('AP_CONTRAST')
% sm.Set_PassedTypeSingle('AP_CONTRAST',10)
% 



