%% Change MAG

TPN = GetMyDir;
if ~exist('FibicsOn','var')
    
    %% Activate ActiveX
    sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
    sm.InitialiseRemoting
    sm.Set_PassedTypeSingle('AP_MAG',25);
    sm.Fibics_Initialise();
    sprintf(' Fibics Initializing, pausing 15 seconds...')
    pause(15)
    FibicsOn = 1;
end
%%
startTime = clock;
Xstart =    50/1000;%  109.9999/1000;
Xstop =    50/1000; %11.6591 / 1000;

FOV = 2000;
W = 2000;
H = 2000;
dwell = 2;
Reps = 3;
WriteBase = [TPN 'test'];
Ystart =      4.4261/1000;%sm.Get_ReturnTypeSingle('AP_STAGE_GOTO_Y')
Ystop =      109.9999/1000;
Xstep = FOV / -1000000;
Xpos = Xstart: Xstep :Xstop;

rowID = 'abcdefghijklmnopqrstuvwxyz';



sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);

for x = 1:length(Xpos)
    
    WriteTo = [WriteBase num2str(x) rowID(1) '.tif'];
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',Xpos(x));
    smwait(sm,'DP_STAGE_IS');
    
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystart);
    smwait(sm,'DP_STAGE_IS');
    'hi'
    sm.Fibics_WriteFOV(FOV)
    'start image'
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    'image going'
    pause(.01)
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',Ystop);
    
    for r = 2:Reps
        while(sm.Fibics_IsBusy),  pause(.2),  end
        %Yend= sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000
        sprintf('imaged strip %d of %d',x,length(Xpos))
    
    pause(.01),'pause to write',r
    WriteTo = [WriteBase num2str(x) rowID(r) '.tif'];
    
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    end
end

while ~strcmp('Idle',sm.Get_ReturnTypeString('DP_STAGE_IS'));
    sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
    pause(.01)
end


sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000

stopTime= clock;

duration = stopTime - startTime
