%%Acquire feducials for snapshot

%%Fibics variables
FOV = stripWidth * 1000;
W = 15000;
H = 15000;
dwell = .1;
TPN = GetMyDir;

    WriteTo = [WriteBase 'F2' '.tif'];
    sm.Set_PassedTypeSingle('AP_STAGE_GOTO_R',0);

    
    %%Start Acquisition
    sm.Fibics_WriteFOV(FOV);
    sm.Fibics_AcquireImage(W,H,dwell,WriteTo);
    
    Xpos = sm.Get_ReturnTypeSingle('AP_STAGE_AT_X')*1000;
    Ypos = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y')*1000;
    
    Fdat.Xpos = Xpos;
    Fdat.Ypos = Ypos;
    Fdat.Rot = 0;
    
    F2dat = Fdat;
    save([TPN 'F2dat.mat'],'F2dat')
    