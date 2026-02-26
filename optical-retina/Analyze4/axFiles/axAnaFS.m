function[]=axAnaFS(TPN,DPN)

%% Check outputs for shifting

    clear ShiftInfo 
    DPNd=[TPN 'data\'];
    
    %get image size
    if exist(DPN) 
        Idir=dir(DPN); Idir=Idir(3:size(Idir,1));
        IdirSize=size(Idir,1);
        ShiftInfo.IdirSize=IdirSize;
    else ShiftInfo.IdirSize='none'; end
    
    %Get skeleton source size
    if exist([TPN 'pics\Irm']) 
        Id=dir([TPN 'pics\Irm']);
        ShiftInfo.IrmSize=size(Id,1)-2;
    elseif exist([TPN  'pics\Irm.mat'])
        load([TPN 'pics\Irm.mat'])
        ShiftInfo.IrmSize=size(Irm,3);
        clear Irm
    elseif exist([TPN 'pics\D.mat'])
        load([TPN 'pics\D.mat'])
        ShiftInfo.IrmSize=size(D,3);
        clear D
    elseif exist([TPN 'data\D.mat'])
        load([TPN 'data\D.mat'])
        ShiftInfo.IrmSize=size(D,3);
        clear D
    else ShiftInfo.IrmSize='none';end
    
    %Get Dot find source size
    if exist([TPN 'pics\BigFilled']) 
        Id=dir([TPN 'pics\BigFilled']);
        ShiftInfo.BFSize=size(Id,1)-2;
    elseif exist([TPN  'pics\BigFilled.mat'])
        load([TPN 'pics\BigFilled.mat'])
        ShiftInfo.BFSize=size(BigFilled,3);
        clear BigFilled
    else ShiftInfo.BFSize='none';end   
 
    ShiftInfo.ShiftDots=ShiftInfo.IdirSize-ShiftInfo.BFSize;
    ShiftInfo.ShiftDend=ShiftInfo.IdirSize-ShiftInfo.IrmSize;
    
    ShiftInfo
    save([TPN 'ShiftInfo.mat'],'ShiftInfo')  
    
end


       


