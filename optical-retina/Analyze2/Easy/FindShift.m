%% There are still a number of files where output is shifted by two 
%%planes


%Get directory name
KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    clear ShiftInfo 
    TPN = [KPN  Kdir(k).name '\'];
    DPN = [TPN 'I\']   
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

%% Gather infor

Head=['IdirSize ' 'IrmSize ' 'BFSize ' 'DotsMin ' 'AllSegMin ' 'AllSegCutMin ']
AllShift=cell(k,7);
for k = 1:size(Kdir,1)
    clear ShiftInfo 
    TPN = [KPN  Kdir(k).name '\'];
    DPN = [TPN 'I\']   ;
    DPNd=[TPN 'data\'];
    
   AllShift(k,1)=cellstr(Kdir(k).name);
   if exist([TPN 'ShiftInfo.mat'])
       load([TPN 'ShiftInfo.mat'])
       
       if isfield(ShiftInfo,'IdirSize')
           AllShift(k,2)=num2cell(ShiftInfo.IdirSize,2);
       end    
       if isfield(ShiftInfo,'IrmSize')
           AllShift(k,3)=num2cell(ShiftInfo.IrmSize,2);
       end       
       if isfield(ShiftInfo,'BFSize')
           AllShift(k,4)=num2cell(ShiftInfo.BFSize,2);
       end       
       if isfield(ShiftInfo,'DotsMin')
           AllShift(k,5)=num2cell(ShiftInfo.DotsMin,2);
       end       
       if isfield(ShiftInfo,'AllSegMin')
           AllShift(k,6)=num2cell(ShiftInfo.AllSegMin,2);
       end
       if isfield(ShiftInfo,'AllSegCutMin')
           AllShift(k,7)=num2cell(ShiftInfo.AllSegCutMin,2);
       end             
   end %if shiftinfo available
       
end

AllShift


