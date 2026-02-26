

KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

UseCells={}
for k = 1:size(Kdir,1)
    TPN = [KPN Kdir(k).name '\'];
    DPN = [TPN 'I\']  
    
    if exist([TPN 'BranchS.mat'])
       UseCells=cat(1,UseCells , TPN) ;
    end
    
end

