

KPN=GetMyDir

Kdir=dir(KPN); Kdir=Kdir(3:size(Kdir,1));

CellName=cell(size(Kdir,1),1);
Status=zeros(size(Kdir,1),3);

for k = 1:size(Kdir,1)
    TPN=([KPN Kdir(k).name '\']);
    N=Kdir(k).name;
    
    CellName(k)={N};
    if exist([TPN 'find\SG.mat'])
        Status(k,1)=1;
    end
    if exist([TPN 'Dots.mat'])
       Status(k,2)=1;
    end
    if exist([TPN 'find\SG.mat'])
        Status(k,3)=1;
    end
    if exist([TPN 'Use.mat'])
        Status(k,4)=1;
    end
    if exist([TPN 'I'])
        Status(k,5)=1;
    end
    if exist([TPN 'images\Combo.tif'])
        Status(k,6)=1;
    end
      
    
end

Incom=sum(Status,2);
Incom=(Incom>0) & (Incom < 5);
Incompletes=CellName(Incom)

