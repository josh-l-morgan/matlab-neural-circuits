
load('UseCells.mat')

for k = 1:size(UseCells,1)
    k
    TPN=char(UseCells(k));
    load([TPN 'Use.mat'])
    DPos=Use.DPos;
    
    DistV10=zeros(size(DPos,1),10);
    DistA10=zeros(size(DPos,1),10);
    clear Near
    
    for d = 1:size(DPos,1)
        
       DistV=dist(DPos,DPos(d,:));
       DistA=sqrt((DPos(:,1)-DPos(d,1)).^2 + (DPos(:,2)-DPos(d,2)).^2);
       
       DistV=sort(DistV);
       DistV10(d,:)=DistV(2:11);
       
       DistA=sort(DistA);
       DistA10(d,:)=DistA(2:11);
        
        
    end
    
    Near.medA=median(DistA10(:,1));
    Near.medV=median(DistV10(:,1));
    save([TPN 'Near.mat'],'Near')
    NearAllA(k,1)=Near.medA;
    NearAllV(k,1)=Near.medV;
    
    
end