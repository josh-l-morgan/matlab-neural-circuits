%%Divide Arbors into pie regions radiating from CB center to look for
%%correlation in P/D, D/A

clear all
TPN = GetMyDir

load([TPN 'Use.mat'])
load([TPN 'CA.mat'])

Territory=CA.Arbor(1).Territory;
CB=Use.Cent(1:2)
DPos=Use.DPos;
Mids=Use.Mids;
NN=Use.NN;


%% Find Rads for Dots and Midpoints
DtRad=zeros(size(DPos,1),1);
DtDist=DtRad;
for i=1:size(DPos,1)
   DtRad(i)=atan2((DPos(i,1)-CB(1)),(DPos(i,2)-CB(2))); 
   DtDist(i)=sqrt((DPos(i,1)-CB(1)).^2+(DPos(i,2)-CB(2)).^2);
end

NNRad=zeros(size(NN,1),1);
NNDist=NNRad;
for i=1:size(NN,1)
    NNRad(i)=atan2((NN(i,1)-CB(1)),(NN(i,2)-CB(2)));
    NNDist(i)=sqrt((NN(i,1)-CB(1)).^2 + (NN(i,2)-CB(2)).^2);
end

DdRad=zeros(size(Mids,1),1);
DdDist=DdRad;
for i=1:size(Mids,1)
    DdRad(i)=atan2((Mids(i,1)-CB(1)),(Mids(i,2)-CB(2)));
    DdDist(i)=sqrt((Mids(i,1)-CB(1)).^2+(Mids(i,2)-CB(2)).^2);
end

[Ty Tx] = find(Territory);
TRad=zeros(size(Ty,1),1);
TDist=TRad;
for i=1:size(Ty,1);
    TRad(i)=atan2((Ty(i)-CB(1)),(Tx(i)-CB(2)));
    TDist(i)=sqrt((Ty(i)-CB(1)).^2 + (Tx(i)-CB(2)).^2);
end    

%% Sample

IDist=10;
ODist=1000;
PieSize=pi/2;

for p = 1 :4
    pLow=(p-1)*(pi/2)-pi;
    pHigh=pLow+PieSize;
    
    SampDt = NNRad>=pLow & NNRad < pHigh & NNDist > IDist & NNDist < ODist;
    UseDd = DdRad>=pLow & DdRad < pHigh & DdDist > IDist & DdDist < ODist;
    SampDd = Use.Length(UseDd);
    SampT = TRad>=pLow & TRad < pHigh & TDist > IDist & TDist < ODist;
    
    Dt(p)=sum(SampDt);
    Dd(p)=sum(SampDd);
    T(p)=sum(SampT);
    
    PA(p)=sum(SampDt)/sum(SampT);
    DA(p)=sum(SampDd)/sum(SampT);
    PD(p)=sum(SampDt)/sum(SampDd);
    
    Quarters.Dt=Dt;
    Quarters.Dd=Dd;
    Quarters.T=T;
    Quarters.PA=PA;
    Quarters.DA=DA;
    Quarters.PD=PD;
    
end

%%  Compare

scatter(DA,PD)

save([TPN 'Quarters.mat'],'Quarters')








