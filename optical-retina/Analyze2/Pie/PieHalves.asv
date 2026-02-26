%%Divide Arbors into pie regions radiating from CB center to look for
%%correlation in P/D, D/A

clear all
load('UseCells.mat')
for k = 1:size(UseCells,1)
    
TPN = char(UseCells(k))

load([TPN 'Use.mat'])
load([TPN 'CA.mat'])

Territory=CA.Arbor(1).Territory;
CB=Use.Cent(1:2);
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

IDist=15;
ODist=1000;
PieSize=pi/2;
c=0;
SampTer=Territory*0;

clear Dt Dd T
for p = 0 : pi/100 : pi
    c=c+1;
    pLow=p-pi;
    pHigh=pLow+pi;
    
    SampDt1 = NNRad>=pLow & NNRad < pHigh & NNDist > IDist ;
    SampDt2 = ~SampDt1 & NNDist > IDist;
    UseDd = DdRad>=pLow & DdRad < pHigh & DdDist > IDist ;
    SampDd1 = Use.Length(UseDd & DdDist > IDist);
    SampDd2 = Use.Length(~UseDd & DdDist > IDist);
    SampT1 = TRad>=pLow & TRad < pHigh & TDist > IDist ;
    SampT2 = ~SampT1 & TDist > IDist; 
    
    Dt(c,:)=[sum(SampDt1)  sum(SampDt2)];
    Dd(c,:)=[sum(SampDd1) sum(SampDd2)];
    T(c,:)=[sum(SampT1) sum(SampT2)];
    
  
end


PA=Dt./T;
DA=Dd./T;
PD=Dt./Dd;

PAd=(PA(:,1)-PA(:,2))./(PA(:,1)+PA(:,2));
DAd=(DA(:,1)-DA(:,2))./(DA(:,1)+DA(:,2));
PDd=(PD(:,1)-PD(:,2))./(PD(:,1)+PD(:,2));


Halves.PA=PA;
Halves.DA=DA;
Halves.PD=PD;
Halves.PAd=PAd;
Halves.DAd=DAd;
Halves.PDd=PDd;

%%  Compare
scatter(PDd,PAd)
scatter(DAd,PDd)
scatter(DAd,PDd),pause(.01)
[R P]=corrcoef(DAd,PDd)

%% Check Peak

Pick=find(abs(DAd)==max(abs(DAd)),1);
Halves.Maxes.maxDAd=DAd(Pick);
Halves.Maxes.DAd2PDd=PDd(Pick);
Halves.Maxes.DAd2PAd=PAd(Pick);


Halves.Maxes
save([TPN 'Halves.mat'],'Halves')

end


%% Look at all cells
clear all
load('UseCells.mat')
for k = 1:size(UseCells,1)

    TPN = char(UseCells(k))
    load([TPN 'Halves.mat'])
    maxDAd(k)=Halves.Maxes.maxDAd;
    DAd2PDd(k)=Halves.Maxes.DAd2PDd;
    
end

scatter(maxDAd,DAd2PDd)
[R P]=corrcoef(maxDAd,DAd2PDd)
