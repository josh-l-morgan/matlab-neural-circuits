%function[] = anaCrit(TPN, DPN)

%% Get File
TPN=GetMyDir
load([TPN 'Dots.mat'])

%% Gather Results
rVol=Dots.Vol;
rCut=Dots.Cut(1:Dots.Num);
rSum=Dots.ItSum;
rDFOf=Dots.DFOf;
rDF=Dots.DF;
rDM=Dots.DistToMask;

%%Default Thresholds
Vol=4;
Cut=2;
ItSum=4;
DFOf=1.5;
DF=5;
DistToMask=1;

%%Apply Constant Threshols
pCut=rCut<Cut;
pDM=rDM>=DistToMask;
eLim=pCut .* pDM;
rVol=rVol .* eLim;
rSum =rSum .* eLim;
rDFOf=rDFOf .* eLim;
rDF=rDF .* eLim;
eLimSum=sum(eLim);

%% Run Thresholds

%%Run Volumes
clear AllSum
pVol=rVol>=Vol;
pSum=rSum>=ItSum;
pDFOf=rDFOf>=DFOf;
pDF=rDF>=DF;

for v=1:30
   pVol=rVol>=v;
   SumVol(v)=sum(pVol); 
   AllSum(v)=sum(pVol & pSum & pDFOf & pDF);
   
end
subplot(4,1,1)
plot(SumVol)
plot(AllSum)


%%Run Volumes
clear AllSum
pVol=rVol>=Vol;
pSum=rSum>=ItSum;
pDFOf=rDFOf>=DFOf;
pDF=rDF>=DF;
for s=1:100
   pSum=rSum>=s;
   SumSum(s)=sum(pSum);
   AllSum(s)=sum(pVol & pSum & pDFOf & pDF);
end
subplot(4,1,2)
plot(SumSum)
plot(AllSum)

%%Run Volumes
clear AllSum
pVol=rVol>=Vol;
pSum=rSum>=ItSum;
pDFOf=rDFOf>=DFOf;
pDF=rDF>=DF;
for df = 1:50
    pDF=rDF>=df;
    SumDF(df)=sum(pDF);
    AllSum(df)=sum(pVol & pSum & pDFOf & pDF);
end
subplot(4,1,3)
plot(SumDF)
plot(AllSum)

%%Run Volumes
clear AllSum
pVol=rVol>=Vol;
pSum=rSum>=ItSum;
pDFOf=rDFOf>=DFOf;
pDF=rDF>=DF;
for df = 1:250
    RecDF(df)=df/10;
    pDFOf=rDFOf>=df/10;
    SumDFOf(df)=sum(pDFOf);   
    AllSum(df)=sum(pVol & pSum & pDFOf & pDF);
end
subplot(4,1,4)
plot(SumDFOf)
plot(AllSum)

%%Apply Thresholds
%{

SumVol=sum(pVol);
SumSum=sum(pSum);
SumDFOf=sum(pDFOf);
SumDF=sum(pDF);

VolAndSum=sum(pVol & pSum);
DFAndDFOfSum=sum(pDFOf & pDF);
SumAndDFOfSum=sum(pSum & pDFOf);
AllSum=sum(pVol & pSum & pDFOf & pDF);
%}
 



