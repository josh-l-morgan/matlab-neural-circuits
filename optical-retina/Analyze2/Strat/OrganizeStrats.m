clear all
load('UseCells.mat')


c=0;

clear Age mDtA mDdA mDDA Devs meanDDs HistDC


NumCells=size(UseCells,1);
Age=zeros(NumCells,1);
mDtA=Age; mDdA=Age; mDDA=Age;
NumNorm=Age; DevsDD=Age; meanDDs=Age;
HistDC=zeros(NumCells,8);
DtoDots=Age; ONs=Age; OFFs=Age; Others=Age; BIs=Age;
Type=cell(NumCells,1);



for k = 1:NumCells
     TPN = char(UseCells(k,:)); 
     
    if exist([TPN 'CASc.mat'])
        c=c+1;
        
        load([TPN 'Cell.mat'])        
        
        Name(c)={Cell.Name};
        Age(c)=str2num(Cell.Age);
        Type{c}=Cell.Type;
        
                  %% read Cell type and mark appropriate matrix
    switch Cell.Type
        case 'ON'
            ONs(c)=1;
        case 'OFF'
            OFFs(c)=1;
        case 'BI'
            BIs(c)=1;
        case 'Other'
            Others(c)=1;
        otherwise
            Unknowns(c)=1;
    end %End Cell type switch
        

       
       load([TPN 'data\DepthDevNoCBmedianSc.mat'])
       NumDots=DepthDev.NumDots;
       NumMids=DepthDev.NumMids;
       DtoDots(c)=DepthDev.DtoDots;
       DtoMids(c)=DepthDev.DtoMids;
       HalfDots(c)=DepthDev.HalfDots;
       HalfMids(c)=DepthDev.HalfMids;
       mAPMean(c)=mean(DepthDev.AllDotsMean(NumDots>10));
       mADMean(c)=mean(DepthDev.AllMidsMean(NumMids>40));
       mAPMed(c)=mean(DepthDev.AllDotsMed(NumDots>10));
       mADMed(c)=mean(DepthDev.AllMidsMed(NumMids>40));       
             
       %Hist
       
       
    end
    k
end
ONs=ONs>0; OFFs=OFFs>0; BIs=BIs>0; Others=Others>0;
Age(Age>30)=35;
Mono=ONs | OFFs | (Age==5);


%% Scale by Resolution
load('Res.mat')


%% Stratification
CompareAges(Age(Mono),DtoMids(Mono),Age(Mono),DtoDots(Mono));
PlotAges(Age(Mono),DtoMids(Mono)'-DtoDots(Mono))
CompareAges(Age(Mono),HalfMids(Mono),Age(Mono),HalfDots(Mono))
CompareAges(Age(Mono),mADMean(Mono),Age(Mono),mAPMean(Mono))
CompareAges(Age(Mono),mADMed(Mono),Age(Mono),mAPMed(Mono))

dTo=(DtoMids'-DtoDots)./(DtoMids'+DtoDots);
dToMed=(mADMed-mAPMed)./(mADMed+mAPMed);
dToMean=(mADMean-mAPMean)./(mADMean+mAPMean);

PlotAges(Age(Mono),dToMean(Mono))
PlotAges(Age(Mono),dToMed(Mono))
PlotAges(Age(Mono),dTo(Mono))


PlotAges(Age(Mono),dToMed(Mono))
CompareAges(Age(Mono),mADMed(Mono),Age(Mono),mAPMed(Mono))
%CompareAges(Half














