
%clear all
'Get Cells'
load(['UseCells.mat'])
[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
GetBIs=find(BIs>0);
BIAge=Age(BIs>0);
'Cells Gotten'
%}



yxum=0.103;
zum=0.3;

c=0;

clear Age mDtA mDdA mDDA Devs meanDDs HistDC



NumCells=size(GetBIs,1);
Age=zeros(NumCells,1);
mDtA=zeros(Age,2); mDdA=zeros(Age,2); mDDA=zeros(Age,2);
NumNorm=Age; DevsDD=Age; meanDDs=Age;
HistDC=zeros(NumCells,8);
DtoDots=Age; 
Type=cell(NumCells,1);



for k = 1:NumCells
     TPN = char(UseCells(GetBIs(k),:)); 
     
    if exist([TPN 'CAbi.mat'])
        c=c+1;
        load([TPN 'CAbi.mat'])
        load([TPN 'Cell.mat'])        
        
        Name(c)={Cell.Name};
        Age(c)=str2num(Cell.Age);
        Type{c}=Cell.Type;
        
                  %% read Cell type and mark appropriate matrix
      
        mDtA(c,:)=CA.mDtA;
        mDdA(c,:)=CA.mDdA;
        mDDA(c,:)=CA.mDDA;    
        
        %record measures for whole cell. 
        mDtAc(c)=sum(CA.mDtA.*CA.Area)/sum(CA.Area);
        mDdAc(c)=sum(CA.mDdA.*CA.Area)/sum(CA.Area);
        mDDAc(c)=sum(CA.mDDA.*CA.Area)/sum(CA.Area);
        

        
        
       load([TPN 'BranchSb.mat'])
       UseBranch=Branch.Length<(Branch.Info.CutBranch+(Branch.Info.CutBranch/10));
       NumNorm(c)=sum(Branch.DotCount(UseBranch)==1)/sum(UseBranch);
       HistDC(c,:)=hist(Branch.DotCount(UseBranch),0:1:7)/sum(UseBranch);
       DevsDD(c)=std(Branch.DD(UseBranch));
       meanDDs(c)=mean(Branch.DD(UseBranch));
       
       load([TPN 'data\DepthDevNoCB.mat'])
       DtoDots(c)=DepthDev.DtoDots;
       DotsA(c)=DepthDev.DotsA;
       DendA(c)=DepthDev.DendA;
       DtoMids(c)=DepthDev.DtoMids;
       clear DepthDev
       
       load([TPN 'GradABi'])


       RDD(c,:)=Grad.RDD;
       ODD(c,:)=Grad.ODD;
       Rdot(c,:)=Grad.Rdots;
       Rdend(c,:)=Grad.Rdend;
       Idot(c,:)=Grad.Idot;
       Odot(c,:)=Grad.Odot;
       Idend(c,:)=Grad.Idend;
       Odend(c,:)=Grad.Odend;
     
       
       %Rad(c)=Grad.Outer;
       
       load([TPN 'Near'])
       NearA(c,1)=Near.medA;
       NearV(c,1)=Near.medV;
             
       clear Grad
       
       load([TPN 'Use.mat'])
       DotNum(c)=size(Use.DPos,1);
       TotLength(c)=sum(Use.Length);
       Area(c,:)=CA.Area;
       
       load([TPN 'ONOFFa.mat'])
        
       DotNumA(c,1)=sum(ONOFF.DotLayer==1);
       DotNumA(c,2)=sum(ONOFF.DotLayer==2);
       TotLengthA(c,1)=sum(Use.Length(ONOFF.MidLayer==1));  
       TotLengthA(c,2)=sum(Use.Length(ONOFF.MidLayer==2));  
       
       % ONOFFS(k,1)=mean(Branch.DD(ONOFF.SegLayer==1));
        %ONOFFS(k,2)=mean(Branch.DD(ONOFF.SegLayer==2));

       
       
       load([TPN 'Neighbors.mat'])
       NNdpos(c)=Neighbors.med.DPos;
       NNnn(c)=Neighbors.med.NN;
       NNrand(c)=Neighbors.med.Rand;
       NNepos(c)=Neighbors.med.EPos;
       
       
    end
    k
end
Age(Age>30)=35;

%% Plot ON and OFF densities
CompareAges(Age,DotNumA(:,1)./Area(:,1),Age,DotNumA(:,2)./Area(:,2))
CompareAges(Age,TotLengthA(:,1)./Area(:,1),Age,TotLengthA(:,2)./Area(:,2))
CompareAges(Age,DotNumA(:,1)./TotLengthA(:,1),Age,DotNumA(:,2)./TotLengthA(:,2))

CompareAges(Age,Rdot(:,1),Age,Rdot(:,2))
CompareAges(Age,Rdend(:,1),Age,Rdend(:,2))
CompareAges(Age,RDD(:,1),Age,RDD(:,2))


%% Plot increase in density over area with age

mDtAd=(mDtA(:,1)-mDtA(:,2))./sum(mDtA,2);
PlotAges(Age,mDtAd+1)
mDdAd=(mDdA(:,1)-mDdA(:,2))./sum(mDdA,2);
PlotAges(Age,mDdAd+1)
mDDAd=(mDDA(:,1)-mDDA(:,2))./sum(mDDA,2);
PlotAges(Age,mDDAd+1)

[R P]=corrcoef(mDdAd, mDDAd) 
[R P]=corrcoef(mDdAd, mDtAd)
[R P]=corrcoef(mDtAd, mDDAd) 


bar([mDdAd(Age==12) mDtAd(Age==12) mDDAd(Age==12)])
ylim([-1.5 1.5])
bar([mDdAd mDtAd mDDAd])
bar([mDdAd(Age==7) mDtAd(Age==7) mDDAd(Age==7)])

%%Raw Scale by area
DdOA=TotLengthA./Area;
DtOA=DotNumA./Area;
DD=DotNumA./TotLengthA;

TLd=(TotLengthA(:,1)-TotLengthA(:,2))./sum(TotLengthA,2);
DNd=(DotNumA(:,1)-DotNumA(:,2))./sum(DotNumA,2);
Ad=(Area(:,1)-Area(:,2))./sum(Area,2);

DdOAd=(DdOA(:,1)-DdOA(:,2))./sum(DdOA,2);
DtOAd=(DtOA(:,1)-DtOA(:,2))./sum(DtOA,2);
DDd=(DD(:,1)-DD(:,2))./sum(DD,2);
%dS=(ONOFFS(:,1)-ONOFFS(:,2))./(ONOFFS(:,1)+ONOFFS(:,2)); %Get dDD from good segs


[R P]=corrcoef(DdOAd, DDd)
[R P]=corrcoef(TLd, DDd)
[R P]=corrcoef(TLd, DNd)
[R P]=corrcoef(DdOAd, DtOAd)
scatter(DdOAd,DtOAd)






PlotAges(Age,DdOAd+1);
PlotAges(Age,DtOAd+1);
bar([DdOAd DtOAd DDd])
scatter(DdOAd, DDd)
bar([ Ad TLd DdOAd DNd DtOAd DDd])

%%Check Out absolute numbers and lengths
DNd=(DotNumA(:,1)-DotNumA(:,2))./mean(DotNumA,2);
TLd=(TotLengthA(:,1)-TotLengthA(:,2))./mean(TotLengthA,2);
DDa=DotNumA./TotLengthA;
DDad=(DDa(:,1)-DDa(:,2))./mean(DDa,2);


PlotAges(Age,DNd+1)
PlotAges(Age,TLd+1)
PlotAges(Age,mean(DDa,2))
PlotAges(Age,DDad+1)
bar([TLd DNd DDad])

PlotAges(Age,DotNum)

%DdOA=TotLength./Area;
DtOA=DotNum./Area;
PlotAges(Age,DdOA);
PlotAges(Age,DtOA);

bar([mDdAd mDtAd DDad])




scatter(DDad,mDDAd)


PlotAges(Age,mDtAc);
pause(.3)
PlotAges(Age,mDdAc);
pause(.3)
PlotAges(Age,mDDAc);

%jmboot(mDtA(Age==12),mDtA(Age==35))


PlotAges(Age,meanDDs)
PlotAges(Age,DevsDD)

%%corrolate density with 
scatter(mDdA(Age==7),mDDA(Age==7))
scatter(mDdA(Age>30),mDDA(Age>30))

%% Show Dots per length scaled
for i = 1:size(HistDC,2),
    i
    PlotAges(Age,HistDC(:,i))
    hold on
end




PlotAges(Age,HistDC(:,3))
%% Stratification as standard deviation from local mean
Mono=ONs | OFFs | (Age==5);
 PlotAges(Age(Mono),DendA(Mono))
 PlotAges(Age(Mono),DotsA(Mono))
 PlotAges(Age(Mono),DtoMids(Mono))
 PlotAges(Age(Mono),DtoDots(Mono))
 PlotAges(Age,DtoDots)
 
 
 CompareAges(Age(Mono),DendA(Mono),Age(BIs),DendA(BIs))
 CompareAges(Age(Mono),DotsA(Mono),Age(BIs),DotsA(BIs))
 CompareAges(Age(Mono),DendA(Mono),Age(Mono),DotsA(Mono));
 
 CompareAges(Age(Mono),DtoMids(Mono),Age(BIs),DtoMids(BIs))
 CompareAges(Age(Mono),DtoDots(Mono),Age(BIs),DtoDots(BIs))
 CompareAges(Age(Mono),DtoMids(Mono),Age(Mono),DtoDots(Mono));
 
 %jmboot(DtoDots(Mono),DtoDots(BIs))
 
%% Plot Gradients
PlotAges(Age,GradV)
PlotAges(Age,RDD)
PlotAges(Age,IDD)
PlotAges(Age,ODD)
PlotAges(Age,Rdot)
PlotAges(Age,Rdend)
PlotAges(Age,RnDD+1)
PlotAges(Age,IDD)
PlotAges(Age,ODD)
PlotAges(Age,RnDot+1)
PlotAges(Age,RnDend+1)



scatter(Age(OFFs),RDD(OFFs),'r')

CompareAges(Age(OFFs),RDD(OFFs),Age(BIs ),RDD(BIs))

%[P R]=jmboot(RDD(Mono),RDD(BIs))

scatter(Rdend,RDD)
scatter(Rdend(Age>10),RDD(Age>10))

[R P]=corrcoef(Rdend,RDD)


%% Radius
PlotAges(Age,Rad)

%% Nearest Neighbor
PlotAges(Age,NearA)
PlotAges(Age,NearV)
PlotAges(Age(Mono),NearA(Mono))
PlotAges(Age(Mono),NearV(Mono))

PlotAges(Age,NNdpos)
PlotAges(Age,NNnn)
PlotAges(Age,NNrand)
PlotAges(Age,NNepos)
PlotAges(Age,NNnn./NNepos)
PlotAges(Age,NNnn./NNrand)

%% Dot Number

PlotAges(Age,TotLength);
PlotAges(Age,DotNum);
PlotAges(Age,DotNum./TotLength);


%% Make Table

%%Table of mDdA mDtA DD   dDend dDot dDD   Age
Arbor=2;
CellNum=size(mDdA,1);
TabL=cell(CellNum,8);
TabL(:,1)=num2cell(DdOA(:,Arbor));
TabL(:,2)=num2cell(DtOA(:,Arbor));
TabL(:,3)=num2cell(DD(:,Arbor));

TabL(:,7)=num2cell(Age);
TabL(:,8)=num2cell(Name);

%%Find Ages
if size(Age,2)>1; Age=Age';end %make Age a column
Ages=[];
for i = 1:size(Age,1)
    if isempty(find(Ages==Age(i)))
        Ages=[Ages;Age(i)];
    end
end
Ages=sort(Ages);

TabLs=[];
for i = 1: size(Ages,1)
   TabLs=[TabLs;TabL(Age==Ages(i),:)]; 
end

Vars=['D/A','P/A','D/P','D/A Grad','P/A Grad','D/P Grad'];

%% Make Table Weighted means

%%Table of mDdA mDtA DD   dDend dDot dDD   Age
CellNum=size(mDdA,1);
TabL=cell(CellNum,5);
TabL(:,1)=num2cell(sum(TotLengthA,2)./sum(Area,2));
TabL(:,2)=num2cell(sum(DotNumA,2)./sum(Area,2));
TabL(:,3)=num2cell(sum(DotNumA,2)./sum(TotLengthA,2));

TabL(:,4)=num2cell(Age);
TabL(:,5)=num2cell(Name);

%%Find Ages
if size(Age,2)>1; Age=Age';end %make Age a column
Ages=[];
for i = 1:size(Age,1)
    if isempty(find(Ages==Age(i)))
        Ages=[Ages;Age(i)];
    end
end
Ages=sort(Ages);

TabLs=[];
for i = 1: size(Ages,1)
   TabLs=[TabLs;TabL(Age==Ages(i),:)]; 
end

Vars=['D/A','P/A','D/P','D/A Grad','P/A Grad','D/P Grad'];













