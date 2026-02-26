clear all
load('UseCells.mat')

yxum=0.103;
zum=0.3;

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
     
    if exist([TPN 'CA.mat'])
        c=c+1;
        load([TPN 'CA.mat'])
        load([TPN 'Cell.mat'])        
        
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
        
        
        
        
        mDtA(c)=CA.mDtA;
        mDdA(c)=CA.mDdA;
        mDDA(c)=CA.mDDA;    
        
       load([TPN 'BranchS.mat'])
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
       
       load([TPN 'Grad'])
       GradV(c)=Grad.Vect;
       IDD(c)=Grad.IDD;
       RDD(c)=Grad.RDD;
       ODD(c)=Grad.ODD;
       Rdot(c)=Grad.Idot/Grad.Odot;
       Rdend(c)=Grad.Idend/Grad.Odend;
       Idot(c)=Grad.Idot;
       Odot(c)=Grad.Odot;
       Idend(c)=Grad.Idend;
       Odend(c)=Grad.Odend;
       
       
       Rad(c)=Grad.Outer;
       
       
       clear Grad
        
    end
    k
end
ONs=ONs>0; OFFs=OFFs>0; BIs=BIs>0; Others=Others>0;
Age(Age>30)=35;



%% Plot increase in density over area with age


PlotAges(Age,mDtA);
pause(.3)
PlotAges(Age,mDdA);
pause(.3)
PlotAges(Age,mDDA);

PlotAges(Age,meanDDs)
PlotAges(Age,DevsDD)

%% Show Dots per length scaled
for i = 1:1%size(HistDC,2),
    PlotAges(Age,HistDC(:,i))
    pause(.1)
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
 
 jmboot(DtoDots(Mono),DtoDots(BIs))
 
%% Plot Gradients
PlotAges(Age,GradV)
PlotAges(Age,RDD)
PlotAges(Age,IDD)
PlotAges(Age,ODD)
PlotAges(Age,Rdot)
PlotAges(Age,Rdend)
PlotAges(Age,Idot)
PlotAges(Age,Odot)
PlotAges(Age,Idend)
PlotAges(Age,Odend)



scatter(Age(OFFs),RDD(OFFs),'r')

CompareAges(Age(OFFs),RDD(OFFs),Age(BIs ),RDD(BIs))

[P R]=jmboot(RDD(Mono),RDD(BIs))


PlotAges(Age,Rad)




