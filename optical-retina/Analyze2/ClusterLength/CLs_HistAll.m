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
        
        
    
    load([TPN 'BranchS.mat'])
       UseBranch=Branch.Length<(Branch.Info.CutBranch+(Branch.Info.CutBranch/10));
       NumNorm(c)=sum(Branch.DotCount(UseBranch)==1)/sum(UseBranch);
       HistDC(c,:)=hist(Branch.DotCount(UseBranch),0:1:7)/sum(UseBranch);
       DevsDD(c)=std(Branch.DD(UseBranch));
       meanDDs(c)=mean(Branch.DD(UseBranch));
       meanNonZeDDs(c)=mean(Branch.DD(UseBranch & Branch.DotCount>0));
       HistDD(c,:)=hist(Branch.DD,0:.1:2); 

    
    
    end
end
      
hold off
for i = 1:size(HistDD,1)
    plot(HistDD(i,:),'b')
    hold on
end
plot(mean(HistDD,1),'r')
hold off



hold off
for i = 1:size(HistDC,1)
    plot([0:7], HistDC(i,:),'k')
    hold on
    xlim([0 7])
    ylim([0 .5])
    pause(.01)
    
end
plot([0:7],mean(HistDC,1),'r')
xlim([0 7])
ylim([0 .5])
hold off



