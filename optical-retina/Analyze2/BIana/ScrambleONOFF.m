%%Scramble ON OFF

clear all
'Get Cells'
load(['UseCells.mat'])
[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
GetBIs=find(BIs>0);
BIAge=Age(BIs>0);
'Cells Gotten'
%}


NumCells=size(GetBIs,1);
Age=zeros(NumCells,1);
mDtA=zeros(Age,2); mDdA=zeros(Age,2); mDDA=zeros(Age,2);
NumNorm=Age; DevsDD=Age; meanDDs=Age;
HistDC=zeros(NumCells,8);
DtoDots=Age; 
Type=cell(NumCells,1);


yxum=0.103;
zum=0.3;

c=0;

clear Age mDtA mDdA mDDA Devs meanDDs HistDC



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

           load([TPN 'Use.mat'])
           DotNum(c)=size(Use.DPos,1);
           TotLength(c)=sum(Use.Length);
           Area(c,:)=CA.Area;

           load([TPN 'ONOFFa.mat'])

           DotNumA(c,1)=sum(ONOFF.DotLayer==1);
           DotNumA(c,2)=sum(ONOFF.DotLayer==2);
           TotLengthA(c,1)=sum(Use.Length(ONOFF.MidLayer==1));  
           TotLengthA(c,2)=sum(Use.Length(ONOFF.MidLayer==2));  


         end
end


%% Calculate real
DdOA=TotLengthA./Area;
DtOA=DotNumA./Area;
DD=DotNumA./TotLengthA;

TLd=(TotLengthA(:,1)-TotLengthA(:,2))./sum(TotLengthA,2);
DNd=(DotNumA(:,1)-DotNumA(:,2))./sum(DotNumA,2);
Ad=(Area(:,1)-Area(:,2))./sum(Area,2);

DdOAd=(DdOA(:,1)-DdOA(:,2))./sum(DdOA,2);
DtOAd=(DtOA(:,1)-DtOA(:,2))./sum(DtOA,2);
DDd=(DD(:,1)-DD(:,2))./sum(DD,2);  


%% Randomize 
reps =10
sampn=size(DotNumA,1);
for r = 1: reps
    picks=fix(rand(sampn,1)*sampn)+1;
    R=DotNumA(picks,:);
    
    
    %% Calculate
    rDdOA=TotLengthA./Area;
    rDtOA=R./Area;
    rDD=R./TotLengthA;

    rDdOAd=(rDdOA(:,1)-rDdOA(:,2))./sum(rDdOA,2);
    rDtOAd=(rDtOA(:,1)-rDtOA(:,2))./sum(rDtOA,2);
    rDDd=(rDD(:,1)-rDD(:,2))./sum(rDD,2);  

    PDtOA=rDD.*rDdOA; %Predict new density
    PDtOAd=(PDtOA(:,1)-PDtOA(:,2))./sum(PDtOA,2);
    
    Collect(r)=mean(abs(PDtOAd));
    
    
end
Real=mean(abs(DtOAd))
Collect

P=sum(Collect<=Real)/reps
    
    
    
    
