
clear all
TPN=GetMyDir
if ~exist([TPN 'find']), mkdir([TPN 'find']); end


'loading Dots'
load([TPN 'Dots.mat'])

%% Define Criteria
DistToMask=1;
Cut=2;

%%Lowest Criteria
Vol=0;
ItSum=0;
ITMax=0;
Con=0;
ConV=0;
DFOf=0;
DofTop=0;
DF=0;
Compact=10;
Long=10;

%%Standard Criteria
% Vol=3;
% ItSum=1;
% ITMax=1;
% Con=0.05;
% ConV=1;
% DFOf=.3;
% DofTop=.3;
% DF=2;
% Compact=10;
% Long=10;




%Record Criteria
SG.StartCrit.Vol=Vol;
SG.StartCrit.Cut=Cut;
SG.StartCrit.ItSum=ItSum;
SG.StartCrit.ITMax=ITMax;
SG.StartCrit.DFOf=DFOf;
SG.StartCrit.DFOfTopHalf=DofTop;
SG.StartCrit.DF=DF;
SG.StartCrit.DistToMask=DistToMask;
SG.StartCrit.Compact=Compact;
SG.StartCrit.Long=Long;


%% Scale ITMax

mB=Dots.MeanBright;
gBG=Dots.Im.GreenBackGround;
mBG=mB-gBG;
It=double(Dots.ITMax);
Contrast=It./mB;%max(1,(mB-(It*2)-gBG));
ContrastV=(It*10)./(Dots.Vol.^(2/3));
%Contrast=It./max(mBG,1);


%% Applly Criteria to Puncta
pVol=Dots.Vol>=Vol;
pCut=Dots.Cut(1:Dots.Num)<Cut;

pItSum=Dots.ItSum>=ItSum;
pITMax=Dots.ITMax>=ITMax;
pContrast=Contrast>=Con;
pContrastV=ContrastV>=ConV;
pDFOf=Dots.DFOf>=DFOf;
pDofTop=Dots.DFOfTopHalf>=DofTop;
pDF=Dots.DF>=DF;
pD2M=Dots.DistToMask<=DistToMask;
pCom=Dots.Round.Compact<=Compact;
pLong=Dots.Round.Long<=Long;

pass=pVol & pCut' & pItSum & pITMax & pContrast & pContrastV & pDFOf & pDofTop & pDF & pD2M & pCom & pLong;

P=find(pass')'; %% list of passing puncta
SG.pass1=pass';
save([TPN 'find\SG.mat'],'SG')

%% Draw image
'Drawing images'
maxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
AllmaxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
maxPassed=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');
YXsize=Dots.ImSize(1)*Dots.ImSize(2);
DisAmAll=maxPassed; DisAm=maxPassed;

%create pic of passed
for i =1:Dots.Num
  maxID(mod(Dots.Vox(i).Ind-1,YXsize)+1)=i; 
  DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)=DisAmAll(mod(Dots.Vox(i).Ind-1,YXsize)+1)+1;
end

for i = 1: size(P,2)
   maxID(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=P(i);  
   maxPassed(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200; %ContrastV(P(i)); 
   DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=DisAm(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1;
end

%ID all
for i = 1:Dots.Num
   AllmaxID(mod(Dots.Vox(i).Ind-1,YXsize)+1)=i;  
end

%load Raw
if ~exist([TPN 'images\maxRaw.mat'])
        maxRaw=imread([TPN 'images\RawMax.tif']);
        save([TPN 'images\maxRaw.mat'],'maxRaw');        
else
    load([TPN 'images\maxRaw.mat'])
end

%%combine and save and image Comparison
TestID=uint16(maxRaw)*2^8;
TestCrit=maxRaw;
TestCrit(:,:,1)=TestCrit(:,:,1);
TestCrit(:,:,2)=TestCrit(:,:,2);
TestID(:,:,3)=maxID;
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
TestCrit(:,:,3)=maxPassed;
imwrite(TestCrit,[TPN 'find\TestCrit.tif'],'Compression','none')
%imwrite(TestID,[TPN 'images\TestID.tif'],'Compression','none')
image(TestCrit),pause(.01)

if ~exist([TPN 'find\yes.tif']) | ~exist([TPN 'find\no.tif'])
    'Please provide yes.tif and no.tif for further analysis'
else

%% Read suggested criteria

YesIM=max(imread([TPN 'find\yes.tif']),[],3);
NoIM=max(imread([TPN 'find\no.tif']),[],3);
%load([TPN 'Dots.tif'])

[YesLab YesNum]=bwlabel(YesIM>0);
[NoLab NoNum]=bwlabel(NoIM>0);

%Find Miss
c=1; %initialize counter
clear Miss
Targ=maxID;
for i = 1:NoNum
   Targ=maxID(find(NoLab==i)); %retreive IDs of dots to be eliminated
   while sum(Targ)
       bTarg=find(Targ,1);
       Miss(c)=Targ(bTarg);
       c=c+1;
       Targ(Targ==Targ(bTarg))=0;
   end %done runing targ
end
SG.manual.Miss=Miss;

%Find Hit
c=1; %initialize counter
clear Hit
for i = 1:YesNum
   FindYes=find(YesLab==i); %retreive IDs of dots to be eliminated
   if ~sum(DisAmAll(FindYes)>1) %if unambiguous
       Targ=AllmaxID(FindYes);
   while sum(Targ)
       gTarg=find(Targ,1);
       Hit(c)=Targ(gTarg); %#ok<AGROW>
       c=c+1;
       Targ(Targ==Targ(gTarg))=0;
   end %done runing targ
   end
end

SG.manual.Hit=Hit;

DHit=maxID*0;
for i =1:size(Hit,2)
  DHit(mod(Dots.Vox(Hit(i)).Ind-1,YXsize)+1)=i; 
end


%% Look at Criteria for yesses and Miss

%Get mins and maxes

MinHitVol=min(Dots.Vol(Hit));
MinHitItSum=min(Dots.ItSum(Hit));
MinHitITMax=min(Dots.ITMax(Hit));
MinHitDFOf=min(Dots.DFOf(Hit));
MinHitDFOfTopHalf=min(Dots.DFOfTopHalf(Hit));
MinHitDF=min(Dots.DF(Hit));
MinHitDistToMask=max(Dots.DistToMask(Hit));
MinHitCompact=min(Dots.Round.Compact(Hit));
MinHitLong=max(Dots.Round.Long(Hit));


%% Scale DFOfTopHalf by Dist2CB
polyD=polyfit(Dots.Dist2CB(Hit),Dots.DFOfTopHalf(Hit)',1)
dic=1:200;
Dic=polyD(1)*dic +polyD(2);
PredDelta=polyD(1)*Dots.Dist2CB +polyD(2);
%DeltaScale=Dots.DFOfTopHalf-PredDelta';

DeltaScale=Dots.DFOfTopHalf./max((polyD(1)*Dots.Dist2CB+polyD(2)),.1)';

showY=Dots.DFOfTopHalf;
scatter(Dots.Dist2CB(Hit),showY(Hit))
hold on
scatter(Dots.Dist2CB(Miss),showY(Miss))
plot(dic,Dic)
hold off

%% Combine Delta and Contrast
DeltaCon=(DeltaScale/var(DeltaScale(Hit))) .* (Contrast/var(Contrast(Hit)));

%% Pick best Thresholds
Step=10;
clear Prop StartSearch StopSearch Interval Range
Prop(1,:)=DeltaScale;
Prop(2,:)=Contrast;
Prop(3,:)=DeltaCon;
Prop(4,:)=Dots.DFOfTopHalf;

Pnum=size(Prop,1);
%%%%%%%%%Replace steps with variance

for rep = 1:5

if rep==1    
    for i = 1:size(Prop,1)
        StartSearch(i)=min(Prop(i,Hit)); %#ok<AGROW>
        StopSearch(i)=max(Prop(i,Miss));
        Interval(i)=(StopSearch(i)-StartSearch(i))/Step;
        Range(i,:)=StartSearch(i):Interval(i):StopSearch(i); %#ok<AGROW>
    end
else
    for i = 1:size(Prop,1)
        StartSearch(i)=min(meanThreshes(i)-Interval(i)); %use previous interval as new range
        StopSearch(i)=max(meanThreshes(i)+Interval(i));
        Interval(i)=(StopSearch(i)-StartSearch(i))/Step;
        Range(i,:)=StartSearch(i):Interval(i):StopSearch(i);
    end    
end


%%Testing Thresholds
'Testing Thresholds'
    c=0;
    Eff=zeros((Step+1)^Pnum,3);
    TList=zeros((Step+1)^Pnum,Pnum);
    for p1=1:size(Range,2)
        for p2=1:size(Range,2)
            for p3=1:size(Range,2)
                for p4=1:size(Range,2)
                    c=c+1;
                    FalseP(1,:)=Prop(1,Miss)>=Range(1,p1);
                    FalseN(1,:)=Prop(1,Hit)<Range(1,p1);
                    FalseP(2,:)=Prop(2,Miss)>=Range(2,p2);
                    FalseN(2,:)=Prop(2,Hit)<Range(2,p2);
                    FalseP(3,:)=Prop(3,Miss)>=Range(3,p3);
                    FalseN(3,:)=Prop(3,Hit)<Range(3,p3);
                    FalseP(4,:)=Prop(4,Miss)>=Range(4,p4);
                    FalseN(4,:)=Prop(4,Hit)<Range(4,p4);
                    

                    FalsePsum=sum(sum(FalseP,1)==size(Prop,1));
                    FalseNsum=sum(sum(FalseN,1)>0);
                    Err=sum(FalsePsum + FalseNsum);
                    %ErrMat(p1,p2,p3,:)=Err;
                    Eff(c,:)=[FalsePsum FalseNsum Err];
                    TList(c,:)=[Range(1,p1) Range(2,p2) Range(3,p3) Range(4,p4)] ; %#ok<AGROW>
                end
            end
        end
    end
    
BestThreshes=find(Eff(:,3)==min(Eff(:,3)));
Balance=Eff(BestThreshes,1)-Eff(BestThreshes,2);
BestBalance=find(abs(Balance)==min(abs(Balance)));
meanThreshes=mean(TList(BestThreshes(BestBalance),:),1);
Cpos=BestThreshes(BestBalance(1));

DeltaThresh=meanThreshes(1)
ConThresh=meanThreshes(2)
DConThresh=meanThreshes(3)
FlatDeltaThresh=meanThreshes(4)
ErrorsMade=Eff(Cpos,3)

end %repeat Threshold finding


SG.FinalCrit.DeltaThresh=DeltaThresh;
SG.FinalCrit.ConThresh=ConThresh;
SG.FinalCrit.DConThresh=DConThresh;
SG.Performance.AllKnown=size(Hit,2)+size(Miss,2);
SG.Performance.KnownHits=size(Hit,2);
SG.Performance.KnownMisses=size(Miss,2);
SG.Performance.FalsePositives=Eff(Cpos,1);
SG.Performance.FalseNegatives=Eff(Cpos,2);
SG.Performance.Errors=Eff(Cpos,3);
SG.Performance.ErrorRate=SG.Performance.Errors/SG.Performance.AllKnown;



%% Applly Criteria to Puncta
pSGdeltaScale=DeltaScale>=DeltaThresh;
pSGcontrast=Contrast>=ConThresh;
pSGdcon=DeltaCon>=DConThresh;
pSGflatDelta=Dots.DFOfTopHalf>=FlatDeltaThresh;

pass2=pSGflatDelta & pSGdcon & pSGdeltaScale & pSGcontrast & pVol & pCut' & pItSum & pITMax & pContrast & pContrastV & pDFOf & pDofTop & pDF & pD2M & pCom & pLong;

%run final pass adding manual Hits and Misses
passF=pass2;
passF(Hit)=1;
passF(Miss)=0;

SG.pass2=pass2';
SG.passF=passF';




%% Save second order properties
SG.SecOrder.DeltaScale=DeltaScale;
SG.SecOrder.Contrast=Contrast;
SG.SecOrder.ContrastV=ContrastV;
SG.SecOrder.DeltaCon=DeltaCon;

save([TPN 'find\SG.mat'],'SG')






P=find(passF')'; %% list of passing puncta



%% Draw image
'Drawing images'

maxPassed2=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');

for i = 1: size(P,2)
   %maxID(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=P(i);  
   maxPassed2(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200; %ContrastV(P(i)); 
end

%combine and save and image Comparison
SGCrit=maxRaw;
SGCrit(:,:,1)=SGCrit(:,:,1);
SGCrit(:,:,2)=SGCrit(:,:,2);
FindEmpty=sum(sum(TestCrit,1),2)==0;
if sum(FindEmpty), Empty=find(FindEmpty,1); else Empty=1; end
SGCrit(:,:,3)=maxPassed2;
imwrite(SGCrit,[TPN 'find\SGCrit.tif'],'Compression','none')
%image(SGCrit), pause(.01)
%%

TPN
SG.Performance

end %find criteria if Yes and No are available
