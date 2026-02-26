
clear all
TPN=GetMyDir
if ~exist([TPN 'find']), mkdir([TPN 'find']); end


'loading Dots'
load([TPN 'Dots.mat'])
'Done Loading'

%% Define Criteria
DistToMask=1;
Cut=2;

Vol=3;
ItSum=1;
ITMax=1;
Con=0.05;
ConV=1;
DFOf=.3;
DofTop=.3;
DF=2;
Compact=10;
Long=10;


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
TestCrit(:,:,1)=TestCrit(:,:,1)*3;
TestCrit(:,:,2)=TestCrit(:,:,2);
TestID(:,:,3)=maxID;
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
DeltaScale=Dots.DFOfTopHalf-PredDelta';
%DeltaScale=Dots.DFOfTopHalf'./(polyD(1)*Dots.Dist2CB+polyD(2));

showY=Dots.DFOfTopHalf;
scatter(Dots.Dist2CB(Hit),showY(Hit))
hold on
scatter(Dots.Dist2CB(Miss),showY(Miss))
plot(dic,Dic)
hold off


%DeltaScale=Dots.DFOfTopHalf;

%% Pick best Thresholds
StartSearchD=min(DeltaScale(Hit));
StopSearchD=max(DeltaScale(Miss));
Drange=StartSearchD:(StopSearchD-StartSearchD)/1000:StopSearchD;
StartSearchC=min(Contrast(Hit));
StopSearchC=max(Contrast(Miss));
Crange=StartSearchC:(StopSearchC-StartSearchC)/100:StopSearchC;


%% Testing Thresholds
'Testing Thresholds'
if StopSearchD<StartSearchD | StopSearchC<StartSearchC
    ConThresh=mean(StopSearchC,StartSearchC);
    DeltaThresh=mean(StopSearchD,StartSearchD);
else
    Eff=zeros(size(Crange,2),size(Drange,2),5);
    for c=1:size(Crange,2);
        for d=1:size(Drange,2)
            FalsePc=Contrast(Miss)>=Crange(c);
            FalseNc=Contrast(Hit)<Crange(c);
            FalsePd=DeltaScale(Miss)>=Drange(d);
            FalseNd=DeltaScale(Hit)<Drange(d);  
            FalseP=sum(FalsePc & FalsePd);
            FalseN=sum(FalseNc | FalseNd);
            Err=sum(FalseP + FalseN);
            Eff(c,d,:)=[Crange(c) Drange(d) FalseP FalseN Err];
        end
    end
end
scatter(Eff(:,1),Eff(:,2))
image(Eff(:,:,5)*3),pause(.01)

BestThreshesIM=Eff(:,:,5)==min(min(Eff(:,:,5))); %create mask of threshes with fewest mistakes
[LabBest lbn]=bwlabel(BestThreshesIM); %label all minima
for lb = 1:lbn; lbSize(lb)=sum(sum(LabBest==lb)); end  %find size of each minima
lbtarg=find(lbSize==max(lbSize)); %target largest patch
[y x]=find(LabBest==lbtarg);
centy=mean(y); centx=mean(x);
Tdist=dist([y x ones(size(y),1)],[centy centx 1]);
UseTarg=find(Tdist==min(Tdist),1);
Cpos=y(UseTarg); %final Contrast position in range
Dpos=x(UseTarg); %final Delta position in range

DeltaThresh=Drange(Dpos);
ConThresh=Crange(Cpos);

SG.FinalCrit.DeltaThresh=DeltaThresh;
SG.FinalCrit.ConThresh=ConThresh;
SG.Performance.AllKnown=size(Hit,2)+size(Miss,2);
SG.Performance.KnownHits=size(Hit,2);
SG.Performance.KnownMisses=size(Miss,2);
SG.Performance.FalsePositives=Eff(Cpos,Dpos,3);
SG.Performance.FalseNegatives=Eff(Cpos,Dpos,4);
SG.Performance.Errors=Eff(Cpos,Dpos,5);
SG.Performance.ErrorRate=SG.Performance.Errors/SG.Performance.AllKnown;



%% Applly Criteria to Puncta
pSGdeltaScale=DeltaScale>=DeltaThresh;
pSGcontrast=Contrast>=ConThresh;

pass2=pSGdeltaScale & pSGcontrast & pVol & pCut' & pItSum & pITMax & pContrast & pContrastV & pDFOf & pDofTop & pDF & pD2M & pCom & pLong;

%run final pass adding manual Hits and Misses
passF=pass2;
passF(Hit)=1;
passF(Miss)=0;

SG.pass2=pass2';
SG.passF=passF';
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
SGCrit(:,:,1)=SGCrit(:,:,1)*3;
SGCrit(:,:,2)=SGCrit(:,:,2);
SGCrit(:,:,3)=maxPassed2;
imwrite(SGCrit,[TPN 'find\SGCrit.tif'],'Compression','none')
%image(SGCrit)
%%

SG.Performance

end %find criteria if Yes and No are available
