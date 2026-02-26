%function[] = anaCrit(TPN, DPN)

clear all
TPN=GetMyDir 

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
Dots.Crit.Vol=Vol;
Dots.Crit.Cut=Cut;
Dots.Crit.ItSum=ItSum;
Dots.Crit.ITMax=ITMax;
Dots.Crit.DFOf=DFOf;
Dots.Crit.DFOfTopHalf=DofTop;
Dots.Crit.DF=DF;
Dots.Crit.DistToMask=DistToMask;
Dots.Crit.Compact=Compact;
Dots.Crit.Long=Long;


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
Dots.OK=pass';
%save([TPN 'Dots.mat'],'Dots')

%% Draw image
maxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
AllmaxID=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
maxPassed=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');
YXsize=Dots.ImSize(1)*Dots.ImSize(2);

%create pic of passed
for i =1:Dots.Num
  maxID(mod(Dots.Vox(i).Ind-1,YXsize)+1)=i;  
end

for i = 1: size(P,2)
   maxID(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=P(i);  
   maxPassed(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200; %ContrastV(P(i));  
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

%combine and save and image Comparison
TestID=uint16(maxRaw)*2^8;
TestCrit=maxRaw;
TestCrit(:,:,1)=TestCrit(:,:,1)*3;
TestCrit(:,:,2)=TestCrit(:,:,2);
TestID(:,:,3)=maxID;
TestCrit(:,:,3)=maxPassed;
imwrite(TestCrit,[TPN 'images\TestCrit.tif'],'Compression','none')
imwrite(TestID,[TPN 'images\TestID.tif'],'Compression','none')
image(TestCrit)


%% Retreive values

[x y] = ginput
targ=maxID(round(y),round(x))

if targ>0
Targ.ItSum=Dots.ItSum(targ);
Targ.ITMax=Dots.ITMax(targ);
Targ.DFOf=Dots.DFOf(targ);
Targ.DofTop=Dots.DFOfTopHalf(targ);
Targ.DF=Dots.DF(targ);
Targ.D2M=Dots.DistToMask(targ);
Targ.Com=Dots.Round.Compact(targ);
Targ.Long=Dots.Round.Long(targ);
Targ.MeanBright=Dots.MeanBright(targ);
Targ.Contrast=Contrast(targ);
Targ.ContrastV=ContrastV(targ);

Targ
else
    'no dot selected'
end
 


