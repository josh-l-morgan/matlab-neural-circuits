%function[] = anaCrit(TPN, DPN)


KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));

for k = 1:size(Kdir,1)
    clear TPN DPN Compare maxRaw Dots
    TPN = [KPN '\' Kdir(k).name '\']     

if exist([TPN 'Dots.mat']) %if there is anything to run
load([TPN 'Dots.mat'])




%% Define Criteria
DistToMask=1;
Cut=2;

Vol=5;
ItSum=20;
ITMax=2;
Con=0.3;
ConV=15;
DFOf=1;
DofTop=1.5;
DF=5;
Compact=10;
Long=5;


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
save([TPN 'Dots.mat'],'Dots')

%% Draw image
maxPassed=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint8');
YXsize=Dots.ImSize(1)*Dots.ImSize(2);

%create pic of passed
for i = 1: size(P,2)
   maxPassed(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=200;     
end

%load Raw
if ~exist([TPN 'images\maxRaw.mat'])
        maxRaw=imread([TPN 'images\RawMax.tif']);
        save([TPN 'images\maxRaw.mat'],'maxRaw');
else
    load([TPN 'images\maxRaw.mat'])
end

%combine and save and image Comparison
Compare=maxRaw;
Compare(:,:,3)=maxPassed;
imwrite(Compare,[TPN 'images\Compare.tif'],'Compression','none')
image(Compare*4)


end %if Dots exist
end %run all cells



 


