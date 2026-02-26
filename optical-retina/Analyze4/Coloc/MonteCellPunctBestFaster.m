
clear all
colormap  gray(256)
TPN=GetMyDir
Reps=10;

%% Get Presynaptic field

'reading Presynaptic Image'
Idir=dir([TPN 'I']); Idir=Idir(3:size(Idir,1));
Ip=imread([TPN 'I\' Idir(1).name]);
[ys xs c] = size(Ip);
zs=size(Idir,1);
I=zeros(ys,xs,zs,'uint8');

%% Extract Presyn channel
for i = 1: zs
   Ip=imread([TPN 'I\' Idir(i).name]);
   I(:,:,i)=Ip(:,:,3);
end
'done reading image'

%% Binarize Presynaptic side
IT = PreSyn(I);
imwriteNp(TPN,IT,'IT.tif')

%%Get Thresholds
for i = 1: zs
   Samp=(I(:,:,max(1,i-2))); 
   MedThresh(i)=median(Samp(:));
   MeanThresh(i)=mean(Samp(:));
end

Imax=max(I,[],3);
image(Imax),pause(.1)


%% Get Post Synaptic Field
'Get Synapse info'

load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])
Pass=SG.passF;
VoxPs=Dots.Vox(Pass);
Pos=Dots.Pos(Pass,:);


%%
SE=strel('ball',3,3);
Cube=ones(7,7,3);
for y=1:7;for x=1:7;for z=1:3;
            Cube(y,x,z)=sqrt((y-4)^2 + (x-4)^2 + (z*3-6)^2);
end, end, end
Cube=Cube<4;
image(sum(Cube,3)*40)

for i = 1: length(VoxPs)
    VoxP=VoxPs(i).Pos;
    mins=min(VoxP,[],1);
    maxes=max(VoxP,[],1);
    DiSpaceV=zeros(maxes-mins+[7 7 3]);
    clear DiVox
    for d = 1:3
        DiVox(:,d)=VoxP(:,d)-mins(d)+1;
    end
    DiVoxInd=sub2ind(size(DiSpaceV),DiVox(:,1)+3,DiVox(:,2)+3,DiVox(:,3)+1);
    DiSpaceV(DiVoxInd)=1;
    DiSpace=DiSpaceV;
    DiSpace=imdilate(DiSpace,Cube);
    %DiSpace=imdilate(DiSpace,SE);
    DiPerim=bwperim(DiSpaceV,26);
    DiSpaceSurround=DiSpace-DiSpaceV+DiPerim;
    %image(sum(DiSpaceSurround,3)*30),pause
    [Sy Sx Sz]=ind2sub(size(DiSpace),find(DiSpaceSurround));
    Sur=[Sy+mins(1)-4,Sx+mins(2)-4,Sz+mins(3)-2];
  
    SurroundA(i)={Sur};
    if mod(i,100)==0; PercentDoneSurrounding=i/length(VoxPs)*100, end
end




%% Get Dendrites
'Get Dendrites'
mdir=dir([TPN 'mask']); mdir=mdir(3:size(mdir,1));
Dend=imread([TPN 'mask\' mdir(1).name]);
[ys xs]=size(Dend); zs=size(mdir,1);
Dend=zeros(ys,xs,zs,'uint8');
for i = 1:zs
    Dend(:,:,i)=imread([TPN 'mask\' mdir(i).name]);
end
Dend=Dend==1;
Dend=bwperim(Dend,6);


sDend=uint8(sum(Dend,3));

%% Measure Real Surrounds

%%Find all in bounds Surrounds
InBounds=ones(1,length(VoxPs));
for i = 1:length(VoxPs)
    Check=cell2mat(SurroundA(i));
    if sum(Check(:)<1) | sum(Check(:,1)>ys) | sum(Check(:,2)>xs) | sum(Check(:,3)>zs)
        InBounds(1,i)=0;
    end
    Check=Check(~((sum(Check<1,2)) |(Check(:,1)>ys) | (Check(:,2)>xs) | (Check(:,3)>zs)),:);
    Surround(i)={Check};
end


InB=find(InBounds);
siz=[ys xs zs];
clear Top Top5 Mean Peak Sum Blues
ShowSurround=I*0;
for i =1: length(VoxPs)
   Plane=round(Pos(i,3)); 
   Check=cell2mat(Surround(i));
   CheckI=sub2ind(siz,Check(:,1),Check(:,2),Check(:,3));
   ShowSurround(CheckI)=ShowSurround(CheckI)+uint8(Pass(i));
   Voxes(i)=size(CheckI,1);
   Blue=I(CheckI);
   Blue=sort(Blue);
   Blues(i)={Blue};
   Sum(i)=sum(Blue);
   Mean(i)=mean(Blue);
   Peak(i)=max(Blue(:));
   Top(i)=mean(Blue(fix(.75*size(Blue,1)):size(Blue,1)));
   Top5(i)=mean(Blue(max(1,size(Blue,1)-10):size(Blue,1)));
   pThresh(i)=sum(Blue>2*MedThresh(Plane));
   BlueT=IT(CheckI);
   sumIT(i)=sum(BlueT);
end
imwriteNp(TPN,ShowSurround,'ShowSurround')

%% Run Monte From Cell

Run=find(InBounds');



[fdy fdx fdz]=ind2sub(siz,find(Dend));
AllDend=[fdy fdx fdz];


mcSum=zeros(size(Run,1),Reps);
mcTop=mcSum; mcTop5=mcSum; mcPeak=mcSum; mcMean=mcSum;
'starting Monte'
New=zeros(siz,'uint8');  
for i = 1:size(Run,1)
fDend=AllDend(AllDend(:,3)>(Pos(Run(i),3))-4 & (AllDend(:,3)<(Pos(Run(i),3))+4),:);    %!!!!!!
for r = 1:Reps
    id = Run(i);
    ID(i)=id;
    
    %%Get relative surround
    rSur=cell2mat(Surround(id));
    for d = 1:3
        rSur(:,d)=rSur(:,d)-round(mean(rSur(:,d)));
    end
    
    %%Restrict Dendrite search area

    while 1 %Search for viable spot
            Pick=fix(rand*size(fDend,1))+1;
            NewSur=[rSur(:,1)+fDend(Pick,1) rSur(:,2)+fDend(Pick,2) rSur(:,3)+fDend(Pick,3)];
            if~(sum(NewSur(:)<1) | sum(NewSur(:,1)>ys) | sum(NewSur(:,2)>xs) | sum(NewSur(:,3)>zs))
                CheckI=sub2ind(siz,NewSur(:,1),NewSur(:,2),NewSur(:,3));
                break
                %if~sum(New(CheckI)), break, end %if space not already occupied
            end
    end
       Plane=round(Pos(i,3)); 
      % New(CheckI)=1;
       Blue=I(CheckI);
       Blue=sort(Blue);
       %Blues(r)={Blue}; %%!!!!!!!!!!!!!!!!!!!! Ditch later or eat Mem
       mcSum(i,r)=sum(Blue);
       mcMean(i,r)=mean(Blue);
       mcPeak(i,r)=max(Blue(:));
       mcTop(i,r)=mean(Blue(fix(.75*size(Blue,1)):size(Blue,1)));
       mcTop5(i,r)=mean(Blue(max(1,size(Blue,1)-10):size(Blue,1)));
       mcpThresh(i,r)=sum(Blue>2*MedThresh(fDend(Pick,3)));
       %%Need a threshold measure
       BlueT=IT(CheckI);
       mcsumIT(i,r)=sum(BlueT); 
 

end
PercentDoneMonte=i/size(Run,1)*100
NewS=max(New,[],3);
%image(sDend*3+NewS*30),pause(.01)

end


%%  Analyze Results
CNum=size(Run,1)
VNum=fix(CNum/20);

Rat=0;
Prop=1;
c=0;
for Rat=0:.1:1
    for Prop=.5:.1:2
        c=c+1;
%%Scale all to puncta mean
%Top5n=Top5(Pass)'./mean(mcTop5,2);

%%Pick Measure 
mcM=mcsumIT * (Rat)+mcpThresh *(1-Rat);
M=sumIT * Rat + pThresh *(1-Rat);


Thresh=median(mcM(:))*Prop;
Colo=M(Run)>Thresh;
Res=sum(Colo);
mcColo=mcM>Thresh;
mcRes=sum(mcColo,1);
P=sum(mcRes>=Res)/Reps;


%%Check Validity
%%Validity is the frequency that you would see any non-coloc in a
%%population where 5% of puncta are random
clear VM
% 
% for r = 1:100
%     for i = 1:size(mcM,2)
%         Rpic=fix(rand(7,1)*CNum)+1;
%         VM(:,i)=mcM(Rpic,i);
%     end
% 
% VColo=VM>Thresh;
% VRes=sum(VColo,1);
% Val(r)=sum(VRes>=VNum)/Reps;
% end
% Validity=mean(Val);

%%Plot possible Synapses
%Cs = probability of colocalization given synapse
Cr=mean(mcRes)/size(Run,1); %Probability of colocalization givin random
Co=Res/size(Run,1); %Observed number of colocalizations
Cs=[0:.01:1];

S=(Co-Cr)./(Cs-Cr);
plot(S);
ylim([0 1]);
minRealC(c)=(Co-Cr)./(1-Cr);
PropC(c)=Prop;
RatC(c)=Rat;
CoC(c)=Co;
CrC(c)=Cr;
    end
end

BestR=find(CrC<.5 & minRealC==max(minRealC(CrC<.5)));
if size(BestR,2)>1,
    BestE=find(CrC(BestR)==min(CrC(BestR)),1); 
    BestR=BestR(BestE);
end


BestCo.minReal=minRealC(BestR);
BestCo.Co=CoC(BestR);
BestCo.Cr=CrC(BestR);
BestCo.Prop=PropC(BestR);
BestCo.Rat=RatC(BestR);

%%Check Validity for Best
mcM=mcsumIT * (BestCo.Rat)+mcpThresh *(1-BestCo.Rat);
M=sumIT * BestCo.Rat + pThresh *(1-BestCo.Rat);

Thresh=median(mcM(:))*BestCo.Prop;
Colo=M(Run)>Thresh;
Res=sum(Colo);
mcColo=mcM>Thresh;
mcRes=sum(mcColo,1);
P=sum(mcRes>=Res)/Reps;

%%Check Validity
%%Validity is the frequency that you would see any non-coloc in a
%%population where 5% of puncta are random
clear VM

for r = 1:100
    for i = 1:size(mcM,2)
        Rpic=fix(rand(7,1)*CNum)+1;
        VM(:,i)=mcM(Rpic,i);
    end

VColo=VM>Thresh;
VRes=sum(VColo,1);
Val(r)=sum(VRes>=VNum)/Reps;
end
Validity=mean(Val);

BestCo.P=P;
BestCo.Val=Validity;
BestCo.NumDots=CNum;
BestCo.Reps=Reps;
BestCo

%%  Analyze Results
CNum=size(Run,1);
VNum=fix(CNum/20);

%%Scale all to puncta mean
%Top5n=Top5(Pass)'./mean(mcTop5,2);
Rat=0;
Prop=1.1;

%%Pick Measure 
mcM=mcsumIT * (Rat)+mcpThresh *(1-Rat);
M=sumIT * Rat + pThresh *(1-Rat);

Thresh=median(mcM(:))*Prop;
Colo=M(Run)>Thresh;
Res=sum(Colo);
mcColo=mcM>Thresh;
mcRes=sum(mcColo,1);
P=sum(mcRes>=Res)/Reps;

%%Check Validity
%%Validity is the frequency that you would see any non-coloc in a
%%population where 5% of puncta are random
clear VM

for r = 1:100
    for i = 1:size(mcM,2)
        Rpic=fix(rand(7,1)*CNum)+1;
        VM(:,i)=mcM(Rpic,i);
    end

VColo=VM>Thresh;
VRes=sum(VColo,1);
Val(r)=sum(VRes>=VNum)/Reps;
end
Validity=mean(Val);

%%Plot possible Synapses
%Cs = probability of colocalization given synapse
Cr=mean(mcRes)/size(Run,1); %Probability of colocalization givin random
Co=Res/size(Run,1); %Observed number of colocalizations
Cs=[0:.01:1];

S=(Co-Cr)./(Cs-Cr);
plot(S);
ylim([0 1]);
minReal=(Co-Cr)./(1-Cr);

StdCo.P=P;
StdCo.minReal=minReal;
StdCo.Co=Co;
StdCo.Cr=Cr;
StdCo.Prop=Prop;
StdCo.Rat=Rat;
StdCo.Val=Validity;
StdCo.NumDots=CNum;
StdCo.Reps=Reps;

StdCo

