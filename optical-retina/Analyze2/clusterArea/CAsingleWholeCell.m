
clear all

TPN=GetMyDir

%% load info
'Load Data'
load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])
load([TPN 'data\AllSegCut.mat'])

%%Extract Dot Positions
OK=SG.passF;
DPos=Dots.Pos(OK,:);
DPos(:,1:2)=DPos(:,1:2)*.103; DPos(:,3)=DPos(:,3)*.3;

%%Extract Dend positions 
Mids=mean(AllSegCut,3);
Length=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 ...
           + (AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 ...
           + (AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);


%% create Nearest Node list (nearest node for each dot and distance)
for i = 1:size(DPos,1)
            Ndist=dist(Mids,DPos(i,:)); %find dist from dot to all nodes
            Near=min(Ndist); %find shortest distance
            Nearest=find(Ndist==Near,1); %get node at that distance
            NN(i,:)=Mids(Nearest,:); %assign that node to NearestNode list for dots
            DotToNN(i,:)=Near; %record that distance for posterity  
end


Mids=round(Mids);
NN=round(NN);

%% Draw maps
'Draw Maps'

%%insure no negs
minDPos=min(NN)-1; minMids=min(Mids)-1;
minBoth=min(minDPos,minMids);
for d = 1:3, 
    if minBoth(d)<0;
        NN(:,d)=NN(:,d)-minBoth(:,d); 
        Mids(:,d)=Mids(:,d)-minBoth(:,d); 
    end
end

Msize=fix(max(max(NN),max(Mids)))+1;

%%Draw Dot map
DotMap=zeros(Msize(1),Msize(2));
for i=1:size(NN,1)
    DotMap(NN(i,1),NN(i,2))=DotMap(NN(i,1),NN(i,2))+1;
end

%%Draw Dend map
DendMap=zeros(Msize(1),Msize(2));
for i=1:size(Mids,1)
    DendMap(Mids(i,1),Mids(i,2))=DendMap(Mids(i,1),Mids(i,2))+Length(i);
end


%% make distance filter

AreaS=20;
Disk=fspecial('disk',AreaS);
DotFilt=imfilter(DotMap,Disk,'same');
DendFilt=imfilter(DendMap,Disk,'same');
BlankFilt=imfilter(double(DendMap>-1000),Disk,'same');
DotFilt=DotFilt./BlankFilt;
DendFilt=DendFilt./BlankFilt;

DendSkel=DendFilt;
DendSkel(DendMap>0)=0;
image(DendSkel*1000)

DDSkel=DotFilt./DendFilt;
DDSkel(DendMap>0)=0;
image(DDSkel*500)

DotSkel=DotFilt;
DotSkel(DendMap>0)=0;
image(DotSkel*2000)


%% Plot Blocks
%{
for s = 20;%1:1 %scaled microns
    clear DotMapS DendMapS
    DPosS=round(NN/s)+1;
    MidsS=round(Mids/s)+1;
    
    Msize=fix(max(max(DPosS),max(MidsS)))+1;

    %%Draw Dot map
    DotMapS=zeros(Msize(1),Msize(2));
    for i=1:size(DPosS,1)
        DotMapS(DPosS(i,1),DPosS(i,2))=DotMapS(DPosS(i,1),DPosS(i,2))+1;
    end

    %%Draw Dend map
    DendMapS=zeros(Msize(1),Msize(2));
    for i=1:size(MidsS,1)
        DendMapS(MidsS(i,1),MidsS(i,2))=DendMapS(MidsS(i,1),MidsS(i,2))+Length(i);
    end

    image(DotMapS*10)
    image(DotMapS./DendMapS*60)
end
%}    

%% Find Ray territory
image(DendMap*100), pause(.3)

Disk2=fspecial('disk',5);
TerFilt=imfilter(DendMap,Disk2,'same');

DFlab=bwlabel(TerFilt);
for i = 1:max(DFlab(:))
    lSize(i)=size(find(DFlab==i),1); %#ok<AGROW>
end
TerFilt=DFlab*0;
TerFilt(DFlab==find(lSize==max(lSize)))=1;
image(TerFilt*200),pause(.3)

%Close
Csize=round(AreaS/2);
SE=strel('disk',Csize);
Buf=Csize*2;
[tys txs] = size(TerFilt);
BufT=zeros(tys+2*Buf,txs+2*Buf);
BufT(Buf+1:Buf+tys,Buf+1:Buf+txs)=TerFilt;
BufC=imclose(BufT,SE);
TerC=BufC(Buf+1:Buf+tys,Buf+1:Buf+txs);
image((TerFilt+TerC)*100), pause(.3)

%Remove holes
TerFill=TerC;
Frame=TerC*0; Frame(1:tys,1)=1; Frame(1:tys,txs)=1; Frame(1,1:txs)=1; Frame(tys,1:txs)=1;
TerHole=bwlabel(~TerC);
for i=1:max(TerHole(:))
    Surround=sum(sum((TerHole==i) .* Frame));
    if ~Surround
        TerFill(TerHole==i)=1;
    end        
end
image(TerFill*100),pause(.01)

%Get Perimeter
TerPerim=bwperim(TerFill,8);
image((TerFill+TerPerim)*100),pause(.3)

%DrawRays
[tpy tpx]=find(TerPerim);
TerRay=TerFill*0;
Cent=Dots.Im.CBpos;
Cent=Cent*.103;
for i = 1: size(tpy,1)
    ydif=Cent(1)-tpy(i); xdif=Cent(2)-tpx(i);
    Long=sqrt(ydif^2 + xdif^2);
    ystep=ydif/Long;
    xstep=xdif/Long;
    
    for l = 1:.5:Long
        TerRay(round(tpy(i)+ystep*l),round(tpx(i)+xstep*l))=1;
    end
end
image(TerRay*100),pause(.01)

DotDist=DotSkel.*TerRay;
image(DotDist*500)
     