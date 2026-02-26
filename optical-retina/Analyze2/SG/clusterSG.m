%% Find clustering
%Concider both absolute and dendritic distances. 
'start'

clear all
%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions

TPN = GetMyDir

load([TPN '\data\D.mat']) %Load dendrite mask
load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])
DPos=Dots.Pos(SG.passF,:);
DPos(:,1:2)=DPos(:,1:2)*.103; DPos(:,3)=DPos(:,3)*.3;
load([TPN 'data\AllSegCut.mat'])
AllSeg=AllSegCut; clear AllSegCut

%% load segments

[Dy Dx Dz]=find3(D); % vectorize by find 3
Dy=Dy*xyum; Dx=Dx*xyum; Dz=Dz*zum; %scale to microns
Dv=[Dy Dx Dz]; clear Dy Dx Dz % Concatinate D find results into vector
clear D Dy Dx Dz

Mids=mean(AllSeg,3);  %find midpoints
for i=1:size(Mids,1)  %get segment lengths
   Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end
Nodes=[AllSeg(:,:,1) ; AllSeg(:,:,2) ; Mids]; %list all nodes (tips and midpoints)


DBins=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200];
%DBins=[5,20];

%% create Nearest Node list (nearest node for each dot and distance)
NN=zeros(size(DPos,1),3);
for i = 1:size(DPos,1)
    Ndist=dist(Mids,DPos(i,:)); %find dist from dot to all nodes
    Near=min(Ndist); %find shortest distance
    Nearest=find(Ndist==Near,1); %get node at that distance
    NN(i,:)=Mids(Nearest,:); %assign that node to NearestNode list for dots
    DotToNN(i,:)=Near; %record that distance for posterity
   
end

%Eliminate extreme displacements
NN=NN(DotToNN<=3,:);

%% run all mid points
'Running Segments'
MidStep=1;
Mid2=Mids;%(1:MidStep:size(Mids,1),:);

for i = 1: size(Mid2,1)  
    %DotDist=dist(Dots,Mid2(i,:));
    MidDist=dist(Mids,Mid2(i,:));
    NNDist=dist(NN,Mid2(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        %LocalDots(i,b)=sum(DotDist<=DBins(b));
        LocalLength(i,b)=sum(Length(MidDist<=DBins(b)));
        LocalN(i,b)=sum(NNDist<=DBins(b));
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Mid2,1)*100,end
end
%LocalDD=LocalDots./LocalLength;
LocalND=LocalN./LocalLength;



%% Display
Mid2=round(Mid2);
LocalDots=LocalND;
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldLD(Mid2(i,1)+1,Mid2(i,2)+1,b)=LocalDots(i,b);
    end        
end
fieldC=fieldLD*0; %make field count
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldC(Mid2(i,1)+1,Mid2(i,2)+1,b)= fieldC(Mid2(i,1)+1,Mid2(i,2)+1,b)+1;
    end        
end
fieldDD=fieldLD*0;
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldDD(Mid2(i,1)+1,Mid2(i,2)+1,b)=fieldDD(Mid2(i,1)+1,Mid2(i,2)+1,b)+LocalND(i,b);
    end        
end
fieldDD(fieldC>0)=fieldDD(fieldC>0)./fieldC(fieldC>0);


%% Make Dot map
fieldDP=zeros(fix(max(Dots(:,1)))+2,fix(max(Dots(:,2))+2));
for i = 1: size(Dots,1)
   fieldDP(round(Dots(i,1))+1,round(Dots(i,2))+1)=fieldDP(round(Dots(i,1))+1,round(Dots(i,2))+1)+1;
end
mask=max(fieldDD,[],3)>0;
xs=max(size(mask,2),size(fieldDP,2));
ys=max(size(mask,1),size(fieldDP,2));
CellPic=zeros(ys,xs,3);
CellPic(1:size(mask,1),1:size(mask,2),1)=mask;
CellPic(1:size(fieldDP,1),1:size(fieldDP,2))=fieldDP;
CellPic(:,:,1)=mask*20;
CellPic(:,:,3)=CellPic(:,:,3)*75;
CellPic=uint8(CellPic);
subplot(6,6,[1:5,7:11,13:17,19:23,25:29])
image(CellPic*10)



%% Play movie
r=1:255;b=255:-1:1;
r=r/255; b=b/255;
rb=r'; rb(:,3)=b';
rb(2:size(rb,1)+1,:)=rb;
rb(1,:)=0;


mask=max(fieldDD,[],3)>0;
image(mask*100)
load('.\temp/cmap.mat')
set(gcf,'Colormap',cmap)
pause(.01)
clear field

meanBright=mean(LocalPlot(:,20));
Times=3;

while 1==1
for b = 1:size(DBins,2)
    Nine=sort(LocalPlot(:,b));
    Ninety=Nine(round(.95*size(Nine,1)));
    scaleDD=255/(meanBright*Times); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    field=fieldDD(:,:,b);
    field=field*scaleDD;
    field(mask)=field(mask)+2;
    subplot(6,6,[1:5,7:11,13:17,19:23,25:29])
    
    image(field)
    subplot(6,6,31:35)
    Diam=([1:size(field,2)]<=DBins(b)*2)*1000;
    Diam(:,:,2)=Diam*255;     Diam(:,:,3)=Diam(:,:,1)*255;
    image(uint8(Diam))
    subplot(6,12,[12,24,36,48])
    Inten=1:256;
    image(Inten')
    htext = uicontrol('Style','text','String',Times);
    DBins(b)
    pause
end
end




'finished'
%%Notes
%{
Could segment each bin
individual variation in density as bins increase
montecarlo for nearest node
run for surface area (and save)
how good is area assumption for multiple diameters

%}

