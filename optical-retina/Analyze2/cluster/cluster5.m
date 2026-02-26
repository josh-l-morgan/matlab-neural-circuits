%% Find clustering
%Concider both absolute and dendritic distances. 
'start'

clear all
%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions


%% Get Folder and Load Critical Data
if exist('.\temp\Last.mat')
     load(['.\temp\Last.mat'])
     if exist(Last)
        TPN=uigetdir(Last)
     else
         TPN=uigetdir
     end
else
    TPN=uigetdir
end

Last=TPN;
if Last>0
save('.\temp\Last.mat','Last')
end

%% load segments
'loading data'

load([TPN '\data\D.mat']) %Load dendrite mask
%sumD=sum(D,3);
%maskD=sumD>0;
%ID=maskD;

[Dy Dx Dz]=find3(D); % vectorize by find 3
Dy=Dy*xyum; Dx=Dx*xyum; Dz=Dz*zum; %scale to microns
Dv=[Dy Dx Dz]; clear Dy Dx Dz % Concatinate D find results into vector

clear D

load([TPN '\data\AllSeg.mat'])
Mids=mean(AllSeg,3);  %find midpoints
for i=1:size(Mids,1)  %get segment lengths
   Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end
Nodes=[AllSeg(:,:,1) ; AllSeg(:,:,2) ; Mids]; %list all nodes (tips and midpoints)
load([TPN '\data\DotStats.mat'])
Dots=DotStats(:,:,3);



DBins=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200];
%DBins=[5,20];

%% create Nearest Node list (nearest node for each dot and distance)
for i = 1:size(Dots,1)
    Ndist=dist(Nodes,Dots(i,:)); %find dist from dot to all nodes
    Near=min(Ndist); %find shortest distance
    Nearest=find(Ndist==Near,1); %get node at that distance
    NN(i,:)=Nodes(Nearest,:); %assign that node to NearestNode list for dots
    DotToNN(i,:)=Near; %record that distance for posterity
    
end



%% run all mid points
'Running Segments'
MidStep=1;
Mid2=Mids(1:MidStep:size(Mids,1),:);

for i = 1: size(Mid2,1)  
    DotDist=dist(Dots,Mid2(i,:));
    MidDist=dist(Mids,Mid2(i,:));
    NNDist=dist(NN,Mid2(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        LocalDots(i,b)=sum(DotDist<=DBins(b));
        LocalLength(i,b)=sum(Length(MidDist<=DBins(b)));
        LocalN(i,b)=sum(NNDist<=DBins(b));
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Mid2,1)*100,end
end
LocalDD=LocalDots./LocalLength;
LocalND=LocalN./LocalLength;
%% run from prespective of dots

'Running Segments'
for i = 1: size(Dots,1)  
    DotToDotDist=dist(Dots,Dots(i,:));
    DotToMidDist=dist(Mids,Dots(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        DotToDots(i,b)=sum(DotToDotDist<=DBins(b));
        DotToLength(i,b)=sum(Length(DotToMidDist<=DBins(b)));
        DotToLength(i,b)=max(DotToLength(i,b),DBins(b)); %set minimum length
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Dots,1)*100,end
end



DotToDD=DotToDots./DotToLength;

%}



%% Calculate Volumes for each Bin

'Checking Volume'
VolumeOn=0;
if VolumeOn == 1

    for i = 1: size(Mid2,1)  
        VolDist=dist(Dv,Mid2(i,:));        
        %% Run all Bins
        for b = 1:size(DBins,2)
            LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
        end
        if mod(i,10)==0, DendPercentDone=i/size(Mid2,1)*100,end
    end



%save('.\temp\LocalVol.mat','LocalVol')

%%Get SurfaceArea
%%assuming cylinders
%a=sqrt(Vol/(L*pi))*L*2*pi;
LocalA=sqrt(LocalVol./(LocalLength*pi)).* (LocalLength * pi * 2);
minArea=LocalLength * pi *2 * .05;  % conciders minimum area for length as 0.1um diam process
LocalA(LocalA<minArea)=minArea(LocalA<minArea); %area is at least minimum for tiny process
LocalDA=LocalDots./LocalA;


%%Run Volume for Dots

    for i = 1: size(Dots,1)  
        VolDist=dist(Dv,Dots(i,:));        
        %% Run all Bins
        for b = 1:size(DBins,2)
            DotToVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
        end
        if mod(i,10)==0, DotPercentDone=i/size(Dots,1)*100,end
    end



%save('.\temp\LocalVol.mat','LocalVol')

%%Get SurfaceArea
%%assuming cylinders
%a=sqrt(Vol/(L*pi))*L*2*pi;
DotToA=sqrt(DotToVol./(DotToLength*pi)).* (DotToLength * pi * 2);
minArea=DotToLength * pi *2 * .05;  % conciders minimum area for length as 0.1um diam process
DotToA(DotToA<minArea)=minArea(DotToA<minArea); %area is at least minimum for tiny process
DotToDA=DotToDots./DotToA;




Local.Mid.LocalVol=LocalVol;
Local.Mid.LocalA=LocalA;
Local.Mid.LocalDA=LocalDA;

Local.Dot.DotToVol=DotToVol;
Local.Dot.DotToA=DotToA;
Local.Dot.DotToDA=DotToA;

end %if Volume on


%% Get some useful numbers
hist(LocalDD(:,3),0:.02:3)
pause(1)
hist(DotToDD(:,3),0:.05:3)
hist(DotToDots(:,3),0:.02:10)

for i=1:size(DBins,2)-1
   diffAll(:,i)=abs(LocalDD(:,i)-LocalDD(:,i+1)); 
end
MeanDiffs=mean(diffAll,1)

Local.DBins=DBins;
Local.Mid.MidStep = MidStep;

Local.Dot.DotToDots=DotToDots;
Local.Dot.DotToLength=DotToLength;
Local.Dot.DotToDD=DotToDD;


Local.Mid.LocalDots=LocalDots;
Local.Mid.LocalLength=LocalLength;
Local.Mid.LocalDD=LocalDD;

%%Nearest Neighbor percent
clear Neighbors

for b= 1: 10
    for n=1:10
        Neighbors(b,n)=sum(DotToDots(:,b)>n)/size(DotToDots,1)*100;
    end
    plot(Neighbors(b,:)),hold on
end
hold off 
Local.Dot.Neighbors=Neighbors;


%save([TPN '\data\Local.mat'],'Local')

%% Display

LocalPlot=LocalDD; %decide who to plot

for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldLD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)=LocalDots(i);
    end        
end
fieldC=fieldLD*0; %make field count
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)= fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)+1;
    end        
end
fieldDD=fieldLD*0;
for i = 1: size(Mid2,1)
    for b = 1 : size(DBins,2)
        fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)=fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b)+LocalPlot(i,b)/fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,b);
    end        
end


%%Make Dot map
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

