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

%% load info from cell
'loading data'

%%Load Skeleton Segments and extract midpoints and Lengths
load([TPN '\data\AllSeg.mat'])
Mids=mean(AllSeg,3);  %find midpoints
for i=1:size(Mids,1)  %get segment lengths
   Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end
Nodes=[AllSeg(:,:,1) ; AllSeg(:,:,2) ; Mids]; %list all nodes (tips and midpoints)

%%Load Dots
load([TPN '\data\DotStats.mat'])
Dots=DotStats(:,:,3);


%%Create List of search radii
DBins=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200];

%%Load Cell information 
load([TPN '\Cell.mat'])
load([TPN '\data\Results.mat'])


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
    %if mod(i,100)==0, PercentDone=i/size(Mid2,1)*100,end
end
LocalDD=LocalDots./LocalLength;
LocalND=LocalN./LocalLength;

%% run from prespective of dots

'Running Dots'
for i = 1: size(Dots,1)  
    DotToDotDist=dist(Dots,Dots(i,:));
    DotToMidDist=dist(Mids,Dots(i,:));
    %% Run all Bins
    for b = 1:size(DBins,2)
        DotToDots(i,b)=sum(DotToDotDist<=DBins(b));
        DotToLength(i,b)=sum(Length(DotToMidDist<=DBins(b)));
        DotToLength(i,b)=max(DotToLength(i,b),DBins(b)); %set minimum length
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    %if mod(i,100)==0, PercentDoneWithDots=i/size(Dots,1)*100,end
end
DotToDD=DotToDots./DotToLength;

%}

%% Run XY map of vis space

%%Run by convolution
MidMap=zeros(round(max(Mids(:,1))),round(max(Mids(:,2))));
for i=1: size(Mids,1)
    x=max(1,round(Mids(i,2)));
    y=max(1,round(Mids(i,1)));
    MidMap(y,x)=MidMap(y,x)+Length(i);
end

DotMap=zeros(round(max(Dots(:,1))),round(max(Mids(:,2))));
for i = 1: size(Dots,1)
    x=max(1,round(Dots(i,2)));
    y=max(1,round(Dots(i,1)));
    DotMap(y,x)=DotMap(y,x)+1;
end




%% Record some useful numbers
Local.Mid.MidMap=MidMap;
Local.Dot.DotMap=DotMap;

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
save([TPN '\data\Local.mat'],'Local')


%% Display

LocalPlot=LocalDD; %decide who to plot
BinsUsed=[5,11];


for i = 1: size(Mid2,1)
    for bb = 1:size(BinsUsed,2)
        b=BinsUsed(bb);
        fieldLD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb)=LocalDots(i);
    end        
end
fieldC=fieldLD*0; %make field count
for i = 1: size(Mid2,1)
    for bb = 1:size(BinsUsed,2)
        b=BinsUsed(bb);
        fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb)= fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb)+1;
    end        
end
fieldDD=fieldLD*0;
for i = 1: size(Mid2,1)
    for bb = 1:size(BinsUsed,2)
        b=BinsUsed(bb);
        fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb)=fieldDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb)+LocalPlot(i,b)/fieldC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,bb);
    end        
end



%% Write Images
image(DotMap*50)
load('.\temp/cmap.mat')
set(gcf,'Colormap',cmap)
clear field

meanBright=mean(LocalPlot(:,2));
Times=3;
Local.Mid.meanBright=meanBright;
Target='\\128.208.64.36\wonglab\Josh\Analyzed\ClusterPics\'
mask=max(fieldDD,[],3)>0;
for b = size(BinsUsed)
    scaleDD=255/(meanBright*Times); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    field=fieldDD(:,:,b);
    field=field*scaleDD;
    field(mask)=field(mask)+2;
    name=[Cell.Name '_r' num2str(DBins(BinsUsed(b))) '.tif']
    imwrite(field,cmap,[Target name])
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
    pause(.01)
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

