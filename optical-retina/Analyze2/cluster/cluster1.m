%% Find clustering
%Concider both absolute and dendritic distances. 
'start'

%clear all
%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions


%% Get Folder and Load Critical Data
if exist('.\temp\Last.mat')
     load(['.\temp\Last.mat'])
     TPN=uigetdir(Last)
else
    TPN=uigetdir
end

Last=TPN;
if Last>0
save('.\temp\Last.mat','Last')
end

%% load segments
'loading data'
load([TPN '\data\AllSeg.mat'])
Mids=mean(AllSeg,3);  %find midpoints
for i=1:size(Mids,1)  %get segment lengths
   Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end
Nodes=[AllSeg(:,:,1) ; AllSeg(:,:,2) ; Mids]; %list all nodes (tips and midpoints)
load([TPN '\data\DotStats.mat'])
Dots=DotStats(:,:,3);
load([TPN '\data\D.mat']) %Load dendrite mask
[Dy Dx Dz]=find3(D); % vectorize by find 3
Dy=Dy*xyum; Dx=Dx*xyum; Dz=Dz*zum; %scale to microns
Dv=[Dy Dx Dz]; clear Dy Dx Dz % Concatinate D find results into vector


DBins=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200];
%DBins=[5,20];

%% run all mid points
'Running Segments'
for i = 1: size(Mids,1)  
    DotDist=dist(Dots,Mids(i,:));
    MidDist=dist(Mids,Mids(i,:));
    %VolDist=dist(Dv,Mids(i,:));
  
        
    %% Run all Bins
    for b = 1:size(DBins,2)
        LocalDots(i,b)=sum(DotDist<=DBins(b));
        LocalLength(i,b)=sum(Length(MidDist<=DBins(b)));
        %LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
    end
    if mod(i,100)==0, PercentDone=i/size(Mids,1)*100,end
end


%% Calculate Volumes for each Bin
VolumeOn=0;
if VolumeOn == 1

    for i = 1: size(Mids,1)  
        VolDist=dist(Dv,Mids(i,:));        
        %% Run all Bins
        for b = 1:size(DBins,2)
            LocalVol(i,b)=sum(VolDist<=DBins(b))*(xyum *xyum*zum);
        end
        if mod(i,100)==0, PercentDone=i/size(Mids,1)*100,end
    end

end %if Volume on


%% Display
for i = 1: size(Mids,1)
    for b = 1 : size(DBins,2)
        fieldLD(round(Mids(i,1))+1,round(Mids(i,2))+1,b)=LocalDots(i);
        fieldDD(round(Mids(i,1))+1,round(Mids(i,2))+1,b)=LocalDots(i,b)/LocalLength(i,b);
    end    
    
end

%% Play movie
colormap gray
while 1==1
for b = 1:size(DBins,2)
    field = zeros(size(fieldDD,1),size(fieldDD,2),3);
    maxDD=max(max(fieldDD(:,:,b)));
    frame=fieldDD(:,:,b);
    meanDD=mean(mean(frame(frame>0)));
    scaleDD=155/meanDD;
    field(:,:,1)=uint8(double(fieldDD(:,:,b))*scaleDD);
    field=uint8(field);
    image(field*3)
    b
    pause
end
end




'finished'
%%Notes
%{
Could segment each bin
%}

