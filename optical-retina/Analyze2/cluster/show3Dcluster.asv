%Make 3D cluster map

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

if exist([TPN '\data\Local.mat'])
    load([TPN '\data\Local.mat'])
else
    'yer shit out of luck. run it yourself'
end


%%Load Skeleton Segments and extract midpoints and Lengths
load([TPN '\data\AllSeg.mat'])
Mid2=mean(AllSeg,3);  %find all  midpoints
for i=1:size(Mid2,1)  %get segment lengths
           Length(i)=dist(AllSeg(i,:,1),AllSeg(i,:,2)); 
end

LocalDD=Local(1).Mid.LocalDD;
LocalDots=Local(1).Mid.LocalDots;

%% Make 3D map
LocalPlot=LocalDD; %decide who to plot

for i = 1: size(Mid2,1)
        field3dLD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1)=LocalDots(i);
 end
field3dC=field3dLD*0; %make field count
for i = 1: size(Mid2,1)
        field3dC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1)= field3dC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1)+1;
    
end
field3dDD=field3dLD*0;
for i = 1: size(Mid2,1)
        field3dDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1)=field3dDD(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1)+LocalPlot(i,5)/field3dC(round(Mid2(i,1))+1,round(Mid2(i,2))+1,round(Mid2(i,3))+1);    
end

meanBright=Local(1).Mid.meanBright;
field3dDD=field3dDD * (255/(meanBright*3));
imwriteNp([TPN '\'],field3dDD,'DD3d')
