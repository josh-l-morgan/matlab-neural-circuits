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

%% Load Neighbor and mean density
load([TPN '\data\Cell.mat'])