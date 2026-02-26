%%Do Everything

'start'
clear all
tic
colormap gray(255) %standard grey colormap

%% Get file names

[DFN,DPN]=uigetfile('*.tif','DialogTitle','Choose first Image of Data Stack')

%Get directory name
f=find(DPN=='\')
f2=f(size(f,2)-1)
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)

if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
if isdir('./historyT')==0, mkdir('./historyT'); end %create directory to store steps

save(['./historyT' TPN(f3:f2-1)],'TPN') %record path in history folder

%%Enter Variables

%%Image Variables
xyum=.103;
zum=.3;
aspect=zum/xyum;% ratio of z to xy dimentions
live=1; %if live tissue then 1, 2 prevents threshold climbing in fixed tissue
save([TPN 'data\ImageName.mat'],'DPN')

global DPN TPN

AnaGui