function[DPN,DFN,TPN]=GetFolder()
global DPN TPN DFN

%% Get file names


[DFN,DPN]=uigetfile('*.tif','DialogTitle','Choose first Image of Data Stack');

%Get directory name
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)

if isdir([TPN 'temp'])==0, mkdir([TPN 'temp']); end %create directory to store steps
if isdir([TPN 'data'])==0, mkdir([TPN 'data']); end %create directory to store steps
if isdir([TPN 'pics'])==0, mkdir([TPN 'pics']); end %create directory to store steps
if isdir('./history')==0, mkdir('./history'); end %create directory to store steps

save(['./history' TPN(f3:f2-1)],'TPN') %record path in history folder
save([TPN 'data\ImageName.mat'],'DPN')

