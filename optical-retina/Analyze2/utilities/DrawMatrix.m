%% Grab a mat to write as stack of 2D tiffs

[DFN DPN] = uigetfile
load([DPN DFN])
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)
Name=DFN(1:size(DFN,2)-4)
%must change name
%imwriteNp(TPN,Name,Name)