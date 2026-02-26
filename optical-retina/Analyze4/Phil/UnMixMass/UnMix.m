%% Seperate signals from two sources that are mixed into two channels. 
%%written by Josh Morgan and Daniel Kerschensteiner


%% Get image
clear all

colormap gray(256)
[DFN DPN] = uigetfile;
TargFile=[DPN DFN];

'reading tiffs'
Iall = tiffread2(TargFile);
for i = 1: size(Iall,2)
I(:,:,i)=Iall(i).data;
end
clear Iall

channel2 = I(:,:,size(I,3)/2+1:size(I,3));
channel1 = I(:,:,1:size(I,3)/2);
clear I 
[ys xs zs]=size(channel2);


%% INPUTDLG GET RATIOS
%  RedSig=[80 116];
%  GreenSig=[33 87];
%  BackSig=[18 7];
% 
'getting inputs', pause(.1)
prompt = {'Signal 1 [Ch1 Ch2] (R G)', 'Signal 2 [Ch1 Ch2] (R G)',...
    'Background [Ch1 Ch2] (R G)', 'Kernel Size [] (1 1 or 2 2)'};
title = 'Input Average Signal Intensities';
nLines = 1;

if exist('.\PreviousSettings.mat')
    load('.\PreviousSettings.mat')
    defAns = {PreviousSettings.rS1, PreviousSettings.rS2, PreviousSettings.bg,PreviousSettings.kS};
else
    defAns = {'', '', '','1 1'};
end
answer = inputdlg(prompt,title,nLines,defAns);

rawSignal1 = str2num(answer{1}); %#ok<ST2NM>
rawSignal2 = str2num(answer{2}); %#ok<ST2NM>
background = str2num(answer{3}); %#ok<ST2NM>
kernelSize = str2num(answer{4}); %#ok<ST2NM>

PreviousSettings.rS1=num2str(rawSignal1);
PreviousSettings.rS2=num2str(rawSignal2);
PreviousSettings.bg=num2str(background);
PreviousSettings.kS=num2str(kernelSize);

save('.\PreviousSettings.mat','PreviousSettings')



%% Median filter
'median filtering',pause(.1)
if kernelSize(1)>1
    for i = 1: size(channel1,3)
       channel1(:,:,i)=medfilt2(channel1(:,:,i),[kernelSize(1) kernelSize(1)]);
    end
end
if kernelSize(2)>1
    for i = 1: size(channel2,3)
       channel2(:,:,i)=medfilt2(channel2(:,:,i),[kernelSize(2) kernelSize(2)]);
    end
end



%% UNMIXING
% 'unmixing'
% signal1 = (rawSignal1 - background)./...
%     (rawSignal1 + rawSignal2 - 2*background);
% 
% signal2 = (rawSignal2 - background)./...
%     (rawSignal1 + rawSignal2 - 2*background);
'unmixing' ,pause(.1)
signal1=rawSignal1-background;
signal2=rawSignal2-background;
signal1=signal1/sum(signal1);
signal2=signal2/sum(signal2);


A = [signal1' signal2'];
b = [channel1(:)'- background(1); channel2(:)'-background(2)];
x = A\single(b);

S1 = zeros(ys, xs, zs);
S1(:)=x(1,:);

S2 = zeros(ys, xs, zs);
S2(:)=x(2,:);
clear x


%% Scale images

S1=S1-median(S1(:));
S1=S1*256/max(S1(:));
S1=uint8(S1);

S2=S2-median(S2(:));
S2=S2*256/max(S2(:));
S2=uint8(S2);

C(:,:,1)=max(S1,[],3);
C(:,:,2)=max(S2,[],3);
C(:,:,3)=C(:,:,2);
image(C)

%% Write image
'writing image', pause(.1)
name=DFN(1:find(DFN=='.')-1);

TPN=[DPN  name ];
if isdir(TPN), rmdir(TPN,'s'), end
mkdir(TPN)

TPN1=[TPN '\Ch1'],mkdir(TPN1)
TPN2=[TPN '\Ch2'],mkdir(TPN2)


NumDig=fix(log10(zs))+1; % get digits
blank=zeros(NumDig,1);
blank=num2str(blank)';

for i=1:zs
    Num=num2str(i);
    Numi=blank;
    Numi(size(Numi,2)-size(Num,2)+1:size(Numi,2))=Num;
    name1=[TPN1 '\' name '_Ch1_' Numi  '.tif'];
    imwrite(S1(:,:,i),name1,'tif','Compression','none')
    name2=[TPN2 '\' name '_Ch2_' Numi  '.tif'];
    imwrite(S2(:,:,i),name2,'tif','Compression','none')
       
    
end








