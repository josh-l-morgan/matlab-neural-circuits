%% Seperate signals from two sources that are mixed into two channels. 
%%written by Josh Morgan and Daniel Kerschensteiner


%% Get image
clear all

colormap gray(256)


'getting image info', pause(.1)
prompt = {'Enter file type [3D tiff = 1, 2D RGBs in folder = 2, OIF = 3]','First unmixed Channel:',...
    'Second unmixed Channel:'};
title = 'Type of image';
nLines = 1;

ImageInfo= inputdlg(prompt,title,nLines,{'2','1','2'});






%% INPUTDLG GET RATIOS
%  RedSig=[80 116];
%  GreenSig=[33 87];
%  BackSig=[18 7];
% 
'getting inputs', pause(.1)
prompt = {'input values? [1 = manual, 2 = gaussian, 3 = extremes, 4 = semiautomated]',...
    'Signal 1 [Ch1 Ch2] (R G)', 'Signal 2 [Ch1 Ch2] (R G)',...
    'Background [Ch1 Ch2] (R G)', 'Kernel Size [] (1 1 or 2 2)'};
title = 'Input Average Signal Intensities';
nLines = 1;

if exist('.\PreviousSettings.mat')
    load('.\PreviousSettings.mat')
    defAns = {'4',PreviousSettings.rS1, PreviousSettings.rS2, PreviousSettings.bg,PreviousSettings.kS};
else
    defAns = {'','', '', '','1 1'};
end
answer = inputdlg(prompt,title,nLines,defAns);

rawSignal1 = str2num(answer{2}); %#ok<ST2NM>
rawSignal2 = str2num(answer{3}); %#ok<ST2NM>
background = str2num(answer{4}); %#ok<ST2NM>
kernelSize = str2num(answer{5}); %#ok<ST2NM>

signal1=rawSignal1-background;
signal2=rawSignal2-background;
signal1=signal1/sum(signal1);
signal2=signal2/sum(signal2);

PreviousSettings.rS1=num2str(rawSignal1);
PreviousSettings.rS2=num2str(rawSignal2);
PreviousSettings.bg=num2str(background);
PreviousSettings.kS=num2str(kernelSize);

save('.\PreviousSettings.mat','PreviousSettings')
pause(.01)

%% Read Image

'reading tiffs'
if str2num(ImageInfo{1}) == 2
    RPN = GetMyDir;
    RPNd=dir(RPN);RPNd=RPNd(3:length(RPNd));
    c=0;
    for i = 1:length(RPNd)
       Nam = RPNd(i).name;
       if sum(Nam(length(Nam)-3:length(Nam))=='.tif')==4
           c=c+1;
           if c==1
              Is=imread([RPN Nam]);
              [ys,xs,cs]=size(Is);
              I=zeros(ys,xs,cs,length(RPNd),'uint16');
              I(:,:,:,c)=Is;
           else
            I(:,:,:,c)=imread([RPN Nam]);
           end
       end
       
    end
    I=I(:,:,:,1:c);
    Isum=find(sum(sum(sum(I,1),2),4));
    if length(Isum)==2
        channel1 = squeeze(I(:,:,Isum(1),:));
        channel2 = squeeze(I(:,:,Isum(2),:));
    else
        channel1 = squeeze(I(:,:,str2num(ImageInfo{2}),:));
        channel2 = squeeze(I(:,:,str2num(ImageInfo{3}),:));

    end
    clear I
    
    TPN=[RPN(1:length(RPN)-1) '_Unmix'];
       if isdir(TPN), rmdir(TPN,'s'), end
    mkdir(TPN)
    name='unmix';
    
elseif str2num(ImageInfo{1}) == 3
    RPN = GetMyDir;
    I=oifread(RPN);
    Isum=find(sum(sum(sum(I,1),2),4));
    if length(Isum)==2
        channel1 = squeeze(I(:,:,Isum(1),:));
        channel2 = squeeze(I(:,:,Isum(2),:));
    else
        channel1 = squeeze(I(:,:,str2num(ImageInfo{2}),:));
        channel2 = squeeze(I(:,:,str2num(ImageInfo{3}),:));

    end
    clear I
    
    TPN=[RPN(1:length(RPN)-1) '_Unmix'];
    if isdir(TPN), rmdir(TPN,'s'), end
    mkdir(TPN)
        name='unmix';
else
    
    [DFN DPN] = uigetfile('.tif');
    TargFile=[DPN DFN];

    Iall = tiffread2(TargFile);
    for i = 1: size(Iall,2)
    I(:,:,i)=Iall(i).data;
    end
    clear Iall

    channel2 = I(:,:,size(I,3)/2+1:size(I,3));
    channel1 = I(:,:,1:size(I,3)/2);
    
    name=DFN(1:find(DFN=='.')-1);
    TPN=[DPN  name ];
    name=DFN(1:find(DFN=='.')-1);


end

channel1=uint16(double(channel1)*255/max(double(channel1(:))));
channel2=uint16(double(channel2)*255/max(double(channel2(:))));

Ch1max=single(max(channel1,[],3));
Ch2max=single(max(channel2,[],3));
Imax=Ch1max*256/max(Ch1max(:)); 
Imax(:,:,2)=Ch2max*256/max(Ch2max(:));
Imax(:,:,3)=Ch2max*256/max(Ch2max(:))*0;
% subplot(2,1,1)
image(uint8(Imax)),pause(.01)
clear I 

[ys xs zs]=size(channel2);
