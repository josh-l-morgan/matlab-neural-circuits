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


%% Get settings automatically
if str2num(answer{1})>1
    
'finding signal ratios',pause(.01)
med1=median(channel1(:));
med2=median(channel2(:));
std1=uint16(std(single(channel1(:))));
std2=uint16(std(single(channel2(:))));

Thresh1=med1+2*std1;
Thresh2=med2+2*std2

Cht1=(channel1>Thresh1) & (channel1 ~= max(channel1(:)));
Cht2=(channel2>Thresh2) & (channel2 ~= max(channel2(:)));
Cht12=(Cht1 | Cht2)& ~(Cht1 & Cht2);
Cht12 = Cht12 & (channel1 < max(channel1(:))-max(channel1(:))/10) & (channel2 < max(channel2(:))-max(channel2(:))/10);

% image(sum(Cht2,3)*30);

BackC1=mean(channel1((channel1<(med1+std1))&(channel2<(med2+std2))));
BackC2=mean(channel2((channel1<(med1+std1))&(channel2<(med2+std2))));


%Rat=single(channel1-BackC1)./single(channel2-BackC2);

Rat=single(channel1-BackC1)./(single(channel2-BackC2)+single(channel1-BackC1));
Rats=Rat(Cht12);

if str2num(answer{1})==2 %% find by fitting gaussians
    Rats=Rats((Rats>0)&(Rats<1));
    options = statset('Display','final');
    RepK=1;
    obj = gmdistribution.fit(Rats,2,'Options',options,'Replicates',RepK);
    RatMeans=obj.mu;
elseif  str2num(answer{1})==3 %find my taking extreames
    SortRat=sort(Rats);
    'start hist'
    [HRat,Xs]=hist(Rats,0:.02:1.02);
    hist(Rats,0:.02:1.02);
    LHrat=length(HRat);
    Xs(find(HRat(1:LHrat/2)==max(HRat(1:LHrat/2)),1))
    Xs(find(HRat(LHrat/2:LHrat)==max(HRat(LHrat/2:LHrat)),1)+LHrat/2)
    LRat=length(SortRat);
    HighRat=median(SortRat(LRat-fix(LRat/3):LRat));
    LowRat=median(SortRat(1:fix(LRat/3)));
    RatMeans=[HighRat;LowRat]; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else

    hist(Rats,0:.02:1);
    while 1
    'Select two peaks and click return'
    [x y]=ginput;
    if length(x)==2,break,end
    
    end
    RatMeans=x; 
end
% Ch2r=1./(RatMeans+1);
% Ch1r=1-Ch2r;
Ch1r=RatMeans;
Ch2r=1-Ch1r;


signal1 = [Ch1r(1) Ch2r(1)];
signal2 = [Ch1r(2) Ch2r(2)];
background = [BackC1 BackC2]; 

end



%% UNMIXING
% 'unmixing'
% signal1 = (rawSignal1 - background)./...
%     (rawSignal1 + rawSignal2 - 2*background);
% 
% signal2 = (rawSignal2 - background)./...
%     (rawSignal1 + rawSignal2 - 2*background);
'unmixing' ,pause(.1)



A = [signal1' signal2'];
b = [channel1(:)'- background(1); channel2(:)'-background(2)];
x = A\single(b);

S1 = zeros(ys, xs, zs);
S1(:)=x(1,:);

S2 = zeros(ys, xs, zs);
S2(:)=x(2,:);
clear x


%% Scale images

S1=S1-mode(single(S1(:)));
S1=S1*255/max(S1(:));
S1=uint8(S1);

S2=S2-mode(single(S2(:)));
S2=S2*255/max(S2(:));
S2=uint8(S2);

clear C

C(:,:,1)=max(S1,[],3);
C(:,:,2)=max(S2,[],3);
C(:,:,3)=C(:,:,2)*0;
% subplot(2,1,2)
image(C)

%% Write image
'writing image', pause(.1)

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




A



