%%%%Find colocalization by shifting channells and comparing correlation
%%%%coefficient


%% Get image
clear all




'getting image info', pause(.1)
prompt = {'Enter file type [3D tiff = 1, 2D RGBs in folder = 2, OIF = 3]','First unmixed Channel:',...
    'Second unmixed Channel:'};
title = 'Type of image';
nLines = 1;

ImageInfo= inputdlg(prompt,title,nLines,{'2','1','2'});



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

        channel1 = squeeze(I(:,:,str2num(ImageInfo{2}),:));
        channel2 = squeeze(I(:,:,str2num(ImageInfo{3}),:));

    clear I
%     
%     TPN=[RPN(1:length(RPN)-1) '_Unmix'];
%        if isdir(TPN), rmdir(TPN,'s'), end
%     mkdir(TPN)
%     name='unmix';
    
elseif str2num(ImageInfo{1}) == 3
    RPN = GetMyDir;
    I=oifread(RPN);
    Isum=find(sum(sum(sum(I,1),2),4));

        channel1 = squeeze(I(:,:,str2num(ImageInfo{2}),:));
        channel2 = squeeze(I(:,:,str2num(ImageInfo{3}),:));

    clear I
    
%     TPN=[RPN(1:length(RPN)-1) '_Unmix'];
%     if isdir(TPN), rmdir(TPN,'s'), end
%     mkdir(TPN)
%         name='unmix';
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
    
%     name=DFN(1:find(DFN=='.')-1);
%     TPN=[DPN  name ];
%     name=DFN(1:find(DFN=='.')-1);


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

%% FIND DOTS Green Channel%%

channel1b=channel1;
channel2b=channel2;

for ch = 1:2
    ChName=['channel' num2str(ch)];
    Ig=eval(ChName);
 
Igm=zeros(ys,xs,zs,'single');
for i=1:zs        
    Igm(:,:,i)=medfilt2(Ig(:,:,i),[3,3]); 
end
clear Ig

MaxDot=7^3;
MinDot=3;

thresholdMap = zeros(ys,xs,zs,'uint8');   %set up matrix to sum passed thresholds
% indVect = 1:ys*xs*zs;
maxIntensity = uint8(max(Igm(:)));
Gmode=mode(single(Igm(:)));
for i = uint16(maxIntensity):-1:uint16(Gmode)
    %run thresholds through all relevant intensities
    clear Igl labels
    [Igl,labels] = bwlabeln(Igm>i,6);%label each area to check

    %reduce bitdepth if possible
    if labels<65536
        Igl=uint16(Igl);
    end
    if labels <= 1
        labels =2;
    end
    nPixel = hist(Igl(Igl>0), 1:labels);
    %run all lables
    for p=1:labels
        % Morphology Filter, Puncta size criteria
        pixelIndex = find(Igl==p);
        if nPixel(p) < MaxDot && nPixel(p) > MinDot

        else
            Igl(pixelIndex)=0;
        end
    end
    %%Add all passing labeled objects to thresholdMap
    thresholdMap(Igl>0)=thresholdMap(Igl>0)+1;
    i
end 
disp('iterative threshold done')
clear Igl peakIndex

Channels(:,:,:,ch)=thresholdMap;
end
channel1=Channels(:,:,:,1);
channel2=Channels(:,:,:,ch);
clear Channels

Chan=max(channel1,[],3)*2;
Chan(:,:,2)=max(channel2,[],3)*2;
Chan(:,:,3)=max(channel2,[],3)*2;
Chan=uint8(Chan);
image(Chan)


%% Shift channels

Range=30;
c=0;
clear CCs Vec
for x = -Range: Range
    for y = -Range : Range
    c=c+1;
    StartX1=max(0,0-x); StopX1=min(0,0-x);
    StartY1=max(0,0-y); StopY1=min(0,0-y);
    StartX2=max(0,0+x); StopX2=min(0,0+x);
    StartY2=max(0,0+y); StopY2=min(0,0+y);
    SampCh1=double(channel1(1+StartY1:ys+StopY1,1+StartX1:ys+StopX1));
    SampCh2=double(channel2(1+StartY2:ys+StopY2,1+StartX2:ys+StopX2));
    CC=corrcoef(SampCh1(:),SampCh2(:));
    CCs(c,1)=CC(1,2);
    Vec(c,:)=[y,x];
    end  %% end x
end %% end y

%% Draw Corrcoef
CCxs=Range*2+1;
CCys=CCxs;

iCC=zeros(CCys,CCxs);
for i = 1: size(Vec,1)
   iCC(Vec(i,1)+Range+1,Vec(i,2)+Range+1)=CCs(i); 
end

subplot(2,1,1)
image(iCC*255/double(max(iCC(:))))

surf(iCC)
colormap(jet)
shading interp


%% Get Data
MaxCorrCoef=max(CCs)
locMax=find(CCs==MaxCorrCoef);
vecMax=Vec(locMax,:);
MaxLocation=mean(vecMax,1)
iML=MaxLocation+Range+1;
HalfiCC=iCC>=(MaxCorrCoef/2);
HalfiCC=bwlabel(HalfiCC);
Targ=HalfiCC(iML(1),iML(2));
HalfiCC=HalfiCC==Targ; %select object including peak
PiCC=bwperim(HalfiCC);
OtherHalf=iCC(~HalfiCC);
meanOtherHalf=mean(OtherHalf)
[PiCCy PiCCx]=find(PiCC);
Dists=sqrt((iML(1)-PiCCy).^2+(iML(2)-PiCCx).^2);
DistToHalf=mean(Dists)


%% Plot as distance
aDists=sqrt((MaxLocation(1)-Vec(:,1)).^2+(MaxLocation(2)-Vec(:,2)).^2);

binwidth=1;
bw=binwidth/2;
c=0;
for i = 0:.1:max(aDists)
    c=c+1;
    ToGet=aDists>=(i-bw) & aDists<=(i+bw);
    bDist(c)=i;
    mCC(c)=mean(CCs(ToGet));
end

subplot(2,1,2)
plot(bDist,mCC)

Negs=find(mCC<=0);
if isempty(Negs)
   ZeroCrossing=0;    
else
    ZeroCrossing=min(bDist(Negs));
end

Data.MaxCorrCoef=MaxCorrCoef;
Data.MaxLocation=MaxLocation;
Data.DistToHalf=DistToHalf;
Data.ZeroCrossing=ZeroCrossing;
Data.MinCorr=min(mCC);

Data


