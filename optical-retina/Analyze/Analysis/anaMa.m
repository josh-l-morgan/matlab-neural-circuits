function[] = anaMa()
% Apply Mask to Data
%%Recieves mask in the form of 0=exterior, 1=dend (to which dots are
%%related), 2=crap (from which sements are removed)

global DFN DPN TPN
colormap gray(255)
'Masking'

%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([TPN 'mask']) %% mask is image in folder 'mask' if exists
    d=dir([TPN 'mask']); %get number of files in directory
    d=d(3:size(d,1));
    
    clear I IM 
    for i=1:size(d,1)
        IM(:,:,i)=imread([TPN 'mask\' d(i).name]);
        PercentRead=i/size(d,1)*100
    end
else %default mask is D
    load([TPN 'data\D.mat']) %load thresholded dendrite arbor
    IM=D>0; %(make sure its in binary)
    clear D
end


yxum=.103; zum=.3;  %% Assign voxel dimensions
clear mask
[y x z]=find3(IM==1); %find all mask
mask=[y x z];
mask(:,1:2)=mask(:,1:2)*yxum;
mask(:,3)=mask(:,3)*zum;



%% Load Dots and Dendrites
if exist([TPN 'dataFix\AllSeg.mat'])
    load([TPN 'dataFix\AllSeg.mat'])
else
    load([TPN 'data\AllSeg.mat'])
end

if exist([TPN 'data\DotStats.mat']),
    load([TPN 'data\DotStats.mat'])
else
    load([TPN 'data\BC.mat'])
    Dots=BC;
end

%%Save Backup for All Seg and Dot Data
if ~exist([TPN 'data\AllSegBackup.mat'])
    save([TPN 'data\AllSegBackup.mat'],'AllSeg')
end
if ~exist([TPN 'data\DotStatsBackup.mat'])
    save([TPN 'data\DotStatsBackup.mat'],'DotStats')
end


%% Check Dots against Mask
%%Find distance between dot and nearest masked pixel
[ys xs zs]=size(IM);

clear Dist mDotToMaskDist
for i =1:size(DotStats,1)
  Dist=dist(mask,DotStats(i,:,3)); %find distance between
  mDotToMaskDist(i)=min(Dist);
  if mod(i,50)==100,PercentDoneWithDots=double(i)/size(DotStats,1)*100,end
  %field=(max(1,    %%%% Dont know what this was supposed to do
end
save([TPN 'mDotToMaskDist.mat'],'mDotToMaskDist')

%change dotdata to those puncta within 2um of mask
DotStats=DotStats(mDotToMaskDist<=2,:,:);
save([TPN 'data\DotStats.mat'],'DotStats')

%% Check Dend against IM Mask pixel by pixel
%%Segment tips rounded to the nearest pixel should be within a Mask Pixel
%%(concider dialating mask)

if max(IM(:))>1 %if there is a crap region (ie, custom mask with 2s)

    'Checking all dendrite Segments'
    
    %convert AllSeg back to pixel form
    AllSegM=AllSeg(:,1:2,:)/yxum;
    AllSegM(:,3,:)=AllSeg(:,3,:)/zum;
    AllSegM=uint16(round(AllSegM));
    
    %%This part is scary and shouldnt be necessary!!!!!!!!!!!!!!
    %% get rid of out of bounds
    AllSegM(:,3,:)=min(AllSegM(:,3,:),zs);
    AllSegM(:,1,:)=min(AllSegM(:,1,:),ys);
    AllSegM(:,2,:)=min(AllSegM(:,2,:),xs);
    
    %%is either end of a segment within a crap region (==2)
    for i=1:size(AllSegM,1)
        SegCheck(i)=~(IM(AllSegM(i,1,1),AllSegM(i,2,1),AllSegM(i,3,1))==2 | IM(AllSegM(i,1,2),AllSegM(i,2,2),AllSegM(i,3,2))==2);
    end
    save([TPN 'data\SegCheck.mat'],'SegCheck')

    %%Save Segments
    AllSeg=AllSeg(SegCheck,:,:);
    save([TPN 'data\AllSeg.mat'],'AllSeg')
end





%% Check Dend against mask (min distance
%%Find distance between dot and nearest masked pixel
%{
for i =1:size(AllSeg,1)
  Dist1=dist(mask,AllSeg(i,:,1)); 
  Dist2=dist(mask,AllSeg(i,:,2));
  Dist=min(Dist1,Dist2);
  
  mDendToMaskDist(i)=min(Dist);
  PercentDoneWithDend=double(i)/size(AllSeg,1)*100
end
save([TPN 'mDendToMaskDist.mat'],'mDendToMaskDist')
%}

[TPN(size(TPN,2)-6:size(TPN,2)-1)]
MaskedAt=uint16(clock)
save([TPN 'data/MaskedAt.mat'],'MaskedAt')

%% Draw Dot and Dend
%{
'Drawing Dots and Dendrites'
xyum=yxum;
Sc=(1/xyum)/2;
DDm=uint8(zeros(round(max(max(AllSeg(:,1,:)))*Sc),round(max(max(AllSeg(:,2,:)))),round(max(max(AllSeg(:,3,:)))*Sc)));

%%Draw Segments

SkelRes=.1;
for i=1:size(AllSeg,1)
        Dist=sqrt((AllSeg(i,1,1)-AllSeg(i,1,2))^2 + (AllSeg(i,2,1)-AllSeg(i,2,2))^2 + (AllSeg(i,3,1)-AllSeg(i,3,2))^2); %find distance
        Length(i)=Dist;
          devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
        for d=1:devs+1
            sy=AllSeg(i,1,1)+((AllSeg(i,1,2)-AllSeg(i,1,1))/devs)*(d-1);
            sx=AllSeg(i,2,1)+((AllSeg(i,2,2)-AllSeg(i,2,1))/devs)*(d-1);
            sz=AllSeg(i,3,1)+((AllSeg(i,3,2)-AllSeg(i,3,1))/devs)*(d-1);
            DDm(round(sy*Sc)+1,round(sx*Sc)+1,round(sz*Sc)+1)=1; %draw Skel
        end
end
clear Dist

%%DrawNodes
for i=1:size(AllSeg,1)
        DDm(round(AllSeg(i,1,1)*Sc)+1,round(AllSeg(i,2,1)*Sc)+1,round(AllSeg(i,3,1)*Sc)+1)=2;
        DDm(round(AllSeg(i,1,2)*Sc)+1,round(AllSeg(i,2,2)*Sc)+1,round(AllSeg(i,3,2)*Sc)+1)=2;
end

%%Draw Dots
for i=1:size(DotStats,1)
    if DotStats(i,3,1)> .5 %Delta F over F filter
        DDm(round(DotStats(i,1,3)*Sc)+1,round(DotStats(i,2,3)*Sc)+1,round(DotStats(i,3,3)*Sc)+1)=3;
    end
end

imwriteNp(TPN,DDm,'DDm')
%}







