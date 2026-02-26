function[]=anaFMo(TPN, DPN)
%% Apply Mask to Data
%%Recieves mask in the form of 0=exterior, 1=dend (to which dots are
%%related), 2=crap (from which sements are removed, 3= Cell body

%{
%Get directory name
DPN=GetMyDir
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)
%}
    


%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'reading image'

if exist([TPN 'mask']) & exist([TPN 'Dots.mat'])%% mask is image in folder 'mask' if exists
        d=dir([TPN 'mask']); %get number of files in directory
        d=d(3:size(d,1));

        clear I IM 
        IM=imread([TPN 'mask\' d(1).name]);
        [ys xs zs] = size(IM);
        IM=zeros(ys, xs, size(d,1),'uint8');
        for im=1:size(d,1)
            IM(:,:,im)=imread([TPN 'mask\' d(im).name]);
            PercentRead=im/size(d,1)*100
        end

    maxMask=max(IM,[],3);
    MaxOfMask=max(maxMask(:));
    maxMask=maxMask*(255/MaxOfMask);
    imwrite(maxMask,[TPN 'images\maxMask.tif'],'Compression','none')
    clear maxMask

%% Load Dots and Dendrites
    load([TPN 'Dots.mat'])
    yxum=.103; zum=.3;  %% Assign voxel dimensions
    DotPos=Dots.Pos;
    DotPos(:,1:2)=DotPos(:,1:2)*yxum;
    DotPos(:,3)=DotPos(:,3)*zum;



%% Check Dend against IM Mask pixel by pixel
    %%Segment tips rounded to the nearest pixel should be within a Mask Pixel
    %%(concider dialating mask)

    load([TPN 'data\AllSeg.mat'])
    load([TPN 'ShiftInfo.mat'])
    AllSeg(:,3,:)=AllSeg(:,3,:)+ShiftInfo.ShiftDend*zum;
    %%Save Segments

    SegCheck=ones(size(AllSeg,1),1);
    if max(IM(:))>1 %if there is a crap region

        'Checking all dendrite Segments'
        
            
        SegMid=mean(AllSeg,3);
        SegMid(:,1:2)=SegMid(:,1:2)/yxum;
        SegMid(:,3)=SegMid(:,3)/zum;
        SegMid=round(SegMid);
        SegMid(SegMid(:,3)>size(IM,3),3)=size(IM,3);
        SegMid(SegMid(:,2)>size(IM,2),1)=size(IM,2);
        SegMid(SegMid(:,1)>size(IM,1),1)=size(IM,1);
        
        
        %%is either end of a segment within a crap region (>1)
        for s=1:size(SegMid,1)
            SegCheck(s)=~(IM(SegMid(s,1),SegMid(s,2),SegMid(s,3))>1);
        end
        save([TPN 'data\SegCheck.mat'],'SegCheck')
    end

    AllSegCut=AllSeg(logical(SegCheck),:,:);
    save([TPN 'data\AllSegCut.mat'],'AllSegCut')


%% Find out what part of the mask 0, 1, 2, 3 the puncta lays on
    clear Cut
    for i = 1: Dots.Num
       Cut(i)=IM(round(Dots.Pos(i,1)),round(Dots.Pos(i,2)),round(Dots.Pos(i,3))); 
    end
    Dots.Cut=Cut;

    save([TPN 'Dots.mat'],'Dots');

    end  %if mask folder exists



