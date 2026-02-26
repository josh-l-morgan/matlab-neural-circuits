%% Automatically draw max and orthoganal projections for any file selected
%%2D images should combine both channels and both Identified Dot and Dend. 


clear all

colormap gray(255)
yxum=.103;  zum=0.3;


%Get directory name
[DFN,DPN]=uigetfile('*.tif','DialogTitle','Pick First Image of Stack');
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)

if ~exist([TPN 'images\arbor']), mkdir([TPN 'images\arbor']); end

load([TPN 'Cell.mat'])
load([TPN 'data\Results.mat'])
    


%% Find Z starts and stops
clear Zstart Zstop
if isfield(Results,'Arbor')
   Zstart(1)=Results.Arbor(1).Top;
   Zstop(1)=Results.Arbor(1).Bottom;
   if size(Results.Arbor,2)>1
       Zstart(2)=Results.Arbor(2).Top;
       Zstop(2)=Results.Arbor(2).Bottom;
    end
end
a=1;  %Which arbor to draw

%% Draw Dend

clear DrawSeg
if exist([TPN 'data\AllSeg.mat'])
    load([TPN 'data\AllSeg.mat'])
    Mids=mean(AllSeg,3);
    Inside=((Mids(:,3)>Zstart(a)) & (Mids(:,3)<Zstop(a)));
    AllSeg=AllSeg(Inside,:,:);
    DrawSeg(:,1:2,:)=AllSeg(:,1:2,:)/yxum;
    DrawSeg(:,3,:)=AllSeg(:,3,:)/zum; 
    DrawSeg=round(DrawSeg); %round to index
    DrawSeg(DrawSeg==0)=1;  %get rid of non indexables
      
    clear AllSeg
    ys=round(max(DrawSeg(:,1,1))*1.1);
    xs=round(max(DrawSeg(:,2,1))*1.1);
    zs=round(max(DrawSeg(:,3,1))*1.1);
    
    'Drawing Segments'
    SegSum=zeros(ys,xs);
    SegOrtho1=zeros(zs,xs);
    SegOrtho2=zeros(xs,zs);
    SkelRes=.5;
    for i=1:size(DrawSeg,1)
        Dist=sqrt((DrawSeg(i,1,1)-DrawSeg(i,1,2))^2 + (DrawSeg(i,2,1)-DrawSeg(i,2,2))^2 + (DrawSeg(i,3,1)-DrawSeg(i,3,2))^2); %find distance
        Length(i)=Dist;
          devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
        for d=1:devs+1
            sy=DrawSeg(i,1,1)+((DrawSeg(i,1,2)-DrawSeg(i,1,1))/devs)*(d-1);
            sx=DrawSeg(i,2,1)+((DrawSeg(i,2,2)-DrawSeg(i,2,1))/devs)*(d-1);
            sz=DrawSeg(i,3,1)+((DrawSeg(i,3,2)-DrawSeg(i,3,1))/devs)*(d-1);
            %%Draw Max
            SegSum(min(round(sy),ys),min(round(sx),xs))=SegSum(min(round(sy),ys),min(round(sx),xs))+1;
            %%Draw Orthos
            %SegOrtho2(min(round(sy),ys),min(round(sz),zs))=SegOrtho2(min(round(sy),ys),min(round(sz),zs))+1;
            %SegOrtho1(min(round(sz),zs),min(round(sx),xs))=SegOrtho1(min(round(sz),zs),min(round(sx),xs))+1;
        end
    end
    clear Dist
    image(SegSum*100)
    %%%%%%%%%%%%SegOrhos too big?
end % if AllSeg exists



%% Draw Dots
if exist([TPN 'data\DotStats.mat'])
    load([TPN 'data\DotStats.mat'])
    Inside=((DotStats(:,3,3)>Zstart(a)) & ( DotStats(:,3,3)<Zstop(a)));
    DotStats=DotStats(Inside,:,:);
    DrawDots=round(DotStats(:,:,2));
    DrawDots(DrawDots<1)=1;
    clear DotStats    
    
    DotSum=zeros(ys, xs);
    DotOrtho1=zeros(zs,xs);
    DotOrtho2=zeros(xs,zs);
    for i=1:size(DrawDots,1)
        DotSum(min(DrawDots(i,1),ys),min(DrawDots(i,2),xs))=DotSum(min(DrawDots(i,1),ys),min(DrawDots(i,2),xs))+1;
      %  DotOrtho2(min(DrawDots(i,1),ys),min(DrawDots(i,3),zs))=DotOrtho2(min(DrawDots(i,1),ys),min(DrawDots(i,3),zs))+1;
       % DotOrtho1(min(DrawDots(i,3),zs),min(DrawDots(i,2),xs))=DotOrtho1(min(DrawDots(i,3),zs),min(DrawDots(i,2),xs))+1;
    end
    
    %%Dialate images
    Shape=ones(5);
    DotSum2=conv2(DotSum,Shape,'same');
    image((DotSum2)*100)
      
end


%% Combine Images and Data

DotDend(:,:,1)=SegSum;
DotDend(:,:,2)=DotSum2;
DotDend(:,:,3)=SegSum*0;
imwrite(DotDend,[TPN 'images\arbor\DotDend' num2str(a) '.tif'],'tif')
        













