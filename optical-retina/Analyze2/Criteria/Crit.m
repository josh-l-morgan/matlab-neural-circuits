%function[] = anaCrit(TPN, DPN)

%Get File
TPN=GetMyDir
  
%Get Image dimensions
DPN = [TPN 'I\']
Idir=dir(DPN);
Idir=Idir(3:size(Idir,1));
Im=imread([DPN Idir(1).name]);
[ys xs c] = size(Im); zs=size(Idir,1);
zum=.3; yxum=.103;



clear DrawSeg DrawDots clear Dots    
%% Draw Segments
    load([TPN 'data\AllSegCut.mat'])
     DrawSeg(:,1:2,:)=AllSegCut(:,1:2,:)/yxum;
     DrawSeg(:,3,:)=AllSegCut(:,3,:)/zum; 
     DrawSeg=round(DrawSeg); %round to index
     DrawSeg(DrawSeg==0)=1;  %get rid of non indexables

     clear AllSeg
     SegSum=zeros(ys,xs);
     SegOrtho1=zeros(zs,xs);
     SegOrtho2=zeros(ys,zs);
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
                        SegOrtho2(min(round(sy),ys),min(round(sz),zs))=SegOrtho2(min(round(sy),ys),min(round(sz),zs))+1;
                        SegOrtho1(min(round(sz),zs),min(round(sx),xs))=SegOrtho1(min(round(sz),zs),min(round(sx),xs))+1;
                    end
                end
                clear Dist AllSEgCut
                %image(SegSum*100),pause(.01)
                %%%%%%%%%%%%SegOrhos too big?


%% Draw Dots
load([TPN 'Dots.mat'])
DrawDots=round(Dots.Pos);
           
OKs=find(Dots.OK+1);
DrawDots(DrawDots<1)=1;

DotSum=zeros(ys, xs);
DotVol=DotSum;DotItSum=DotSum;DotDFOf=DotSum;DotDF=DotSum;DotD2M=DotSum;DotCut=DotSum;
DotOrtho1=zeros(zs,xs);
DotOrtho2=zeros(ys,zs);
            
for i=1:size(OKs,1)
                DotSum(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=DotSum(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))+1;
                DotOrtho2(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),3),zs))=DotOrtho2(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),3),zs))+1;
                DotOrtho1(min(DrawDots(OKs(i),3),zs),min(DrawDots(OKs(i),2),xs))=DotOrtho1(min(DrawDots(OKs(i),3),zs),min(DrawDots(OKs(i),2),xs))+1;
                DotVol(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.Vol(OKs(i));
                DotItSum(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.ItSum(OKs(i));
                DotDFOf(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.DFOf(OKs(i));
                DotDF(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.DF(OKs(i));
                DotD2M(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.DistToMask(OKs(i));
                DotCut(min(DrawDots(OKs(i),1),ys),min(DrawDots(OKs(i),2),xs))=Dots.Cut(OKs(i));
end

            
%% Define Criteria
%%Make filter
fSize=5;
Gh=fspecial('gaussian',fSize*2+1,fSize);
Gh(Gh<Gh(1,fSize))=0;
Gh=Gh*1/max(Gh(:));

DotSumf=imfilter(DotSum*255,Gh,'same');
SegSumf=imfilter(SegSum*255,Gh,'same');

clear IMf
b=DotSumf+SegSumf;
IMf(:,:,1)=(SegSumf>0)*100;
IMf(:,:,2)=(DotSumf*0);
IMf(:,:,3)=(DotSumf>0)*100;
IMf=uint8(IMf);
image(IMf)
pause(.01)

%% Run New Thresholds
Vol=1
%Cut=2
ItSum=2
DFOf=4
DF=1
DistToMask=1       

%%Draw Pass
Ipass=(DotVol>=Vol) & (DotItSum>=ItSum) & (DotDFOf>=DFOf) & (DotDF>=3) & (DotCut<=1) & (DotD2M<=DistToMask);
Ipassf=imfilter(Ipass*255,Gh,'same');
IMf(:,:,2)=Ipassf*255;
image(IMf);pause(.01)




%{
Name=[KPN 'Combos\' Kdir(k).name '.tif'];
imwrite(IMf,Name,'Compression','none')
if ~exist([TPN 'images']), mkdir([TPN 'images']);end
imwrite(IMf,[TPN 'images\IMf.tif'],'Compression','none')
%}
    



%% Select Puncta
clear pass
for i = 1: Dots.Num
   v=Dots.Vol(i)>=Vol;
   c=Dots.Cut(i)<2;
   itsum=Dots.ItSum(i)>=ItSum;
   dfof=Dots.DFOf(i)>=DFOf;
   df=Dots.DF(i)>=DF;
   dm=Dots.DistToMask(i)<=DistToMask;
   pass(i)=v  & c & itsum & dm & df & dfof ;
end

P=find(pass')'; %% list of passing puncta
Dots.OK=pass';
%save([TPN 'Dots.mat'],'Dots')

%% Draw puncta matrix

%{

%%ID puncta

IDs=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
for i = 1 : size(P,1)
   for v=1:size(Dots.Vox(P(i)).Pos,1)
        IDs(Dots.Vox(P(i)).Pos(v,1),Dots.Vox(P(i)).Pos(v,2))=P(i); 
   end
end
imwrite(IDs,[TPN 'images\IDs.tif'],'Compression','none')



Passed = zeros(Dots.ImSize,'uint8');
for i = 1: size(Dots.Vox,2)
    Passed(Dots.Vox(i).Ind)=50;    
end

for i = 1: size(P,1)
   Passed(Dots.Vox(P(i)).Ind)=200;     
end

colormap gray(255)
maxPassed=max(Passed,[],3);
image(maxPassed)
%{
imwriteNp(TPN,Passed,'Passed')
imwrite(maxPassed,[TPN 'images\maxPassed.tif'],'Compression','none')
%}

%}








