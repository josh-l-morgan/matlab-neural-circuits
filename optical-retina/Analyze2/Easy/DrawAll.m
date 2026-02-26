KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 65:65%size(Kdir,1)
    TPN = [KPN '\' Kdir(k).name '\'];
    DPN = [TPN 'I\']   
    
    clear DrawSeg DrawDots clear Dots
    if exist([TPN 'data\AllSegCut.mat']) & exist([TPN 'Dots.mat'])%& ~exist([TPN 'images\IMf.tif'])
        
            Idir=dir(DPN);
    Idir=Idir(3:size(Idir,1));
    Im=imread([DPN Idir(1).name]);
    [ys xs c] = size(Im); zs=size(Idir,1);
        
        
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
                clear Dist
                %image(SegSum*100),pause(.01)
                %%%%%%%%%%%%SegOrhos too big?




        %% Draw Dots
            load([TPN 'Dots.mat'])
            DrawDots=round(Dots.Pos);
            DrawDots=DrawDots(Dots.OK,:);
            DrawDots(DrawDots<1)=1;

            DotSum=zeros(ys, xs);
            DotOrtho1=zeros(zs,xs);
            DotOrtho2=zeros(ys,zs);
            for i=1:size(DrawDots,1)
                DotSum(min(DrawDots(i,1),ys),min(DrawDots(i,2),xs))=DotSum(min(DrawDots(i,1),ys),min(DrawDots(i,2),xs))+1;
                DotOrtho2(min(DrawDots(i,1),ys),min(DrawDots(i,3),zs))=DotOrtho2(min(DrawDots(i,1),ys),min(DrawDots(i,3),zs))+1;
                DotOrtho1(min(DrawDots(i,3),zs),min(DrawDots(i,2),xs))=DotOrtho1(min(DrawDots(i,3),zs),min(DrawDots(i,2),xs))+1;
            end


        %% dilate images
        fSize=5;
        Gh=fspecial('gaussian',fSize*2+1,fSize);
        Gh(Gh<Gh(1,fSize))=0;
        Gh=Gh*1/max(Gh(:));

        DotSumf=imfilter(DotSum*255,Gh,'same');
        SegSumf=imfilter(SegSum*5,Gh,'same');
        
        clear IMf
        b=DotSumf+SegSumf;
        IMf(:,:,1)=SegSumf;
        IMf(:,:,2)=DotSumf;
        IMf(:,:,3)=b*0;
        IMf=uint8(IMf);
        image(IMf)
        pause(.01)

    else
        IMf = zeros(100);
    end

if ~exist([TPN 'images']), mkdir([TPN 'images']);end
imwrite(IMf,[TPN 'images\IMf.tif'],'Compression','none')    


%Create image index
IMf2=imresize(IMf,yxum*2)*2;
clear Cell
if exist([TPN 'Cell.mat'])
    load([TPN 'Cell.mat'])
    Name=[KPN 'Combos\' 'P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif'];
else
    Name=[KPN 'Combos\' Kdir(k).name '.tif'];
end
if isempty(find(Name=='?'));
    imwrite(IMf2,Name,'Compression','none')
end




end %run all cells