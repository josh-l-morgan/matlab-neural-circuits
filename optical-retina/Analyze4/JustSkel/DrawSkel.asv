function[]=DrawSkel(TPN)

colormap gray(255)

load([TPN 'D.mat'])
load([TPN 'data\AllSeg.mat'])
[ys xs zs] = size(D);
load([TPN 'ImageInfo.mat'])
xyum=ImageInfo{1};
zum=ImageInfo{2};

AllSegV=AllSeg/xyum;
AllSegV(:,3,:)=AllSeg(:,3,:)/zum;


Skel= zeros([ys xs zs],'uint8');
Skel=uint8(D);
Slope = AllSegV(:,:,2)-AllSegV(:,:,1);
Length = sqrt(Slope(:,1).^2 + Slope(:,2).^2 + Slope(:,3).^2);

Divs=2;
Sc=0:Divs;

for l = 1: fix(max(Length))+1
    Get=Length<l;
    Divs=l*10;
    for i = 0: Divs
        Marks=round(AllSegV(Get,:,1)+(Slope(Get,:)/Divs)*i);
        Skel(sub2ind([ys xs zs],Marks(:,1),Marks(:,2),Marks(:,3)))=1000;
    end
end

image(max(Skel,[],3)*100)
imwriteNp(TPN,Skel,'Skel')
