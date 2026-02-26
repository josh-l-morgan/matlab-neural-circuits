function[aI] = edgeVessel(I,s1,thresh)

kSize = [100 100];

if ~exist('s1','var')
    s1 =3;
end

if ~exist('thresh','var')
    thresh = .2;
end
lowthresh = thresh/4;

sI = edge(I,'Canny',[lowthresh thresh],s1);
se = strel('disk',1);
dI = imdilate(sI,se);
lI = bwlabel(~dI,8);

props = regionprops(lI,'Area','PixelIdxList');
areas = [props.Area];
pass = find((areas<200) & (areas>2));
middles = cat(1,props(pass).PixelIdxList);

aI = lI * 0;
aI(middles) = 255;







