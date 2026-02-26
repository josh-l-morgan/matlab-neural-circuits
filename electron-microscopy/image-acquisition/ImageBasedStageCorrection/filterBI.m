% 
[BFN BPN] = GetMyFile;
[IFN IPN] = GetMyFile;

%%
colormap gray(256) 
B = double(imread([BPN BFN]));
I = double(imread([IPN IFN]));

%% process I
kSize = [100 100];
s1 =3;
s2 = 5;

kRad = (kSize + 1)/2;
kern = zeros(kSize);

[y x z] = ind2sub(kSize,find(kern==0));
dists = sqrt(((y-kRad(1))).^2 + ((x - kRad(2))).^2);

cKern = 1 * exp(-.5 * (dists/s1).^2);
cKern = cKern/sum(cKern(:));
sKern = 1 * exp(-.5 * (dists/s2).^2);
sKern = sKern/sum(sKern(:));
kern(:) = cKern - sKern;

subplot(2,1,1)
plot(kern(round(kRad(1)),:))

%%Convolve

Imf = fastCon(I,kern);
Bmf = fastCon(B,kern);
Imf = Imf(pixClip +1:end-pixClip,pixClip +1:end-pixClip);
Bmf = Bmf(pixClip + 1:end - pixClip,pixClip +1:end-pixClip);
Imf(Imf<0)=0;
Bmf(Bmf<0)=0;

%%image
subplot(2,2,1)
image(B)
subplot(2,2,3)
image(Bmf * 10)

subplot(2,2,2)
image(I)
subplot(2,2,4)
image(Imf * 10)

%% find objects
%%{
subplot(1,2,1)
image(I)
subplot(1,2,2) 
sI = edge(Imf,'Canny',[.01 .3],s1/2);
se = strel('disk',1);
dI = imdilate(sI,se);
lI = bwlabel(~dI,8);
% gKern = gaus3d([6 6 1],2);
% sI = fastCon(sI,gKern);
props = regionprops(lI,'Area','PixelIdxList');
areas = [props.Area];
pass = find((areas<200) & (areas>2));
middles = cat(1,props(pass).PixelIdxList);
aI = lI * 0;
aI(middles) = 255;
image(aI)
%}













