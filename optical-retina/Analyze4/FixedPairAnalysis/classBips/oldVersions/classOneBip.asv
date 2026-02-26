%% Read data
clear all
colormap gray(256)
TPN = GetMyDir
'reading data'
TPNd = dir(TPN); TPNd = TPNd(3:length(TPNd));
TiffNames={};
for i = 1: size(TPNd,1)
    siz=length(TPNd(i).name);
    if TPNd(i).name(siz-2:siz)== 'tif'
        if TPNd(i).name(siz-8:siz-7)~='-R'
            TiffNames(length(TiffNames)+1,1)={TPNd(i).name};
        end
    end
end

for i = length(TiffNames):-1:1
    I(:,:,i) = imread([TPN TiffNames{i}]);
end

[yi xi zi] = ind2sub(size(I),find(I>0));
subI = [min(yi) max(yi);min(xi) max(xi);min(zi) max(zi)];
I = double(I(subI(1,1):subI(1,2),subI(2,1):subI(2,2),...
    subI(3,1):subI(3,2)));
image(sum(I,3))

D = I>0;
[ys xs zs] = size(D);
[Dys Dxs Dzs] = ind2sub(size(D),find(D>0));
D= D(min(Dys):max(Dys),min(Dxs):max(Dxs),min(Dzs):max(Dzs));
[ys xs zs] = size(D);


%% Find shaft (currently not removing shaft)
'find shaft'
sumD = sum(D,3);
subplot(2,2,1)
image(sumD * 255/max(sumD(sumD>0)))
% [maxYs maxXs] = find(sumD == max(sumD(:)));
% 
% cdRad = 30;
% shaft = [mean(maxYs) mean(maxXs)];
% cutDisk = fspecial('disk',cdRad);
% [nSy nSx] = ind2sub(size(cutDisk),find(cutDisk>0));
% nearShaft = sub2ind([ys xs],round(nSy + shaft(1)-cdRad), round(nSx + shaft(2)-cdRad));
% for z = 1:size(D,3)
%    Dplane = D(:,:,z);
%    Dplane(nearShaft) = 0;
%    Dc(:,:,z) = Dplane;
% end
% sumDc = sum(Dc,3)
% subplot(2,2,2)
% image(sumDc * 255/max(sumDc(sumDc>0)))
% % shift orthogonally
% 'shift to princomp'
% 
% [Yp Xp Zp] = ind2sub([ys xs zs], find(D>0));
% [coef, score, latent] = princomp([Yp Xp Zp]);
% 
% score(:,1) = score(:,1) - min(score(:,1)) + 1;
% score(:,2) = score(:,2) - min(score(:,2)) + 1;
% score(:,3) = score(:,3) - min(score(:,3)) + 1;
% 
% nY = fix(max(score(:,1)))+1;
% nX = fix(max(score(:,2)))+1;
% nZ = fix(max(score(:,3)))+1;
% 
% P = zeros(nY, nX, nZ);
% score = round(score);
% newInd = sub2ind([nY nX nZ],score(:,1),score(:,2),score(:,3));
% P(newInd) = 1;
% 
% for i = 1: size(score,1)
%    P(score(i,1),score(i,2),score(i,3)) = 1; 
% end
% xD = squeeze(sum(D,1));
% yD = squeeze(sum(D,2));
% xP = squeeze(sum(P,1));
% yP = squeeze(sum(P,2));
% subplot(2,2,1),image(xD * 100/median(xD(xD>0)));
% subplot(2,2,2),image(yD * 100/median(yD(yD>0)));
% subplot(2,2,3),image(xP * 100/median(xP(xP>0)));
% subplot(2,2,4),image(yP * 100/median(yP(yP>0)));
% pause(.01)
%% Render mask
P = D;
'render cell'
subplot(1,1,1)
p = patch(isosurface(P,.1));
%isonormals(xs,ys,zs,D,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1/.06 1/.06 1/.2])
view(3); axis tight
camlight 
lighting gouraud
pause
%% Characterize morphology
%% Depth profile
dP = squeeze(sum(sum(P,1),2));
bar(dP)
whm = find(dP>=(max(dP)/2));
tempDP = dP;
tempDP(whm) = 0;
hold on
bar(tempDP,'b')
hold off

topA = min(whm):size(P,3);

%% Polygone
M = P;
sumTop = sum(M(:,:,topA),3);
image(sumTop * 100/median(sumTop(sumTop>0)));


[y x] = find(sumTop>0);
numArea = length(x);
k = convhull(y,x);
scatter(y,x,'.')
hold on
plot(y(k),x(k),'r-')
hold off
polyArea = polyarea(y(k),x(k)) ;

%% Load Skeleton

load([TPN 'data/AllSeg.mat'])
scaleSeg = AllSeg;
scaleSeg(:,1,:) = scaleSeg(:,1,:)/0.06 + min(y);
scaleSeg(:,2,:) = scaleSeg(:,2,:)/0.06 + min(x);
scaleSeg(:,3,:) = scaleSeg(:,3,:)/0.2;
hold on
for i = 1:size(scaleSeg,1)
    plot([scaleSeg(i,1,1) scaleSeg(i,1,2)],[scaleSeg(i,2,1) scaleSeg(i,2,2)],'g')
end
hold off

xyLengths = sqrt((scaleSeg(:,1,1)-scaleSeg(:,1,2)).^2+(scaleSeg(:,2,1)-scaleSeg(:,2,2)).^2);
totLength = sum(xyLengths);

pause 

%% Find Butons
subplot(2,1,1)
image(sumTop * 100/median(sumTop(sumTop>0)));
fftI = abs(fft2(sumTop));
subplot(2,1,2)
image(fftI)

%% Band Pass filter
freak = [0 .1 .2 .3 .4 .5 .6 .7 .8 .8 1];
amp = [1 0 0 0 0 0 0 0 0 0 0];
myBP = fir2(10, freak, amp);
myBPr = ftrans2(myBP);
image(myBPr * 100/max(myBPr(:))+ 100)
%freqz2(myBPr)
filtTop = conv2(sumTop,myBPr,'same');
image(filtTop * 255/max(filtTop(:)))


%%Ratio different filts

%% 3D erosion
kern = ones(3,3,3)/27;
F  = zeros(10,10,10);
F(1:3,1:3,1:3) = 1;
F = P;
colF = uint8(sumTop * 20);
colF(:,:,2) = sumTop * 0;
colF(:,:,3) = sumTop * 20;
% for r = 1:30
%    
%     F = convn(F,kern,'same');
%     F = F>=.9;
%      F = convn(F,kern,'same');
%     F = F>=.1;
%     colF(:,:,1) = sum(F,3) * 20;
%     image(colF),pause(.10)
%     
% end

F = P;
F = convn(F,kern,'same');
F = F>=.5;
for r = 1:10
    r
    perF = bwperim(F);
    F(perF>0)=0;
    %image(sum(F,3)*100);pause
    F = convn(F,kern,'same');
    F = F>=.15;
    
    colF(:,:,1) = sum(F,3) * 20;
    image(colF*30),pause(.1)
end

roundish = sum(F(:));  %% record roundish volume. 

%% find buotns by watershed
tic
dI =bwdist(I==0,'euclidean');
toc

for i = 1:size(dI,3)
    image((dI(:,:,i)*20)),pause(.1)
end
%% fft convolve
kSize = 5;
fKern = zeros(kSize,kSize,kSize);
fKern(round(kSize/2),round(kSize/2),round(kSize/2)) = 1;
dfKern = bwdist(fKern);
dfKern = kSize/2-dfKern;
dfKern = dfKern/max(dfKern(:));

gKern = 1 * exp(-.5 * ((dfKern-kSize+1)/(kSize/3)).^2);
for i = 1:size(gKern,3)
   image(fitH(gKern,i)),pause(.1); 
end

cI = fastCon(dI,gKern);

for i = 1:size(cI,3)
    i
    subplot(2,1,1)
    image(dI(:,:,i)*20),pause(.1)
    subplot(2,1,2)
   image(fitH(cI,i)),pause(.1)
end

%%
dI2 = max(cI(:))-cI;
dI2 = dI2 * 100/max(dI2(:));
dI2(~I) = -Inf;
dI2 = imhmin(dI2,5); %
%dI2(dI2<-10) = -10;
wI = watershed(dI2);
for i = 1:size(wI,3)
    imshow(label2rgb(wI(:,:,i),'jet','w')),pause
end


%%  Collect data
%{
maskArea
polyArea
maskVol
%}

