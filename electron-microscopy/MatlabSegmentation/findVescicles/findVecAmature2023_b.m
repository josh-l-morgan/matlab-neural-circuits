colormap gray(256)

%[TFN TPN] =  GetMyFile('.tif')

TFN = 'Tile_r2-c1_w8_2_sec006_580_10820.tif';
TPN = 'Z:\Active\mahsa\containers\DeblurGAN\deblur-gan\train\B\'


Iraw = imread([TPN TFN]);
%I = Iraw(200:750,1:350);
I = double(Iraw);

colormap gray(256)
subplot(2,2,1)
image(255-I)


Im = medfilt2(double(I), [3 3]);

[ys, xs] = size(I);

%% Create hollow gaussian with negative surround

R = 50;

fsize = R * 2 + 1;
distMap = zeros(fsize,fsize);
for x = 1:fsize,
    for y = 1:fsize
        distMap(y,x) = sqrt((y-R-1)^2 + (x-R-1)^2);
    end
end
useK = find(distMap<=R);
blankK = find(distMap>R);
% 
% memA = 10
% memSD = 5
% memK = memA * exp(-.5 * ((distMap)/memSD)^2);
% memK(useK) = memK(useK)-mean(memK(useK));
% memK = memK /max(memK(:));
% memK(blankK) = 0;
% image(memK*100+100)


centA = 1
vRad = 3;
centSD = vRad

reMap = distMap;
reMap = reMap-vRad;
reMap(reMap<0) = 0;
center = centA * exp(-.5 * ((reMap)/centSD).^2);
if 0
    center = center /max(center(:))* 1.7;
    %center(center > mean([mean(center(:)) max(center(:))])) = mean(center(:));

    surA = 1
    surSD = 30;
    surround = surA * exp(-.5 * ((distMap)/surSD)^2);
    surround = surround/max(surround(:))*2;

    sur2A = 1
    sur2SD = 50;
    surround2 = sur2A * exp(-.5 * ((distMap)/sur2SD)^2);
    surround2 = surround2/max(surround2(:))*.4;

    kRat = mean(center(useK))/mean(surround(useK));
    kern =  surround - surround2 - center;
    mean(kern(:))
    kern(useK) = kern(useK)-mean(kern(useK));
    kern(blankK) = 0;

    %%Ring kern
    ring = distMap*0;
    r = [3 ,5,20];
    ring(distMap<r(3)) = -1;
    ring(distMap<r(2)) = 1;
    ring(distMap<r(1)) = -1;
    ring(ring>0)=ring(ring>0)/sum(ring(ring>0));
    ring(ring<0)=ring(ring<0)/abs(sum(ring(ring<0)));
    kern = ring;
else
    kern = center;
    kern = kern / sum(kern(:));
end



%image(kern * 155/max(kern(:)) + 100 )
subplot(2,2,2)
plot(kern(R+1,:))
%ylim([-1 2])


%% Convolve
subplot(2,2,3)
Imf = fastCon(Im,kern);

sdI = std(Imf(:)); 
Imf = Imf - mean(Imf(:));
Imf = Imf * 50/sdI;
Imf = Imf + 120;
image(256-Imf)

col(:,:,1) = Im;
col(:,:,3) = Imf;
subplot(2,2,4) 
image(uint8(256-col))
drawnow





%% find peaks

BW = imregionalmax(Imf);
[py px] = find(BW>0);

image(BW*10000)

Icol = cat(3,I,I,I);
Icol(:,:,1) = BW*1000;
subplot(2,2,4) 
image(Icol)


%% shape filters

%%Define shapes
lookDist = 20;
aBinWidth = .2;
aRange = 0:.1:.9;
toriR = [0 3; 4 6; 8 10];


[y x] = find(ones(lookDist*2));
y = y-lookDist;
x = x - lookDist;

a = atan2(y,x);
a = a - min(a);
a = a/max(a);
d = sqrt(y.^2 + x.^2);



%%Make pie slices
clear pie
for s = 1:length(aRange);

    isOr = 0;
    aStart = aRange(s)-aBinWidth;
    if aStart<0
        aStart = 1+aStart;
        isOr = 1;
    end
    aStop = aRange(s) + aBinWidth;
    if aStop>1
        aStop = aStop - 1;
        isOr = 1;
    end

    if isOr
            isIn = ((a>=aStart) | (a <= aStop)) & (d<=lookDist);

    else
    isIn = (a>=aStart) & (a <= aStop) & (d<=lookDist);
    end
     
    pie(s).isIn = isIn;
     pie(s).x = x(isIn);
     pie(s).y = y(isIn);
     pie(s).d = d(isIn);


end

%%Make tori
clear tori
for t = 1:size(toriR,1);

        isIn = (d>=toriR(t,1)) & (d<=toriR(t,2));
        tori(t).isIn = isIn;
        tori(t).x = x(isIn);
        tori(t).y = y(isIn);
        tori(t).d = d(isIn);

end

%%Make checkNear
isIn = (d<=5);
checkNear.isIn = isIn;
checkNear.x = x(isIn);
checkNear.y = y(isIn);
checkNear.d = d(isIn);


%%MakePT
for s = 1:length(pie);
    for t = 1:length(tori);

        pt(s,t).isIn = pie(s).isIn & tori(t).isIn;
        pt(s,t).x = x(pt(s,t).isIn);
        pt(s,t).y = y(pt(s,t).isIn);
        pt(s,t).d = d(pt(s,t).isIn);

    end
end

%% Search peaks


usePeaks = find((py>lookDist) & (py<(size(I,1)-lookDist)) & ...
    (px>lookDist) & (px<(size(I,2)-lookDist)));


ptMed = zeros(length(py),size(pt,1),size(pt,2));
showLook = I*0;
for p = 1:length(py)
    if sum(usePeaks==p)
    meds = zeros(size(pt,1),size(pt,2));
    showLook = zeros(size(I));

    sY = py(p) + checkNear(1).y;
    sX = px(p) + checkNear(1).x;
    sInd = sub2ind(size(I),sY,sX);
    nearVals = I(sInd);
    lowest = find(nearVals == min(nearVals),1);

    my = sY(lowest);
    mx = sX(lowest);


    clf
    %%Get middle
    sY = my + tori(1).y;
    sX = mx + tori(1).x;
    sInd = sub2ind(size(I),sY,sX);
    pv(p).centV = mean(Im(sInd));
    showLook(sInd) = 1000;

    sY = my + tori(3).y;
    sX = mx + tori(3).x;
    sInd = sub2ind(size(I),sY,sX);
    pv(p).outerV = median(Im(sInd));

    subplot(1,3,1)
    Icol(:,:,3) = showLook*1000;
    image(uint8(Icol))

    pieMeds = zeros(size(pt,1),2);
    for s = 1:size(pt,1) %search pie slices

        for t = 1:2
            sY = my + pt(s,t+1).y;
            sX = mx + pt(s,t+1).x;

            goodLook = 1;
            sInd = sub2ind(size(I),sY,sX);

            pieVals = Im(sInd);
            pieMeds(s,t) = median(pieVals);
        end

    end
    pv(p).pieMeds = pieMeds;

    pieDif = pieMeds(:,1)-pieMeds(:,2);
    centDif = pieMeds(:,1)-pv(p).centV ;
    minDif = min(pieDif,centDif);

    subplot(1,3,2)
    cla
    hold on
    plot([1:10]*0+ pv(p).centV  ,'r')
    plot(pv(p).pieMeds(:,2) ,'b')
    plot(pv(p).pieMeds(:,1),'g')
    ylim([0 256])

    subplot(1,3,3)
    plot(minDif,'k')
    pv(p).minDif = minDif;
    pv(p).aveBin = mean(minDif>0);
    pv(p).aveDif = mean(minDif);
    else
        pv(p).aveBin = 0;
    end

end



aveBin = [pv.aveBin];

Icol = cat(3,I,I,I);
for i = 0:.1:1
    i
    clf
    isV = find(aveBin>=i);
    showLook = zeros(size(I));
    for p = 1:length(isV)
        sY = py(isV(p)) + checkNear.y;
        sX = px(isV(p)) + checkNear.x;
        sInd = sub2ind(size(showLook),sY,sX);
        showLook(sInd) = 1;
    end
    Icol(:,:,3) = showLook*1000;
    image(uint8(Icol))
    pause(1)
end


return
%% Threshold

It = Imf < (median(Imf(:)) * .5);
image(It*1000);

Il = bwlabel(It,4);

Ivec = It * 0; 
Ilarge = Ivec;
vPos = [];
Icent = It * 0;
Ilong = It * 0;
Ismall = It * 0;
for i = 1: max(Il(:))
    lPos = find(Il == i);
    if length(lPos)>200
        Ilarge(lPos) = 1;
    elseif length(lPos)<20
        Ismall(lPos) = 1;
    else
        [y x ] = ind2sub([ys xs],lPos);
        [a pc c d] = pca([y x]);
        lengthWidth = max(pc,[],1) - min(pc,[],1);
        if lengthWidth(1)/lengthWidth(2) > 2;
            Ilong(lPos) = 1;
        else
            
            Ivec(lPos) = 1;
            [y x] = ind2sub([ys xs],lPos);
            vPos(size(vPos,1)+1,:) = [mean(y) mean(x)];
            Icent(round(mean(y)), round(mean(x))) = 1;
        end    
    end
end

image(Icent * 10000)
Imem = Ilarge | Ilong;

%% group

R = 20;

fsize = R * 2 + 1;
distMap = zeros(fsize,fsize);
for x = 1:fsize,
    for y = 1:fsize
        distMap(y,x) = sqrt((y-R-1)^2 + (x-R-1)^2);
    end
end
useK = find(distMap<=R);
blankK = find(distMap>R);

% groupA = R^2
% groupSD = R
% groupK = groupA * exp(-.5 * ((distMap)/groupSD)^2);
% groupK = groupK - groupK(1,R+1);

groupK = R - distMap;
groupK(groupK<0) = 0;
groupK = groupK/max(groupK(:));

Igroup = conv2(Icent, groupK,'same');
Igroup(Igroup < .6) = 0;
Igroup(Ilarge>0) = 0;

SE = strel('disk', 5);
Iopen = imopen(Igroup,SE);
Iopen(Ilarge>0) = 0;

[IgLab IglN] = bwlabel(Iopen,4);

Ipass = Igroup * 0;
for i = 1: IglN
   groupPos = find(IgLab == i);
%    gVal(i) = sum(Igroup(groupPos));
   if length(groupPos) > 1000
       numVec(i) = sum(Icent(groupPos));
       Ipass(groupPos) = 50 + numVec(i)* 5;
   end
   
end

image(Ipass * 1000)


image(Igroup* 100/mean(Igroup(:)))

% 
% SE = strel('disk',15);
% Igroup = imfilter(Ivec,fspecial('disk',20));
% 
% Iclose = imclose(Ivec,SE);
% Iclose(Ilarge>0) = 0;
% Iclust = Iclose * 0;
% 
% 
% 
% [Icl gN] = bwlabel(Iclose,4);
% for i = 1: gN
%    groupPos = find(Icl == i);
%    if length(groupPos) > 1000
%        vecAr = sum(Ivec(groupPos));
%        if vecAr/length(groupPos) > .3
%             Iclust(groupPos) = 1;
%        end
%    end
%     
% end



%% Show
clf
Ic = Ipass;
Ic(:,:,2) = I;
Ic(:,:,3) = (Ivec * 1000);

image(uint8(Ic))
imwrite(I,[TPN 'sample.tif'],'compression','none')
imwrite(uint8(Ic),[TPN 'found.tif'],'compression','none')
