


%% Get image Names
SPN = '\\storage1.ris.wustl.edu\jlmorgan\Active\morganLab\DATA\ImageProcessing\MitochondriaSegmentation\model-200\'
dSPN = dir([SPN '*.png']);
iNams = {dSPN.name}

%% Parse images
oCount = 0;
rCount = 0;
clear outName outSec rawName rawSec
for i = 1 :length(iNams);

    nam = iNams{i};

    if strcmp(nam(1:6),'output')
        oCount = oCount+1;
        sec = str2num(nam(8:end-4));
        outName{oCount} = nam;
        outSec(oCount) = sec;


    elseif strcmp(nam(1:3),'raw')
        rCount = rCount+1;
        sec = str2num(nam(5:end-4));
        rawName{rCount} = nam;
        rawSec(rCount) = sec;

    end


end

secs = unique(intersect(outSec,rawSec));
numSecs = length(secs);


%% Show start
clear stackR stackO
for i = 1:length(secs)

    subplot(1,2,1)
    targR = find(rawSec==secs(i));
    R = imread([SPN rawName{targR}]);
    image(R)
    subplot(1,2,2)
    targO = find(outSec==secs(i));
    O = imread([SPN outName{targO}]);
    image(O)
    drawnow
    stackR(:,:,i) = mean(R,3);
    stackO(:,:,i) = mean(O,3);

end

%% Merge stacks

[rsY rsX rsZ] = size(stackR);
[osY osX osZ] = size(stackO);
bufY = (rsY-osY)/2
bufX = (rsX-osX)/2;

clear stackM
stackM(:,:,:,2) = stackR;
stackM(bufY+1:bufY+osY,bufX+1:bufX+osX,:,3) = stackO;
stackM(:,:,:,1) = stackR * .3;
stackM(bufY+1:bufY+osY,bufX+1:bufX+osX,:,1) = stackR(bufY+1:bufY+osY,bufX+1:bufX+osX,:);
stackRC = stackR(bufY+1:bufY+osY,bufX+1:bufX+osX,:);

colormap gray(255)
for i = 1:size(stackM,3)
    subplot(1,2,1)
    % image(uint8(squeeze(stackM(:,:,i,:))))
    image(stackM(:,:,i,1))
    subplot(1,2,2)
    image(stackM(:,:,i,3))
    pause
end


%% Watershed

minWatArea = 500 * 3;
stackW = stackO * 0;



sp1 = subplot(1,3,1);
Ir = squeeze(mean(stackRC(:,:,:,1),3));
colormap(sp1,gray(256))
image(Ir)

sp2 = subplot(1,3,2)

I = stackO(:,:,:);

Iw = 256 - I;
Iw = imhmin(Iw,1,26);
[W ] =  watershed(Iw,26);
lNum = max(W(:));
image(mean(I,3))
colormap(sp2,gray(256))

sp3 = subplot(1,3,3)
image(max(W,[],3))
colormap(sp3,hsv(double(lNum)))
pause(.01)

props = regionprops(W,I,'maxintensity','area','BoundingBox');
A = cat(1,props.Area);
MaxI = cat(1,props.MaxIntensity);
medI = median(I(:));
propV = MaxI-medI;
bbs = cat(1,props.BoundingBox);
bbs(:,1:3) = ceil(bbs(:,1:3));
bbs(:,4:6) = bbs(:,4:6)-1;

bI = I *0;
bigWat = find(A>=minWatArea);
for w = 1:length(bigWat)
    bb = bbs(bigWat(w),:);
    Ibb = I(bb(2):bb(2)+bb(5),bb(1):bb(1)+bb(4), bb(3):bb(3)+bb(6));
    Wbb = W(bb(2):bb(2)+bb(5),bb(1):bb(1)+bb(4), bb(3):bb(3)+bb(6));
    Ibb(Wbb ~=bigWat(w)) = 0;
    Iw = Ibb >= (propV(bigWat(w))/2 + medI);
    Ibbt = Ibb.*Iw;
    Ibbt(Ibbt>0) = sum(Ibbt(:));
    bI(bb(2):bb(2)+bb(5),bb(1):bb(1)+bb(4), bb(3):bb(3)+bb(6)) = ...
        bI(bb(2):bb(2)+bb(5),bb(1):bb(1)+bb(4), bb(3):bb(3)+bb(6)) +Ibbt;

end
image(sum(bI,3)*.001)
pause(.01)

stackW = bI;
sumW = sum(stackW,3);
image(sumW * 256/max(sumW(:)))




