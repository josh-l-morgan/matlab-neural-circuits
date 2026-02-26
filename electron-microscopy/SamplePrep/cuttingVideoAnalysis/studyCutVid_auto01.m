
colormap gray(256)
subplot(1,1,1)

[TFN TPN] = GetMyFile;  %Get file location
vidObj = mmreader([TPN TFN]);  %Get movie information as object
%vidVals = get(vidObj); %to look at available movie information
nFrames = vidObj.NumberOfFrames 
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

manual = 0;

if manual

%% Mark Section
firstI = read(vidObj, 1);  %read first few images
image(firstI),pause(.01)
cF = gcf;
sprintf('click on edges of section and press return')
[gx gy] = ginput;

else 
    
%% Find Section
rawI = double(read(vidObj, [1 100]));  %read first few images
I = squeeze(mean(rawI,3));  % convert to black/white
difI = I(:,:,2:end)-I(:,:,1:end-1);
image(sum(difI,3))

for i = 1:size(difI,3)
   image(difI(:,:,i)*3),pause
end

meanI = mean(I,3); %take average of stack
subI = I - repmat(meanI,[1 1 size(I,3)]);
% for i = 1:size(subI,3)
%     sI = subI(:,:,i);
%     sI = sI/max(sI(:));
%     subI(:,:,i) = sI;
%     image(subI(:,:,i)*256),pause
% end
difI = mean(abs(subI),3);
image(difI)

yMean = mean(difI,2);
yBar = find(yMean == max(yMean),1);
xMean = mean(difI(yBar-5:yBar+5,:),1);
plot(xMean)


%%
cKern = gaussmf([-100 : 100],[5 0])
cKern = cKern/sum(cKern);
filtX = filtfilt(cKern, 1, xMean);
plot(filtX)

thresh = max(filtX) - mean(filtX);
passes = find(filtX> (mean(filtX) + thresh/2));

gx = passes;
gy = yBar;

end



%%  Create image from line scans
ysamp = 1;
yave = 1;

lx = round(min(gx));
rx = round(max(gx));
by = round(max(gy));
ty = by - yave+1;

strip = zeros(nFrames * ysamp,rx-lx+1,3);
for i = 1: nFrames
    frame = read(vidObj, i);  %read first few images
    sy = (i-1) * ysamp + 1;
    strip(sy :(sy + ysamp - 1),:,:) = mean(frame(ty:by,lx:rx,:),1);
end

image(uint8(strip))
%%
mStrip = squeeze(mean(strip,2));
plot(mStrip(:,1),'r')
hold on
plot(mStrip(:,2),'g')
plot(mStrip(:,3),'b')

ratStrip = mean(mStrip(:,1:2),2)./mStrip(:,3);
plot(ratStrip*100,'k')
hold off


%% Filter strip
x = -10 : 10;
sig = 3;
cKern = gaussmf([-100 : 100],[5 0])
cKern = cKern/sum(cKern);

sKern = gaussmf([-100 : 100],[30 0])
sKern = sKern/sum(sKern);
mKern = cKern - sKern;

plot(mKern)

fStrip = filtfilt(mKern,1, mean(mStrip(:,1:2),2));
plot((fStrip))
hold on
mmStrip = mean(mStrip(:,1:2),2);
%plot(mmStrip - mean(mmStrip),'k')
hold off

%% Test frame size
testf = [1:300];
L = length(fStrip);
n = fix(L/max(testf));

for f = testf
   sfram = reshape(fStrip(1:f*n),f,n);
   prodf(f) = mean(prod(sfram,2));
   image(sfram * 5),pause
end
subplot(2,1,1)
plot(prodf)

[pos vals] = getPeaks(prodf); %get peaks
fsize = pos(find(abs(vals) > std(abs(vals))*2,1)) %get first big peak

fnum = fix(L/max(fsize));
subplot(2,1,2)
mframs = reshape(fStrip(1:fsize*fnum),fsize,fnum);
image(mframs* 5)

template = mean(mframs(:,2:11),2);
plot(template)

%% find frames
findF = filter(template,1,fStrip);
plot(findF)
[pos vals] = getPeaks(findF);

[svals idx] = sort(vals,'descend');
nufVals = svals(fnum);
spos = pos(idx);

testVals = svals(1:fnum);
std(testVals)

sweepPos = sort(spos(1:fnum),'ascend');

scatter(spos(1:fnum),svals(1:fnum),'r');
hold on
scatter(spos(fnum+1:end),svals(fnum+1:end),'b')
hold off

%% Collect sweeps
clear sweeps

for i = 1:length(sweepPos)-1
    sweeps(:,:,:,i) = strip(sweepPos(i):sweepPos(i)+fsize,:,:);
end
meanSweep = mean(mean(mean(sweeps,3),4),2);
plot(meanSweep)
offset = find(meanSweep == min(meanSweep));

sweepPos2 = sweepPos + offset;
for i = 1:length(sweepPos2)-2
    sweeps(:,:,:,i) = strip(sweepPos2(i):sweepPos2(i)+fsize,:,:);
end
meanSweep = mean(mean(mean(sweeps,3),4),2);


%% find edge
clear eStack
for i = 1: size(sweeps,4)
    F = sweeps(:,:,:,i);
    F = mean(F,3);
    eF = edge(F,'canny',.1,.7);
    subplot(2,1,1)
    image(F)
    subplot(2,1,2)
    image(eF*1000)
    pause(.1)
    eStack(:,:,i) = eF;
end

sumE = sum(eStack,3);
image(sum(eStack,3)*10)
plot(sum(sumE,2))

dE = imfilter(sumE,fspecial('disk',4));
%dE = dE * 10/mode(dE(:));
watMin = std(dE(:));
minC2 = imhmin(dE,watMin);
wI2 = watershed(minC2,8);
image(wI2*10)
pause(.01)

[ys xs] = size(wI2);
snag = wI2(1:round(ys/2), round(xs/3):end-round(xs/3));
sec = mode(double(snag(:)));

[y x] = find(wI2 == sec);
tip = max(y);


sweepPos2 = sweepPos + offset + tip;
clear secs
for i = 1:length(sweepPos2)-2
    secs(:,:,:,i) = strip(sweepPos2(i):sweepPos2(i)+fsize,:,:);
end
meanSecs = mean(mean(mean(secs,3),4),2);


for i = 1:size(secs,4)
    image(uint8(secs(:,:,:,i))),pause(.01)
end

image(uint8(mean(secs,4)))

%% Analyze sections

mSecs = squeeze(mean(mean(secs,1),2))';
for c = 1:3
    mSecs(:,c) = mSecs(:,c)/mean(mSecs(:,c));
end
plot(mSecs)

ratSec = mean(mSecs(:,1:2),2) ./ mSecs(:,3);
plot(ratSec)
meanSec = mean(secs,4);
%image(uint8(meanSec))

difSec = secs - repmat(meanSec,[1 1 1 size(secs,4)]);




%%




%% read movie


imwrite(uint8(strip),[TPN 'strip.tif'],'Compression','none')
