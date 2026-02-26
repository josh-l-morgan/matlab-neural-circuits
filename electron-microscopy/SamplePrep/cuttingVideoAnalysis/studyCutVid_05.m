%clear all
colormap gray(256)


[TFN TPN] = GetMyFile;  %Get file location
% vidObj = mmreader([TPN TFN]);  %Get movie information as object
% %vidVals = get(vidObj); %to look at available movie information
% nFrames = vidObj.NumberOfFrames 
% vidHeight = vidObj.Height;
% vidWidth = vidObj.Width;

Info = imfinfo([TPN TFN]);
nFrames = length(Info);
vidHeight = Info(1).Height;
vidWidth = Info(1).Width;

manual = 1;
sampFrames = 1;
subplot(1,1,1)
firstFrames = 10000;
 
%% Find Section
%rawI = double(read(vidObj, [1 min(nFrames,600)]));  %read first few images

rawI = zeros(vidHeight,vidWidth,3,min(nFrames,600));
for i = 1: min(nFrames,600)
    rawI(:,:,:,i) = imread([TPN TFN],'index',i);
end
firstI = uint8(mean(rawI,4));

%% Mark Section
image(firstI),pause(.01)
cF = gcf;
sprintf('click on edges of section and press return')
[gx gy] = ginput;
meany = round(mean(gy));

slice = rawI(meany,:,:,:);
slice = permute(slice,([4 2 3 1]));
image(uint8(slice));

sprintf('mark front of section')
[gxSec gySec] = ginput;
secy = round(gySec(1)/sampFrames);

%%  Create image from line scans
yave = 1;

lx = round(min(gx));
rx = round(max(gx));

width = round(max(gx)-min(gx));
lw = max(1,lx - width);
rw = min(vidWidth,rx+width);
by = meany+fix(yave/2);
ty = by-yave+1;

readF = 1:fix(nFrames/sampFrames);

nFrames = 1000;
strip = zeros(nFrames,rw-lw+1,3);


tic
for i = 1:nFrames
     if ~mod(i,100)
        sprintf('Reading frame %d of %d.', i,nFrames)
    end
   lscan = imread([TPN TFN],'Index',i,'PixelRegion',{[ty by],[lw rw]});
   sy = (i-1)+ 1;
   strip(sy :(sy ),:,:) = lscan;
end

toc
% 
% 
% for i = 1: firstFrames%nFrames
%     if ~mod(i,100)
%         sprintf('Reading frame %d of %d.', i,nFrames)
%     end
%     frame = read(vidObj, i);  %read first few images
%     sy = (i-1)+ 1;
%     strip(sy :(sy ),:,:) = mean(frame(ty:by,lw:rw,:),1);
% end

lx = lx-lw; %rescale boarders to cliped strip
rx = rx-lw;

save([TPN 'strip.mat'],'strip')
image(uint8(strip))

%% average x
mStrip = squeeze(mean(strip,2));
L = squeeze(mean(strip(:,1:lx,:,:),2));
S = squeeze(mean(strip(:,lx:rx,:,:),2));
R = squeeze(mean(strip(:,rx:end,:,:),2));
plot(L(:,1),'r')
hold on
plot(S(:,1),'g')
plot(R(:,1),'b')
hold off

B = mean([mean(L,2) mean(R,2)],2);


%% Filter strip
x = -10 : 10;
sig = 3;
cKern = gaussmf([-100 : 100],[3 0])
cKern = cKern/sum(cKern);

sKern = gaussmf([-100 : 100],[30 0])
sKern = sKern/sum(sKern);
mKern = cKern - sKern;

plot(mKern)

fStrip = filtfilt(mKern,1,B);
plot((fStrip))

[pos vals] = getPeaks(fStrip);

scatter(pos,vals)

lows = pos(vals<mean(vals));
gaps = lows(2:end)-lows(1:end-1);
fsize = mode(gaps)

%% Test frame size
testf = [1:1000];
L = length(fStrip);
n = fix(L/max(testf));

for f = testf
   sfram = reshape(fStrip(1:f*n),f,n);
   prodf(f) = mean(prod(sfram,2));
   %image(sfram * 5),pause
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

[svals idx] = sort(vals,'ascend');
nufVals = svals(fnum);
spos = pos(idx);

testVals = svals(1:fnum);
std(testVals)

sweepPos = sort(spos(1:fnum),'ascend');

scatter(spos(1:fnum),svals(1:fnum),'r');
hold on
scatter(spos(fnum+1:end),svals(fnum+1:end),'b')
hold off
lows = spos(1:fnum)
%% Collect Sections

rel = find((lows - secy)<0);
offset =  secy - lows(rel(end));

sweeps = zeros(fsize+1,size(strip,2),3,length(lows));
for i = 1:length(lows)
    s1 = max(lows(i) - fsize + offset,1);
    s2 = min(lows(i) + offset,size(strip,1));
    s2-s1
    sweeps(1:(s2-s1)+1,:,:,i) = strip(s1:s2,:,:);
end

for i = 1:size(sweeps,4)
    image(uint8(sweeps(:,:,:,i))),pause
end

save([TPN 'sweeps.mat'],'sweeps');

%% Average sections
Lm = squeeze(mean(mean(sweeps(:,1:lx,:,:),1),2));
Sm = squeeze(mean(mean(sweeps(:,lx:rx,:,:),1),2));
Rm = squeeze(mean(mean(sweeps(:,rx:end,:,:),1),2));

Lrat = mean(Lm(1:2,:),1)./Lm(3,:);
Rrat = mean(Rm(1:2,:),1)./Rm(3,:);
Srat = mean(Sm(1:2,:),1)./Sm(3,:);

Lg = mean(Lm(1:2,:),1);
Rg = mean(Rm(1:2,:),1);
Sg = mean(Sm(1:2,:),1);

backBright = mean([median(Lg) median(Rg)]);
backRat = mean([median(Lrat) median(Rrat)]);
secBright = median(Sg);
secRat = median(Srat);

Lg = (Lg- backBright) / (secBright - backBright);
Rg = (Rg- backBright) / (secBright - backBright);
Sg = (Sg- backBright) / (secBright - backBright);


Lrat = (Lrat- backRat) / (secRat - backRat);
Rrat = (Rrat- backRat) / (secRat - backRat);
Srat = (Srat- backRat) / (secRat - backRat);

subplot(2,1,1)
hold off
scatter(1:length(Lg),Lg,'c','.')
plot(Lg,'c')
hold on
scatter(1:length(Rg),Rg,'b','.')
plot(Rg,'b')
scatter(1:length(Sg),Sg,'k','.')
plot(Sg,'k')
scatter(1:length(Srat),Srat,'r','.')
plot(Srat,'r')

%%
 

%imwrite(uint8(strip),[TPN 'strip.tif'],'Compression','none')
