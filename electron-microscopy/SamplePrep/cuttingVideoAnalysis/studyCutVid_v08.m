clear all
colormap gray(256)

%% input variables
maxFrames = 1000000000; %debugging limit
sampFrames = 1;

[TFN TPN] = GetMyFile;  %Get file location
vidObj = mmreader([TPN TFN]);  %Get movie information as object
%vidVals = get(vidObj); %to look at available movie information
nFrames = vidObj.NumberOfFrames
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

nam = TFN(1:end-4);

manual = 1;
subplot(1,1,1)
nFrames = min(maxFrames,nFrames);

%% Find Section
startRead = max(1,nFrames/2-300);
rawI = double(read(vidObj, [startRead min(nFrames/2+300,nFrames)]));  %read first few images
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

sprintf('mark front and back of section')
[gxSec gySec] = ginput;
secy = round(gySec/sampFrames)+startRead ;

%%  Create image from line scans
yave = 1;

lx = round(min(gx));
rx = round(max(gx));

width = rx-lx;
lw = max(1,lx - width);
rw = min(vidWidth,rx+width);
by = meany+fix(yave/2);
ty = by-yave+1;

readF = 1:fix(nFrames/sampFrames);
tic
strip = zeros(nFrames,rw-lw+1,3);
for i = 1: nFrames
    if ~mod(i,100)
        sprintf('Reading frame %d of %d.', i,nFrames)
    end
    frame = read(vidObj, i);  %read first few images
    sy = (i-1)+ 1;
    strip(sy :(sy ),:,:) = mean(frame(ty:by,lw:rw,:),1);
end
toc
lx = lx-lw; %rescale boarders to cliped strip
rx = rx-lw;

%image(uint8(strip)) %better to just save.
save([TPN nam '_stripSpace.mat'])
imwrite(uint8(strip),[TPN nam '_strip.tif'],'Compression','none')

%% find sections

exSec = strip(min(secy):max(secy),:,:);
sSize = max(secy)-min(secy);
image(uint8(exSec)),pause(.1)

% 
% tic
% %matches = conv2(strip,exSec,'shape','valid');
% matches = fastCon(mean(strip,3),mean(exSec,3));
% image(matches/max(matches(:))*255)
% plot(mean(matches,2))
% toc

strip1d = squeeze(mean(strip,2));
for c = 1:3
    sStrip = Strip(:,:,c)-


exSec1d = mean(mean(exSec,3),2);
exSec1d = abs(exSec1d-max(exSec1d));
strip1d = mean(mean(strip,3),2);
strip1d = abs(strip1d-max(strip1d));
findEx = filter(exSec1d,1,strip1d)

fStrip = filtfilt(exSec1d,1,strip1d);
%fStrip = filter(exSec1d,1,strip1d);
%fStrip = conv(strip1d,exSec1d)

plot(strip1d/max(strip1d),'r')
hold on
plot(fStrip/max(fStrip),'g')
hold off

%% get peaks
cent = fspecial('gaussian',3,10);
x=-sSize:1:sSize;
cent=gaussmf(x,[1 0]);
sur = gaussmf(x,[sSize/2 0]);
kern = cent/sum(cent)-sur/sum(sur);
plot(x,kern)
xlabel('gaussmf, P=[2 5]')

hStrip = filtfilt(kern,1,fStrip);
plot(hStrip)

%% find peaks top down
[pos vals] = getPeaks(hStrip);
pos = pos(vals>0);
vals = vals(vals>0);
cStrip = hStrip;
for i = 1:length(cStrip)
    targ = find(cStrip == max(cStrip),1);
    loc(i) = targ;
    cStrip(max(1,targ-sSize*.7):min(length(cStrip),targ+sSize * .7)) = 0;
    plot(cStrip),pause(.01)
    if ~sum(cStrip>0),break,end
    
end
loc = sort(loc,'ascend')

%% Collect Sections

refID = find(loc>min(secy) & loc<max(secy),1)
ref = loc(refID);
if isempty(ref)
    ref = min(secy);
end
back = max(secy)-ref;
front = min(secy)-ref;


sweeps = zeros(sSize+1,size(strip,2),3,length(loc));
for i = 1:length(loc)
    s1 = max(loc(i) + front,1);
    s2 = min(loc(i) + back,size(strip,1));
    s2-s1
    sweeps(1:(s2-s1)+1,:,:,i) = strip(s1:s2,:,:);
end

for i = 1:size(sweeps,4)
    image(uint8(sweeps(:,:,:,i))),pause(.1)
end

save([TPN 'sweeps.mat'],'sweeps');

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

kdir = 1;
f = 1;
subplot(2,1,2)
while 1
    clear ink
    ink = lower(input('move: ','s'))
    pause(.1)
    if isempty(ink)
    elseif ink == 'a'
        kdir = -1
        
    elseif ink == 'd'
        kdir = 1
    elseif str2num(ink)
        f = str2num(ink)-kdir;
    end
    f = f+ kdir
    if f> size(sweeps,4), f = size(sweeps,4); end
    if f< 1, f = 1; end
    image(uint8(sweeps(:,:,:,f))),pause(.1)
    
end

