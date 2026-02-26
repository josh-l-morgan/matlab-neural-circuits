subplot(1,1,1)
colormap gray(256)


[TFN TPN] = GetMyFile;  %Get file location
vidObj = mmreader([TPN TFN]);  %Get movie information as object
%vidVals = get(vidObj); %to look at available movie information
nFrames = vidObj.NumberOfFrames 
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

manual = 1;


%% Find Section
rawI = double(read(vidObj, [1 min(nFrames,200)]));  %read first few images
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
secy = round(gySec(1));

%%  Create image from line scans
yave = 1;

lx = round(min(gx));
rx = round(max(gx));

width = round(max(gx)-min(gx));
lw = max(1,lx - width);
rw = min(vidWidth,rx+width);
by = meany+fix(yave/2);
ty = by-yave+1;

strip = zeros(nFrames * ysamp,rw-lw+1,3);
for i = 1: nFrames
    frame = read(vidObj, i);  %read first few images
    sy = (i-1) * ysamp + 1;
    strip(sy :(sy + ysamp - 1),:,:) = mean(frame(ty:by,lw:rw,:),1);
end

lx = lx-lw; %rescale boarders to cliped strip
rx = rx-lw;

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

ratL = mean(L(:,1:2),2)./L(:,3);
ratS = mean(S(:,1:2),2)./S(:,3);
ratR = mean(R(:,1:2),2)./R(:,3);
plot(ratL(:,1),'r')
hold on
plot(ratS(:,1),'g')
plot(ratR(:,1),'b')
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
std(vals)/3
hist(vals,std(vals/3))
lows = pos(vals<mean(vals));
gaps = lows(2:end)-lows(1:end-1);
fsize = mode(gaps)

%% Collect Sections

rel = find((lows - secy)<0);
offset =  secy - lows(rel(end));

clear sweeps
for i = 1:length(sweepPos)
    s1 = max(lows(i) - fsize + offset,1);
    s2 = min(lows(i) + offset,size(strip,1));
    sweeps(:,:,:,i) = strip(s1:s2,:,:);
end

for i = 1:size(sweeps,4)
    image(uint8(sweeps(:,:,:,i))),pause(.01)
end

%% Average sections
Lm = squeeze(mean(mean(sweeps(:,1:lx,:,:),1),2));
Sm = squeeze(mean(mean(sweeps(:,lx:rx,:,:),1),2));
Rm = squeeze(mean(mean(sweeps(:,rx:end,:,:),1),2));

Lrat = mean(Lm(1:2,:),1)./Lm(3,:);
Rrat = mean(Rm(1:2,:),1)./Rm(3,:);
Srat = mean(Sm(1:2,:),1)./Sm(3,:);

Lg = mean(Lm,1);
Rg = mean(Rm,1);
Sg = mean(Sm,1);

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
plot(Lrat,'r')
hold on
plot(Rrat,'g')
plot(Srat,'k')

subplot(2,1,2)
plot(Lg,'r')
hold on
plot(Rg,'g')
plot(Sg,'k')



%imwrite(uint8(strip),[TPN 'strip.tif'],'Compression','none')
