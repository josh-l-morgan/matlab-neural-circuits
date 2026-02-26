

clear all
SPN = 'E:\Bassnet\ReflectionAlpha\speckleBoundary.oif.files\'
shouldWrite = 0;
voxWidth = 0.31;
voxDepth = 1;
aspect = voxDepth/voxWidth;

dSPN = dir([SPN '*.tif']);
nams = {dSPN.name};
iNum = length(nams);

%% Parse color z
c = zeros(iNum,1);
z = zeros(iNum,1);
for i = 1:iNum
    nam = nams{i};
    C = regexp(nam,'C');
    C = C(end);
    Z = regexp(nam,'Z');
    Z = Z(end);
    dot = regexp(nam,'.tif');
    c(i) = str2num(nam(C+1:Z-1));
    z(i) = str2num(nam(Z+1:dot-1));
end

zNum = max(z);
cNum = max(c);
zs = zNum;

%% Read in images
imfo = imfinfo([SPN nams{1}]);
ys = imfo.Height;
xs = imfo.Width;
Ir = zeros(ys,xs,zNum,cNum,'double');

for i = 1:iNum
    I = double(imread([SPN nams{i}]));
    Ir(:,:,z(i),c(i)) = I;
end


%% Generate max
Imax = squeeze(max(Ir,[],3));

Idif = Ir(:,:,:,2)./Ir(:,:,:,3);

Fmax = max(Idif,[],3);
Fmean = mean(Idif,3);
image(Fmax*60)
image(Fmean*100)
Icol = cat(3,Fmean*200, Imax(:,:,1)*.07, Imax(:,:,3)*.07);

subplot(5,5,[1:4 6:9 11:14 16:19])

Icol(:,:,1) = Fmean*50/median(Fmean(:));
if 0
Ic1 = Fmean;
Ic2 = Imax(:,:,1);
Ic3 = Imax(:,:,3);
Ic1 = Ic1 * 150/median(Ic1(:));
Ic2 = Ic2 * 100/median(Ic2(:));
Ic3 = Ic3 * 30/median(Ic3(:));


else
    Ic1 = Imax(:,:,3);
    Ic2 = Imax(:,:,2);
    Ic3 = Imax(:,:,1);
    Ic1 = Ic1 * 150/median(Ic1(:));
    Ic2 = Ic2 * 100/median(Ic2(:));
    Ic3 = Ic3 * 30/median(Ic3(:));
end


Icol(:,:,1) = Ic1;
Icol(:,:,2) = Ic2;
Icol(:,:,3) = Ic3;

image(uint8(Icol))




%% median filter
Im = Ir;
for i = 1:zNum
    fprintf('median filtering %d of %d\n',i,zNum)
    for c = 1:3
        Im(:,:,i,c)   = medfilt2(Ir(:,:,i,c),[3 3]);
    end
end

subplot(5,5,[5 10 15 20])
if 0
for i = 1:xs
    Icol = squeeze(Im(:,i,:,:));
    Ic1 = Icol(:,:,1);
    Ic2 = Icol(:,:,2);
    Ic3 =Icol(:,:,3);
    Ic1 = Ic1 * 50/median(Ic1(:));
    Ic2 = Ic2 * 50/median(Ic2(:));
    Ic3 = Ic3 * 50/median(Ic3(:));

    Icol(:,:,1) = Ic3;
    Icol(:,:,2) = Ic2;
    Icol(:,:,3) = Ic1;
    image(uint8(Icol));
    drawnow
end
end

ImedX = squeeze(median(Im,2));
image(uint8(ImedX*.4))

%% median of max plot
Imax = max(Im,[],3);
Imed = squeeze(median(Imax,2));
subplot(5,1,5)
cla
hold on
cols = {'r' 'g' 'b'};
for c = 1:size(Imed,2)
    p = Imed(:,c);
    p = p-min(p);
    p = p/max(p);
    Imed(:,c) = p;
    plot(Imed(:,c),cols{c})
end


%% find surface
Ig1 = imgaussfilt3(Im(:,:,:,1),3);
Ig2 = imgaussfilt3(Im(:,:,:,2),3);
Ig3 = imgaussfilt3(Im(:,:,:,3),3);
IgD = Ig2./Ig3;
IgDt = IgD-min(IgD(:));
IgDt = IgDt * 300/max(IgDt(:));


Ig1max = max(Ig1,[],3);
Ig2max = max(Ig2,[],3);
Ig3max = max(Ig3,[],3);

surfV = Ig2;
surfVmax = max(surfV,[],3);
maxInd = find(surfV == repmat(surfVmax,[1 1 zNum]));
[y x z] = ind2sub(size(surfV),maxInd);
yxInd = sub2ind([ys xs],y,x);
uInd = unique(yxInd);
cInd = hist(yxInd,uInd);
mult = find(cInd>1);
sing = find(cInd==1);
maxZ = zeros(length(uInd),1);
maxZ(yxInd) = z; %insert zs. will work for single results  
%%insert z for locations where multiple zs were at max intensity
for i = 1:length(mult)
    hit = find(yxInd==(cInd(mult(i))));
    maxZ(mult(i)) = mean(z(hit));
end
Iz = zeros(ys,xs);
Iz(:) = maxZ;
Izt = Iz-min(Iz(:));
Izt = Izt * 256/max(Izt(:));


zp = median(Izt,2);
zp = zp-min(zp);
zp = zp/max(zp);

subplot(5,1,5)
hold on
plot(zp,'k')
cols = {'r' 'g' 'b'};
for c = 1:size(Imed,2)
    plot(Imed(:,c),cols{c})
end


%% Writing
if shouldWrite
    TPN = 'E:\Bassnet\ReflectionAlpha\writeOutput2\';
    if ~exist(TPN,'dir'),mkdir(TPN),end
end

%% Write side view series
if 0
    binRad = 5;
    sampFreq = 10;
    sampPlane = [1:sampFreq:xs];
    for i = 1:length(sampPlane)
        p = sampPlane(i);
        s1 = max(1,p-binRad);
        s2 = min(xs,p+binRad);
        iNam = sprintf('sec%03.0f.tif',i);
        I1 = squeeze(mean(Ir(:,s1:s2,:,1),2)) * .2;
        I2 = squeeze(mean(Ir(:,s1:s2,:,2),2)) * .1;
        I3 = squeeze(mean(IgDt(:,s1:s2,:),2)) * 2;


        Ic = cat(3,I3,I1,I2);
        Ic2 = imresize(Ic,[ys zNum * aspect]);
        imshow(uint8(Ic2))
        drawnow
        imwrite(uint8(Ic2),[TPN iNam]);

    end
end


%% Write top Down view series
if 0
    binRad = 0;
    sampFreq = 1;
    sampPlane = [1:sampFreq:zs];
    for i = 1:length(sampPlane)
        p = sampPlane(i);
        s1 = max(1,p-binRad);
        s2 = min(xs,p+binRad);
        iNam = sprintf('sec%03.0f.tif',i);
        I1 = squeeze(mean(Ir(:,:,s1:s2,1),3)) * .2;
        I2 = squeeze(mean(Ir(:,:,s1:s2,2),3)) * .1;
        I3 = squeeze(mean(IgDt(:,:,s1:s2),3)) * 2;

        Ic = cat(3,I3,I1,I2);
        Ic2 = Ic;%imresize(Ic,[ys zNum * aspect]);
        imshow(uint8(Ic2))
        drawnow
        imwrite(uint8(Ic2),[TPN iNam]);

    end
end






