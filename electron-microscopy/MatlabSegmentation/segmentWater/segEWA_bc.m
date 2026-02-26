%function[lI] = segEWA(I);
%% Matlab 2d segmentor
%%-Josh Morgan
%% Get pics
TPN = GetMyDir;
inams = getPics(TPN);  %find all tifs

if ~exist([TPN 'labeled'])
    mkdir([TPN 'labeled']);
end
profile on
info = imfinfo([TPN inams{1}]);
i = 10;
I = imread([TPN inams{1}]);
%I = I(900:1200,900:1200);

%% Variables

%%Watershed variables
gKern2 = gaus3d([3 3 1],1); %df = ([ 20 20 1], 3)
watMin =30; %df = 10, minimum intensity change for watershed

%%Edge varibles
gKern1 = gaus3d([3 3 1],1); %kernal for edge convolution
threshs= [ .5 .01];%[.8:.1:.8]; % df = [.3:.1:.5], thresolds to check for edges.  later scaled
sizes = [ 2 5];  % sizes to check for edges
minEdge = 40; % default 50, how small can an edge be

%%Trim Watershed variables
searchRadius = 20; %df = 20
middleSize = 3; % size of sampling area at wattershed border
grabNearest = 30; %number of edge values to grab
borderSize =3; %dilattion of border to fill small holes

%% Read Data
colormap gray(256)
tic
[ys, xs] = size(I);

I3 = cat(3,I,I,I); %Color image

subplot(1,2,1)
image(uint8(I))
subplot(1,2,2)

%% Compare newBorders to BreadCrumbs

load([TPN 'binCrumb.mat'])

obs = binCrumb.obs;
clear binCrumb

%% Transform crumbs
% 150 - z, x * 2, y * 2

mins = [ 10 ^10 10^10 10^10];
maxs = [1 1 1];
for o = 1:length(obs)
    ob = obs{o};
    ob(:,1) = 151 - ob(:,1);
    ob(:,2:3) = ob(:,2:3) * 2;
    obs{o} = ob;

    mins = min(min(ob,[],1),mins);
    maxs = max(max(ob,[],1),maxs);
end

%% turn into planes

planes = cell(maxs(1),1);

for o = 1:length(obs)
    ob = obs{o};
    for i = 1:size(ob,1)
        plane = planes{ob(i,1)};
        plane = [plane; o ob(i,2:3)];
        planes{ob(i,1)} = plane;
    end
end



%% Watershed
subplot(1,2,2)
c2 = fastCon2d(I,gKern2);
c2 = max(c2(:)) - c2;
minC2 = imhmin(c2,watMin);

wI2 = watershed(minC2,8);
image((wI2==0)*1000 + c2 * .8);
pause(.01)


%% Edge
subplot(1,2,2)

c1 = fastCon2d(I,gKern1);
%image(fitH(c1))
%eI = sumEdge(c1,threshs,sizes,0); %edge detection
eI = c1* 0;
for i = 1:length(threshs);
    eI = eI + edge(c1,'Canny',[.001 threshs(i)],sizes(i)); %edge detection
end

image(fitH(c1)*.5 + eI*1000)
pause(.01)

%%Filter Edges
subplot(1,2,2)
eI = bwlabel(eI,8);
eProps = regionprops(eI,'MinorAxisLength','MajorAxisLength');
minAx = [eProps.MinorAxisLength];
majAx = [eProps.MajorAxisLength];

use = (majAx > minEdge);
use = [0 use];
uI = use(eI + 1);

image(fitH(c1)*.5 + uI*1000)
pause(.01)




%% Trim Water2


[wy wx] = find(wI2==0);  %grab bordrs
sRad = searchRadius; % define search radius

%%Make distance Kernal
kSize = sRad * 2+1;
dKern = zeros(kSize);
[y x] = ind2sub(size(dKern),find(dKern==0));
dKern(:) = sqrt(((y-sRad-1)).^2 + ((x - sRad -1)).^2);

%Make sub Kernals
[y x ] = find(dKern<3);
my = y  - sRad-1;
mx = x - sRad -1;
[y x] = find(dKern<=sRad);
dList = dKern(dKern<=sRad);
ny = y - sRad -1;
nx = x - sRad -1;


%% Pad
uI = pad(uI,sRad);
c2 = pad(c2,sRad);
wy = wy + sRad;
wx = wx + sRad;
newBorders = double(uI * 0);
pys = size(uI,1);
for i = 1:length(wy)
    %transform search kernal
    dy = wy(i) + ny;
    dx = wx(i) + nx;
    sInd = dy + (dx-1) * pys;

    sEdge = uI(sInd);
    dVal = dList(sEdge>0); %Grab distances
    if ~isempty(dVal) %if there are edges

        [dVal idx] = sort(dVal);
        nidix = idx(dVal <= dVal(min(grabNearest,length(dVal))));
        sVal = c2(sInd); %get brightness
        eVal = sVal(sEdge>0); %get edge brightness
        useVal = eVal(nidix);

        %Get middle values
        dy = wy(i) + my;
        dx = wx(i) + mx;
        mVal = c2(dy + (dx-1) * pys);
        %newBorders(wy(i),wx(i)) = max(mVal) - median(sVal);

        newBorders(wy(i),wx(i)) = median(mVal)-median(useVal);
    end
end

%%unPad
uI = unpad(uI,sRad);
c2 = unpad(c2,sRad);
newBorders = unpad(newBorders,sRad);



%% Use breadCrumbs
subplot(1,1,1)
plane = planes{1}; pid = plane(:,1);
pind = plane(:,3) + (plane(:,2)-1) * 2048;
Ip = zeros(2048,2048);
Ip(pind) = pid;
Ip = imdilate(Ip,strel('disk',2));
%Ip = Ip(900:1200,900:1200);
image(I+uint8(Ip)*1000)

pI = newBorders;
pI(~wI2) = pI(~wI2)-min(pI(~wI2))+1;


se = strel('disk',1);
M= zeros(ys,xs);
allM = zeros(ys,xs,10);
[Lab numLab] = bwlabel(~pI,4);
for i = 1:100
    finished = 1;
    for o = 1: numLab
        if ~mod(o,100)
            sprintf('Checking object %d of %d',o,numLab)
        end
        L = Lab ==o;
        vals = Ip(L);
        ids = vals(vals>0);
        if isempty(ids)
            finished = 0;
            D = imdilate(L,se);

            E = imerode(D,se);
            Did = find(D & ~L);
            pI(Did(pI(Did) == min(pI(Did)))) = 0;
            pI(E) = 0;
        end

    end
    [Lab numLab] = bwlabel(~pI,4);
    if finished, break,end
end

lI = Lab * 0;
conflict = 0;
for i = 1:numLab
    L = Lab ==i;
        vals = Ip(L);
        ids = vals(vals>0);
        lI(L) = mode(ids);
        if length(ids)>1
            conflict = conflict + 1;
        end
        
end


%% Label

subplot(2,2,1)
image(I)
subplot(2,2,2)
image(uint8(cat(3,wI2*1000,c2,uI*1000)))
subplot(2,2,3)
for i = 1:1
    myCol = hsv(max(lI(:))+1)*40;

    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
    image(I3 * .7 + colO + uint8(cat(3,Ip,Ip,Ip)) * 1000),pause(.01)
end
autoTime = toc
profile off

ShowMess = I3 * .7 + colO + uint8(cat(3,Ip,Ip,Ip)) * 1000;
mkdir([TPN 'pics']);
imwrite(ShowMess,[TPN 'pics\Showmess.tif'],'Compression','none')