function[lI] = segEWA(I);
%% Matlab 2d segmentor 
%%-Josh Morgan
% 
%% Get pics
TPN = GetMyDir;
inams = getPics(TPN);  %find all tifs

if ~exist([TPN 'labeled'])
    mkdir([TPN 'labeled']);
end
%profile on
info = imfinfo([TPN inams{2}]);
inams{1}
 i = 10;
 I = imread([TPN inams{2}]); 
 I = 255-I;
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

subplot(1,1,1)

%newBorders = fastCon(newBorders,gaus3d([5 5 1],3))>.5;
newBorders = imdilate(newBorders,strel('disk',borderSize));
image(newBorders * 100)
%% Decide on final image

pI = newBorders;
%pI(blankW>0) = 0;

pause(.01)
%% Label
subplot(1,1,1)


pI(1:10,:) = 1;
pI(end -10:end,:) = 1;
pI(:,1:10) = 1;
pI(:,end-10:end) = 1;

lI = bwlabel(~pI,8);
%subplot(2,2,1)
image(I)
%subplot(2,2,2)
image(uint8(cat(3,wI2*1000,c2,uI*1000)))
%subplot(2,2,3)
for i = 1:1
    myCol = hsv(max(lI(:))+1)*40;

    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
    image(I3 * .7 + colO),pause(.01)
end
autoTime = toc
profile off
