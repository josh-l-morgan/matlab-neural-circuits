function[lI] = segEWA(I);
%% Matlab 2d segmentor 
%%-Josh Morgan


%% Variables

%%Edge varibles
gKern1 = gaus3d([5 5 1],.5); %kernal for edge convolution
threshs= [.3:.1:.5]; % df = [.3:.1:.5], thresolds to check for edges.  later scaled
sizes = 1:1:5;  % sizes to check for edges

%%Watershed variables
gKern2 = gaus3d([15 15 1],2); %df = ([ 20 20 1], 3)
watMin = 10; %df = 10, minimum intensity change for watershed

%%Trim Watershed variables
searchRadius = 20; %df = 20

%% Read Data
colormap gray(256)
tic
%I = I(5500:6000,6500:7000);
[ys, xs] = size(I);

I3 = cat(3,I,I,I); %Color image

subplot(1,2,1)
image(I)
subplot(1,2,2)

%% Edge
subplot(1,2,1)

c1 = fastCon(I,gKern1);
%image(fitH(c1))
eI = sumEdge(c1,threshs,sizes,0); % robustish edge detection
image(fitH(c1)*.5 + eI*1000)
pause(.01)
%%Filter Edges
subplot(1,2,2)
eI = bwlabel(eI,8);
eProps = regionprops(eI,'MinorAxisLength','MajorAxisLength');
minAx = [eProps.MinorAxisLength];
majAx = [eProps.MajorAxisLength];

use = (majAx > 50);
use = [0 use];
uI = use(eI + 1);

image(fitH(c1)*.5 + uI*1000)
pause(.01)

%% Watershed sum Edges
%{
%theres got to be a use for this
subplot(1,2,1);

%%Low pass to find wells
dI1 = fastCon(eI,gaus3d([6 6 1],1));%dI1 = imhmin(dI,.04);
%dI1 = imhmin(dI1,.000005);
image(fitH(dI1)),pause(.01)

wI = watershed(dI1,8);
wProps = regionProps(wI,I,'PixelIdxList','MeanIntensity');
mi = [0 wProps.MeanIntensity];
image(fitH(mi(wI+1))),pause(.01)

%% Sort Wat1 

subplot(1,2,2)
[wy wx] = find(wI==0);
near = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 1; 1 -1];
allCon = zeros(length(wy) * 3,2);
c = 0;
for i = 1: length(wy)
    nyx = wall([wy(i) + near(:,1) wx(i) + near(:,2)],size(I));
    cons = unique(wI(sub2ind(size(wI),nyx(:,1),nyx(:,2))));
    cons = cons(cons>0);
    conp = combntns(cons,2);
    allCon(c+1:c+size(conp,1),:) = conp;
    c = c + size(conp,1);
end
allCon = allCon(1:c,:);
allCon = unique(allCon,'rows');

fR = sort(unique(allCon));
%%
passReg = I * 0;
pass = zeros(length(fR)+1,1);
for i = 1: length(fR)
   [fy fx] = find(allCon == fR(i));
   rVal = wProps(fR(i)).MeanIntensity;
   surR = allCon(fy,:);
   surR = setdiff(surR(:),fR(i));
   sVal = [wProps(surR).MeanIntensity]; 
   pass(i+1) = ( rVal - min(sVal)) > (max(sVal)-rVal);
   %pass(i+1) = (sum(sVal>rVal)/length(sVal))>=.5;
end

image(pass(wI+1)*50)
blankW = pass(wI + 1);

%% Watershed Dark
subplot(1,1,1)
gKern2 = gaus3d([20 20 1],1);
c2 = fastCon(I,gKern2);

image(fitH(c2))
c2 = c2;
image(fitH(c2))
minC2 = imhmin(c2,10);

wI2 = watershed(minC2,8);
image((wI2==0)*1000 + c2 * .8);
pause(.01)

%w2Props = regionProps(wI2,I,'MeanIntensity');
%}
%% Watershed
subplot(1,2,2)

c2 = fastCon(I,gKern2);

image(fitH(c2))
c2 = max(c2(:)) - c2;
image(fitH(c2))
minC2 = imhmin(c2,watMin);

wI2 = watershed(minC2,8);
image((wI2==0)*1000 + c2 * .8);
pause(.01)

%w2Props = regionProps(wI2,I,'MeanIntensity');
%}

%% Trim Water2
[wy wx] = find(wI2==0);
sRad = searchRadius;
sur = gaus3d([sRad * 2 + 1 sRad * 2 + 1 1],3);
mid = sur>=sur(sRad,sRad-2);
sur = (sur>=sur(sRad,1)) & sur<sur(sRad,3);
[y x ] = find(mid);
myx = [y x] -sRad;
[y x] = find(sur);
yx = [y x] - sRad;
newBorders = double(I * 0);
for i = 1:length(wy)
    dy = wy(i) + yx(:,1);
    dx = wx(i) + yx(:,2);
    dyx = wall([dy dx],size(I));
    sInd = sub2ind(size(I),dyx(:,1),dyx(:,2));
    sVal = c2(sInd);
    sEdge = uI(sInd);
    eVal = sVal(sEdge>0);
    
    dy = wy(i) + myx(:,1);
    dx = wx(i) + myx(:,2);
    dyx = wall([dy dx],size(I));
    mVal = c2(sub2ind(size(I),dyx(:,1),dyx(:,2)));
    
    
    %newBorders(wy(i),wx(i)) = max(mVal) - median(sVal);
    if ~isempty(eVal)
        newBorders(wy(i),wx(i)) = max(mVal)-median(eVal); 
        %change to median(mVal) to reduce oversegmentation
    end
end
subplot(1,1,1)

%newBorders = fastCon(newBorders,gaus3d([5 5 1],3))>.5;
newBorders = imdilate(newBorders,strel('disk',3));
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

for i = 1:10
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
%% Get mouse inputs
%{
%%Left click fuse, right click cut

%%MakeCuttingdisk
sRad = 3;
diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],2);
diskE = diskE>=diskE(sRad,1);
pDisk = bwperim(diskE);
diskE(pDisk) = 0;
[dy dx] = find(diskE);
[py px] = find(pDisk);

flipCol = 0;

while 1
   
    [mx my button] = ginput;
    if length(button) >1
    myx = wall([my mx], size(I));
    mx = myx(:,2); my = myx(:,1);
    lid = sub2ind(size(I),round(my),round(mx));
    labs = sort(lI(lid));
    labs = labs(labs>0);
    if (numel(button) >0 ) 
        if button(1) == 1
            for i = 2:length(labs)
                lI(lI == labs(i)) = labs(1);
            end
        elseif (button(1) == 2)
            'you hit two'
%             fr = imfreehand(gca,'Closed',0)
%             position = wait(fr)
%             pos = getPosition(fr)
%             delete(fr)
% 
% 
%             mx = pos(:,2);
%             my = pos(:,1);
%             button = 1;
        elseif (button(1) ==3 )

            %toOpen = imopen(toOpen,strel('disk',1));
            clearInd = []; grabInd = [];
            for s = 1:length(mx)-1;
                L = sqrt((mx(s+1)-mx(s)).^2 + (my(s+1)-my(s)).^2);
                step = (mx(s+1)-mx(s))/L * .3;
                lSlope = (my(s+1)-my(s))/(mx(s+1)-mx(s));
                lx = mx(s) :step: mx(s+1);
                ly = (lx-mx(s)) * lSlope + my(s);
                for l = 1: length(ly)
                    lyx = wall(round([dy+ly(l) dx + lx(l)]),size(I));
                    clearInd = cat(1,clearInd,sub2ind(size(I),lyx(:,1),lyx(:,2)));
                    gyx = wall(round([py+ly(l) px + lx(l)]),size(I));
                    grabInd = cat(1,grabInd,sub2ind(size(I),gyx(:,1),gyx(:,2)));

                end
            end

            %%Break and label
            toOpen = lI * 0;

            effected = unique(lI(grabInd));
            effected = effected(effected>0);

            toOpen(ismember(lI,effected)) = 1;
            image(fitH(toOpen))
            toOpen(clearInd) = 0;
            toOpen = bwlabel(toOpen,4);
            newObs = unique(toOpen(grabInd));
            newObs = newObs(newObs>0);
            lI(clearInd) = 0;
            for nO = 2:length(newObs)
                lI(toOpen == newObs(nO)) = max(lI(:))+1;
            end
        end

    end %if button
    end
    
    flipCol = ~flipCol;
    myCol = hsv(max(lI(:))+1) * 40 * flipCol;


    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
    image(I3 * .7 + colO)
%}
end

