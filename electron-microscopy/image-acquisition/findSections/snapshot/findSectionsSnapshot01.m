%%Find some things in a 2D gray scale

clear all
colormap gray(255)
[TFN TPN] = GetMyFile;
tic

%% Set Variables
estimateWidth = 25;
eW = estimateWidth;

%%  Get Image
Iraw = imread([TPN TFN]);

%%

I =  imresize(Iraw,[size(Iraw,1)*.05,size(Iraw,2)*.05]);

%I = Iraw(4700:5300,2200:3200,:);
[ys xs] = size(I);

I = mean(I,3);
I = max(I(:)) - I;

subplot(1,2,1)
image(fitH(I)),pause(.01)
subplot(1,2,2)


%% difference of gaussians.   for some reason.
% gKern = gaus3d([10 10 1],1,[1 1 1]);
% smI = fastCon(I,gKern);
% gKern = gaus3d([200 200 1],40,[1 1 1]);
% bgI = fastCon(I,gKern);
% difI = (smI - bgI);
% image(fitH(difI))

%% sum edges
threshs= [.01:.01:.03]; % thresolds to check for edges.  later scaled
sizes = .5:.5:1.5;  % sizes to check for edges
sI = sumEdge(I,threshs,sizes,0); % robustish edge detection
image(fitH(sI))



%% Watershed sum Edges

%%Low pass to find wells
dI = fastCon(sI,gaus3d([16 16 1],3));
dI1 = imhmin(sI,.04);
image(fitH(dI1)),pause(.01)

wI = watershed(dI1,4);
cmap = rand(max(wI(:))+1,1)*255+1; %random colormap
image(cmap(wI+1)),pause(.01)

%% Get region properties
xProps = regionprops(wI,dI,'MinorAxisLength','MajorAxisLength',...
    'MeanIntensity','Centroid','MinIntensity','Orientation');
minAx = [xProps.MinorAxisLength];
majAx = [xProps.MajorAxisLength];
meanIntensity = [xProps.MeanIntensity];
ori = [xProps.Orientation];

%% Grab tape
subplot(1,1,1)
tape = (minAx < 120 ) & (minAx > 50) ;
tape = tape & (majAx >200);
tape = tape & ((90 -abs(ori))<10);
tape = find(tape);
tI = I * 0;
for i = 1:length(tape)
  tI(wI == tape(i)) = 1; 
end
image((wI == 0) * 1000 + tI * 100),pause(.01)

%% Edit Tape

emptyReturns = 0;
'Hit Return when tape is labeled'
while emptyReturns < 1
   [y x button] = ginput(1);
              ID = wI(round(x),round(y));
  if button ==1 
           tape = union(tape,ID);
  elseif button == 3
       tape = setdiff(tape,ID);
  elseif isempty(button)
      emptyReturns = emptyReturns + 1;
  end
  
tI = I * 0;
for i = 1:length(tape)
  tI(wI == tape(i)) = 1; 
end
image((wI == 0) * 1000 + tI * 100),pause(.01)

end
tape = tape(tape>0);
%% filter regions

xProps = regionprops(wI,dI,'MinorAxisLength','MajorAxisLength',...
    'MeanIntensity','Centroid','MinIntensity');
minAx = [xProps.MinorAxisLength];
majAx = [xProps.MajorAxisLength];
meanIntensity = [xProps.MeanIntensity];

pI = I * 0;
for i = 1: length(xProps)
    pI(wI == i) = xProps(i).MeanIntensity;
end
image(fitH(pI))

%% Grab sections
sections = (minAx > 10) & (majAx < 50);
secMap = [0 sections];
sections = find(sections);

bigSec = secMap(wI+1);
se = strel('disk',3);
bigSec = imdilate(bigSec,se);
image(bigSec*1000)
bigSec = bwlabel(bigSec);


inTape = [];
for i = 1:max(bigSec(:))
   wTouch = unique(wI(bigSec ==i));
   sTouch = intersect(wTouch,sections);
   tTouch = intersect(wTouch,tape);
   if (length(tTouch) == 1 ) %& (length(sTouch) == 1)
       inTape = [inTape sTouch];
   end
end

sections = intersect(sections,inTape);
showI = I * 0;
for i = 1:length(sections)
   showI(wI == sections(i)) = 1000; 
end
image(showI),pause(.01)

toc

%% Color result
subplot(1,1,1)
gI = I * 0;
cents = cat(1,xProps(sections).Centroid);

gDisk = gaus3d([11 11 1],3);
gDisk = gDisk>gDisk(1,6);
[gy gx] = find(gDisk);
gy = gy - 6; gx = gx-6;
for i = 1:size(cents,1)
    dy = max(1,round(gy + cents(i,2)));
    dx = max(1,round(gx + cents(i,1)));
    dy(dy>ys) = ys; dx(dx>xs) = xs;
    gI(sub2ind(size(gI),dy,dx)) = 1;
end
cgI = I;
cgI(gI>0) = 0;
colI = uint8(cat(3,cgI,cgI + gI * 255,cgI));
image(colI)
   


%% Edit
emptyReturns = 0;
'Hit Return to quit and save'
while emptyReturns < 1
   [y x button] = ginput(1);
   
  if button ==1 
           cents = cat(1,cents,round([y x]));
% elseif button == 2
%            break
  elseif button == 3
           dists = sqrt((cents(:,1)-y).^2 + (cents(:,2)-x).^2);
           if min(dists) < 50;
            keep = ~(dists == min(dists));
            cents = cents(keep,:);
           end
  elseif isempty(button)
      emptyReturns = emptyReturns + 1;
  end

gI = gI * 0;   
for i = 1:size(cents,1)
    dy = max(1,round(gy + cents(i,2)));
    dx = max(1,round(gx + cents(i,1)));
    gI(sub2ind(size(gI),dy,dx)) = 10000;
end
cgI = I;
cgI(gI>0) = 0;
colI = uint8(cat(3,cgI,cgI + gI * 255, cgI));
image(colI),pause(.01)
   
end

%% save labeled

newTFN = ['labeled_' TFN];
%imwrite(colI,[TPN newTFN],'Compression','none');