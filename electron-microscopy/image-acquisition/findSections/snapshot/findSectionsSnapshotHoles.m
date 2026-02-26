%%Find some things in a 2D gray scale

clear all
colormap gray(255)
[TFN TPN] = GetMyFile;
%% Set Variables
tapeMin = 225;
tapeMax = 275;
majThresh = 60;
minThresh = 10;

%%  Get Image
Iraw = imread([TPN TFN]);
tic

%%
scl = 3000/max(size(Iraw));
I =  imresize(double(Iraw),[size(Iraw,1)*scl,size(Iraw,2)*scl]);
[ys xs] = size(I);
I = max(I(:)) - I;

subplot(1,2,1)
image(I),pause(.01)
subplot(1,2,2)

%% sum edges
threshs= [.05:.05:.05]; % thresolds to check for edges.  later scaled
sizes = 1:1:1;  % sizes to check for edges
sI = sumEdge(I,threshs,sizes,2); % robustish edge detection
image(fitH(sI))

%% Watershed sum Edges

dI1 = imhmin(sI,.1);
image(fitH(dI1)),pause(.01)

wI = watershed(dI1,8);
%edgeVal = max(wI(:))+1;
wI(wI == 0) = 1;
cmap = rand(max(wI(:))+1,1)*255+1; %random colormap
image(cmap(wI+1)),pause(.01)

%% Get region properties
xProps = regionprops(wI,I,'MinorAxisLength','MajorAxisLength',...
    'MeanIntensity','Centroid','MinIntensity','Orientation',...
    'PixelIdxList','BoundingBox');
minAx = [xProps.MinorAxisLength];
majAx = [xProps.MajorAxisLength];
meanIntensity = [xProps.MeanIntensity];
ori = [xProps.Orientation];

%% Grab tape
subplot(1,1,1)
tape = (minAx < tapeMax ) & (minAx > tapeMin) ;
tape = tape & (majAx >tapeMax);
tape = tape & ((90 -abs(ori))<20);
tape = find(tape);
tI = I * 0;
tapeIds = cat(1,xProps(tape).PixelIdxList);
tI(tapeIds) = 1;
image((wI == 0) * 1000 + tI * 100),pause(.01)
tapeVal = median(meanIntensity(tape));
tapeHoles = bwlabel(tI==0,4)>1;

%% Edit Tape

% emptyReturns = 0;
% 'Hit Return when tape is labeled'
% while emptyReturns < 1
%    [y x button] = ginput(1);
%               ID = wI(round(x),round(y));
%   if button ==1 
%            tape = union(tape,ID);
%   elseif button == 3
%        tape = setdiff(tape,ID);
%   elseif isempty(button)
%       emptyReturns = emptyReturns + 1;
%   end
%   
% tI = I * 0;
% for i = 1:length(tape)
%   tI(wI == tape(i)) = 1; 
% end
% image((wI == 0) * 1000 + tI * 100),pause(.01)
% 
% end
% tape = tape(tape>0);
%% filter regions

sizeSec =(minAx> minThresh) & (majAx< majThresh);
potSec = find(sizeSec & (meanIntensity>tapeVal));
secPix = cat(1,xProps(potSec).PixelIdxList);
secPixID = wI(secPix);
secPixHole = tapeHoles(secPix);
secInHole = unique(secPixID(secPixHole));

secPix = cat(1,xProps(secInHole).PixelIdxList);
pI = I * 0;
pI(:) = meanIntensity(wI(:));
pI(secPix) = 0;
image(pI)

sections = secInHole;

%% Color result
subplot(1,1,1)
gI = I * 0;
cents = cat(1,xProps(sections).Centroid);

gDisk = gaus3d([21 21 1],3);
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
% emptyReturns = 0;
% 'Hit Return to quit and save'
% while emptyReturns < 1
%    [y x button] = ginput(1);
%    
%   if button ==1 
%            cents = cat(1,cents,round([y x]));
% % elseif button == 2
% %            break
%   elseif button == 3
%            dists = sqrt((cents(:,1)-y).^2 + (cents(:,2)-x).^2);
%            if min(dists) < 50;
%             keep = ~(dists == min(dists));
%             cents = cents(keep,:);
%            end
%   elseif isempty(button)
%       emptyReturns = emptyReturns + 1;
%   end
% 
% gI = gI * 0;   
% for i = 1:size(cents,1)
%     dy = max(1,round(gy + cents(i,2)));
%     dx = max(1,round(gx + cents(i,1)));
%     gI(sub2ind(size(gI),dy,dx)) = 10000;
% end
% cgI = I;
% cgI(gI>0) = 0;
% colI = uint8(cat(3,cgI,cgI + gI * 255, cgI));
% image(colI),pause(.01)
%    
% end

%% save labeled

newTFN = ['labeled_' TFN];
%imwrite(colI,[TPN newTFN],'Compression','none');

%% Grab High res

for i = 1:length(sections)

    
    bBox = xProps(sections(i)).BoundingBox;
    bBox(1:2) = bBox(1:2)-5;
    bBox(3:4) = bBox(3:4)+10;
    bBox = round(bBox/scl);
    sPic = Iraw(bBox(2):bBox(2)+ bBox(4),bBox(1):bBox(1)+bBox(3));
    image(sPic),pause(.01)
    sPics{i} = sPic;
end


toc


