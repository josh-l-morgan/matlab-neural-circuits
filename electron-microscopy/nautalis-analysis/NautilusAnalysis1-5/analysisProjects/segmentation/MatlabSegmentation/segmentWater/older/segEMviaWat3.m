colormap gray(256)
tic
[TFN TPN] = GetMyFile;


Iraw = imread([TPN TFN]);
I = Iraw;
I = I(500:1200,600:1500);
I3 = cat(3,I,I,I);

subplot(1,2,1)

image(I)
subplot(1,2,2)
%I = 255-I;
[ys, xs] = size(I);

%% Filter for membranes
gKern1 = gaus3d([15 15 1],3);
c1 = fastCon(I,gKern1);
image(fitH(c1))

gKern2 = gaus3d([200 200 1],30);
c2 = fastCon(I,gKern2);
image(fitH(c2))

difI = c1 - c2;
image(fitH(difI))
pause(.01)
%% Edge
gKern1 = gaus3d([15 15 1],2);
c1 = fastCon(I,gKern1);
image(fitH(c1))
eI = edge(c1,'Canny',[.000005 .3],3);
image(fitH(eI*1000+c1))
pause(.01)
%% close dark
% eI(1:10,:) = 0;
% eI(end -10:end,:) = 0;
% eI(:,1:10) = 0;
% eI(:,end-10:end) = 0;
cI = eI;
sRad = 10;
for r = 1:1
    pI = bwperim(cI);
    edgeIDs = find(eI);
    diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],3);
    diskE = diskE>=diskE(sRad,1);
    [dy dx] = find(diskE);
    yx = [dy dx]-sRad;
    for i = 1:length(edgeIDs)
        dVal = c1(edgeIDs(i));
        [cy cx] = ind2sub(size(I),edgeIDs(i));
        y = yx(:,1) + cy; x = yx(:,2) + cx;
        uyx = wall([y x],size(I));
        iyx = sub2ind(size(I),uyx(:,1),uyx(:,2));
        sVal = c1(iyx);  % get values of surround
        cI(iyx(sVal<dVal)) = 1;   %set dark surround to 1
        %image(fitH(cI + eI)),pause
    end
    image(fitH(cI + eI)),pause(1)
end
%% Watershed
gKern2 = gaus3d([10 10 1],3);
c2 = fastCon(I,gKern2);

image(fitH(c2))
c2 = max(c2(:)) - c2;
image(fitH(c2))
minC2 = imhmin(c2,20);

wI = watershed(minC2,8);
image((wI==0)*1000 + c2 * .8);
pause(.01)

%%
catCol =uint8(cat(3,eI * 1000, wI * 1000, cI * 1000));
image(catCol)
pause(.01)
%%

pI = I * 0;
pI = wI==0;
pI = imdilate(pI,strel('disk',4));
%pI(eI) = 0;
pI(cI==0) = 0;
image(pI * 100)
%Arbitrary trimming
% image((wI==0) * 100 + (xI)*100)
% xI(wI==0) = 0;
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
%%Left click fuse, right click cut

%%MakeCuttingdisk
sRad = 3;
diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],2);
diskE = diskE>=diskE(sRad,1);
pDisk = bwperim(diskE);
diskE(pDisk) = 0;
[dy dx] = find(diskE);
[py px] = find(pDisk);

while 1
[mx my button] = ginput;
lid = sub2ind(size(I),round(my),round(mx));
labs = sort(lI(lid));
labs = labs(labs>0);
if numel(button) >0 
if button(1) == 1
for i = 2:length(labs)
    lI(lI == labs(i)) = labs(1);
end
elseif button(1) ==3
 
    %toOpen = imopen(toOpen,strel('disk',1));
    clearInd = []; grabInd = [];
    for s = 1:length(mx)-1;
       step = mx(s+1)-mx(s);
       step = abs(step)/(step * 3);
        lx = mx(s) :step : mx(s+1);
       lSlope = (my(s+1)-my(s))./(mx(s+1)-mx(s));
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
    lI(toOpen == newObs(1)) = max(lI(:))+1;
end

end %if button

myCol = hsv(max(lI(:))+1)*40;

[r rix] = sort(rand(size(myCol,1)),1);
myCol= myCol(rix,:);
myCol = cat(1,[0 0 0],myCol);
red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
skipCol = 5 + round(rand*10);

colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
image(I3 * .7 + colO)

end
%%
