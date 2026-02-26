%% Matlab 2d segmentor 
%%-Josh Morgan



%% Read Data
colormap gray(256)
tic
[TFN TPN] = GetMyFile;

Iraw = imread([TPN TFN]);
I = Iraw;
%I = I(5500:6000,6500:7000);
[ys, xs] = size(I);

I3 = cat(3,I,I,I); %Color image

subplot(1,2,1)
image(I)
subplot(1,2,2)

%% Edge
gKern1 = gaus3d([5 5 1],1);
c1 = fastCon(I,gKern1);
image(fitH(c1))
eI = edge(I,'Canny',[.000005 .2],2);
image(fitH(eI*1000+c1))
pause(.01)
%% close dark
tic
cI = eI;
sRad = 10;
for r = 1:1
    pI = bwperim(cI);
    edgeIDs = find(eI);
    diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],3);
    diskE = (diskE>=diskE(sRad,1)) & (diskE<= diskE(sRad,round(sRad/2)));
    [dy dx] = find(diskE);
    yx = [dy dx]-sRad;
    for i = 1:length(edgeIDs)
        dVal = c1(edgeIDs(i));
        [cy cx] = ind2sub(size(I),edgeIDs(i));
        y = yx(:,1) + cy; x = yx(:,2) + cx;
        uyx = wall([y x],size(I));
        iyx = sub2ind(size(I),uyx(:,1),uyx(:,2));
        sVal = c1(iyx);  % get values of surround
        cI(iyx(sVal<(dVal*.8))) = 1;   %set dark surround to 1
        %image(fitH(cI + eI)),pause
    end
    %cI2 = imerode(cI,strel('disk',sRad * .5));
    image(fitH(cI + eI)),pause(1)
end
toc
%% connect edges
% tic
% 
% 
% sRad = 20;
%  diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],3);
%   threshE = (diskE>=diskE(sRad,1)) & (diskE<= diskE(sRad,2));
%    threshE = threshE + (diskE>=diskE(sRad,fix(sRad/2))) & (diskE<= diskE(sRad,fix(sRad/2)+1));
%   [dy dx] = find(threshE);
%  yx = [dy dx]-sRad; %neighbors
% [ey ex] = find(eI);
% links = zeros(length(ey)*length(dy),2);
% lPos = 1;
% for i = 1:length(ey)
%     y = yx(:,1) + ey(i); x = yx(:,2) + ex(i);
%     uyx = wall([y x],size(I));
%     iyx = sub2ind(size(I),uyx(:,1),uyx(:,2));
%     sVal = eI(iyx);
%     conInd = iyx(sVal>0);
%     numCon = length(conInd);
%     links(lPos:lPos+numCon-1,:) = [ones(numCon,1)*sub2ind(size(I),ey(i),ex(i)) conInd];
%     %links = cat(1,links,[ones(length(conInd),1)*sub2ind(size(I),ey(i),ex(i)) conInd]);
%     lPos = lPos + length(conInd);
% end
% links = links(links(:,1)>0,:);
% toc
%% Analyze links
% tic
% linkI = I * 0;
% pools = zeros(size(links,1),1);
% for l = 1:size(links,1)
%     [my mx] = ind2sub(size(I),links(l,:));
%     lSlope = (my(2)-my(1))/(mx(2)-mx(1));
%     if lSlope <= 1  & ~isinf(lSlope)
%         lx = mx(1) :abs(mx(2)-mx(1))/(mx(2)-mx(1)): mx(2);
%         ly = round((lx-mx(1)) * lSlope + my(1));
%     else
%         lSlope = (mx(2)-mx(1))/(my(2)-my(1));
%         ly = my(1) :abs(my(2)-my(1))/(my(2)-my(1)): my(2);
%         lx = round((ly-my(2)) * lSlope + mx(2));
%     end
%     linkID = sub2ind(size(I),ly,lx);
%     vCurv = c1(linkID);
%     plot(vCurv)
%     rim = min(c1(links(l,:)));
%     wells = vCurv(vCurv<rim);
%     pool = (rim - mean(wells))/rim;
%     %pause
%     pools(l) = pool;
%     if pool>.3
%         linkI(sub2ind(size(I),ly,lx)) =   linkI(sub2ind(size(I),ly,lx))+1;
%     end
%     
% end
% 
% image(double(linkI) * 10 + eI * 1000)
% toc 


%% Watershed
gKern2 = gaus3d([20 20 1],3);
c2 = fastCon(I,gKern2);

image(fitH(c2))
c2 = max(c2(:)) - c2;
image(fitH(c2))
minC2 = imhmin(c2,10);

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

pI(~cI) = 0;
%pI(linkI==0) = 0;
pI = imdilate(pI,strel('disk',4));
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

end
%%
