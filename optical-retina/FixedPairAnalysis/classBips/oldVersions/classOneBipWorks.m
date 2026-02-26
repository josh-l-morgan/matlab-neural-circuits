%% Read data

clear all
tic
colormap gray(256)
%TPN = GetMyDir
TPN = 'C:\Users\joshm\Documents\MyWork\myData\bipolarData\bipDump\2\'
'reading data'
TPNd = dir(TPN); TPNd = TPNd(3:length(TPNd));
TiffNames={};
for i = 1: size(TPNd,1)
    siz=length(TPNd(i).name);
    if TPNd(i).name(siz-2:siz)== 'tif'
        if TPNd(i).name(siz-8:siz-7)~='-R'
            TiffNames(length(TiffNames)+1,1)={TPNd(i).name};
        end
    end
end

for i = length(TiffNames):-1:1
    I(:,:,i) = imread([TPN TiffNames{i}]);
end

[rys rxs rzs] = size(I);
[yi xi zi] = ind2sub(size(I),find(I>0));
subI = [min(yi) max(yi);min(xi) max(xi);min(zi) max(zi)];
I = double(I(subI(1,1):subI(1,2),subI(2,1):subI(2,2),...
    subI(3,1):subI(3,2)));
image(sum(I,3))

D = I>0;
[ys xs zs] = size(D);
[Dys Dxs Dzs] = ind2sub(size(D),find(D>0));
D= D(min(Dys):max(Dys),min(Dxs):max(Dxs),min(Dzs):max(Dzs));
[ys xs zs] = size(D);


%% Find shaft (currently not removing shaft)
'find shaft'
sumD = sum(D,3);
subplot(2,2,1)
image(sumD * 255/max(sumD(sumD>0)))
% [maxYs maxXs] = find(sumD == max(sumD(:)));
%
% cdRad = 30;
% shaft = [mean(maxYs) mean(maxXs)];
% cutDisk = fspecial('disk',cdRad);
% [nSy nSx] = ind2sub(size(cutDisk),find(cutDisk>0));
% nearShaft = sub2ind([ys xs],round(nSy + shaft(1)-cdRad), round(nSx + shaft(2)-cdRad));
% for z = 1:size(D,3)
%    Dplane = D(:,:,z);
%    Dplane(nearShaft) = 0;
%    Dc(:,:,z) = Dplane;
% end
% sumDc = sum(Dc,3)
% subplot(2,2,2)
% image(sumDc * 255/max(sumDc(sumDc>0)))
% % shift orthogonally
% 'shift to princomp'
%
% [Yp Xp Zp] = ind2sub([ys xs zs], find(D>0));
% [coef, score, latent] = princomp([Yp Xp Zp]);
%
% score(:,1) = score(:,1) - min(score(:,1)) + 1;
% score(:,2) = score(:,2) - min(score(:,2)) + 1;
% score(:,3) = score(:,3) - min(score(:,3)) + 1;
%
% nY = fix(max(score(:,1)))+1;
% nX = fix(max(score(:,2)))+1;
% nZ = fix(max(score(:,3)))+1;
%
% P = zeros(nY, nX, nZ);
% score = round(score);
% newInd = sub2ind([nY nX nZ],score(:,1),score(:,2),score(:,3));
% P(newInd) = 1;
%
% for i = 1: size(score,1)
%    P(score(i,1),score(i,2),score(i,3)) = 1;
% end
% xD = squeeze(sum(D,1));
% yD = squeeze(sum(D,2));
% xP = squeeze(sum(P,1));
% yP = squeeze(sum(P,2));
% subplot(2,2,1),image(xD * 100/median(xD(xD>0)));
% subplot(2,2,2),image(yD * 100/median(yD(yD>0)));
% subplot(2,2,3),image(xP * 100/median(xP(xP>0)));
% subplot(2,2,4),image(yP * 100/median(yP(yP>0)));
% pause(.01)

%% Depth profile
dP = squeeze(sum(sum(P,1),2));
bar(dP)
whm = find(dP>=(max(dP)/2));
tempDP = dP;
tempDP(whm) = 0;
hold on
bar(tempDP,'b')
hold off

topA = min(whm):size(P,3);pause(.1)

%% Polygone
M = P;
sumTop = sum(M(:,:,topA),3);
image(sumTop * 100/median(sumTop(sumTop>0)));


[y x] = find(sumTop>0);
numArea = length(x);
k = convhull(y,x);
scatter(y,x,'.')
hold on
plot(y(k),x(k),'r-')
hold off
polyArea = polyarea(y(k),x(k)) ;

%% Load Skeleton

load([TPN 'data/AllSeg.mat'])
scaleSeg = AllSeg;
scaleSeg(:,1,:) = scaleSeg(:,1,:)/0.06 + min(y);
scaleSeg(:,2,:) = scaleSeg(:,2,:)/0.06 + min(x);
scaleSeg(:,3,:) = scaleSeg(:,3,:)/0.2;
hold on
for i = 1:size(scaleSeg,1)
    plot([scaleSeg(i,1,1) scaleSeg(i,1,2)],[scaleSeg(i,2,1) scaleSeg(i,2,2)],'g')
end
hold off

xyLengths = sqrt((scaleSeg(:,1,1)-scaleSeg(:,1,2)).^2+(scaleSeg(:,2,1)-scaleSeg(:,2,2)).^2);
totLength = sum(xyLengths);

pause(.01)



%% find buotns by watershed
gKern = gaus3d([7 7 3],3,[1 1 3]);
gKern = gKern/sum(gKern(:));
cI = fastCon(I,gKern);
cI = cI>.5;

pI = bwperim(~I);
pI(:,:,1) = 0; pI(:,:,zs) = 0;
pI(:,1,:) = 0; pI(:,xs,:) = 0;
pI(1,:,:) = 0; pI(ys,:,:) = 0;
dI =bwdistsc(~I & ~pI,[.06 .06 .2]);

%% fft convolve
%%Make gaussian kernal
gKern = gaus3d([15 15 5],2,[1 1 3]);
cI = fastCon(dI,gKern);


%% watershed
dI2 = max(cI(:))-cI;
%dI2 = dI2 * 100/median(dI2(dI2>0));
dI2(~I) = -Inf;%min(dI2(:));
%gKern = gaus3d([15 15 5],2,[1 1 3]);
%cI = fastCon(dI2,gKern);
%cI = imhmin(cI,.1); %
%dI2(dI2<-10) = -10;
dI2(pI) = max(dI2(:));
wI = watershed(dI2,26);
numWI = max(wI(:));
subplot(1,1,1)
% for i = 1:size(wI,3)
%     imshow(label2rgb(wI(:,:,i),'jet','w')),pause
% end


%% Analyze subregions
rProps = regionprops(wI,cI,'WeightedCentroid',...
    'PixelIdxList','PixelValues','MaxIntensity','Area');
rP = wI * 0;
for i = 1:length(rProps)
    num(i) = rProps(i).Area * prod([.06 .06 .2]);
    maxD(i) = rProps(i).MaxIntensity;
    eqVol(i) = 4/3 * pi * maxD(i)^3;
    roundPart(i) = eqVol(i)/num(i);
    rP(rProps(i).PixelIdxList) = roundPart(i)* 255;
    cents2(i,:) = rProps(i).WeightedCentroid;
end
cents(:,1) = cents2(:,2);
cents(:,2) = cents2(:,1);
cents(:,3) = cents2(:,3);
% for i = 1:size(rP,3)
%    image(rP(:,:,i)),pause
% end

%% Connect Watershed
clear nears
gapI = I & ~wI;

near = [1 0 0;-1 0 0; 0 1 0 ; 0 -1 0; 0 0 1; 0 0 -1];
[ny nx nz] = ind2sub([3 3 3],find(ones(3,3,3)));
near = [ny-2 nx-2 nz-2];
wI2 = wI;
wI2(I & (wI2 ==1)) = 0;


for r = 1:16
    gaps = find(I & ~wI2);
    showI2 = wI2;
    c = 0;
    colPairs = {};
    for i = 1:length(gaps)
        [y x z] = ind2sub(size(wI2),gaps(i));
        nears(:,1) = near(:,1) + y;
        nears(:,2) = near(:,2) + x;
        nears(:,3) = near(:,3) + z;
        nears = wall(nears,size(wI2));
        nearby = wI2(sub2ind(size(wI2),nears(:,1),nears(:,2),nears(:,3)));
        things = unique(nearby(nearby>1));
        if length(things) == 1
            wI2(y,x,z) = things;
            showI2(y,x,z) = things*20;
        end
        if length(things) >1
            c = c+1;
            colPairs{c,1} = things;
            pairPos(c) = gaps(i);
        end
    end % search all gaps
end % repeat spreading
for i = 1:size(gapI,3)
    image(showI2(:,:,i)),pause(.1)
end
gapI = I & ~wI2;
toc
for i = 1:size(gapI,3)
    subplot(2,1,1)
    image(wI(:,:,i)*30)
    subplot(2,1,2)
    image((wI2(:,:,i)>1)* 50 )%+ gapI(:,:,i)*1000),
    pause(.01)
end

%% findPairs
%%Grab pairs out of potential triples, quads
pairList = [];
listPos = [];
rdUp = 10^length(num2str(numWI));
for i = 1:length(colPairs)
    pairing = sort(colPairs{i});
    allComs = combntns(pairing,2);
    for p = 1: size(allComs,1)
        pairID = allComs(p,1) * rdUp + allComs(p,2);
        pairList(length(pairList)+1) = pairID;
        listPos(length(listPos)+1) = pairPos(i);
    end
end

%%All unique Pairs
uPairs = unique(pairList);

%%Find positions of joints
[lGap nGap] = bwlabeln(gapI,8);
showI = gapI * 0;
for i = 1:length(uPairs)
%     showI(listPos(pairList == uPairs(i)))=100;
%     gapLab = unique(lGap(listPos(pairList == uPairs(i))));
%     gID = [];
%     for l = 1:length(gapLab)
%         gID = cat(1,gID,find(lGap == gapLab(l)));
%     end
%     [jPy jPx jPz] = ind2sub(size(wI),gID);
    [jPy jPx jPz] = ind2sub(size(wI),(listPos(pairList == uPairs(i))));
    jointPos(i,:) = [mean(jPy) mean(jPx) mean(jPz)];
    jointMem(i,:) = [fix(uPairs(i)/rdUp) mod(uPairs(i),rdUp)];
    jointVol(i) = length(jPy);
end



%% Make Skel
%more conservative solution to internal skeleton

tips = [];tipPos = [];
for i = 2: numWI
    [participates fb] = find(jointMem == i);
    nodes = jointPos(participates,:);
    nodes = cat(1,cents(i,:),nodes);
    jIDs = cat(2,i,uPairs(participates));
    if size(nodes,1) == 1  % no joints
        [y x z] = ind2sub(size(wI),find(wI == i));
        [COEFF,Score] = princomp([y x z]);
        minPos = find(Score == min(Score(:,1)),1);
        maxPos = find(Score == max(Score(:,1)),1);
        nodes = [cents(i,:);y(minPos) x(minPos) z(minPos);y(maxPos) x(maxPos) z(maxPos)];
        nPair = [1 2; 1 3];
        jIDs = [jIDs -i -i-rdUp];
        tipPos(size(tipPos,1)+1:size(tipPos,1)+2,:) = nodes(2:3,:);
        tips(length(tips)+1:length(tips)+2) = jIDs(2:3);
    elseif size(nodes,1) == 2  % 1 joint (end point)
        [y x z] = ind2sub(size(wI),find(wI == i));
        dists = dist([y x z],nodes(2,:));
        furthest = find(dists == max(dists),1);
        nodes(3,:) = [y(furthest) x(furthest) z(furthest)];
        nPair = [1 2;1 3];
        jIDs = [jIDs -i];
          tipPos(size(tipPos,1)+1,:) = nodes(end,:);
        tips(length(tips)+1) = jIDs(end);
    else % muliple joints
        nPair = zeros(size(nodes,1),2);
        for n = 1:size(nodes,1)  %find first pairings (nearest neighbors)
            dists = dist(nodes,nodes(n,:));
            nPair(n,:) = sort([n find(dists == min(dists(dists>0)))]);
        end
        nPair = unique(nPair,'rows');
        numPair = size(nPair,1);
        %%Find second pairings
        passed1 = 0;
        while passed1 < numPair
            passed1 = 0; %Count how many pairs pass being connected
            for n = 1:size(nPair,1)
                oNodes = paired(nPair,n);
                if isempty(find(oNodes==1))
                    otherNodes = setdiff(nPair(:),oNodes);
                    minDist = zeros(length(oNodes),1);
                    minPair = minDist;
                    for o = 1:length(oNodes)
                        dists = dist(nodes(otherNodes,:),nodes(oNodes(o),:));
                        minDist(o) = min(dists);
                        minPair(o) = otherNodes(find(dists == min(dists),1));
                    end
                    minnest = find(minDist==min(minDist));
                    newPair = [oNodes(minnest) minPair(minnest)];
                    nPair = cat(1,nPair,newPair);
                else
                    passed1 = passed1+1;
                end  %end if connected to 1
            end % end run all pairs
        end % end repeat until all connected to center (1)
    
    end % if multiple joints

skel(i).nodes = nodes;
skel(i).ids = jIDs(nPair);
skel(i).pairs = nPair;
end

%% Break Loops
%%Using -> uPairs,jointVol,links
links = []; allPos = []; allIds = [];
for i = 1: length(skel)
    links = cat(1,links,skel(i).ids);
end
allIds = cat(2,1:size(cents,1), uPairs, tips);
allPos = cat(1,cents,jointPos, tipPos);


links = breakLoop(links,uPairs,jointVol);

%% PlotSkel
%%Make all Seg
subplot(1,2,1)
hold off
plot(1,1)
whitebg('k')
P = D;
'render cell'

p = patch(isosurface(I>0,.1));
%isonormals(xs,ys,zs,D,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1/.06 1/.06 1/.2])
alpha(.5)
view(3); axis tight
camlight
lighting gouraud
pause(.01)
hold on

lCol = 'rgbcmwy';
lCol = 'bg';
[paths, showP] = getPaths(links);
c = 0;
% [y x z] = ind2sub(size(I),find(I));
% sub = mod(1:length(y),100) == 1;
% scatter3(y(sub), x(sub) , z(sub)/3,'.')
% whitebg
% hold on
for p = 1:length(paths)
    path = paths{p};
    if length(path)>1
        for n = 2:length(path)
            c = c+1;
            pos(1,:) = allPos(allIds ==path(n-1),:);
            pos(2,:) = allPos(allIds == path(n),:);
            useCol = lCol(mod(c,length(lCol))+1);
            plot3(pos(:,2),pos(:,1),pos(:,3),useCol,'LineWidth',3)
            hold on
        end
    end
end
hold off
daspect([1/.06 1/.06 1/.2])
pause(1)

%%
subplot(1,2,2)
AllSeg = [];

for i = 1:size(links,1)
    AllSeg(size(AllSeg,1)+1,:,:) = cat(3,allPos(allIds ==links(i,1),:), allPos(allIds ==links(i,2),:));
    
end

sumwI = max(wI,[],3);
image(sumwI),
colormap colorcube(255)
[ploty plotx] = ind2sub(size(sumwI),find(sumwI>0));
%scatter(plotx,ploty,'.','b')
hold on
for i = 1:size(AllSeg,1)
    plot([AllSeg(i,2,1) AllSeg(i,2,2)],[AllSeg(i,1,1) AllSeg(i,1,2)],'r',...
        'LineWidth',2)
 
end
hold off


for i = 1:3, AllSeg(:,i,:) = AllSeg(:,i,:) + subI(i,1); end
AllSeg(:,3,:) = (AllSeg(:,3,:)-rzs) ;
AllSeg(:,1,:) = AllSeg(:,1,:) - rys;
AllSegb = AllSeg;
AllSegb(:,1,:) = AllSeg(:,2,:);
AllSegb(:,2,:) = AllSeg(:,1,:);
AllSeg = AllSegb;
save([TPN 'data\AllSegW.mat'],'AllSeg')

%%  Collect data
%{
maskArea
polyArea
maskVol
%}
%% Render mask
% for i = 1: max(wI(:))
% 'render cell'
% subplot(1,1,1)
% p = patch(isosurface(wI == i,.1));
% %isonormals(xs,ys,zs,D,p)
% set(p,'FaceColor','red','EdgeColor','none');
% daspect([1/.06 1/.06 1/.2])
% view(3); axis tight
% camlight
% lighting gouraud
% pause
% end


toc