

for i = 1:length(sms)
    smsVgcs(i) = sms(i).sm.cid;
end

filterPostCids = [];

filterPostCids = [1107 1184 3103 3110 3111 3159 6033 6061 ]; %onlyl check triads with 4i rgcs


sameVG3 = 1;
vgcCids = [2 3 4 5 10 11 13 14 20 21 25 6720 6735 6738];
pre = tis.syn.pre;
post = tis.syn.post;
isBip = tis.syn.preClass == 7;
isRGC = tis.syn.postClass == 1;

%% find triad cellular motifs 
clear tri
pre2vgc = [];
post2vgc = [];
tri.vgcCid = [];
tri.bipCid = [];
tri.rgcCid = [];
for v = 1:length(vgcCids)

    pre2vgc =  pre((post== vgcCids(v)) & isBip);
    post2vgc = post((pre== vgcCids(v)) & isRGC);

    if ~isempty(filterPostCids)
        post2vgc = intersect(post2vgc,filterPostCids);
    end


    pre2vgc = unique(pre2vgc);
    post2vgc = unique(post2vgc);
    triMat = zeros(length(pre2vgc),length(post2vgc));

    for i = 1:length(pre2vgc);
        for j = 1:length(post2vgc);
            triMat(i,j) = sum((pre==pre2vgc(i)) & (post == post2vgc(j)));
        end
    end

    [ii jj]  = find(triMat>0);

    L1 = length(tri.vgcCid)+1;
    L2 = length(tri.vgcCid) + length(ii);
    if ~isempty(ii)
        tri.vgcCid(L1:L2,1) = ii*0+vgcCids(v);
        tri.bipCid(L1:L2,1) = pre2vgc(ii);
        tri.rgcCid(L1:L2,1) = post2vgc(jj);
    end
    tri
end

%% Get synapse ids for triads

%%IDs for tis
tNum = length(tri.bipCid);
for t = 1:tNum
    tri.id(t).b2rId = find((pre==tri.bipCid(t)) & ( post == tri.rgcCid(t)));
    tri.id(t).v2rId = find((pre==tri.vgcCid(t)) & ( post == tri.rgcCid(t)));
    tri.id(t).b2vId = find((post==tri.vgcCid(t)) & ( pre == tri.bipCid(t)));
end

%%IDs for sm
for t = 1:tNum
    curVid = find(smsVgcs==tri.vgcCid(t));
    syn = sms(curVid).sm.syn;

    tri.vid(t).v2rId = find((syn.pre==tri.vgcCid(t)) & ( syn.post == tri.rgcCid(t)));
    tri.vid(t).b2vId = find((syn.post==tri.vgcCid(t)) & ( syn.pre == tri.bipCid(t)));
end

%% get distances for each bip to rgc direct

pos = tis.syn.pos;
goodTri = ones(length(tri),1);
for t = 1:tNum
    c = 0;
    dEuc = [];
    dLin = [];
    isSameVG3 = [];

    curVid = find(smsVgcs==tri.vgcCid(t));
    d = sms(curVid).sm.skel2skel.linDist;
    vpos = sms(curVid).sm.syn.pos;

    for i = 1:length(tri.id(t).b2rId)
        b2rPos = pos(tri.id(t).b2rId(i),:);
        for v2 = 1:length(tri.vid(t).v2rId);
            v2rPos = vpos(    (tri.vid(t).v2rId(v2)),:);

            d1 = sqrt((b2rPos(1)-v2rPos(1)).^2 + (b2rPos(2)-v2rPos(2)).^2 + ...
                (b2rPos(3)-v2rPos(3)).^2);

            for b2 = 1:length(tri.vid(t).b2vId);

                b2vPos = vpos(tri.vid(t).b2vId(b2),:);


                d2 = sqrt((b2rPos(1)-b2vPos(1)).^2 + (b2rPos(2)-b2vPos(2)).^2 + ...
                    (b2rPos(3)-b2vPos(3)).^2);

                d3 = sqrt((v2rPos(1)-b2vPos(1)).^2 + (v2rPos(2)-b2vPos(2)).^2 + ...
                    (v2rPos(3)-b2vPos(3)).^2);

                c = c+1;
                dEuc(c,:) = [d1 d2 d3];

                dLin(c) = d(tri.vid(t).b2vId(b2),tri.vid(t).v2rId(v2));

            end
        end
    end

    if isempty(dEuc)
        dEuc = [inf inf inf];
        goodTri(t) = 0;
        dLin = inf;
    end


    tri.dEuc{t} = dEuc;
    tri.dLin{t} = dLin;
    sumEuc = sum(dEuc,2);
    sumEuc(~isSameVG3) = inf;
    minEuc = min(sumEuc);
    minEucId = find(sumEuc==minEuc,1);
    minLin = min(dLin);
    minLinId = find(dLin == minLin,1);

    tri.dEucSmallest(t,:) = dEuc(minLinId,:);
    tri.minEuc(t) = minEuc;
    tri.minEucId(t) = minEucId;

    tri.dLinSmallest(t,1) = dLin(minLinId);
    tri.minLin(t) = minLin;
    tri.minLinId(t) = minLinId;


    otherSynNum =  length(tri.vid(t).v2rId);
    hasSharedRibbon = sum(tri.dEuc{t}(:,2)==0) / otherSynNum;
    tri.b2vNum(t) = length(tri.vid(t).b2vId);
    tri.hasSharedRibbon(t) = hasSharedRibbon;

end


%% display
dEuc = tri.dEucSmallest;
dLin = tri.dLinSmallest;
isOnSame = tri.minEuc<inf;

%maxEuc = sum(dEuc(:,[1 3]),2);
%maxEuc = dEuc(:,3);
useDist = dLin;


fromDyad = find((tri.dEucSmallest(:,2) ==0) & isOnSame');
notDyad = find((tri.dEucSmallest(:,2) >0) & isOnSame');
closeNotDyad = find((tri.dEucSmallest(:,2) >0) & isOnSame' & (dEuc(:,2) < 5));


hRange = [0:1:150];
histMinEuc = histcounts(useDist(isOnSame),hRange);
histMinEucDyad = histcounts(useDist(fromDyad ),hRange);
histMinEucNotDyad = histcounts(useDist(notDyad ),hRange);
histMinEucCloseNotDyad = histcounts(useDist(closeNotDyad ),hRange);


maxH = max(histMinEuc);
subplot(4,1,1)
bar(hRange(1:end-1),histMinEuc,'k','edgecolor','none')
ylim([0 maxH])
subplot(4,1,2)
bar(hRange(1:end-1),histMinEucDyad,'k','edgecolor','none')
ylim([0 maxH])
subplot(4,1,3)
bar(hRange(1:end-1),histMinEucNotDyad,'k','edgecolor','none')
ylim([0 maxH])
subplot(4,1,4)
bar(hRange(1:end-1),histMinEucCloseNotDyad,'k','edgecolor','none')
ylim([0 maxH])


median(useDist(fromDyad))
median(useDist(notDyad))
median(useDist(closeNotDyad))


mean(useDist()<1)

mean(useDist(fromDyad)<1)
mean(useDist(notDyad)<1)
mean(useDist(closeNotDyad)<1)

ranksum((useDist(fromDyad)),(useDist(notDyad)))

%sum(histMinEuc(1:3))/sum(histMinEuc(:))


%% influence

lc = 16;
useDs = useDist(isOnSame)
d = median(useDs)
W = exp(-d/lc) % Apply length constant
mean(W)
median(W)

mean(useDist(isOnSame)<=2)


%% dyad

closeNonDyad = useDist(intersect( notDyad,find(dEuc(:,2) < 2)))

median(closeNonDyad)


sum(tri.hasSharedRibbon(isOnSame))/sum(tri.b2vNum(isOnSame))


mean(tri.hasSharedRibbon(isOnSame)>0)
sum(tri.hasSharedRibbon(isOnSame)>0)



