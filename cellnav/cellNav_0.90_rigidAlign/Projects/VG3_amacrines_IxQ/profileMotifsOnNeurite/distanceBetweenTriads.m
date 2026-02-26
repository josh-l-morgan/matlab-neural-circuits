

sameVG3 = 1;
vgcCids = [2 3 4 5 10 11 13 14 20 21 25 6720 6735 6738];
pre = tis.syn.pre;
post = tis.syn.post;
isBip = tis.syn.preClass == 7;
isRGC = tis.syn.postClass == 1;


pre2vgc = [];
post2vgc = [];
for v = 1:length(vgcCids)
    pre2vgc = [pre2vgc; pre((post== vgcCids(v)) & isBip)];
    post2vgc = [post2vgc; post((pre== vgcCids(v)) & isRGC)];
end
pre2vgc = unique(pre2vgc);
post2vgc = unique(post2vgc);
post2vgc = [2002 3051 3119]; %!!!!!!!!!!!

triMat = zeros(length(pre2vgc),length(post2vgc));
for i = 1:length(pre2vgc);
    for j = 1:length(post2vgc);
        triMat(i,j) = sum((pre==pre2vgc(i)) & (post == post2vgc(j)));
    end
end

[ii jj] = find(triMat>0);

clear tri
tri.bipCid = pre2vgc(ii);
tri.rgcCid = post2vgc(jj);

tNum = length(tri.bipCid);
for t = 1:tNum
    tri.id(t).b2rId = find((pre==tri.bipCid(t)) & ( post == tri.rgcCid(t)));

    v2rId = [];
    v2rCid = [];
    for v = 1:length(vgcCids)
        ids = find((pre==vgcCids(v)) & ( post == tri.rgcCid(t)));
        v2rId = [v2rId; ids];
        v2rCid = [v2rCid; ids * 0+vgcCids(v)];
    end
    tri.id(t).v2rId = v2rId;
    tri.id(t).v2rCid = v2rCid;

    b2vId = [];
    b2vCid = [];
    for v = 1:length(vgcCids)
        ids = find((post==vgcCids(v)) & ( pre == tri.bipCid(t)));
        b2vId = [b2vId; ids];
        b2vCid = [b2vCid; ids * 0+vgcCids(v)];
    end
    tri.id(t).b2vId = b2vId;
    tri.id(t).b2vCid = b2vCid;

end

%% get distances for each bip to rgc direct

pos = tis.syn.pos;
for t = 1:tNum
    c = 0;
    dEuc = []; 
    isSameVG3 = [];
    for i = 1:length(tri.id(t).b2rId)
        b2rPos = pos(tri.id(t).b2rId(i),:);
        for v2 = 1:length(tri.id(t).v2rId);
            v2rPos = pos(    (tri.id(t).v2rId(v2)),:);

            d1 = sqrt((b2rPos(1)-v2rPos(1)).^2 + (b2rPos(2)-v2rPos(2)).^2 + ...
                (b2rPos(3)-v2rPos(3)).^2);

            for b2 = 1:length(tri.id(t).b2vId);
                b2vPos = pos(tri.id(t).b2vId(b2),:);

                d2 = sqrt((b2rPos(1)-b2vPos(1)).^2 + (b2rPos(2)-b2vPos(2)).^2 + ...
                    (b2rPos(3)-b2vPos(3)).^2);

                d3 = sqrt((v2rPos(1)-b2vPos(1)).^2 + (v2rPos(2)-b2vPos(2)).^2 + ...
                    (v2rPos(3)-b2vPos(3)).^2);

                c = c+1;
                isSameVG3(c) = tri.id(t).v2rCid(v2) == tri.id(t).b2vCid(b2);
                dEuc(c,:) = [d1 d2 d3];

            end
        end
    end
    tri.dEuc{t} = dEuc;
    tri.isSameVG3{t} = isSameVG3;
    sumEuc = sum(dEuc,2);
    sumEuc(~isSameVG3) = inf;
    minEuc = min(sumEuc);
    minEucId = find(sumEuc==minEuc,1);

    tri.dEucSmallest(t,:) = dEuc(minEucId,:);
    tri.minEuc(t) = minEuc;
    tri.minEucId(t) = minEucId;
end

%% display
dEuc = tri.dEucSmallest;
isOnSame = tri.minEuc<inf;

%maxEuc = sum(dEuc(:,[1 3]),2);
%maxEuc = dEuc(:,3);
maxEuc = dEuc(:,3);


fromDyad = find((tri.dEucSmallest(:,2) ==0) & isOnSame');
notDyad = find((tri.dEucSmallest(:,2) >0) & isOnSame');

hRange = [0:1:50];
histMinEuc = histcounts(maxEuc(isOnSame),hRange);
histMinEucDyad = histcounts(maxEuc(fromDyad ),hRange);
histMinEucNotDyad = histcounts(maxEuc(notDyad ),hRange);


maxH = max(histMinEuc); 
subplot(3,1,1)
bar(hRange(1:end-1),histMinEuc,'k','edgecolor','none')
ylim([0 maxH])

subplot(3,1,2)
bar(hRange(1:end-1),histMinEucDyad,'k','edgecolor','none')
ylim([0 maxH])
subplot(3,1,3)
bar(hRange(1:end-1),histMinEucNotDyad,'k','edgecolor','none')
ylim([0 maxH])

mean(maxEuc(fromDyad))
mean(maxEuc(notDyad))

sum(histMinEuc(1:3))/sum(histMinEuc(:))








