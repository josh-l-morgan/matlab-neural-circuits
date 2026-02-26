function[cellCounts] = getSynapseNumberWithLocalExclusion(sms,minDist,checkTypeID)

%%Generate number of rgcsubtypes or bipolar cell subtypes innervated by 
%%VG3s in a way that reduces the influence of local clustering

global tis

if ~exist('minDist','var')
  minDist = 0;
end  
if ~exist('checkTypeID','var')
    checkTypeID = 1;
end

%count synapses between pairs of VG3s and RGCs with minDist exclusion length

for v = 1:length(sms)
    sm = sms(v).sm;
    sm.syn.pre;
    uPost = unique(sm.syn.pre);
    uPost = setdiff(uPost,[sm.cid 0]);
    cellsHit{v} = uPost;
    D = sm.syn2Skel.syn2SynDist;
    goodHits{v} = zeros(length(uPost),1);
    for p = 1:length(uPost)
        hits = find(sm.syn.pre==(uPost(p)));
        left = hits;
        use = [];
        while ~isempty(left)
            grab = ceil(rand*length(left));
            new = left(grab);
            use = [use new];
            d = D(new,left);
            d(grab) = -100;
            ok = d>minDist;
            left = left(ok);
        end
            goodHits{v}(p) = length(use);
    end
end

%% Fuse VG3s
posts = [];
hitNums = [];
for v = 1:length(cellsHit)
    if ~isempty(cellsHit{v});
        posts = cat(1,posts,cellsHit{v});
        hitNums = cat(1,hitNums,goodHits{v});
    end
end
uPosts = unique(posts);
hits = uPosts * 0;
for p = 1:length(uPosts);
    hits(p) = sum(hitNums(posts==uPosts(p)));
end

sum(hits)

%% Get types

types = uPosts * 0;
subTypes = uPosts * 0;

for i = 1:length(uPosts)
    targ = find(tis.cids == uPosts(i));
    if isempty(targ)
        disp(sprintf('help, I cant find cid %d',uPosts(i)));
    else
        types(i) = tis.cells.type.typeID(targ);
        subTypes(i) = tis.cells.type.subTypeID(targ);
    end
end

%% trim to RGCs

uPosts = uPosts(types == checkTypeID);
hits = hits(types == checkTypeID);
subTypes = subTypes(types == checkTypeID);

hitSub = unique(subTypes);
hitSub = setdiff(hitSub,0);

subHit = hitSub * 0;
for i = 1:length(hitSub)
    subHit(i) = sum(hits(subTypes==hitSub(i)));
end

%% Names
clear cellCounts
for s = 1:length(hitSub)
    cellCounts{s,1} = tis.cells.type.subTypeNames{checkTypeID}{hitSub(s)};
    cellCounts{s,2} = subHit(s);
end


