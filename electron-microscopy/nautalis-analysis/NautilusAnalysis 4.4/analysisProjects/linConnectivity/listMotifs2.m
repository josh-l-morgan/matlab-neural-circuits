
load('MPN.mat');
load([MPN 'obI.mat']);
targCell = 125;
maxDist = 5;

disp('combining spreadsheet and segmentation synapses')
sm = addDatToSynMat(obI)
disp('finding topological distance between synapses on target cell within max range')
sm  = getTopoEucDistBetweenSyn(sm);
sm = getTopoEucDistBetweenSkelAndSyn(sm)

%% analyze Motifs


isTarg =(sm.pre == targCell) | (sm.post == targCell);
isRGC = sm.preClass == 1;
isTC = sm.postClass == 2;
isPreLIN = sm.preClass == 3;
isPostLIN = sm.postClass == 3;
isPreUNK = sm.preClass == 4;


%search all RGC input, check for non RGC TC
%search all LIN input, check for isolated LIN output
%search for all UNK input






%% find checkIn relative motifs


[checkRGC] = getMotifsAroundSyn(sm,find(isRGC),maxDist,targCell)
[checkPreLIN] = getMotifsAroundSyn(sm,find(isPreLIN),maxDist,targCell)
[checkPreUNK] = getMotifsAroundSyn(sm,find(isPreUNK),maxDist,targCell)

checkS = find(isRGC & isTarg);

check.checkS = checkS;
check.targCell = targCell;
check.maxDist = maxDist;

for i = 1: length(checkS)
    
    checkPos = sm.pos(checkS(i),:);
    preID = sm.pre(checkS(i));
    postID = sm.post(checkS(i));
    preIDClass = sm.preClass(checkS(i));
    postIDClass = sm.postClass(checkS(i));
    
    check.pos(i,:) = checkPos;
    
    close = find(sm.dist(checkS(i),:)<=maxDist);
    synNum(i) = length(close);
    
    preClass = sm.preClass(close);
    postClass = sm.postClass(close);
    pre = sm.pre(close);
    post = sm.post(close);
    
    %identify all pre and post
    pre2pre = pre(post==preID);
    pre2preClass = preClass(post==preID);
    post2pre = post(pre == preID);
    post2preClass = postClass(pre == preID);
    pre2post = pre(post == postID);
    pre2postClass = preClass(post == postID);
    post2post = post(pre == postID);
    post2postClass = postClass(pre == postID);
    
    %% identify triads
    [recipPre ia ib] = intersect(pre2pre,post2pre);
    [recipPost ia ib] = intersect(pre2post,post2post);
    [io ioa iob] = intersect(pre2pre,post2post);
    ioClass = pre2preClass(ioa);
    [ii iia iib] = intersect(pre2pre,pre2post);
    iiClass = pre2preClass(iia);
    [oo ooa oob] = intersect(post2pre,post2post);
    ooClass = post2preClass(ooa);
    [oi oia oib] = intersect(post2pre,pre2post);
    oiClass = post2preClass(oia);
    
    %%check for these motifs
    %     m.preLab = length(checkS);
    %     m.preTc = 0;
    %     m.preTcDoub = 0;
    %     m.preLin = 0;
    %     m.preLinTc = 0;
    %     m.preLinLinIn = 0;
    %     m.preLinLinOut = 0;
    %     m.preLinLinDoub = 0;
    %     m.preOneOut = 0;
    %     m.preLinTcDi = 0;
    %     m.preLinLinDi = 0; %chain of two
    %     m.recip = 0;
    %     m.autapse = 0;
    m.preClass = preIDClass;
    m.postClass = postIDClass;
    
    m.preTcTri = 0;
    m.preLinTri =0;
    m.preOnly = 0;
    m.preTcDi = 0;
    m.preLinDi = 0;
    m.preLinInTri = 0;
    m.postOnly =0;

    %% check pre motif by mapping against post2post

   
    if isempty(post2post) %local dead end
        m.preOnly = 1; % RGC has no other inputs
    elseif preID % not zero
        zeroClass = post2postClass(post2post == 0);
        m.preTcDi = m.preTcDi + sum(zeroClass==2);
        m.preLinDi = m.preLinDi + sum(zeroClass==3);
        checkPost = find(post2post>0);
        for o = 1:length(labPost) % check against all targets of post
            checkPostClass = post2postClass(checkPost(o));
            checkPostID = post2post(checkPost(o));
            if sum(oo == checkPostID); %if preRGC and post LIN converge on target
                if checkPostClass == 2
                    m.preTcTri = m.preTcTri + 1;
                elseif checkPostClass == 3
                    m.preLinTri = m.preLinTri + 1;
                end
            else % no triad
                if checkPostClass == 2
                    m.preTcDi =  m.preTcDi  + 1;
                elseif checkPostClass == 3
                    m.preLinDi = m.preLinDi  + 1;
                end
            end
        end
    end
    
    %%Find cell that both receives from pre and innervates post
    otherMid = setdiff(oi,0);
    m.preLinInTri = sum(otherMid);
    
    %% determine if post is isolated (otherwise motif will be described from pre perspective)
    if isempty(pre2pre)
        m.postOnly = 1;
    end
    
%         zeroClass = pre2preClass(pre2pre == 0);
%         m.postRgcDi = m.postRgcDi + sum(zeroClass==1);
%         m.postLinDi = m.postLinDi + sum(zeroClass==3);
%         checkPost = find(pre2pre>0);
%         for o = 1:length(labPost) % check against all targets of post
%             checkPostClass = post2postClass(checkPost(o));
%             checkPostID = post2post(checkPost(o));
%             if sum(oo == checkPostID); %if preRGC and post LIN converge on target
%                 if checkPostClass == 2
%                     m.preTcTri = m.preTcTri + 1;
%                 elseif checkPostClass == 3
%                     m.preLinTri = m.preLinTri + 1;
%                 end
%             else % no triad
%                 if checkPostClass == 2
%                     m.preTcDi =  m.preTcDi  + 1;
%                 elseif checkPostClass == 3
%                     m.preLinDi = m.preLinDi  + 1;
%                 end
%             end
%         end
        
        
    end
    
    
end



end














        
        
        %% Check for LIN-LIN-LIN motifs
        if sum(linIn) & sum(linOut)
            preLin = pre(linIn);
            postLin = post(linOut);
            recip = intersect(preLin,postLin)
            m.recipNum = length(recip);
            m.autapse = sum(recip==targCell);
            
            preLin = setdiff(preLin,[0 targCell]);
            postLin = setdiff(postLin,[0 targCell]);
            prePre = 0; prePost = 0; postPre = 0; postPost = 0;
            for s1 = 1:length(preLin)
                prePre = sum((pre==preLin(s1)) & (post==preLin(s1)));
                for s2 = 1:length(postLin)
                    prePost = sum((pre==preLin(s1)) & (post==postLin(s2)));
                    postPre = sum((pre==postLin(s2)) & (post==preLin(s1)));
                end
            end
            for s2 = 1:length(postLin)
                postPost = sum((pre==postLin(s2)) & (post==postLin(s2)));
            end
            m.linTriType = [prePost postPre prePre postPost];
            m.linTri = sum(m.linTriType);
        end
        
        
    end
    
    
    
    skelM(i).m = m;
    
    
    
end












return

colormap(jet(100))
scatter3(sm.skelNodes(:,1),sm.skelNodes(:,2),sm.skelNodes(:,3),50,synNum*1000,'o','filled')



%% find RGC motifs
rgcIN = find(isTarg & isRGC);
maxDist = 5;
for i = 1:length(rgcIN)
    
    close = find(sm.dist(rgcIN(i),:)<=maxDist);
    
    
    
end












%% Get data
mot = getMotifs(obI);
synMat = getSynMat(obI);
synStruct = getSynMat(obI);
synPos = getSynPos(1);
targ = 125;

%% Match synapses
%%Find all synapses in synPos that are not in synMat.

onTarg = ((synMat.pre == targ) | (synMat.post == targ))  &...
    (sum(synMat.synPos,2)>0);
oPos = synMat.synPos(onTarg,:);
oPreClass = synMat.preClass(onTarg);
oPostClass = synMat.postClass(onTarg);
oPre = synMat.pre(onTarg);
oPost = synMat.post(onTarg);
oUse = ones(size(oPre));

res = obI.em.res;
dSamp =  (res .* [4 4 1])./1000;

%%Pool labels
pos = synPos.postRGCPos;
pos = pos(:,[2 1 3]);
pos(:,1) = pos(:,1)*dSamp(1);
pos(:,2) = pos(:,2)*dSamp(2);
pos(:,3) = pos(:,3)*dSamp(3);
L = size(pos,1);
nPos = pos;
preClass = zeros(L,1) + 1;
postClass = zeros(L,1) + 3;
pre = zeros(L,1);
post = zeros(L,1) + targ;


pos = synPos.postUnkPos;
pos = pos(:,[2 1 3]);
pos(:,1) = pos(:,1)*dSamp(1);
pos(:,2) = pos(:,2)*dSamp(2);
pos(:,3) = pos(:,3)*dSamp(3);
L = size(pos,1);
nPos = [nPos; pos];
preClass = [preClass; zeros(L,1) + 4];
postClass = [postClass; zeros(L,1) + 3];
pre = [pre; zeros(L,1)];
post = [post; zeros(L,1) + targ];


pos = synPos.postLinPos;
pos = pos(:,[2 1 3]);
pos(:,1) = pos(:,1)*dSamp(1);
pos(:,2) = pos(:,2)*dSamp(2);
pos(:,3) = pos(:,3)*dSamp(3);
L = size(pos,1);
nPos = [nPos; pos];
preClass = [preClass; zeros(L,1) + 3];
postClass = [postClass; zeros(L,1) + 3];
pre = [pre; zeros(L,1)];
post = [post; zeros(L,1) + targ];


dif = cat(3,oPos(:,1)-nPos(:,1)',oPos(:,2)-nPos(:,2)',oPos(:,3)-nPos(:,3)');
dist = sqrt(sum(dif.^2,3));
minDists = min(dist,[],1);


scatter3(oPos(:,1),oPos(:,2),oPos(:,3),'.','k')
hold on
scatter3(nPos(:,1),nPos(:,2),nPos(:,3),'.','r')
hold off

maxDist = 2;
nUse = ones(size(nPos,1),1);
for i = 1:size(pos,1)
    isType = (oPreClass == preClass(i)) &  (oPostClass == postClass(i));
    minDist = min(dist(isType & oUse,i));
    hit = find(isType & oUse & dist(:,i) == minDist,1);
    if minDist<=maxDist
        oUse(hit) = 0;
        nUse(i) = 0;
    end
end







%% use synmat to get post tcr
is125pre = synMat.pre == 125;
post125class = synMat.postClass(is125pre)

synMat.typeNames

[sum(post125class==1) ...
    sum(post125class==2) ...
    sum(post125class==3) ...
    sum(post125class==4)]

is125post = synMat.post == 125;
pre125class = synMat.preClass(is125post);

[sum(pre125class==1) ...
    sum(pre125class==2) ...
    sum(pre125class==3) ...
    sum(pre125class==4)]


%% use syn pos to get post tcr

fromRGC_dat = dsAnchors(synPos.postRGCPos,obI,[2 1 3]);
fromLIN_dat = dsAnchors(synPos.postLinPos,obI,[2 1 3]);


%% motifs

%anaMot = analyzeMotifs(obI);

rgcs = mot.cel.types.rgcs;
tcrs = mot.cel.types.tcrs;
lins = mot.cel.types.lins;

%% syns
syns = synMat;
from125 = synMat.pre == 125;
to125 = synMat.post == 125;
toRGC  = from125 & (synMat.postClass == 1);
toTCR  = from125 &(synMat.postClass== 2);
toLIN = from125 & (synMat.postClass == 3);
toUNK = from125 & (synMat.postClass == 4);
fromRGC = to125 & (synMat.preClass == 1);
fromLIN = to125 & (synMat.preClass == 3);
fromUNK = to125 & (synMat.preClass == 4);

sum(toRGC)
sum(toTCR)
sum(toLIN)
sum(toUNK)

sum(fromRGC)
sum(fromLIN)
sum(fromUNK)

synTarg = synMat.pre(to125);
synTargRGC = intersect(synTarg,rgcs);
synTargTCR = intersect(synTarg,tcrs);
synTargLIN = intersect(synTarg,lins);

synPos125toTCR = cat(1,synMat.synPos(toTCR,:));
synPos125toLIN = cat(1,synMat.synPos(toLIN,:));

synPos125fromRGC = cat(1,mot.syn.synPos{fromRGC,:});
synPos125fromLIN = cat(1,mot.syn.synPos{fromLIN,:});
synPos125fromUNK = cat(1,mot.syn.synPos{fromUNK,:});

%% triads
tri = mot.tri;
triCell = tri.triCell;
triCell = triCell(triCell(:,1) == 125,:);
triTarg = triCell(:,3);
triTargRGC = intersect(triTarg,rgcs);
triTargTCR = intersect(triTarg,tcrs);
triTargLIN = intersect(triTarg,lins);


prim125 = find(tri.triCell(:,1)==125);
useTri = prim125;
%for u = 1:length(prim125);%1:size(tri.triCell,1);

triPoints = drawTriads(tri,useTri);

%synCol = cat(1,triPoints.lineCol,triPoints.ballCol);
% synRad = cat(1,triPoints.lineRad,triPoints.ballRad);
% synType = cat(1,triPoints.lineRad*0+1,triPoints.ballRad*0+2);
% useSyn = cat(2,triPoints.lineGroup, triPoints.ballGroup);
primCell = unique(triPoints.cellGroup(:,1));
secCell = unique(triPoints.cellGroup(:,2));
tertCell = unique(triPoints.cellGroup(:,3));



%% diads

diads = mot.di.diCell;
diads = diads(diads(:,1) == 125,:);
diTarg = diads(:,3);
diTarg = setdiff(diTarg,triTarg);

diTargRGC = intersect(diTarg,mot.cel.types.rgcs);
diTargTCR = intersect(diTarg,mot.cel.types.tcrs);
diTargLIN = intersect(diTarg,lins);


%% color relationship to 124

postTCR = intersect(tcrs,isPost);
triTarg = intersect(tertCell,tcrs);
diTarg = setdiff(diTarg,triTarg);



