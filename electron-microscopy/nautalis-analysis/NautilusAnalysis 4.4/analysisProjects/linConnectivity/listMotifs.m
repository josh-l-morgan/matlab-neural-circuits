
load('MPN.mat');
load([MPN 'obI.mat']);
targCell = 125;

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


%% find all combinations

targSyns = find(isTarg);
maxDist = 5;

for i = 1:size(sm.skelNodes,1)
    disp(sprintf('checking motif %d of %d',i,size(sm.skelNodes,1)))
    close = find(sm.skelDist(i,:)<=maxDist);
    synNum(i) = length(close);
    
    if isempty(close)
        mot{i} = 0;
    else
        %%
        preClass = sm.preClass(close);
        postClass = sm.postClass(close);
        pre = sm.pre(close);
        post = sm.post(close);
        
        
        rgcIn = (preClass == 1) & (post == targCell);
        tcOut = (postClass == 2) & (pre == targCell);
        linIn = (preClass == 3) & (post == targCell);
        linOut = (postClass == 3) & (pre == targCell);
        unkIn = (preClass == 4) & (post == targCell);
        
        
        %% if no inputs
        if ~sum(rgcIn | linIn | unkIn);
           m.outOnlyTC = sum(postClass == 2);
           m.outOnlyLin = sum(postClass == 3);
        else
            m.outOnlyTC = 0;
            m.outOnlyLin = 0;
        end
        
        %% check RGC inputs
        m.rgcIn = sum(rgcIn);
        
        checkRGC = unique(pre(rgcIn));
        checkRGC = setdiff(checkRGC,0);
        m.rgcLab = length(checkRGC);        
        m.rgcTc = 0;
        m.rgcTcDoub = 0;
        m.rgcLin = 0;
        m.rgcLinTc = 0;
        m.rgcLinLinIn = 0;
        m.rgcLinLinOut = 0;
        m.rgcLinLinDoub = 0;
        m.rgcOneOut = 0;
        m.rgcLinTcDi = 0;
        m.rgcLinLinDi = 0; %chain of two
        
        for r = 1:length(checkRGC)
            rgcID = checkRGC(r);
            isID = pre==rgcID;
            rPostID = post(isID);
            rPostClass = postClass(isID);
            
            if isempty(setdiff(rPostID,targCell))
                m.rgcOneOut = m.rgcOneOut + 1;
            end
            
            m.rgcTc = m.rgcTc + sum(rPostClass == 2);
            m.rgcLin = m.rgcLin + sum(rPostClass == 3);

            checkPost = setdiff(unique(rPostID),[0 targCell]);
            m.rgcTcDoub = m.rgcTcDoub + sum((rPostClass==2) & (rPostID==0));
            m.rgcLinLinDoub = m.rgcLinLinDoub + sum((rPostClass==3) & (rPostID==0));
            
            %% find rgc triads
            for p = 1:length(checkPost)
                
                if sum((post == checkPost(p)) & (postClass == 2)) % if post is TC
                    isPost = (pre==targCell) & (post == checkPost(p)) ;
                    if sum(isPost)
                        m.rgcLinTc = m.rgcLinTc + 1;
                    else
                        m.rgcTcDoub = m.rgcTcDoub + 1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3)) %if post is LIN
                    isPre = (pre == checkPost(p)) & (post == targCell)& (preClass == 3);
                    if sum(isPre)
                        m.rgcLinLinIn = m.rgcLinLinIn + 1;
                    else
                        m.rgcLinLinDoub = m.rgcLinLinDoub +1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3))%if post is LIN
                    isPost = (pre==targCell) & (post == checkPost(p)) & (postClass == 3);
                    if sum(isPost)
                        m.rgcLinLinOut = m.rgcLinTc + 1;
                    else
                        m.rgcLinLinDoub = m.rgcLinLinDoub +1;
                    end
                end
                
            end
            
            %%Find RGC LIN TC diads
           postTarg = post((pre == targCell) & (postClass == 2));
           m.rgcLinTcDi = m.rgcLinTcDi + sum(postTarg == 0);
           m.rgcLinTcDi = m.rgcLinTcDi + length(setdiff(postTarg,[0; checkPost]));
           
            %%Find RGC LIN LIN diads
           postTarg = post((pre == targCell) & (postClass == 3));
           m.rgcLinLinDi = m.rgcLinLinDi + sum(postTarg == 0);
           m.rgcLinLinDi = m.rgcLinLinDi + length(setdiff(postTarg,[0; checkPost]));
            
        end
        
        %% check UNK inputs
        m.unkIn = sum(unkIn);
        
        checkUNK = unique(pre(unkIn));
        checkUNK = setdiff(checkUNK,0);
        m.unkLab = length(checkUNK);        
        m.unkTc = 0;
        m.unkTcDoub = 0;
        m.unkLin = 0;
        m.unkLinTc = 0;
        m.unkLinLinIn = 0;
        m.unkLinLinOut = 0;
        m.unkLinLinDoub = 0;
        m.unkOneOut = 0;
        m.unkLinTcDi = 0;
        m.unkLinLinDi = 0; %chain of two
        
        for r = 1:length(checkUNK)
            unkID = checkUNK(r);
            isID = pre==unkID;
            rPostID = post(isID);
            rPostClass = postClass(isID);
            
            if isempty(setdiff(rPostID,targCell))
                m.unkOneOut = m.unkOneOut + 1;
            end
            
            m.unkTc = m.unkTc + sum(rPostClass == 2);
            m.unkLin = m.unkLin + sum(rPostClass == 3);

            checkPost = setdiff(unique(rPostID),[0 targCell]);
            m.unkTcDoub = m.unkTcDoub + sum((rPostClass==2) & (rPostID==0));
            m.unkLinLinDoub = m.unkLinLinDoub + sum((rPostClass==3) & (rPostID==0));
            
            %% find lin triads
            for p = 1:length(checkPost)
                
                if sum((post == checkPost(p)) & (postClass == 2)) % if post is TC
                    isPost = (pre==targCell) & (post == checkPost(p)) ;
                    if sum(isPost)
                        m.unkLinTc = m.unkLinTc + 1;
                    else
                        m.unkTcDoub = m.unkTcDoub + 1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3)) %if post is LIN
                    isPre = (pre == checkPost(p)) & (post == targCell)& (preClass == 3);
                    if sum(isPre)
                        m.unkLinLinIn = m.unkLinLinIn + 1;
                    else
                        m.unkLinLinDoub = m.unkLinLinDoub +1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3))%if post is LIN
                    isPost = (pre==targCell) & (post == checkPost(p)) & (postClass == 3);
                    if sum(isPost)
                        m.unkLinLinOut = m.unkLinTc + 1;
                    else
                        m.unkLinLinDoub = m.unkLinLinDoub +1;
                    end
                end
                
            end
            
            %%Find RGC LIN TC diads
           postTarg = post((pre == targCell) & (postClass == 2));
           m.unkLinTcDi = m.unkLinTcDi + sum(postTarg == 0);
           m.unkLinTcDi = m.unkLinTcDi + length(setdiff(postTarg,[0; checkPost]));
           
            %%Find RGC LIN LIN diads
           postTarg = post((pre == targCell) & (postClass == 3));
           m.unkLinLinDi = m.unkLinLinDi + sum(postTarg == 0);
           m.unkLinLinDi = m.unkLinLinDi + length(setdiff(postTarg,[0; checkPost]));
            
        end
        
         %% check LIN inputs
        m.linIn = sum(linIn);
        
        checkLIN = unique(pre(linIn));
        checkLIN = setdiff(checkLIN,0);
        m.linInLab = length(checkLIN);        
        m.linInTc = 0;
        m.linInTcDoub = 0;
        m.linInLin = 0;
        m.linInLinTc = 0;
        m.linInLinLinIn = 0;
        m.linInLinLinOut = 0;
        m.linInLinLinDoub = 0;
        m.linInOneOut = 0;
        m.linInLinTcDi = 0;
        m.linInLinLinDi = 0; %chain of two
        
        for r = 1:length(checkLIN)
            linID = checkLIN(r);
            isID = pre==linID;
            rPostID = post(isID);
            rPostClass = postClass(isID);
            
            if isempty(setdiff(rPostID,targCell))
                m.linInOneOut = m.linInOneOut + 1;
            end
            
            m.linInTc = m.linInTc + sum(rPostClass == 2);
            m.linInLin = m.linInLin + sum(rPostClass == 3);

            checkPost = setdiff(unique(rPostID),[0 targCell]);
            m.linInTcDoub = m.linInTcDoub + sum((rPostClass==2) & (rPostID==0));
            m.linInLinLinDoub = m.linInLinLinDoub + sum((rPostClass==3) & (rPostID==0));
            
            %% find lin triads
            for p = 1:length(checkPost)
                
                if sum((post == checkPost(p)) & (postClass == 2)) % if post is TC
                    isPost = (pre==targCell) & (post == checkPost(p)) ;
                    if sum(isPost)
                        m.linInLinTc = m.linInLinTc + 1;
                    else
                        m.linInTcDoub = m.linInTcDoub + 1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3)) %if post is LIN
                    isPre = (pre == checkPost(p)) & (post == targCell)& (preClass == 3);
                    if sum(isPre)
                        m.linInLinLinIn = m.linInLinLinIn + 1;
                    else
                        m.linInLinLinDoub = m.linInLinLinDoub +1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3))%if post is LIN
                    isPost = (pre==targCell) & (post == checkPost(p)) & (postClass == 3);
                    if sum(isPost)
                        m.linInLinLinOut = m.linInLinTc + 1;
                    else
                        m.linInLinLinDoub = m.linInLinLinDoub +1;
                    end
                end
                
            end
            
            %%Find RGC LIN TC diads
           postTarg = post((pre == targCell) & (postClass == 2));
           m.linInLinTcDi = m.linInLinTcDi + sum(postTarg == 0);
           m.linInLinTcDi = m.linInLinTcDi + length(setdiff(postTarg,[0; checkPost]));
           
            %%Find RGC LIN LIN diads
           postTarg = post((pre == targCell) & (postClass == 3));
           m.linInLinLinDi = m.linInLinLinDi + sum(postTarg == 0);
           m.linInLinLinDi = m.linInLinLinDi + length(setdiff(postTarg,[0; checkPost]));
            
        end
        
        %% 
         %% check LIN inputs
        m.targOut = sum(pre==targCell);
        
        checkLIN = targCell;
        m.targTc = 0;
        m.targTcDoub = 0;
        m.targLin = 0;
        m.targLinTc = 0;
        m.targLinLinIn = 0;
        m.targLinLinOut = 0;
        m.targLinLinDoub = 0;
        m.targOneTc = 0;
        m.targOneLin = 0;
        m.targLinTcDi = 0;
        m.targLinLinDi = 0; %chain of two
        
        for r = 1:length(checkLIN)
            linID = checkLIN(r);
            isID = pre==linID;
            rPostID = post(isID);
            rPostClass = postClass(isID);
            
%             if isempty(setdiff(rPostID,targCell))
%                 m.targOneOut = m.targOneOut + 1;
%             end
            
            m.targTc = m.targTc + sum(rPostClass == 2);
            m.targLin = m.targLin + sum(rPostClass == 3);

            checkPost = setdiff(unique(rPostID),[0 ]);
            m.targTcDoub = m.targTcDoub + sum((rPostClass==2) & (rPostID==0));
            m.targLinLinDoub = m.targLinLinDoub + sum((rPostClass==3) & (rPostID==0));
            
            %% find lin triads
            for p = 1:length(checkPost)
                
                if sum((post == checkPost(p)) & (postClass == 2)) % if post is TC
                    isPost = (pre==targCell) & (post == checkPost(p)) ;
                    if sum(isPost)
                        m.targLinTc = m.targLinTc + 1;
                    else
                        m.targTcDoub = m.targTcDoub + 1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3)) %if post is LIN
                    isPre = (pre == checkPost(p)) & (post == targCell)& (preClass == 3);
                    if sum(isPre)
                        m.targLinLinIn = m.targLinLinIn + 1;
                    else
                        m.targLinLinDoub = m.targLinLinDoub +1;
                    end
                end
                
                if sum((post == checkPost(p)) & (postClass == 3))%if post is LIN
                    isPost = (pre==targCell) & (post == checkPost(p)) & (postClass == 3);
                    if sum(isPost)
                        m.targLinLinOut = m.targLinTc + 1;
                    else
                        m.targLinLinDoub = m.targLinLinDoub +1;
                    end
                end
                
            end
            
            %%Find RGC LIN TC diads
           postTarg = post((pre == targCell) & (postClass == 2));
           m.targLinTcDi = m.linInLinTcDi + sum(postTarg == 0);
           m.targLinTcDi = m.linInLinTcDi + length(setdiff(postTarg,[0; checkPost]));
           
            %%Find RGC LIN LIN diads
           postTarg = post((pre == targCell) & (postClass == 3));
           m.targLinLinDi = m.targLinLinDi + sum(postTarg == 0);
           m.targLinLinDi = m.targLinLinDi + length(setdiff(postTarg,[0; checkPost]));
            
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



