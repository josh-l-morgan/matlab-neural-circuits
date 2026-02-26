
%%Test the probability of axons with same seed cell identity being linked
%%together by thalamocortical cells.

clear all
load('MPN.mat')
load([MPN 'obI.mat'])
allEdges = obI.nameProps.edges(:,[2 1]);
seedList = [108 201 903 907];
useCells = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);


preList = useCells.preList;
postList = setdiff(useCells.postList,seedList);
clear preSeed
for p = 1:length(preList)
    isPre = allEdges(:,1) == preList(p);
    isPost = allEdges(isPre,2);
    [isSeed seedNum] = intersect(seedList,isPost);
    preSeed{p} = seedNum;
end

con = zeros(length(preList),length(postList));
for y = 1:length(preList)
    for x = 1:length(postList)
        con(y,x) = sum((allEdges(:,1) == preList(y)) & (allEdges(:,2) == postList(x)));
    end
end

[y x] = find(con>0);
edges = [y x];

seedCon = zeros(length(preList));
for y = 1:length(preList)
    for x = 1:length(preList)
        seedCon(y,x) = length(intersect(preSeed{y},preSeed{x}))>0;
    end
end
subplot(3,1,1)
image(seedCon*100)


tic

preCon = zeros(length(preList));
for y = 1:length(preList)
    
    crossMat = repmat(con(y,:),[size(con,1),1]);
    addMat = crossMat & con;
    preCon(y,:) = sum(addMat,2);
end
subplot(3,1,2)
image(preCon*100)
toc

%% real measures

inSeed = seedCon & preCon;
realInRat = sum(inSeed(:))/sum(preCon(:));

%% rand measures


reps = 1000;

for r = 1:reps
    r
    pick = randperm(size(edges,1));
    newEdges = [edges(pick,1) edges(:,2)];
    
    randCon = con * 0;
    %randCon(sub2ind(size(con),edges(:,1),edges(:,2))) = 1;
        randCon(sub2ind(size(con),newEdges(:,1),newEdges(:,2))) = 1;

    newCon = zeros(length(preList));
    for y = 1:length(preList)
        
        crossMat = repmat(randCon(y,:),[size(con,1),1]);
        addMat = crossMat & randCon;
        newCon(y,:) = sum(addMat,2);
    end
    
    if ~mod(r,100)
        subplot(3,1,3)
        image(newCon*100),pause(.1)
    end
    
    inSeedRand = seedCon & newCon;
    randInRat(r) = sum(inSeedRand(:))/sum(preCon(:));
    
    
end

%%
realInRat
rangeX(randInRat)
median(randInRat)
P = sum(randInRat>= realInRat)/length(randInRat)































