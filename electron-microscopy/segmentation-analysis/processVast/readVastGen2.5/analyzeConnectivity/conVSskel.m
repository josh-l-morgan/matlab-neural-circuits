
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];
load([TPN 'cellArbors'])
load([MPN 'obI.mat'])

%% Get syn data and select synapses
targCell = 108;
synapses = obI.nameProps.edges;
preWithTarg = unique(synapses(synapses(:,1)==targCell,2));
preWithTarg = preWithTarg(preWithTarg>0);
synWithPre = [];
for i = 1:length(preWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==preWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
synPos = ceil(obI.colStruc.anchors(synObj,:)/8);

conArbor.targCell  = targCell;

%% find distance between pre and post at synapses
minBranchLength = 0;
% syn 3 pos should be  904   537   413

collectMids = [];
goodCheck = synObj*0 + 1;
showPlots = 0;
nearestPost = synObj * 0;
nearestPre = synObj * 0;
procDist = synObj * 0;
for i = 1 : length(synObj)
    
    if showPlots
        scatter(synPos(:,2),synPos(:,3),'b','.')
        hold on
    end
    
    checkPos = synPos(i,:);
    %%Correct for bad sub conversion
    checkPos(1) = checkPos(1) - 1024/8;
    checkPos(2) = checkPos(2) - 1024/8;
    
    preCell = find(cellArbors.cellName == synPre(i));
    postCell = find(cellArbors.cellName == synPost(i));
    
    %% Check pre
    checkArbor = cellArbors.arbors(preCell);
    if length(checkArbor)>0
        useBranch = checkArbor.branchLengths>minBranchLength;
        allMids = cat(1,checkArbor.branches(useBranch).edgeMids);
        collectMids = cat(1,collectMids,allMids);
        if showPlots
            scatter(allMids(:,2),allMids(:,3),'m','.')
        end
        
        %allLengths = cat(2,checkArbor.branches(useBranch).edgeLengths);
        midDist = sqrt((allMids(:,1)-checkPos(1)).^2 + ...
            (allMids(:,2)-checkPos(2)).^2 +(allMids(:,3)-checkPos(3)).^2);
        prePos = allMids(find(midDist== min(midDist),1),:);
        nearestPre(i) = min(midDist);
    else
        goodCheck(i) = 0;
    end
    
    %% Check post
    checkArbor = cellArbors.arbors(postCell);
    if length(checkArbor)>0
        useBranch = checkArbor.branchLengths>minBranchLength;
        allMids = cat(1,checkArbor.branches(useBranch).edgeMids);
        midDist = sqrt((allMids(:,1)-checkPos(1)).^2 + ...
            (allMids(:,2)-checkPos(2)).^2 +(allMids(:,3)-checkPos(3)).^2);
        nearestPost(i) = min(midDist);
        postPos = allMids(find(midDist== min(midDist),1),:);
        
        
        if showPlots
            scatter(allMids(:,2),allMids(:,3),'r','.')
            scatter(checkPos(:,2),checkPos(:,3),'g')
            hold off
            pause
        end
        
        %disp(sprintf('%d and %d',synPre(i),synPost(i)));
        %disp(sprintf('min dist %.1f %.1f',nearestPre(i),nearestPost(i)));
        
        procDist(i) = sqrt((postPos(1)-prePos(1))^2 + (postPos(2)-prePos(2))^2 + (postPos(3)-prePos(3))^2);
        
    else
        goodCheck(i) = 0;
    end
end

betterCheck = goodCheck & (nearestPost<100) & (nearestPre<100);

usePost = nearestPost(betterCheck>0);
usePre = nearestPre(betterCheck>0);
useProc = procDist(betterCheck>0);

distAtSyn.syn2post = usePost;
distAtSyn.syn2pre = usePre;
distAtSyn.pre2post = useProc;

barPre = hist(usePre,0:1:10);
barPost = hist(usePost,0:1:10);
barProc = hist(useProc,0:1:10);

bar([barPre;barPost; barProc]')

conArbor.distAtSyn = distAtSyn;

%% find closest point and total length within 2x dist of closest point

preList = unique(synPre);
preList = preList(preList>0);
postList = unique(synPost);
postList = postList(postList>0);
lookDist = 30;

goodCheck = synObj*0 + 1;
showPlots = 0;
synMat = zeros(length(preList),length(postList));
lengthMat = zeros(length(preList),length(postList),lookDist);

for pre = 1 : length(preList)
    
    disp(sprintf('running pre %d of %d',pre,length(preList)))
    preCell = find(cellArbors.cellName == preList(pre));
    arbor1 = cellArbors.arbors(preCell);
    
    parfor post = 1: length(postList)
        synMat(pre,post) = sum(synPre==preList(pre) & synPost == postList(post));
        
        postCell = find(cellArbors.cellName == postList(post));
        arbor2 = cellArbors.arbors(postCell);
        
        lengthOverlaps = threshOverlap(arbor1,arbor2,lookDist);
        lengthMat(pre,post,:) = lengthOverlaps;
        
        %     lengthOverlaps = postAccess(arbor1,arbor2,lookDist);
        %     lengthMat(pre,post,:) = lengthOverlaps;
        %
        %
        %     lengthOverlaps = postAccess(arbor2,arbor1,lookDist);
        %     lengthMat(pre,post,:) = lengthOverlaps;
        
    end
end

%%
useValsAtDist = 5;
grabMat = lengthMat(:,:,useValsAtDist);
useVals = grabMat(:)>0;


for i = 1:size(lengthMat,3)
    grabMat = lengthMat(:,:,i);
    scatter(grabMat(:),synMat(:),'.','b')
    Xdat = grabMat(useVals);
    Ydat = synMat(useVals);
    
    hold on
    scatter(Xdat,Ydat,'.','r')
    i
    [fit1 gof] = fit(Xdat,Ydat,'poly1');
    maxX = max(Xdat);
    X = [0 maxX];
    Y = X* fit1.p1 + fit1.p2;
    line(X ,Y);
    [rho pc] = corr(Xdat,Ydat);
    text(0,10, num2str(rho));
    hold off
    %xlim([0 mean(grabMat(:))])
    
    trackRho(i) = rho;
    trackSSE(i) = gof.sse;
    pause(.01)
end

bestFit = find(trackRho == max(trackRho(:)));

grabMat = lengthMat(:,:,bestFit);
Xdat = grabMat(useVals);
Ydat = synMat(useVals);
[rho pc] = corr(Xdat,Ydat)

scatter(Xdat,Ydat,'.','r')
hold on

maxX = max(Xdat);
X = [0 maxX];
[fit1 gof] = fit(Xdat,Ydat,'poly1');
Y = X* fit1.p1 + fit1.p2;
line(X ,Y);
hold off

fitDat.useValsAtDist = useValsAtDist;
fitDat.Xdat = Xdat;
fitDat.Ydat = Ydat;
fitDat.fit1 = fit1;
fitDat.gof = gof;
fitDat.rho = rho;
fitDat.corrP = pc;

conArbor.fitDat = fitDat;

mats.preList = preList;
mats.postList = postList;
mats.synMat = synMat;
mats.lengthMat = lengthMat(:,:,bestFit);
mats.fullLengths.lengthMat = lengthMat;
mats.fullLengths.bestFit = bestFit;
mats.fullLengths.trackRho = trackRho;
mats.fullLengths.trackSSE = trackSSE;
mats.fullLengths.readme = 'the full lengthMat includes every calculation of the length overlap of arbors from 1 : the third dimention of the mat.  The length number is in the arbor voxel space';

conArbor.mats = mats;

%% Show syn vs length

colMat = synMat*256/max(synMat(:));
 conArbor.mats.lengthMat;
colMat(:,:,2) = conArbor.mats.lengthMat * 256/max(conArbor.mats.lengthMat(:));
colMat(:,:,3) = conArbor.mats.lengthMat * 0;

image(uint8(colMat))

%%
save([TPN 'conArbor.mat'],'conArbor')


