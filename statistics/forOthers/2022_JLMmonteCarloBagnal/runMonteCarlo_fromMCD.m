



clear all
datFolder = 'E:\Bagnal\2022\';
load([datFolder 'mcd.mat'])
f = figure;



%% Assign Variables

useAffSyn = 1; %Use proximity to afferent synapses (1) or any afferent skeleton node (0)
distThresh = 5; % maximum distance counted in model of afferent proximity
useLinDist = 1; %% compute linear distance between synapses along neurite topology
binaryDistance = 0; %% count all distances under dist Thresh as 1 (not 0-1 by proximity)
binarizeProbabilities = 0; %% instead of weighted probabilites of synapse formation, assign everything a 1 or 0;


bin = 10; %%Range of distances used for modeling proximity
reps = 10; %number of model repeats

clf

%% Unpack mcd
skels = mcd.skels;
affs = mcd.affs;

if 0
    useID = mcd.skelIDs([skels.skelType]>0);
    useAffID = mcd.affIDs([affs.skelType]>0);
else
    
    affIDs = unique(mcd.con.conSource);
    skelIDs = unique(mcd.con.conTarget);
    useID = skelIDs;
    useAffID = affIDs;
end


%% Find synaptic probability for nodes based on proximity to afferents

%%Collect
allAffPos = [];
allAffSynPos = [];
for i = 1:length(useAffID)
    targ = find(affIDs == useAffID(i),1);
    s = affs(targ);
    pos = s.pos;
    sPos = pos(s.sNid,:);
    allAffPos = cat(1,allAffPos,pos);
    allAffSynPos = cat(1,allAffSynPos,sPos);
end

if useAffSyn
    prePos = allAffSynPos;
else
    prePos = allAffPos;
end


%%Assign synapse probabilities to all nodes

for i = 1:length(skels)
    disp(sprintf('measuring distances for %d of %d',i,length(skels)))
    s = skels(i);
    pos = s.pos;
    
    %%Find real proximity presynaptic syn nodes to all post synaptic cell nodes
    dif1 = prePos(:,1) - pos(:,1)';
    dif2 = prePos(:,2) - pos(:,2)';
    dif3 = prePos(:,3) - pos(:,3)';
    distsAff2Post = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
    
    
    threshDistsAff2Syn = distThresh-distsAff2Post;
    if binaryDistance %% binarize distances (1 for everything within distance threshold)
        threshDistsAff2Syn =threshDistsAff2Syn>0;
    else
        threshDistsAff2Syn(threshDistsAff2Syn<0) = 0;
    end
    
    
    
    sumWeightAll = sum(threshDistsAff2Syn,1);
    %%Get rid of closest synapse
    maxPower = max(threshDistsAff2Syn,[],1);
    sumWeightAll = sumWeightAll-maxPower;
    
    skels(i).nodeWeights = sumWeightAll;
    skels(i).synWeights = sumWeightAll(s.sNid);
    
end

%%Compare all weights
allSynWeights = [];
allNodeWeights = [];
for i = 1:length(skels)
    allSynWeights = cat(2,allSynWeights,skels(i).synWeights);
    allNodeWeights = cat(2,allNodeWeights,skels(i).nodeWeights);
end


%%normalize weights
maxWeight = max(allNodeWeights(:));
sumWeightSynN = allSynWeights/ maxWeight * 100;
sumWeightAllN = allNodeWeights/ maxWeight * 100;

%%Find fractions of synapses
checkRange = [1:max(sumWeightAllN(:))];

numAll = zeros(length(checkRange),1);
numSyn = numAll;
for i = 1:length(checkRange)
    numAll(i) = sum((sumWeightAllN> (checkRange(i)-bin)) & (sumWeightAllN<=checkRange(i)+bin));
    numSyn(i) = sum((sumWeightSynN> (checkRange(i)-bin)) & (sumWeightSynN<=checkRange(i)+bin));
end
plot(numAll/max(numAll),'k');
hold on
plot(numSyn/max(numSyn),'g');
plot(numSyn./numAll,'r')
hold off
pause(2)
synProb = [0; numSyn./numAll]; %shift one to the right to allow for easy 0 values


%%Assign synapse probabilities to all nodes

for i = 1:length(skels)

    if binarizeProbabilities %everything within distThresh of afferents gets a 1
        skels(i).nodeSynProb = skels(i).nodeWeights>0;
    else
        s = skels(i);
        sumWeightAllN = skels(i).nodeWeights/ maxWeight * 100;
        weightInd = ceil(sumWeightAllN)+1;
        %weightInd(weightInd>100) = 100; %!!!!!!!!!!!
        skels(i).nodeSynProb = synProb(weightInd);
    end

end

mcd.skels = skels;
save([datFolder 'mcd.mat'],'mcd');

sumWeightAllN = skels(i).nodeWeights/ maxWeight * 100;
        weightInd = ceil(sumWeightAllN)+1;
        %weightInd(weightInd>100) = 100; %!!!!!!!!!!!
        skels(i).nodeSynProb = synProb(weightInd);


%% Get real values

meanDistReal = zeros(length(useID),1);
synNums = meanDistReal;
for i = 1:length(useID)
    targ = find(skelIDs==useID(i),1);
    s = skels(targ);
    useNodes = 1:length(s.nid);
    nodeNum = length(useNodes);
    synNum = length(s.sNid);
    synNums(i) = synNum;
    
    sPos = s.pos(s.sNid,:);
    
    if ~useLinDist %euclidian or linear distances
        %%Get distances
        dif1 = sPos(:,1) - sPos(:,1)';
        dif2 = sPos(:,2) - sPos(:,2)';
        dif3 = sPos(:,3) - sPos(:,3)';
        dists = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
        meanDist = sum(dists(:))/(synNum^2-synNum);
        meanDistReal(i) = meanDist;
        
    else
        %%Get linear distance
        d = s.skel2skelLinDist;
        sNid = s.sNid;
        sD = d(sNid,sNid);
        [y x] = find(sD+1);
        sDists = sD(y<x);
        meanDistReal(i) = mean(sDists)/1000;
    end
    
    if 1
        clf
        scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
        hold on
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
        hold off
        title(sprintf('cell %d, %d synapses, mean distance = %0.1f',useID(i),synNum, meanDistReal(i)))
        drawnow
        
    end
    
    
end

mcd.res.meanDistReal = meanDistReal;
save([datFolder 'mcd.mat'],'mcd');

%% Run randomization, do each cell separately for speed
tic
meanDistRand = zeros(length(useID),reps);
clear saveRedistributions
for i = 1:length(useID)
    disp(sprintf('running randomization for %d, %d of %d',useID(i),i,length(useID)))
    targ = find(skelIDs==useID(i),1);
    s = skels(targ);
    useNodes = 1:size(s.pos,1);
    useWeights = s.nodeSynProb(useNodes);
    nodeNum = length(useNodes);
    synNum = length(s.sNid);
    saveRedistributions(i).sPos = [];
    if 1
        clf
        sPos = s.pos(s.sNid,:);
        
        scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
        hold on
        scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
        scatter3(allAffPos(:,1),allAffPos(:,2),allAffPos(:,3),'b','.')
        scatter3(allAffSynPos(:,1),allAffSynPos(:,2),allAffSynPos(:,3),'m','o')
        drawnow
    end
    
    for r = 1:reps
        
        %sNidR = useNodes(randperm(nodeNum,synNum));
        
        sNidR = randsample(useNodes,synNum,true,useWeights);
        sPos = s.pos(sNidR,:);
        
        
        if r<10
            scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','g')
            saveRedistributions(i).sPos = cat(1,saveRedistributions(i).sPos,sPos);
            drawnow
        end
        
        if 0
            scatter3(s.pos(:,1),s.pos(:,2),s.pos(:,3),'.','k')
            hold on
            scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','r')
            hold off
            pause(.1)
            
        end
        
        if ~useLinDist
            %%Get distances
            dif1 = sPos(:,1) - sPos(:,1)';
            dif2 = sPos(:,2) - sPos(:,2)';
            dif3 = sPos(:,3) - sPos(:,3)';
            dists = sqrt(dif1.^2 + dif2.^2 + dif3.^2)/1000;
            meanDist = sum(dists(:))/(synNum^2-synNum);
            meanDistRand(i,r) = meanDist;
            
        else
            %%Get linear distance
            d = skels(i).skel2skelLinDist;
            sNid = sNidR;
            sD = d(sNid,sNid);
            [y x] = find(sD+1);
            sDists = sD(y<x);
             meanDistRand(i,r) = mean(sDists)/1000;
        end
        
        
        
    end
    pause(.1)
end

toc

mcd.res.meanDistRand = meanDistRand;
save([datFolder 'mcd.mat'],'mcd');

%% analyze results

mMeanDistReal = mean(meanDistReal)
mMeanDistRandom = mean(meanDistRand,1);

realRes = mMeanDistReal;
sVals = sort(mMeanDistRandom,'ascend');

P = mean(sVals<=realRes)
ci95 = [sVals(max(1,round(0.025*reps))) sVals(round(0.975*reps))]

mcd.res.P = P;
mcd.res.ci95 = ci95;
mcd.res.sVals = sVals;
mcd.res.mMeanDistReal = mMeanDistReal;
mcd.res.testStat = realRes;

save([datFolder 'mcd.mat'],'mcd');

clf
subplot(2,1,1)
hRange = [0:.1:max(sVals)];
hRand = hist(sVals,hRange);
hRand = hRand/reps;
bar(hRange,hRand,'edgecolor','b','facecolor','b')
hold on
halfHeight = max(hRand(:))/2;
plot([ci95(1) ci95(1)],[halfHeight-.1*halfHeight halfHeight+.1*halfHeight],'k')
plot([ci95(2) ci95(2)],[halfHeight-.1*halfHeight halfHeight+.1*halfHeight],'k')
plot([ci95(1) ci95(2)],[halfHeight halfHeight],'k')

plot([realRes realRes],[0 halfHeight],'r')
scatter([realRes],[halfHeight],'r','o','filled')

title(sprintf('toSyn %d, distThresh %d, linear dist %d\n, binaryDist %d, bin %d, reps %d' ,...
    useAffSyn, distThresh, useLinDist, binaryDistance, bin, reps))
    

%% Display randomizations
clf
axis 'equal'
axis 'off'
f.Color = [1 1 1]
set(gca,'clipping','off')
view(-4,19)
hold on
try 
    aNscat.delete
    aSscat.delete
end
dSize = 5;
aNscat = scatter3(allAffPos(:,1),allAffPos(:,2),allAffPos(:,3),2*dSize,'g','o','filled');
aSscat = scatter3(allAffSynPos(:,1),allAffSynPos(:,2),allAffSynPos(:,3),10*dSize,'m','o','filled');
aNscat.MarkerEdgeAlpha = 0;
aNscat.MarkerFaceAlpha = .1;
aSscat.MarkerEdgeAlpha = 0;
aSscat.MarkerFaceAlpha = 1;

for i = 1:length(skels)
    sPos = saveRedistributions(i).sPos;
    pos = skels(i).pos;
    scatter3(pos(:,1),pos(:,2),pos(:,3),2*dSize,[.4 .2 0],'filled',...
        'markeredgealpha',0,'markerfacealpha',.3)
    scatter3(sPos(:,1),sPos(:,2),sPos(:,3),10*dSize,'o','k','filled',...
        'markeredgealpha',0,'markerfacealpha',.6)
    drawnow
end

%% Print figure
if 0
    %fDir = uigetdir;
    filename = ['D:\WorkDocs\Publications\Bagnal\2022_revisions\pics' ...
        '\' 'unbiased']
    set(gcf,'renderer','Painters')
    print('-depsc','-tiff','-r300', '-painters',[filename,'.eps'])
    
end






