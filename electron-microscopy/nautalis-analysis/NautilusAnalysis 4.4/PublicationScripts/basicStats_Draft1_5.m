
load('MPN.mat')
load([MPN 'obI.mat'])


%% Seed Cell Stats
targSeed = 201;

%%Pick Data
useList = obI2cellList_seedInput_RGC_TCR(obI,targSeed);
con = useList.con;
allEdges = obI.nameProps.edges(:,[2 1]);

%%Pre and Post
preTarg = preTo(allEdges,targSeed)
postTarg = postTo(allEdges,targSeed)

axNum = size(preTarg,1)
synNum  = sum(preTarg(:,2))

postList = useList.postList
allSynSum = sum(con(:))
offSeedSynSum = allSynSum-synNum

offSeedSynSum/(length(postList)-1)

%% axon divergence

seedList =  [108 201]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
con = useList.con;

diverge = sum(con>0,2)
axSyn = sum(con,2)
scatter(axSyn,diverge)

bigEnough = diverge(axSyn>30)
mean(bigEnough)
length(bigEnough)
std(bigEnough)
rangeX(bigEnough,.95)


%% dLGN volum/ cell number

TCinVolume = 3000; %~
emVolume = .380 * .300 * .280;
brainAtlasSections = [
    120 400 200
    120 500 250
    120 500 250
    120 600 400
    120 600 400
    120 500 400
    120 500 400
    120 500 400
    120 500 300
    120 400 300
    120 400 150
    ]
bAVolume = sum(prod(brainAtlasSections/1000,2));
averageVol = 0.26;%mm^3
averageNeuronNumber = 16900;
emFrac = emVolume/averageVol    
estimatedTCnum = TCinVolume/emFrac
    
    
    
    