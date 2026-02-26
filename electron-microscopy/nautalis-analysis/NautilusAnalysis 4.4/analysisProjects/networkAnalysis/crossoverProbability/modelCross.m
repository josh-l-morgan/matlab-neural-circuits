
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])
allEdges = obI.nameProps.edges(:,[2 1]);

useList = obI2cellList_seedInput_RGC_TCR(obI,[108 201]);

targ = find(useList.preList == 1121);
notTarg = setdiff([1:length(useList.preList)],targ);
useList.preList = useList.preList(notTarg);
useList.con = useList.con(notTarg,:);

conPref = seedPreferences([108 201],useList);

crossoverAxons = [2032	2033	2034	2035]



%%
cellNum = length(conPref.cellList);
shared = conPref.sharedAx

targA = find(useList.postList == 108);
Arat = shared(:,targA);

targB = find(useList.postList == 201);
Brat = shared(:,targB);
useShared = setdiff(1:length(shared),[targA targB]);
shared = shared(:,useShared);
shared = shared(:,sum(shared,1)>0);

summed = sum(shared,1)
crossProb = .05;
mixed = min(shared,[],1)./max(shared,[],1)

sortReal = sort(mixed);
bar(sortReal)

reps = 1000;
randMixed = zeros(reps,length(summed));
for r = 1:reps
    for t = 1: length(summed)
    crossed = sum(rand(summed(t),1)<=crossProb);
    mix = [crossed summed(t)-crossed];
    rMixed(t) = min(mix)/max(mix);
    end
    randMixed(r,:) = sort(rMixed);
    
end
meanRandMixed = mean(randMixed,1);
 
subplot(2,1,1)
bar(sortReal)
subplot(2,1,2)
bar(meanRandMixed)
subplot(1,1,1)
bar([sortReal;meanRandMixed]')
