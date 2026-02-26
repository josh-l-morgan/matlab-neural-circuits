%may9

analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\';
mergeDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\';
load([analDir 'tis.mat']);
load([mergeDir 'dsObj.mat']);
curTis=tis;

%clean the tis file
cleaning=1;
if cleaning
    curTis=findBadSynapse(curTis);
end

preTypes=cid2type(curTis.syn.edges(:,2),curTis);
postTypes=cid2type(curTis.syn.edges(:,1),curTis);

%Get all of the synapses that are a ribbon and have a vg3 post
allBPC2VGCidx = find(postTypes{1}==8 & postTypes{3}==1 & curTis.syn.synType'==2);

bpcSubsAll=preTypes{3}(allBPC2VGCidx)';
bpcSubsAlllkup=bpcSubsAll;
bpcSubsAlllkup(bpcSubsAlllkup==0)=24;
bpcSubsAlllkup(bpcSubsAlllkup>26)=24;
%nameRA=curTis.cells.type.subTypeNames{7}(:);
bpcSubNamesAll=curTis.cells.type.subTypeNames{7}(bpcSubsAlllkup)';
obIDsAll=curTis.syn.obID(allBPC2VGCidx);
oldNamesAll=curTis.obI.nameProps.oldNames(obIDsAll);
newNamesAll=curTis.obI.nameProps.names(obIDsAll);
edgesAll=curTis.syn.edges(allBPC2VGCidx,:);
ribPosAll=curTis.syn.pos(allBPC2VGCidx,:);
ribPosAllV=ribPosAll.*[250 250 25];






%% 05/24/22
vgcCidList=[2 3 4 5 10 11 13 14 20];
inputIDs=find(ismember(curTis.syn.edges(:,1),vgcCidList));
inputTypes=preTypes{1}(inputIDs);

outputIDs=find(ismember(curTis.syn.edges(:,2),vgcCidList));
outputTypes=postTypes{1}(outputIDs);

figure(); histogram(outputTypes);

