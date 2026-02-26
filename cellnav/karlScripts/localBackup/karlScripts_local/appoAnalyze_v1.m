%% check the syn to appo dists
allRGCcids=type2cid({'rgc'},{'all'},curTis);
vgcCidList=[2 3 4 5 10 11 13 14 20];
allRGCcids=allRGCcids{1}';
rgcIDs=find(ismember(curTis.syn.edges(:,1),allRGCcids));
vgcIDs=find(ismember(curTis.syn.edges(:,2),vgcCidList));

synIDs=intersect(rgcIDs,vgcIDs);
synLocs=curTis.syn.pos(synIDs,:);
%synLocs=synLocs.*[250 250 25];
[gcl,inl,synDepths]=getIPLdepth(synLocs(:,3),synLocs(:,2),synLocs(:,1),[],[]);
%distMat=pdist2(appoLocs,synLocs);
%figure(); histogram(distMat(:));
%appoSynMinDist=min(distMat,[],2);
%figure(); histogram(appoSynMinDist,115);



%% get all outputs for a vgc
vgcCidList=2;
for i=1:length(vgcCidList)
    curVGCcid=vgcCidList(i);
    curV2RIDS=intersect(synIDs,find(curTis.syn.edges(:,2)==curVGCcid));
    curVtargCids=curTis.syn.edges(curV2RIDS,1);
    targetCounts=arrayfun(@(x)length(find(curVtargCids==x)),unique(curVtargCids),'Uniform',false);
    targetCountMat=cell2mat(targetCounts);
    results=horzcat(unique(curVtargCids),targetCountMat);
    resultsSorted=sortrows(results,2,'descend');
    
end

%% as a test, bpc appositions close to the cid2-cid2002 synapses.
synList=find(curTis.syn.edges(:,1)==2002 & curTis.syn.edges(:,2)==2);
synLocs=curTis.syn.pos(synList,:);
outSynList=find(curTis.syn.edges(:,2)==2);
randSynList=randsample(outSynList,length(synList));
randSynLocs=curTis.syn.pos(randSynList,:);
% manually created the apposition lists from cellNav. - working on the
% getting the function to work. Uses the globs from the GUI, so I'll have
% to hack those to use the function independently.

distSynApp=pdist2(bpc_cid2_appo,synLocs);
distRandApp=pdist2(bpc_cid2_appo,randSynLocs);

histcounts(distRandApp(:));

%allBPCcids=cell2mat(type2cid({'bpc'},{'all'},curTis));
allBPCcids=cid2;
all4owcids=cell2mat(type2cid({'rgc'},{'4ow'},curTis));

appoArray=cell(3,1);
appoCids={2,2002,allBPCcids};

appoArray{1}=getAppos(appoCids{1},appoCids{2});
appoArray{2}=getAppos(appoCids{1},appoCids{3});
appoArray{3}=getAppos(appoCids{2},appoCids{3});

synArray=curTis.syn.pos(curTis.syn.edges(:,1)==2002&curTis.syn.edges(:,2)==2,[3 1 2]);

BVRmat=pdist2(synArray,appoArray{3});
distFilt=min(BVRmat,[],1);

%% test
aaa=appoArray{3}(:,[3 2 1]).*[250 250 25];

bpcAppoNearVGC=aaa(distFilt<10,:);

for i=length(bpcAppoNearVGC(:,1)):-1:1
curLoc=round(bpcAppoNearVGC(i,:));
clipboard('copy',curLoc);
pause();
end

%% Feb 21st. Do VGCs form RGC outputs near BPC inputs?
%starting with a couple structures from earlier. The vgcCidList and the
%allBPCcids
%get the IDs of all the BPC inputs to VGCs
allBVsynIDs=find(ismember(curTis.syn.edges(:,1),vgcCidList)&ismember(curTis.syn.edges(:,2),allBPCcids));
%get all the positions of the BPC-VGC synapses
syn_BPC_VGC=curTis.syn.pos(allBVsynIDs,:);

%get all the appositions of VGCs to OFFa
%appo_VGC_OFFa = [];
%get all the appositions of BPCs to VGCs
%appo_BPC_VGC = [];

%get all the OFFa cids
offaCidList=cell2mat(type2cid({'rgc'},{'4ow'},curTis));
%get the W3 cids
w3CidList=cell2mat(type2cid({'rgc'},{'51'},curTis));
%get the amc cids
amcCidList=cell2mat(type2cid({'amc'},{'all'},curTis));
%get all the synapses for VGCs to OFFa
allVOsynIDs=find(ismember(curTis.syn.edges(:,1),offaCidList)&ismember(curTis.syn.edges(:,2),vgcCidList));
%get all the synapses for VGCs to W3
allVWsynIDs=find(ismember(curTis.syn.edges(:,1),w3CidList)&ismember(curTis.syn.edges(:,2),vgcCidList));
%get the positions of those syns
syn_VGC_AMCa=curTis.syn.pos(curTis.syn.edges(:,1)==0&ismember(curTis.syn.edges(:,2),vgcCidList),:);
syn_VGC_AMCb=curTis.syn.pos(ismember(curTis.syn.edges(:,1),amcCidList)&ismember(curTis.syn.edges(:,2),vgcCidList),:);
syn_VGC_AMC=vertcat(syn_VGC_AMCa,syn_VGC_AMCb);

syn_VGC_OFFa=curTis.syn.pos(allVOsynIDs,:);
syn_VGC_W3=curTis.syn.pos(allVWsynIDs,:);
figz=0;
if figz
plotPts={appo_VGC_OFFa,appo_BPC_VGC,syn_BPC_VGC,syn_VGC_OFFa};
colors=[[0 1 1];[1 0 1];[1 1 0];[0 1 0]];
figure();
hold on
for i=3:length(plotPts)
    curDat=plotPts{i};
    scatter3(curDat(:,1),curDat(:,2),curDat(:,3),25,colors(i,:),'o','filled');
end
end

%set up for the distance calculations between the different classes of
%synapses.
distMat_BPC_VGC_OFFa=pdist2(syn_BPC_VGC,syn_VGC_OFFa);
distMat_BPC_VGC_AMC=pdist2(syn_BPC_VGC,syn_VGC_AMC);
distMat_BPC_VGC_W3=pdist2(syn_BPC_VGC,syn_VGC_W3);

distMin_BPC_VGC_OFFa=min(distMat_BPC_VGC_OFFa,[],1);
distMin_BPC_VGC_OFFa_inv=min(distMat_BPC_VGC_OFFa,[],2);
distMin_BPC_VGC_W3=min(distMat_BPC_VGC_W3,[],1);
distMin_BPC_VGC_W3_inv=min(distMat_BPC_VGC_W3,[],2);
%this won't work because there are a lot more synapses to amcs than to RGCs
%need to rethink what I'm doing here.

distMin_BPC_VGC_AMC=min(distMat_BPC_VGC_AMC,[],1);
distMin_BPC_VGC_AMC_inv=min(distMat_BPC_VGC_AMC,[],2);

%distMat_BPC_VGC_AMC_short=randsample(distMin_BPC_VGC_AMC,length(distMin_BPC_VGC_OFFa));
histDat_BPC_VGC_OFFa=histcounts(distMin_BPC_VGC_OFFa,[0:0.2:5]);
%histDat_BPC_VGC_AMCsh=histcounts(distMat_BPC_VGC_AMC_short,[0:0.2:5]);
figA=figure();
hold on
scatter(1:length(distMin_BPC_VGC_OFFa_inv),distMin_BPC_VGC_OFFa_inv,25,'c.');
scatter(1:length(distMin_BPC_VGC_AMC_inv),distMin_BPC_VGC_AMC_inv,25,'m.');
scatter(1:length(distMin_BPC_VGC_W3_inv),distMin_BPC_VGC_W3_inv,25,'g.');




%plot([0:0.1:4.9],histDat_BPC_VGC_OFFa,10,'m.');
for i=1:100
    distMat_BPC_VGC_AMC_short=randsample(distMin_BPC_VGC_AMC,length(distMin_BPC_VGC_OFFa));
    %distMat_BPC_VGC_AMC_short_inv=randsample(distMin_BPC_VGC_AMC_inv,length(distMin_BPC_VGC_OFFa_inv));
    histDat_BPC_VGC_AMCsh=histcounts(distMat_BPC_VGC_AMC_short,[0:0.2:5]);
    scatter([0:0.2:4.8]+rand(1)/25,histDat_BPC_VGC_AMCsh,5,'c');
    %scatter(1:478,distMin_BPC_VGC_OFFa_inv,5,'m.');
end
scatter([0:0.2:4.8],histDat_BPC_VGC_OFFa,25,'mo','filled');

%histogram(distMin_BPC_VGC_OFFa,50);

%% Tuesday evening

figure();
hold on
histDat_BPC_VGC_OFFa_inv=histcounts(distMin_BPC_VGC_OFFa_inv,[0:0.5:10]);
histDat_BPC_VGC_AMC_inv=histcounts(distMin_BPC_VGC_AMC_inv,[0:0.5:10]);
plot([0:0.5:9.5],histDat_BPC_VGC_OFFa_inv,'m');
plot([0:0.5:9.5],histDat_BPC_VGC_AMC_inv,'c');

