%% EI balance by RGC target subtype

%This script is for determining whether different types of RGC targets are
%receiving information with different balances of excitation and
%inhibition.
%It is known that VGCs receive strong OFF inhibition from surrounding AMCs,
%but how this input is combined with center BPC excitation and formed into
%output signals to different classes of RGC is unknown.

%load skeletons if necessary
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
vgcCidList=[2 3 4 5 13 14];
loadSkel=0;
if loadSkel
    allSkels=cell(1,6);
    for i=1:length(vgcCidList)
        
        curCid=vgcCidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end

%% Main
%% quick check on the types of rgc targets and frequencies
targetTypeDat=cid2type(curTis.syn.edges(ismember(curTis.syn.edges(:,2),vgcCidList),1),curTis);
rgcTargSubs=targetTypeDat{3}(targetTypeDat{1}==1);
histDat=histcounts(rgcTargSubs,[1:length(curTis.cells.type.subTypeNames{1})]);
targSubIDs=find(histDat>10);
targSubsNames=curTis.cells.type.subTypeNames{1}(targSubIDs);
f0=figure();
hold on
bar(histDat(targSubIDs));
xticks([1:length(targSubIDs)]);
xticklabels(targSubsNames);
ylabel('# of synapses from VG3');
xlabel('RGC subtype');
title('RGC types receiving >10 VG3 synapses');

%% Go through skels
lc=13;
% 4ow=23; 51=24;
rgcSubTypeList=[21 23 24 29 51];
bpcOFFsubs=[1:5,15,18,19,24];
bpcONsubs=[6:12,14,17,20,21];

f1=figure();
hold on
colMapA=hsv(8);

%f2=figure();
if fig3
scf3=figure();
title(['Input near VG3->RGC outputs by RGC subtype; lc=' num2str(lc)])
hold on
end

for i=1:length(allSkels)
    curCid=vgcCidList(i);
    curSM=allSkels{i}.sm;
    curSyns=curSM.syn;
    [n,m,curSyns.depth]=getIPLdepth(curSyns.pos(:,3),curSyns.pos(:,1),curSyns.pos(:,2),[],[]);
    inputs=curSyns.post==curCid;
    outputs=curSyns.pre==curCid;
    preTypeDat=cid2type(curSyns.pre,curTis);
    postTypeDat=cid2type(curSyns.post,curTis);
    %super sketch, but I'm just saying that input&~bpc = amc
    bpcSynIDs=find(preTypeDat{1}==7 & inputs'==1);
    amcSynIDs=find(inputs'==1 & preTypeDat{1}~=7);
    onBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcONsubs));
    offBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcOFFsubs));
    curDistMat=curSM.syn2Skel.syn2SynDist;
    vMat=zeros(size(curSM.syn2Skel.syn2SynDist));
    vMat(bpcSynIDs,:)=1;
    vMat(amcSynIDs,:)=-1;
    infMat=vMat.*exp(-curDistMat./lc);
    EIbalMat=sum(infMat,1);
    
    figure(f1);
    scatter(curSyns.depth(onBpcSynIDs),1+i/20,25,colMapA(i,:),'filled');
    scatter(curSyns.depth(offBpcSynIDs),2+i/20,25,colMapA(i,:),'filled');
    scatter(curSyns.depth(amcSynIDs),3+i/20,25,colMapA(i,:),'filled');
    
    %histogram(curSyns.depth);
    if sanityCheck2
        scf1=figure();
        histogram(EIbalMat,50);
        title(['cid' num2str(curCid) ': E/I Bal; all syns; lc=' num2str(lc)]);
    end
    %curSM.syn2Skel.syn2SynDist;
    
    colMapB=hsv(256);
    
    for j=2%1:length(rgcSubTypeList)
        curRGCsubType=rgcSubTypeList(j);
        curSynIDs=find(postTypeDat{1}==1 & postTypeDat{3}==curRGCsubType & outputs'==1);
        figure(f1);
        scatter(curSyns.depth(curSynIDs),3+i/20+j,25,colMapA(i,:),'filled');
        %        scatter(curSyns.depth(curSynIDs),3+i/20+j,25,colMapB(round(j/length(rgcSubTypeList)*100+35),:),'filled');
        if sanityCheck
            if fig2
                scf2=figure();
                title(['cid' num2str(curCid) ': E/I bal to ' curTis.cells.type.subTypeNames{1}{curRGCsubType} ...
                    ' x VGC Input synapse'])
                hold on
                for k=1:length(curSynIDs)
                    bpcIn=curDistMat(bpcSynIDs,curSynIDs(k));
                    bpcInNear=bpcIn(bpcIn<30);
                    amcIn=curDistMat(amcSynIDs,curSynIDs(k));
                    amcInNear=amcIn(amcIn<30);
                    scatter(repmat(k,length(bpcInNear)),bpcInNear,25,'cs','filled');
                    scatter(repmat(k+0.25,length(amcInNear)),amcInNear,25,'md','filled');
                end
                xlabel('synapse #');
                ylabel('distance from synapse (um)');
            end
            if fig3
                figure(scf3);
                for k=1:length(curSynIDs)
                    bpcIn=curDistMat(bpcSynIDs,curSynIDs(k));
                    bpcOffIn=curDistMat(offBpcSynIDs,curSynIDs(k));
                    bpcOnIn=curDistMat(onBpcSynIDs,curSynIDs(k));
                    bpcInNear=bpcIn(bpcIn>0.05);
                    bpcOffInNear=bpcOffIn(bpcOffIn>0.05);
                    bpcOnInNear=bpcOnIn(bpcOnIn>0.05);
                    amcIn=abs(curDistMat(amcSynIDs,curSynIDs(k)));
                    amcInNear=amcIn(amcIn>0.05);
                    %scatter(repmat(k,length(bpcInNear)),bpcInNear,25,'cs','filled');
                    scatter(repmat(j*10+i,length(bpcOffInNear)),bpcOffInNear,35,'md','filled');
                    scatter(repmat(j*10+i+.1,length(bpcOnInNear)),bpcOnInNear,35,'cs','filled');
                    scatter(repmat(j*10+i+0.25,length(amcInNear)),amcInNear,15,'go');
                end
%                 xlabel('synapse #');
%                 ylabel('influence');
%                 set(gca, 'YDir','reverse');
%                 xticks([0:length(rgcSubTypeList)]+.2);
%                 xticklabels({'4i','4ow','51','63','37'});
                
            end
        end
    end
    
    figure(scf3);
    %xlabel('synapse #');
    ylabel('influence');
    set(gca, 'YDir','reverse');
    xticks([14:10:length(rgcSubTypeList)*10+12]);
    xticklabels({'4i','4ow','51','63','37'});
    xlim([5,length(rgcSubTypeList)*10+5]);
    %     figure();
    %     hold on
    %     scatter3(curSM.arbor.nodes.pos(:,1),curSM.arbor.nodes.pos(:,2), ...
    %         curSM.arbor.nodes.pos(:,3),curSM.arbor.nodes.rad*5,'.k');
    %     scatter3(curSyns.pos(:,1)*10,curSyns.pos(:,2)*10,curSyns.pos(:,3)*10, ...
    %         25,synCol,'o','filled');
    %     title('Average ratio of topological to euclidean distance of synapse to all other synapses?');
end

figure(f1);
xline(0.45);
xlim([0.1,0.8]);
yticks([1:length(rgcSubTypeList)+3]+.2);
yticklabels({'ON BPC in','OFF BPC in','AMC in','4i','4ow','51','63','37'});
xlabel('IPL depth');
ylabel('target type');
title('IPL depth of Synapses');

%%
%figure out what is going on with the ON-BPC inputs in the OFF plexus
atoneLocs=[];
atoneDepths=[];
atoneEdges=[];
for i=1:length(allSkels)
    curCid=vgcCidList(i);
    curSM=allSkels{i}.sm;
    curSyns=curSM.syn;
    [n,m,curSyns.depth]=getIPLdepth(curSyns.pos(:,3),curSyns.pos(:,1),curSyns.pos(:,2),[],[]);
    inputs=curSyns.post==curCid;
    outputs=curSyns.pre==curCid;
    preTypeDat=cid2type(curSyns.pre,curTis);
    postTypeDat=cid2type(curSyns.post,curTis);
    %super sketch, but I'm just saying that input&~bpc = amc
    bpcSynIDs=find(preTypeDat{1}==7 & inputs'==1);
    amcSynIDs=find(inputs'==1 & preTypeDat{1}~=7);
    onBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcONsubs));
    offBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcOFFsubs));
    
    ONlayerSynIDs=find(curSyns.depth>0.45);
    OFFlayerSynIDs=find(curSyns.depth<0.45);
    curList=intersect(onBpcSynIDs,OFFlayerSynIDs);
    curListOFF=intersect(offBpcSynIDs,ONlayerSynIDs);
    atoneLocs=vertcat(atoneLocs,curSyns.pos(curListOFF,:));
    atoneDepths=vertcat(atoneDepths,curSyns.depth(curListOFF));
    atoneEdges=vertcat(atoneEdges,curSyns.edges(curListOFF,:));
end

if 0
    for i=1:size(atoneLocs,1)
        i
        atoneEdges(i,1:2)
        atoneDepths(i)
        curLoc=atoneLocs(i,:).*[250 250 25];
        clipboard('copy',curLoc([2 1 3]));
        pause();
    end
end



%% extra
if 0
    
    for i=1:length(vgcCidList)
        curSM=allSkels{i}.sm;
        synDistMat=curSM.syn2Skel.syn2SynDist;
        skelDistMat=curSM.syn2Skel.syn2SkelDist;
        
        
        
        
    end
    
end