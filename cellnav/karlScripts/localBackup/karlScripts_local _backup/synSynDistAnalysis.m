%% Find out if the ON and OFF bipolar inputs are as
% close on the arbor as they are in space
if ~exist('allSkels')
load('c:\work\allSkels.mat');
end

f1=figure();
tlo1=tiledlayout(1,1);
colors=colormap(lines);
cmod=1;

hoodDat={};
%f2=figure();
thicc=[0 0 0 0 0];
if 0
    [q,w,testDeps1]=getIPLdepth([20,30],[mean(pos(:,1)) mean(pos(:,1))], [mean(pos(:,2)) mean(pos(:,2))],[],[]);
    thicc(1)=testDeps1(1)-testDeps1(2);
    prtl=[10 10];
    [q,w,testDeps2]=getIPLdepth([20,30],[prctile(pos(:,1),prtl(1)) prctile(pos(:,1),prtl(1))], [prctile(pos(:,1),prtl(2)) prctile(pos(:,1),prtl(2))],[],[]);
    thicc(2)=testDeps2(1)-testDeps2(2);
    prtl=[10 90];
    [q,w,testDeps3]=getIPLdepth([20,30],[prctile(pos(:,1),prtl(1)) prctile(pos(:,1),prtl(1))], [prctile(pos(:,1),prtl(2)) prctile(pos(:,1),prtl(2))],[],[]);
    thicc(3)=testDeps3(1)-testDeps3(2);
    prtl=[90 90];
    [q,w,testDeps4]=getIPLdepth([20,30],[prctile(pos(:,1),prtl(1)) prctile(pos(:,1),prtl(1))], [prctile(pos(:,1),prtl(2)) prctile(pos(:,1),prtl(2))],[],[]);
    thicc(4)=testDeps4(1)-testDeps4(2);
    prtl=[90 10];
    [q,w,testDeps5]=getIPLdepth([20,30],[prctile(pos(:,1),prtl(1)) prctile(pos(:,1),prtl(1))], [prctile(pos(:,1),prtl(2)) prctile(pos(:,1),prtl(2))],[],[]);
    thicc(5)=testDeps5(1)-testDeps5(2);
end

for i=1:length(allSkels)
    %% load skels, etc.
    curSM=allSkels{i}.sm;
    curCid=curSM.cid;
    %get the types of all the synapse partners
    preTypes=cid2type(curSM.syn.edges(:,2),curTis);
    postTypes=cid2type(curSM.syn.edges(:,1),curTis);
    %get the fv library object
    %curCidFv=load([fvDir num2str(curCid) '.mat']);
    curCidFv=curSM.nep.fv;
    % get depth info for nodes/syns/etc
    [inl,gcl,synDepths]=getIPLdepth(curSM.syn.pos(:,3),curSM.syn.pos(:,1),curSM.syn.pos(:,2),[],[]);
    [inl,gcl,nepDepths]=getIPLdepth(curSM.nep.pos(:,3),curSM.nep.pos(:,1),curSM.nep.pos(:,2),[],[]);
    
    %find bpc inputs and get the partner cids
    bpcInIDs2=find(curSM.syn.preClass==7);
    bpcInIDs=find(preTypes{1}'==7);
    %find the IDs of those synapses
    
    bpcInObIDs=curSM.syn.obID(bpcInIDs);
    bpcInsynIDs=curSM.syn.synID(bpcInIDs);
    
    offTypes=[1:5 15 18 19];
    onTypes=[6:12 14 17 20 21 24];
    
    bpcInIDs=find(preTypes{1}'==7);
    NOTbpcInIDs=find(preTypes{1}'~=7);
    bpcOFFInIDs=find(preTypes{1}'==7 & ismember(preTypes{3}',offTypes));
    bpcONInIDs=find(preTypes{1}'==7 & ismember(preTypes{3}',onTypes));
    
    
    
    syn2synLin=curSM.syn2Skel.syn2SynDist;
    for r=1:size(syn2synLin,1)
        syn2synLin(r,r)=inf;
    end
    
    
    acrossSynDists=curSM.syn2Skel.syn2SynDist(bpcOFFInIDs,bpcONInIDs);
    OFFSynDists=curSM.syn2Skel.syn2SynDist(bpcOFFInIDs,bpcOFFInIDs);
    ONSynDists=curSM.syn2Skel.syn2SynDist(bpcONInIDs,bpcONInIDs);
    minXdist=zeros(1,size(syn2synLin,1));
    minXdist(bpcOFFInIDs)=min(acrossSynDists,[],2);
    minXdist(bpcONInIDs)=min(acrossSynDists,[],1);
    
    if 0
    figure(); hold on;
    onoffDists=histcounts(acrossSynDists,[0:5:150]);
    offDists=histcounts(OFFSynDists,[0:5:150]);
    onDists=histcounts(ONSynDists,[0:5:150]);
    plot([1:length(onoffDists)],smooth(onoffDists));
    plot([1:length(offDists)],smooth(offDists));
    plot([1:length(onDists)],smooth(onDists));
    legend({'ON-OFF','OFF-OFF','ON-ON'})
    end
    %get the euc dists for all the synapses
    pos = curSM.syn.pos;
    [gcl,inl,synDepths]=getIPLdepth(pos(:,3),pos(:,1),pos(:,2),[],[]);
    off2onMinDists=min(acrossSynDists,[],2);
    offDepths=synDepths(bpcOFFInIDs);
    %figure(); hold on
    figure(f1);
    nexttile(1)
    hold on
    scatter(offDepths,off2onMinDists,25,colors(i*cmod,:),'x');
    xlim([0.2 0.75]);
    ylim([0 100]);
    on2offMinDists=min(acrossSynDists,[],1);
    onDepths=synDepths(bpcONInIDs);
    nexttile(1)
    hold on
    scatter(onDepths,on2offMinDists,25,colors(i*cmod,:),'o');
    xlim([0.2 0.75]);
    ylim([0 100]);
    % This is showing that some are significantly farther in euclidean
    % space than in linear space. no bueno
    dif1 = pos(:,1)-pos(:,1)';
    dif2 = pos(:,2)-pos(:,2)';
    dif3 = pos(:,3)-pos(:,3)';
    syn2synEuc = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
    if 0
    figure(); hold on;
    scatter(syn2synEuc(:),syn2synLin(:));
    xlim([0 140])
    ylim([0 140])
    end
    badDist=find(syn2synLin(:)<syn2synEuc(:));
    difDist=syn2synLin(:)-syn2synEuc(:);
    linD=syn2synLin(badDist);
    [badx,bady]=ind2sub(size(syn2synEuc),badDist);
    
    % how polar is the neighborhood of RGC outputs?
    routInds=find(postTypes{1}'==1);
    routSubs=postTypes{3}(routInds)';
    unSubs=unique(routSubs);
    for j=1:length(unSubs)
        curSub=unSubs(j);
        if curSub>0
            %get the sminds of the rgc outputs for this subtype
            curRoutSynInds=routInds(find(routSubs==curSub));
            %find the closest bpc input
            curSynDepths=synDepths(curRoutSynInds);
            curDistMat=syn2synLin(curRoutSynInds,:);
            curDistMat(:,NOTbpcInIDs)=inf;
            [minVals,minInds]=mink(curDistMat,3,2);
            minSubXdist=minXdist(minInds);
            meanSubXdist=mean(minSubXdist,2);
            hoodDat{i,curSub,1}=meanSubXdist;
            hoodDat{i,curSub,2}=curSynDepths;
        end
    end
        
    %figure(f2)
    
    
    
    
end
IPLmid=0.45;
xline(IPLmid);
plot([IPLmid IPLmid+.2],[0 20],'k--')
plot([IPLmid-0.2 IPLmid],[20 0],'k--')
legend({'off','on'});
ylabel('distance from opposing input');
xlabel('IPL depth');
ylim([0 60])

sgtitle('Distance to nearest opposite polarity BPC input')
    

%% same plot, but split by RGC type, not the VG3 cid
figure();
hold on;
for i=1:size(hoodDat,2)


    
end
%%
cumSubMinDist=cell(1,length(curTis.cells.type.subTypeNames{1}));
for i=1:length(allSkels)
    for j=1:length(curTis.cells.type.subTypeNames{1})
        curSubDat=hoodDat{i,j};
        cumSubMinDist{j}=[cumSubMinDist{j}; curSubDat];
    end
end

figure();
hold on
plotNames={};
plotIt=1;
for j=1:length(curTis.cells.type.subTypeNames{1})
    curPlotDat=cumSubMinDist{j};
    if length(curPlotDat)>5
        plotNames={plotNames{:} curTis.cells.type.subTypeNames{1}{j}};
        sw{plotIt}=swarmchart(repmat(plotIt,length(curPlotDat),1),curPlotDat);
        sw{plotIt}.XJitterWidth=0.6
        plotIt=plotIt+1;
    end
end
xticks([1:plotIt]);
xticklabels(plotNames);
ylabel('distance (um)');
xlabel('rgc subtype');
title('Mean distance to BPC input of opposite polarity from 3 nearest BPC inputs to RGC outputs');
    
    
    