bpc2offaIDs=find(allPreTypes{1}==7 & allPostTypes{1}==1 & allPostTypes{3}==23);
preOffaBpcSubtypes=allPreTypes{3}(bpc2offaIDs);
preOffaKnownSubtypes=find(preOffaBpcSubtypes>0);
preOffaBpcSubtypeNames=curTis.cells.type.subTypeNames{7}(preOffaBpcSubtypes(preOffaKnownSubtypes));

%% Make a graph of the different types of BPCs and what kinds of RGCs they connect to
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
bpc2rgc=getClassConn('bpc','all','rgc','all',curTis);
bpcTypes=unique(allPreTypes{3}(bpc2rgc));
rgcTypes=unique(allPostTypes{3}(bpc2rgc));
results=zeros(length(bpcTypes),length(rgcTypes));
for i=1:length(bpc2rgc)
    curID=bpc2rgc(i);
    curPreSub=allPreTypes{3}(curID);
    curPostSub=allPostTypes{3}(curID);
    results(find(bpcTypes==curPreSub),find(rgcTypes==curPostSub))= ...
    results(find(bpcTypes==curPreSub),find(rgcTypes==curPostSub))+1;
end


%%
f0=figure();
hold on
%ax0=axes;
%set(ax0, 'XAxisLocation', 'top','YAxisLocation','Left');
displayResults=results*(floor(500/max(results(:))));
displayResults(displayResults>205)=205;
synFactor=floor(500/max(results(:)));
sumUnk=1;
if sumUnk
keepCols=[1 3:14 16];
keepRows=[1:13 15];
sumCols=[1 2 15];
sumRows=[1 14];
displayResults(1,:)=sum(displayResults(sumRows,:),1);
displayResults(:,1)=sum(displayResults(:,sumCols),2);
displayResults2=displayResults(keepRows,keepCols);
else
    displayResults2=displayResults;
end
image(displayResults2)
%color = get(f0,'Color');
%set(ax0,'XColor',color,'YColor',color)
cmp=colormap(turbo);
cmp=cmp([1:205],:);
cb=colorbar();
cb.Limits=[0 205];
cb.Ticks=[0:5*synFactor:200];
colorbarLabs=[0:5:200/synFactor];
finalLabs=num2cell(colorbarLabs);
cb.TickLabels=finalLabs;
cb.Label.String='# synapses';
xticks([1:length(keepCols)]);
yticks([1:length(keepRows)]);
xtickNames=curTis.cells.type.subTypeNames{1}(rgcTypes(rgcTypes>0));
xtickNames={'unk' xtickNames{keepCols(2:end)-1}};
ytickNames=curTis.cells.type.subTypeNames{7}(bpcTypes(bpcTypes>0));
ytickNames={'unk' ytickNames{keepRows(2:end)-1}};
xticklabels(xtickNames);
yticklabels(ytickNames);
%ax0.XTick=[1:length(rgcTypes)];
%ax0.XTickLabel=xtickNames;
%%
f1=figure();
hold on
plotDat=sum(displayResults2,1);
b1=bar(plotDat);


f2=figure();
hold on
plotDat=sum(displayResults2,2);
b1=barh(plotDat);

%% Shortest path from bipolar input to rgc output
resultTotal={};
nearestRGCsynTotal={}; %zeros(length(curTis.cells.type.subTypeNames{7}),1);
for bpcTypeIt=1:length(curTis.cells.type.subTypeNames{7})
    curBPC_vgInds=find(allPreTypes{1}==7 & allPreTypes{3}==bpcTypeIt & ismember(curTis.syn.edges(:,1),allVGcids)');
    curBPC_vgIDs=curTis.syn.edges(curBPC_vgInds,3);
    results={};
    nearestRGCsyn={};
    if ~isempty(curBPC_vgIDs)
        for i=1:length(curBPC_vgIDs)
            curID=curBPC_vgIDs(i);
            curInd=curBPC_vgInds(i);
            vgcTarg=curTis.syn.edges(curInd,1);
            if ismember(vgcTarg,vgcCidList)
                curSM=allSkels{vgcCidList==vgcTarg}.sm;
                smInd=find(curSM.syn.edges(:,3)==curID);
                if ~isempty(smInd)
                    results{i}=getLocalHood(curSM,smInd,30);
                end
            end
        end
    end
    
    for k=1:length(results)
        curDat=results{k};
        if ~isempty(curDat)
            postTypeDat=cid2type(curDat(:,3),curTis);
            preTypeDat=cid2type(curDat(:,4),curTis);
            vg_rgc_IDs=find(postTypeDat{1}==1);
            if ~isempty(vg_rgc_IDs)
                nearestRGCout=find(postTypeDat{1}==1,1);
                targType=postTypeDat{3}(nearestRGCout);
                nearestRGCsyn{k}=targType;
            end
        end
    end
    
    resultTotal{bpcTypeIt}=results;
    nearestRGCsynTotal{bpcTypeIt}=cell2mat(nearestRGCsyn);
end

f3=figure();
hold on
TL1=tiledlayout('flow');
histbins=[-0.5:1:length(curTis.cells.type.subTypeNames{1})+0.5];
nameList=curTis.cells.type.subTypeNames{1};
nameList={'none' nameList{:}};
for i=1:length(nearestRGCsynTotal)
    curDat=nearestRGCsynTotal{i};
    if ~isempty(curDat)
    nexttile
    title([curTis.cells.type.subTypeNames{7}(i) num2str(length(curDat))]);
    hold on
    curPieDat=histcounts(curDat,histbins);
    plotIDs=find(curPieDat>0);
    bar(curPieDat(plotIDs));
    names=nameList(plotIDs);
    xticks([1:length(plotIDs)]);
    xticklabels(names);
    %pie(curPieDat,nameList);
    end
end
%%


allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
vgcCidList=[2 3 4 5 13 14];
bpcPre=find(allPreTypes{1}==7);
amcPre=find(allPreTypes{1}==0|allPreTypes{1}==8);
vgcPost=find(ismember(curTis.syn.edges(:,1),vgcCidList)');
v2rIDs=find(ismember(curTis.syn.edges(:,1),vgcCidList)'&allPostTypes{1}==1);
bpcPreIDs=curTis.syn.edges(bpcPre,3);
amcPreIDs=curTis.syn.edges(amcPre,3);
vgcPostIDs=curTis.syn.edges(vgcPost,3);
b2vIDs=intersect(bpcPreIDs,vgcPostIDs);
a2vIDs=intersect(amcPreIDs,vgcPostIDs);
% find the RGCs that share ribbons
b2vEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),b2vIDs),:);
b2vEdges=b2vEdges(~ismember(b2vEdges(:,1),vgcCidList),:);
b2vEdgeTypes=cid2type(b2vEdges(:,1),curTis);

testList=curTis.syn.edges(find(allPreTypes{1}'==7&ismember(curTis.syn.edges(:,1),vgcCidList)),3);

targRibLocMat=struct();
% get the ribbons shared between rgcs and vgc
targType=[1 23];
if length(targType)==1
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1),3);
else
    targRibIds=b2vEdges(b2vEdgeTypes{1}==targType(1)&b2vEdgeTypes{3}==targType(2),3);
end
targRibNames=curTis.obI.nameProps.names(targRibIds)';
targRibEdges=curTis.syn.edges(ismember(curTis.syn.edges(:,3),targRibIds),:);
targRibLocs=curTis.obI.colStruc.anchors(targRibIds,:);
targRibLocMat(1).name='all rgc';
targRibLocMat(1).targIDs=targRibIds;
targRibLocMat(1).targNames=targRibNames;
targRibLocMat(1).targEdges=targRibEdges;
targRibLocMat(1).targLocs=targRibLocs;


figure(); histogram(b2vEdgeTypes{3},[-0.5:1:59.5]);
xticks([1:59])
xticklabels(curTis.cells.type.subTypeNames{1})
hold on
v2rTargType=cid2type(curTis.syn.edges(find(allPostTypes{1}'==1&ismember(curTis.syn.edges(:,2),allVGcids)),1),curTis);
histogram(v2rTargType{3},[-0.5:1:59.5]);


%%
