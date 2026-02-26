bpc2offaIDs=find(allPreTypes{1}==7 & allPostTypes{1}==1 & allPostTypes{3}==23);
preOffaBpcSubtypes=allPreTypes{3}(bpc2offaIDs);
preOffaKnownSubtypes=find(preOffaBpcSubtypes>0);
preOffaBpcSubtypeNames=curTis.cells.type.subTypeNames{7}(preOffaBpcSubtypes(preOffaKnownSubtypes));

%% Make a graph of the different types of BPCs and what kinds of RGCs they connect to

vgcAllCids=[2 3 4 5 6 10 11 13 14 20];

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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same thing, but for the influence
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
v2rInds=getClassConn('amc','vgc','rgc','all',curTis);
v2rIDs=curTis.syn.edges(v2rInds,3);

infDat=cell(length(curTis.cells.type.subTypeNames{7}),length(curTis.cells.type.subTypeNames{1}));

for i=1:length(v2rIDs)
    curID=v2rIDs(i);
    curInd=v2rInds(i);
    curPreVG=curTis.syn.edges(curInd,2);
    curPostSub=allPostTypes{3}(curInd);
    if curPostSub==0
        curPostSub=1;
    end
    if ~isempty(curPreVG) & ismember(curPreVG,vgcCidList)
        curSM=allSkels{vgcCidList==curPreVG}.sm;
        curSMpreTypes=cid2type(curSM.syn.edges(:,2),curTis);
        curSMpostTypes=cid2type(curSM.syn.edges(:,1),curTis);
        bpcPresmIDs=find(curSMpreTypes{1}==7);
        smInd=find(curSM.syn.edges(:,3)==curID);
        curEdges=curSM.syn.edges(smInd,:);
        if ~isempty(smInd)
            curDist=curSM.syn2Skel.syn2SynDist(:,smInd);
            curDist(smInd)=200;
            infMat=exp(-curDist./lc);
            for bIt=1:length(curTis.cells.type.subTypeNames{7})
                curF=infMat(curSMpreTypes{1}==7&curSMpreTypes{3}==bIt);
                curDat=infDat{bIt,curPostSub};
                curDat=[curDat sum(curF)];
                infDat{bIt,curPostSub}=curDat;
            end
        end
    end
end

totInf=zeros(length(curTis.cells.type.subTypeNames{7}),length(curTis.cells.type.subTypeNames{1}));
for i=1:size(totInf,1)
    for j=1:size(totInf,2)
        if length(infDat{i,j})>0
        totInf(i,j)=sum(infDat{i,j})./length(infDat{i,j});
        end
            
    end
end
horSum=sum(totInf,2);
verSum=sum(totInf,1);
grows=find(horSum>0);
gcols=find(verSum>0);
bnames=curTis.cells.type.subTypeNames{7}(grows);
rnames=curTis.cells.type.subTypeNames{1}(gcols);
shoInf=totInf(grows,gcols);
%% %% fig
f10=figure();
hold on
displayResults=shoInf;
%displayResults(displayResults>205)=205;
%synFactor=floor(500/max(results(:)));
displayResults2=displayResults*80;
image(displayResults2);
cmp=colormap(turbo);
%cmp=cmp([1:205],:);
cb=colorbar();
%cb.Limits=[0 205];
%cb.Ticks=[0:5*synFactor:200];
%colorbarLabs=[0:5:200/synFactor];
%finalLabs=num2cell(colorbarLabs);
%cb.TickLabels=finalLabs;
%cb.Label.String='# synapses';
xticks([1:length(rnames)]);
yticks([1:length(bnames)]);
%xtickNames=curTis.cells.type.subTypeNames{1}(rgcTypes(rgcTypes>0));
%xtickNames={'unk' xtickNames{keepCols(2:end)-1}};
%ytickNames=curTis.cells.type.subTypeNames{7}(bpcTypes(bpcTypes>0));
%ytickNames={'unk' ytickNames{keepRows(2:end)-1}};
xticklabels(rnames);
yticklabels(bnames);
%ax0.XTick=[1:length(rgcTypes)];
%ax0.XTickLabel=xtickNames;
%% Pie Time!
pieFig=figure();
hold on
tlo=tiledlayout('flow');
for m=1:length(rnames)
    if length(infDat{1,gcols(m)})>0
    curPieDat=zeros(size(infDat,1),length(infDat{1,gcols(m)}));
    for n=1:length(bnames)
        curPieDat(n,:)=infDat{n,gcols(m)};
    end
    pieSum=sum(curPieDat,2);
    goodSlices=find(pieSum>0);
    sliceNames=curTis.cells.type.subTypeNames{7}(goodSlices);
    sliceDat=pieSum(goodSlices);
    curPie=pie(sliceDat,sliceNames);
    title([rnames{m} ' n=' num2str(length(infDat{1,m}))]);
    nexttile
    end
end


%%
f1=figure();
hold on
plotDat=sum(displayResults2,1);
b1=bar(plotDat);

f2=figure();
hold on
plotDat=sum(displayResults2,2);
b1=barh(plotDat);









%% Shortest path from SHARED bipolar input to rgc output
% THIS IS USELESS
while 0
%get a list of shared ribbons
vpost=curTis.syn.edges(ismember(curTis.syn.edges(:,1),vgcAllCids),3);
rpost=curTis.syn.edges(find(allPostTypes{1}==1),3);
sharedRibIDs=intersect(vpost,rpost);

nearestRGCsynType=zeros(length(sharedRibIDs),1);
nearestRGCsynDist=zeros(length(sharedRibIDs),1);
resultTotal={};
ribParts={};
missing=0;
missingRibVGC=[];
missingRibID=[0 0 0 0];
for ribIt=1:length(sharedRibIDs)
    results=[];
    curRibID=sharedRibIDs(ribIt);
    curEdges=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,:);
    vgcTarg=intersect(curEdges(:,1),vgcCidList);
    otherTarg=curEdges(~ismember(curEdges(:,1),vgcAllCids),1);
    ribParts{ribIt}=otherTarg;
    if ~isempty(vgcTarg)
        for vIt=1%:length(vgcTarg)
            curSM=allSkels{vgcCidList==vgcTarg(vIt)}.sm;
            smInd=find(curSM.syn.edges(:,3)==curRibID);
            if isempty(smInd) & length(vgcTarg)>1
                curSM=allSkels{vgcCidList==vgcTarg(2)}.sm;
                smInd=find(curSM.syn.edges(:,3)==curRibID);
            end
            if ~isempty(smInd)
                results=getLocalHood(curSM,smInd,30);
            else
                missing=missing+1
                missingRibVGC=[missingRibVGC vgcTarg(1)];
                missingRibID=[missingRibID; curEdges];
            end
            curDat=results;
            if ~isempty(curDat)
                postTypeDat=cid2type(curDat(:,3),curTis);
                preTypeDat=cid2type(curDat(:,4),curTis);
                vg_rgc_IDs=find(postTypeDat{1}==1);
                if ~isempty(vg_rgc_IDs)
                    nearestRGCout=find(postTypeDat{1}==1,1);
                    targType=postTypeDat{3}(nearestRGCout);
                    targDist=curDat(nearestRGCout,2);
                    %nearestRGCsyn=targType;
                end
            end
            
            if ~isempty(results)
                resultTotal{ribIt}=results;
                nearestRGCsynType(ribIt)=targType;
                nearestRGCsynDist(ribIt)=targDist;
            end
            
        end
    end
end
%% show the results
figure();
hold on
histbins=[-0.5:1:length(curTis.cells.type.subTypeNames{1})+0.5];

nameList=curTis.cells.type.subTypeNames{1};
nameList={'none' nameList{:}};
curPieDat=histcounts(nearestRGCsynType,histbins);
displayResults=curPieDat;
plotIDs=find(curPieDat>0);
displayResults=displayResults(plotIDs);
keepCols=[1 3:14 16];
sumCols=[1 2 15];
displayResults(1)=sum(displayResults(sumCols),2);
bar(displayResults(keepCols));
names=nameList(plotIDs);
names=names(keepCols);
xticks([1:length(plotIDs)]);
xticklabels(names);

title('RGC target subtype for synapse closest to shared RGC/VG3 ribbon')

end
%%
%% What RGC outputs are within X microns of a shared ribbon?
%get a list of shared ribbons
% THIS IS ALSO A DUMB IDEA, PROBABLY
while 0
hoodSize=20;

vpost=curTis.syn.edges(ismember(curTis.syn.edges(:,1),vgcAllCids),3);
rpost=curTis.syn.edges(find(allPostTypes{1}==1),3);
sharedRibIDs=intersect(vpost,rpost);

nearestRGCsynType={}; %zeros(length(sharedRibIDs),1);
nearestRGCsynDist={}; %zeros(length(sharedRibIDs),1);
resultTotal={};
ribParts={};
missing=0;
missingRibVGC=[];
missingRibID=[0 0 0 0];
for ribIt=1:length(sharedRibIDs)
    results=[];
    curRibID=sharedRibIDs(ribIt);
    curEdges=curTis.syn.edges(curTis.syn.edges(:,3)==curRibID,:);
    vgcTarg=intersect(curEdges(:,1),vgcCidList);
    otherTarg=curEdges(~ismember(curEdges(:,1),vgcAllCids),1);
    ribParts{ribIt}=otherTarg;
    if ~isempty(vgcTarg)
        for vIt=1%:length(vgcTarg)
            curSM=allSkels{vgcCidList==vgcTarg(vIt)}.sm;
            smInd=find(curSM.syn.edges(:,3)==curRibID);
            if isempty(smInd) & length(vgcTarg)>1
                curSM=allSkels{vgcCidList==vgcTarg(2)}.sm;
                smInd=find(curSM.syn.edges(:,3)==curRibID);
            end
            if ~isempty(smInd)
                results=getLocalHood(curSM,smInd,hoodSize);
            else
                missing=missing+1
                missingRibVGC=[missingRibVGC vgcTarg(1)];
                missingRibID=[missingRibID; curEdges];
            end
            curDat=results;
            if ~isempty(curDat)
                postTypeDat=cid2type(curDat(:,3),curTis);
                preTypeDat=cid2type(curDat(:,4),curTis);
                vg_rgc_IDs=find(postTypeDat{1}==1);
                if ~isempty(vg_rgc_IDs)
                    nearestRGCout=find(postTypeDat{1}==1,1);
                    targType=postTypeDat{3}(vg_rgc_IDs);
                    targDist=curDat(vg_rgc_IDs,2);
                    %nearestRGCsyn=targType;
                end
            end
            
            if ~isempty(results)
                resultTotal{ribIt}=results;
                nearestRGCsynType{ribIt}=targType;
                nearestRGCsynDist{ribIt}=targDist;
            end
            
        end
    end
end
allNearRGCtypes=[];
for m=1:length(nearestRGCsynType)
    allNearRGCtypes=[allNearRGCtypes nearestRGCsynType{m}];
end

figure();
hold on
histbins=[-0.5:1:length(curTis.cells.type.subTypeNames{1})+0.5];

nameList=curTis.cells.type.subTypeNames{1};
nameList={'none' nameList{:}};
curPieDat=histcounts(allNearRGCtypes,histbins);
displayResults=curPieDat;
plotIDs=find(curPieDat>0);
names=nameList(plotIDs);
displayResults=displayResults(plotIDs);
keepCols=[1:length(plotIDs)];
sumCols=contains(names,{' ','unk'});
keepCols=keepCols(~sumCols);
names=names(keepCols);
displayResults(1)=sum(displayResults([1;find(sumCols(:))]),2);
bar(displayResults(keepCols));
%names=names(keepCols);
xticks([1:length(plotIDs)]);
xticklabels(names);
title('target RGC types within 20um of a shared ribbon input');
end
%%
%Need to figure out how many bpc inputs there are to each RGC for the paper
test=getClassConn('bpc','all','rgc','all',curTis);
b2rCids=curTis.syn.edges(test,1);
b2rTypes=cid2type(b2rCids,curTis);
rgCidList = unique(b2rCids);
out = [rgCidList,histc(b2rCids(:),rgCidList)];
[a,b] = sort(out(:,2),'descend');
output=horzcat(rgCidList(b),a);
test2=cid2type(output(1:10,1),curTis);
%%

f3=figure();
hold on
TL1=tiledlayout('flow');
histbins=[-0.5:1:length(curTis.cells.type.subTypeNames{1})+0.5];
nameList=curTis.cells.type.subTypeNames{1};
nameList={'none' nameList{:}};
for i=1:length(nearestRGCsynType)
    curDat=nearestRGCsynType{i};
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
v2rTargType=cid2type(curTis.syn.edges(find(allPostTypes{1}'==1&ismember(curTis.syn.edges(:,2),allVG3cids)),1),curTis);
histogram(v2rTargType{3},[-0.5:1:59.5]);


%%
