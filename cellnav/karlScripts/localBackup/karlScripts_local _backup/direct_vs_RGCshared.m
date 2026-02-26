%% Make a graph of the different types of BPCs and what kinds of RGCs they connect to
allPreTypes=cid2type(curTis.syn.edges(:,2),curTis);
allPostTypes=cid2type(curTis.syn.edges(:,1),curTis);
bpc2rgc=getClassConn('bpc','all','rgc','all',curTis);
bpcTypes=unique(allPreTypes{3}(bpc2rgc));
rgcTypes=unique(allPostTypes{3}(bpc2rgc));
results=zeros(length(bpcTypes),length(rgcTypes));
for i=1:length(bpc2rgc)Sub))= ...
        results(find(bpcTypes==curPreSub),find(rgcTypes==curPostSub))+1;
end

%%
yOrder=[1 2 3 4 5 11 6 7 8 9 10 13 12 14];
xOrder=[1 2 4 5 6 7 8 9 10 3 12 11 13 14];

    curID=bpc2rgc(i);
    curPreSub=allPreTypes{3}(curID);
    curPostSub=allPostTypes{3}(curID);
    results(find(bpcTypes==curPreSub),find(rgcTypes==curPost
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
% get the orders of the things right.
displayResults2=displayResults2(yOrder,xOrder);

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
xticklabels(xtickNames(xOrder));
yticklabels(ytickNames(yOrder));
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
yOrder=[1 2 3 8 4 5 6 7 9 10 12];
xOrder=[1 3 5 6 7 8 9 4 10 11 12 14 13 15 17];
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
    nexttile
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
    title([rnames{m} ' n=' num2str(length(infDat{1,gcols(m)}))]);
    
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
%get a list of shared ribbons
vpost=curTis.syn.edges(ismember(curTis.syn.edges(:,1),vgcAllCids),3);
rpost=curTis.syn.edges(find(allPostTypes{1}==1),3);
bpre=curTis.syn.edges(find(allPreTypes{1}==7),3);
apost=curTis.syn.edges(find(allPostTypes{1}==8 & ~ismember(curTis.syn.edges(:,2)',vgcAllCids) & ~ismember(curTis.syn.edges(:,1)',vgcAllCids)),3);
zeroPost=curTis.syn.edges(find(curTis.syn.edges(:,1)'==0 & ~ismember(curTis.syn.edges(:,2)',vgcAllCids)),3);
sharedRibIDs=intersect(vpost,rpost);
amcPartRibIDs=intersect(vpost,[apost;zeroPost]);
nearestRGCsynType=zeros(length(sharedRibIDs),1);
nearestRGCsynDist=zeros(length(sharedRibIDs),1);
nearestRGCsynCid=zeros(length(sharedRibIDs),1);
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
    ribPartTypes=cid2type(otherTarg,curTis);
    rgcPartSubtype=ribPartTypes{3}(ribPartTypes{1}==1);
    rgcPartCid=otherTarg(ribPartTypes{1}==1);
    if ~isempty(vgcTarg)
        for vIt=1%:length(vgcTarg)
            curSM=allSkels{vgcCidList==vgcTarg(vIt)}.sm;
            smInd=find(curSM.syn.edges(:,3)==curRibID);
            if isempty(smInd) & length(vgcTarg)>1
                curSM=allSkels{vgcCidList==vgcTarg(2)}.sm;
                smInd=find(curSM.syn.edges(:,3)==curRibID);
            end
            if ~isempty(smInd)
                results=getLocalHood(curSM,smInd,50);
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
                    nearestRGCout=find(ismember(curDat(:,3),rgcPartCid),1);
                    %nearestRGCout=find(postTypeDat{1}==1 & ismember(postTypeDat{3},rgcPartSubtype),1);
                    if ~isempty(nearestRGCout)
                        targCid=curDat(nearestRGCout,3);
                        targType=postTypeDat{3}(nearestRGCout);
                        targDist=curDat(nearestRGCout,2);
                        %nearestRGCsyn=targType;
                    end
                end
            end
            
            if ~isempty(results) & ~isempty(nearestRGCout)
                resultTotal{ribIt}=results;
                nearestRGCsynCid(ribIt)=targCid;
                nearestRGCsynType(ribIt)=targType;
                nearestRGCsynDist(ribIt)=targDist;
            end
            
        end
    end
end
%% Find how many of the RGCs sharing ribbons are innervated by VG3
sharedRibInds=find(ismember(curTis.syn.edges(:,3),sharedRibIDs))    ;
postVgcRGCcids=curTis.syn.edges(ismember(curTis.syn.edges(:,2),vgcAllCids)&allPostTypes{1}'==1,1);

ribPartRGCcids=curTis.syn.edges(sharedRibInds,1);
ribPartRGCcids=ribPartRGCcids(~ismember(ribPartRGCcids,vgcAllCids));

ribPartRGCknown=ismember(ribPartRGCcids,postVgcRGCcids);
'percentage of ribbon partner RGCs which have a VG3 input somewhere'
mean(ribPartRGCknown)

ribPartAMCinds=find(ismember(curTis.syn.edges(:,3),amcPartRibIDs));
ribPartAMCcids=curTis.syn.edges(ribPartAMCinds,1);
ribPartAMCcids=ribPartAMCcids(~ismember(ribPartAMCcids,vgcAllCids));

pieDat=[length(ribPartAMCcids),length(ribPartRGCcids(ribPartRGCknown)),length(ribPartRGCcids(~ribPartRGCknown))];
sliceLabs={'350','128','37'};
f11=figure();
pie1=pie(pieDat,sliceLabs);
patchHand = findobj(pie1, 'Type', 'Patch'); 
patchHand(1).FaceColor = [0.8 0.8 0.8];
patchHand(3).FaceColor = [.8 0 .8];
patchHand(2).FaceColor = [0 .8 .8];

ribIDs=curTis.syn.edges(find(ismember(curTis.syn.edges(:,1),vgcAllCids) & allPreTypes{1}'==7),3);
ribIDs=unique(ribIDs);
%ribIDcell=num2str(ribIDs);
%histDat=curTis.syn.edges(ismember(curTis.syn.edges(:,3),ribIDs),3);
ribFreq=zeros(size(ribIDs));
for i=1:length(ribIDs)
    ribFreq(i)=sum(curTis.syn.edges(:,3)==ribIDs(i));
end
ribFreq=ribFreq';
freqResults=horzcat(ribIDs,ribFreq');
%% figure out the ribPartner types and compare to the nearest types

ribPartRGCcidTrm=zeros(size(nearestRGCsynType));
ribPartSubtypes=zeros(size(nearestRGCsynType));
ffInd=zeros(size(nearestRGCsynType));
sameTypeBool=zeros(size(nearestRGCsynType));
for n=1:length(ribParts)
    curPartList=ribParts{n};
    ribPartRGCcidTrm(n)=curPartList(1);
    curPartTypeDat=cid2type(curPartList,curTis);
    curPartSubType=curPartTypeDat{3};
    if length(curPartSubType)>1
    curPartSubType(curPartSubType==0)=[];
    end
    if isempty(curPartSubType)
        curPartSubType=0;
    end
    if length(curPartSubType)>1
    curPartSubType=curPartSubType(1);
    end
    ribPartSubtypes(n)=curPartSubType;
    
    if ismember(nearestRGCsynCid(n),curPartList) & nearestRGCsynCid(n)>0
        ffInd(n)=1;
    end
end

shortDist=nearestRGCsynDist(ribPartSubtypes>0 & nearestRGCsynType>0);

partNearCompare=horzcat(ribPartSubtypes,nearestRGCsynType);
sameTypeBool=ribPartSubtypes==nearestRGCsynType;
IDdonly=sameTypeBool(ribPartSubtypes>0 & nearestRGCsynType>0);
'how often is the subtype matching?'
mean(IDdonly)
sameTypeNotCid=sameTypeBool-ffInd;
sameTypeNotCidIDd=sameTypeNotCid(ribPartSubtypes>0 & nearestRGCsynType>0);
'how often FF to other cid?'
mean(sameTypeNotCidIDd)
FFid=ffInd(ribPartSubtypes>0 & nearestRGCsynType>0);
'how often FF?'
mean(FFid)

shortDist(find(FFid))

'average dist to direct FF'
mean(shortDist(find(FFid)))
'average dist to subtype FF'
mean(shortDist(find(sameTypeNotCidIDd)))

check1=horzcat(nearestRGCsynType,nearestRGCsynCid,ribPartSubtypes,ribPartRGCcidTrm,ffInd,nearestRGCsynDist); 

%% post josh meeting

%tlo=tiledlayout('flow');

ribPartSubtypesComb=ribPartSubtypes;
ribPartSubtypesComb(ribPartSubtypesComb==0)=1;
ribPartSubtypesComb(ribPartSubtypesComb==57)=1;
%check3=horzcat(cell2mat(nearestRGCsynType'),nearestRGCsynCid,ribPartSubtypesComb',ribPartRGCcidTrm',ffInd',cell2mat(nearestRGCsynDist')); 

partTypeList=unique(ribPartSubtypesComb);

swarmDat={};
swarmDatNames={};
histBins=[0:2:30];
histDat=zeros(length(partTypeList),length(histBins)-1);
for i=1:length(partTypeList)
    curSub=partTypeList(i);
    curSwat=nearestRGCsynDist(ribPartSubtypesComb==curSub & nearestRGCsynType==curSub & curSub~=1);
    swarmDat{i}=horzcat(repmat(i,size(curSwat)),curSwat);
    swarmDatNames(i)=curTis.cells.type.subTypeNames{1}(partTypeList(i));
    %[swarmDatNames{i} ' n=' num2str(length(swarmDat{i}))]
    curHist=histcounts(curSwat,histBins);
    histDat(i,:)=curHist;
end

f11=figure();
%tlo = tiledlayout(2,1,'TileSpacing','Compact');
tlo = tiledlayout('Flow');
% tlo.TileSpacing = 'compact';
% tlo.Padding = 'compact';
%trying to make the histograms.
for i=[1:12]%1:length(partTypeList)
    nexttile(tlo);
    curBar=bar(histDat(i,:));
    ylim([0 10])
    title([swarmDatNames{i} ' n=' num2str(length(swarmDat{i}))]);
    xlabel('Distance (μm)','interpreter', 'tex');
    ylabel('# of synapses');
    yticks([2:2:8]);
    yticklabels([2:2:8]);
    xticklabels([2:2:30]);
end
sgtitle({'Distance from shared ribbon inputs','to nearest feed-forward RGC synapse'});

%%
v2rInds=getClassConn({'amc'},{'vgc'},{'rgc'},{'all'},curTis);
rgcTargTypes=cid2type(curTis.syn.edges(v2rInds,1),curTis);
rgcTargSubs=rgcTargTypes{3};
targTab=tabulate(rgcTargSubs);
targTab(find(targTab(:,1)==1),2)=targTab(find(targTab(:,1)==0),2) ...
    + targTab(find(targTab(:,1)==1),2) + targTab(find(targTab(:,1)==57),2);
targTab(find(targTab(:,1)==0),2)=0;
targTab(find(targTab(:,1)==57),2)=0;
[sDat, sOrd]=sort(targTab(:,2),'descend');
histGraphDat=targTab(sOrd,2);
graphSubIndsPlus=targTab(sOrd,1)+1;
nams={'none' curTis.cells.type.subTypeNames{1}{:}};
nams{2}='unknown';
%nams='unknown';
histGraphNames=nams(targTab(sOrd,1)+1);
histCutOff=max(find(histGraphDat>4))+0.5;
targTypeFig=figure();
hold on
bar(histGraphDat);
xticks([1:length(histGraphNames)])
xticklabels(histGraphNames);
xlim([0.5 histCutOff]);
ylabel('# of synapses');
xlabel('RGC subtype');
title('Number of VG3-RGC synapses by target subtype');


%% Last one. the midrange. 
% Same individuals, but not limited to the same ribbon
for i=1:length(b2vIDs)
    curID=b2vIDs(i);
    curRibTargs=curTis.syn.edges(curTis.syn.edges(:,3)==curID,1);
    
    
end

%%
figure();
hold on
partTable=tabulate(ribPartSubtypesComb');
partTable2=targTab;
partTable2(1,2)=0;
for j=2:size(partTable2,1)
    if ~isempty(find(partTable(:,1)==partTable2(j,1)))
    partTable2(j,2)=partTable(find(partTable(:,1)==partTable2(j,1)),2);
    else
        partTable2(j,2)=0;
    end
end
histGraphDat=partTable2(sOrd,2);
bar(histGraphDat);
xticks([1:length(histGraphNames)])
xticklabels(histGraphNames);
xlim([0.5 histCutOff]);
ylabel('# of synapses');
xlabel('RGC subtype');
title('Number VGC-RGC shared ribbons by partner subtype');

%%
figure();
hold on
swit=1;
s1={};
for i=1:length(swarmDat)
    curDat=swarmDat{i};
    if length(curDat)>2
        swnames{swit}=swarmDatNames{i};
        s1{i}=swarmchart(repmat(swit,length(curDat(:,2)),1),curDat(:,2),15,HCcolmap(swit,:),'filled');
        %s1{i}.
        s1{i}.YJitter='none';
        s1{i}.XJitter='density';
        s1{i}.XJitterWidth=0.5;
        swit=swit+1;
    end
end

%yticks([1 2]);
%yticklabels({'Feed-Forward','non-Feed-Forward'})
%xlim([0 20])
xticks([1:swit-1]);
xticklabels(swnames)
ylabel('Distance (μm)','interpreter', 'tex');
title('Distance from shared ribbon inputs to nearest feed-forward RGC synapse');

%%
ffBool=(ffInd|((ribPartSubtypes==nearestRGCsynType)&ribPartSubtypes>0));
nffBool=(~ffInd&(ribPartSubtypes~=nearestRGCsynType));
check2=horzcat(ffBool,nffBool);
mean(sum(check2,2));

%ffIDs=find(ffInd&nearestRGCsynType>0);
%nffIDs=find(~ffInd&ribPartSubtypes~=nearestRGCsynType);

%ffDists=nearestRGCsynDist(find(ffInd&nearestRGCsynType>0));
%nffDists=nearestRGCsynDist(find(~ffInd&ribPartSubtypes~=nearestRGCsynType));
