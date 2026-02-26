%%synfluence creator
%% params to switch often
loadSkels=1;
debug=1;
lotsaplots=1;
allGraphs=0;

%% model params, taken from previous script
resistivity = 2000;
Rm = 2000;
radius = .01;
Vo = 1;
a=radius/10000;
%x=d/10000;
ra=resistivity/(pi*a^2);
rm=Rm/(2*pi*a);
lc=sqrt(rm/ra);

%% set up locations
baseDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\';
WPN=[baseDir 'Analysis\'];
MPN=[baseDir 'Merge\'];
SMDir=[baseDir 'Analysis\SMs\'];
load([MPN 'obI.mat']);
load([MPN 'dsObj.mat']);
load([WPN 'tis.mat']);

%% load skeletons if necessary
if loadSkels==1
vgcCidList=[2 3 4 5 13];
SMs={};
uniqueVGCs=unique(vgcCidList);
for curSMIt=1:length(uniqueVGCs)
    fileStr=['sm_cid' + string(uniqueVGCs(curSMIt)) + '.mat'];
    if isfile([SMDir + fileStr])
        load([SMDir + fileStr]);
        SMs{curSMIt}=sm;
    end
end
end

%% get the skeleton data into the supermatrix
allEdges=[];
allPos=[];
allSynType=[];
axisObIDList=[];

for curVGC=1:length(SMs)
    curSM=SMs{curVGC};
    distMats{curVGC}=curSM.syn2Skel.syn2SynDist;
    axisObIDList=[axisObIDList;curSM.syn.obID];
    allEdges=[allEdges; curSM.syn.edges];
    allPos=[allPos; curSM.syn.pos];
    allSynType=[allSynType; curSM.syn.synType];
    %allSynType
end

allPreTypes=cid2type(allEdges(:,2),tis);
allPostTypes=cid2type(allEdges(:,1),tis);

bigMat=blkdiag(distMats{1},distMats{2},distMats{3},distMats{4},distMats{5});
infMat=Vo*exp(-bigMat/lc/10000);
infMat(bigMat==0)=0;
typeMat=allSynType*allSynType';

%% diagPlots
if lotsaplots==1
    figure();
    colormap(turbo);
    image(bigMat);
    title('Distances between synapses');
    
    figure();
    colormap(turbo);
    image(infMat.*255+typeMat.*2);
    title('Synapse Influence');
    
    figure();
    histogram(bigMat(bigMat>0));
    title('Distance between synapses histo');
    
    figure();
    colormap(turbo);
    image((typeMat-1)*250);
    title('Synapse Type');
end

%% cool plots
%bipolar influence on rgc output synapses.

bpc2rgcMat=infMat(allPreTypes{1}==7,allPostTypes{1}==1);
figure();
colormap(turbo);
image(bpc2rgcMat.*255);
title('Synapse BPC -> RGC');

%% dig out the data for a particular cell
rgcCidList=type2cid({'rgc', 'rgc', 'rgc', 'rgc', 'rgc'},{'4ow','37','51','63','85'},tis);
bpcList=type2cid({'bpc'},{'all'},tis);
bpcListTypes = cid2type(bpcList{1},tis);
%get an ordered bpcList for making pretty grphs
bpcOrderedSynIDs=[];
bpcOrderedTicks=[0];
bpcOrderedTickLabels={'0'};
bpcSubTypeSynIDLists={};
k=2;
for j=1:length(tis.cells.type.subTypeNames{7})
        curSubStr=tis.cells.type.subTypeNames{7}(j);
        curBPCSubtypeList=type2cid({'bpc'},curSubStr,tis);
        BPCSubtypeIDs=find(ismember(allEdges(:,2),curBPCSubtypeList{:}));
        bpcSubTypeSynIDLists{j}=BPCSubtypeIDs;
        tot=bpcOrderedTicks(end);
        if length(BPCSubtypeIDs)>0
        [ord,idx]=sort(allEdges(BPCSubtypeIDs,2));
        bpcOrderedSynIDs=[bpcOrderedSynIDs;BPCSubtypeIDs(idx)];
        bpcOrderedTicks=[bpcOrderedTicks tot+length(BPCSubtypeIDs)];
        bpcOrderedTickLabels{k}=curSubStr{:};
        k=k+1;
        end
end

%get an ordered bpcList for making pretty grphs
rgcOrderedSynIDs=[];
rgcOrderedTicks=[0];
rgcOrderedTickLabels={'0'};
rgcSubTypeSynIDLists={};
rgcSubtypeList={};
k=2;
for j=1:length(tis.cells.type.subTypeNames{1})
        curSubStr=tis.cells.type.subTypeNames{1}(j);
        curRGCSubtypeList=type2cid({'rgc'},curSubStr,tis);
        rgcSubtypeList{j} = curRGCSubtypeList;
        RGCSubtypeIDs=find(ismember(allEdges(:,1),curRGCSubtypeList{:}));
        rgcSubTypeSynIDLists{j}=RGCSubtypeIDs;
        tot=rgcOrderedTicks(end);
        if length(RGCSubtypeIDs)>0
        [ord,idx]=sort(allEdges(RGCSubtypeIDs,1));
        rgcOrderedSynIDs=[rgcOrderedSynIDs;RGCSubtypeIDs(idx)];
        rgcOrderedTicks=[rgcOrderedTicks tot+length(RGCSubtypeIDs)];
        rgcOrderedTickLabels{k}=curSubStr{:};
        k=k+1;
        end
end


%loop
for curSubIt=1:length(rgcCidList)
    subTypeCidList=rgcCidList{curSubIt};
    rgcSynIDs=find(ismember(allEdges(:,1),subTypeCidList));
    for j=1:length(tis.cells.type.subTypeNames{7})
        curSubStr=tis.cells.type.subTypeNames{7}(j);
        curBPCSubtypeList=type2cid({'bpc'},curSubStr,tis);
        BPCSubtypeIDs=find(ismember(allEdges(:,2),curBPCSubtypeList{:}));
        subGraph=infMat(rgcSynIDs,BPCSubtypeIDs);
        if allGraphs==1
            figure();
            hold on
            title(curSubStr);
            colormap(turbo);
            image(subGraph.*255);
        end
    end
end

%% Figure showing all bpc types onto RGC types.
figure();
hold on
title('BPC type influence RGC type synapses');
colormap(turbo);
subG=infMat(bpcOrderedSynIDs,rgcOrderedSynIDs);
subG2=infMat(bpcOrderedSynIDs,:);
subG2=subG2(:,rgcOrderedSynIDs);
image(subG.*255);

% clean up labels
xticks(rgcOrderedTicks(2:end));
yticks(bpcOrderedTicks(2:end));
xticklabels(rgcOrderedTickLabels(2:end));
yticklabels(bpcOrderedTickLabels(2:end));

figure();
hold on

scatter3(allPos(:,1),allPos(:,2),allPos(:,3),'b.');
scatter3(allPos(isInfluenced,1),allPos(isInfluenced,2),allPos(isInfluenced,3),'ro');

bpc1=getCidVox(1096,2,dsObj,tis);
vox=bpc1{:};
vox=vox./10;
scatter3(vox(:,1),vox(:,2),vox(:,3),'g.')

finalResults=zeros(length(bpcSubTypeSynIDLists),length(rgcSubTypeSynIDLists));
% 
% for y=1:length(bpcSubTypeSynIDLists)
%     for x=1:length(rgcSubTypeSynIDLists)
%         subGraph=infMat(bpcSubTypeSynIDLists{y},:);
%         subGraph=subGraph(:,rgcSubTypeSynIDLists{x});
%         curResult=sum(subGraph(:));
%         
%         stdSub(y,x,:) = std(subGraph,2)
%         finalResults(y,x)=curResult;
%         synNs(y,x,:) = size(subGraph);
%     end
%     
%     
% end

%%Sum For to measure each RGC
sumBipMat = [];
for y= 1:length(bpcSubTypeSynIDLists)
    subGraph=infMat(bpcSubTypeSynIDLists{y},:);
    subBipMap(y,:) = sum(subGraph,1);
end
stdSub = [];
for x=1:length(rgcSubTypeSynIDLists)
    subGraph=subBipMap(:,rgcSubTypeSynIDLists{x});
    rgcSynSubGraphs{x} = subGraph;
    curResult=sum(subGraph,2);
    stdSub(:,x) = std(subGraph')';
    finalResults(:,x)=curResult;
    
end
stdSub(isnan(stdSub))=0;


isData = find(sum(finalResults,1)>0);
numIsDat = length(isData);
% try,close isd,end
% isd = figure;

for i = 1:length(isData)
    subplot(ceil(numIsDat/3),3,i)
    
    subG = rgcSynSubGraphs{isData(i)}
    norm = max(subG(:));
    subG = subG/norm;
    
 
    cProf = (finalResults(:,isData(i)));
    cProf = cProf/norm;
    bar(cProf);
    hold on
    xPos =  repmat( [1:size(subG,1)]',[1 size(subG,2)]);
    xPos = xPos + randn(size(xPos))/10;
    scatter(xPos(:),subG(:),5,'r')
    title(tis.cells.type.subTypeNames{1}(isData(i)));
    xticklabels(tis.cells.type.subTypeNames{7})
    rgcCellNums = length(rgcSubtypeList{isData(i)}{:});
    
    rgcIns = length(rgcSubTypeSynIDLists{isData(i)});
    text(0,-.4,sprintf('%d syn, %d cells',rgcIns,rgcCellNums ),'HorizontalAlignment','left','VerticalAlignment','top')
    ylim([0 norm * 1.1])
    hold off
end
%% 
if 0
subTypeCidList=rgcCidList{1};
figure();
hold on
title('BPC type influence on VGC->OFF-A RGC synapses');
colormap(turbo);
rgcSynIDs=find(ismember(allEdges(:,1),subTypeCidList));
rgcSynIDs=find(allEdges(:,1)==1110);
subG=infMat(rgcSynIDs,bpcOrderedSynIDs);
image(subG.*255);
end



%% junk
if 0
% set up main loop
%tis.syn.obID has 3363 syns (2850 of those are unique)
%skeleton structures also have obIDs. So that should be how we track them.
allSyn=length(tis.syn.obID);
emptyMat=zeros(length(tis.syn.obID),length(tis.syn.obID));

vgcCidList=[2 3 4 5 13];

for i=1:allSyn
    curObID=tis.syn.obID(i);
    curEdges=tis.syn.edges(i,:);
    w=ismember(curEdges(1:2),vgcCidList);
    if sum(w)>0
        curVGCCid=curEdges(w);
        
        
    end
end
end