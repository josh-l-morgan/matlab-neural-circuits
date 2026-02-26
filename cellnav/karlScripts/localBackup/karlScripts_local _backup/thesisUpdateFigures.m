% thesis update figures
analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\';
mergeDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\';
dataDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\';
fvLib='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
load('Y:\MATLAB\cellNav\karlScripts\localBackup\karlScripts_local\aliasList.mat');
load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\obI.mat');

load([analDir 'tis.mat']);
load([mergeDir 'dsObj.mat']);
load('Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\COI.mat');
curTis=tis;

%clean the tis file
cleaning=0;
if cleaning
    curTis=findBadSynapse(curTis);
end

%newCols=HCcolmap;


aliasEdges=aliasChecker(curTis.syn.edges(:,[1 2]),[]);

preTypes=cid2type(curTis.syn.edges(:,2),curTis);
postTypes=cid2type(curTis.syn.edges(:,1),curTis);

%% Table stuff for update
vgcCidList=[2 3 4 5 10 11 13 14 20];
inputIDs=find(ismember(curTis.syn.edges(:,1),vgcCidList));
outputIDs=find(ismember(curTis.syn.edges(:,2),vgcCidList));
rgcPost=find(postTypes{1}==1);
bpcPre=find(preTypes{1}==7);
acPre=find(preTypes{1}==8|preTypes{1}==0);
acPost=find(postTypes{1}==8|postTypes{1}==0);
bpcPost=find(postTypes{1}==7);
b2v=length(intersect(inputIDs,bpcPre));
v2r=length(intersect(outputIDs,rgcPost));
a2v=length(intersect(inputIDs,acPre));
v2a=length(intersect(outputIDs,acPost));
v2b=length(intersect(outputIDs,bpcPost));
%curTis.cells.type.subTypeNames{7}
offSubTypes=[1 2 3 4 5 15 18 19];
onSubTypes=[6:12 14 17 21 24];
offPre=find(preTypes{1}==7&ismember(preTypes{3},offSubTypes));
onPre=find(preTypes{1}==7&ismember(preTypes{3},onSubTypes));
off2v=length(intersect(inputIDs,offPre));
on2v=length(intersect(inputIDs,onPre));

% amke a quick fig for showing the locations of the synapses
figure();
hold on
vgFVs={};
for o=2:3
    curCid=vgcCidList(o);
    curFVname=[num2str(curCid) '.mat'];
    curFV=load([fvLib curFVname]);
    curPatch=patch(curFV.fv);
    curPatch.EdgeAlpha=0;
    curPatch.FaceAlpha=0.3;
    curPatch.FaceColor=HCcolmap(o,:);
    
end
xlim([15 35])
view(90,0)


%preTypeTest=cid2type(aliasEdges(:,2),curTis);
%postTypeTest=cid2type(aliasEdges(:,1),curTis);

%% show a VG3 with the different synapses using the FV display
fig1=1;

figCidList=[2 3 4 5 13 14];

if fig1
    fover=figure();
    %tlo=tiledlayout('flow');
    %tlo=tiledlayout(2,3,'TileSpacing','compact') ;
    for curCid=1:length(figCidList)
        %nexttile
        %figure();
        
    figCid=figCidList(curCid);
    
    hold on
    vgc13vox=getCidVox(figCid,10,dsObj,curTis);
    vgc13vox=vgc13vox{1};
    vgc13fv=load([fvLib num2str(figCid) '.mat']);
    
    vg13patch=patch(vgc13fv.fv);
    vg13patch.EdgeColor='none';
    vg13patch.FaceColor=[.1 .1 .1];
    vg13patch.FaceAlpha=.1;
    hold on
    
    allSyns=find(curTis.syn.edges(:,1)==figCid | curTis.syn.edges(:,2)==figCid);
    inputs=find(curTis.syn.edges(:,1)==figCid);
    outputs=find(curTis.syn.edges(:,2)==figCid);
    bpcPre=find(preTypes{1}==7);
    amcPre=find(curTis.syn.edges(:,2)'==0 | preTypes{1}==8);
    rgcPost=find(postTypes{1}==1);
    amcPost=find(curTis.syn.edges(:,1)'==0 | postTypes{1}==8);
    
    bpcInputs=intersect(bpcPre,inputs);
    amcInputs=intersect(amcPre,inputs);
    rgcOutputs=intersect(rgcPost,outputs);
    amcOutputs=intersect(amcPost,outputs);
    
    synMat={bpcInputs amcInputs rgcOutputs amcOutputs};
    sizeNum=35;
    sizes=[sizeNum sizeNum sizeNum sizeNum];
    colors=[HCcolmap(8,:);HCcolmap(1,:);HCcolmap(6,:);HCcolmap(2,:)];
    symbols={'v' 'v' '^' '^'};
    for i=1:length(synMat)
        curSyns=synMat{i};
        curSynPos=curTis.syn.pos(curSyns,:);
        scatter3(curSynPos(:,3),curSynPos(:,1),curSynPos(:,2), ...
            sizes(i),colors(i,:),symbols{i},'filled');
    end
    view(-90,-90);
    set(gca,'visible','off')
    axis tight
    axis equal
    %leglabs={'vgc arbor' 'BPC in' 'AMC in' 'RGC out' 'AMC out'};
    %leggy=legend(leglabs,'Location','northeast');
    title(num2str(figCid));
    scabar=plot3([45 45],[110 120],[125 125],'LineWidth',3,'Color',[.5 .5 .5]);
    
    %pause();
    view(-90,0)
    %pause();
    scabar.Visible='off'
    disp 'ready for next'
    clf
    end
    %leglabs={'vgc arbor' 'BPC in' 'AMC in' 'RGC out' 'AMC out'};
    %leggy=legend(leglabs,'Location','northeast');
     
     
     %run the following as needed
     %for cid3
     %plot3([40 40],[110 120],[125 125],'LineWidth',3,'Color',[.5 .5 .5])
     %for cid13
     %plot3([40 40],[130 120],[125 125],'LineWidth',3,'Color',[.5 .5 .5])
     
     %view(-90,0)
end

%% fig for the depths of the synapses
vgcCidList=[2 3 4 5 10 11 13 14 20];
allSynMat={[] [] [] []};
allVox=[0 0 0];
for i=1:length(vgcCidList)
    figCid=vgcCidList(i);
    allSyns=find(curTis.syn.edges(:,1)==figCid | curTis.syn.edges(:,2)==figCid);
    inputs=find(curTis.syn.edges(:,1)==figCid);
    outputs=find(curTis.syn.edges(:,2)==figCid);
    bpcPre=find(preTypes{1}==7);
    amcPre=find(curTis.syn.edges(:,2)'==0 | preTypes{1}==8);
    rgcPost=find(postTypes{1}==1);
    amcPost=find(curTis.syn.edges(:,1)'==0 | postTypes{1}==8);
    
    bpcInputs=intersect(bpcPre,inputs);
    amcInputs=intersect(amcPre,inputs);
    rgcOutputs=intersect(rgcPost,outputs);
    amcOutputs=intersect(amcPost,outputs);
    synMat={bpcInputs amcInputs rgcOutputs amcOutputs};
    allSynMat={[allSynMat{1};bpcInputs] [allSynMat{2};amcInputs] ...
        [allSynMat{3};rgcOutputs] [allSynMat{4};amcOutputs]};
    vgcVox=getCidVox(figCid,5,dsObj,curTis);
    vgcVox=vgcVox{1};
    if figCid~=2
        allVox=vertcat(allVox,vgcVox);
    end
end

%get depths for everything
bins=[0:.03:1];
allVox=allVox/10;
[a,b,voxDepths]=getIPLdepth(allVox(:,3),allVox(:,1),allVox(:,2),[],[]);
voxHist=histcounts(voxDepths,bins);
allHists=zeros(size(voxHist).*[length(allSynMat)+1 1]);
allHists(length(allSynMat)+1,:)=voxHist;
for i=1:length(allSynMat)
    curMat=allSynMat{i};
    curPos=curTis.syn.pos(curMat,:);
    [a,b,curDepths]=getIPLdepth(curPos(:,3),curPos(:,1),curPos(:,2),[],[]);
    curHist=histcounts(curDepths,bins);
    allHists(i,:)=curHist;
end
%figure(); plot(voxHist);
%%
smoothBool=0;
figure();
%t=tiledlayout(1,1);
hold on;
title('Arbor Density and Depth of Synapses in IPL');
%normalize input synapse data.
plotLabels={'arbor Density','BPC inputs','AMC inputs','RGC outputs','AMC outputs','OFF/ON boundary'};
plotColors=[colors;0.8 0.8 0.8];
set(gca, 'YDir','reverse');
bins=[0:.03:1];
depthAxis=bins(1:end-1).*100;
order=[5 1 2 3 4];
for i=1:size(allHists,1)
    j=order(i);
    curPlotDat=allHists(j,:);
    if smoothBool
        curPlotDat=smooth(curPlotDat);
    end
    if j==5
        %ax1=axes(t);
        %plot(curPlotDat(:,2),curPlotDat(:,1),'Color',plotColors(i,:), ...
        %    'LineWidth',2);
        area(depthAxis,curPlotDat,'FaceColor',plotColors(j,:));
    else
        %ax2=axes(t);
        plot(depthAxis,curPlotDat*200,'Color',plotColors(j,:), ...
            'LineWidth',2);
    end
end
%ax2.XAxisLocation = 'top';
%ax2.Color = 'none';
xline(48);
xlabel('IPL depth (%)');
legend(plotLabels(1:5));
yticks([])
ylabel('density a.u.')
view(0,-90);

%% Make a better figure for the depth of the different input synapses and the density of the arbor.


%% plot the relative density of AMC and BPC inputs
smoothBool=1;
figure();
%t=tiledlayout(1,1);
hold on;
title('Inhibitory to Excitatory Synapse Ratio by IPL Depth');
%normalize input synapse data.
plotLabels={'arbor Density','Inhib/Excit ratio','BPC inputs','AMC outputs'}; %,'AMC outputs','OFF/ON boundary'};
plotColors=[0 1 1; 1 0 1; 0 1 0; 1 1 0;0.8 0.8 0.8;0 0 1];
set(gca, 'YDir','reverse');
depthAxis=[0:2:98];
allHists(6,:)=allHists(2,:)./allHists(1,:);
allHists(find(isinf(allHists)))=0;
allHists(find(isnan(allHists)))=0;
order=[5 6 1 2];
for i=1:4
    j=order(i);
    curPlotDat=allHists(j,:);
    if smoothBool
        curPlotDat=smooth(curPlotDat);
    end
    if j==5
        %ax1=axes(t);
        %plot(curPlotDat(:,2),curPlotDat(:,1),'Color',plotColors(i,:), ...
        %    'LineWidth',2);
        area(depthAxis,curPlotDat/3000,'FaceColor',plotColors(j,:));
    elseif j==6
        %ax2=axes(t);
        plot(depthAxis,curPlotDat,'Color',plotColors(j,:), ...
            'LineWidth',2);
    else
        %ax2=axes(t);
        plot(depthAxis,curPlotDat/100,'Color',plotColors(j,:), ...
            'LineWidth',2);
    end
end
%ax2.XAxisLocation = 'top';
%ax2.Color = 'none';
xline(47.5);
xlabel('IPL depth (%)');
legend(plotLabels(1:4));

view(-90,-90);



%% loading skellingtons


% Make sure that you have cellNav open and on the AprilMerge volume so that
% the tis.mat is loaded correctly. Alternatively, the tis.mat could be
% loaded separately.

%number of random iterations
randLength=100;
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
vgcCidList=[2 3 4 5 13 14];
loadSkel=0;
if ~exist('allSkels')
    allSkels=cell(1,6);
    for i=1:length(vgcCidList)
        curCid=vgcCidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end
end
if ~exist('curTis')
    global glob tis
    curTis=tis;
end
loadAll=0;
if loadAll
    for i=1:length(allSkels)
        curSkel=allSkels{i};
        curSWC=nep2swc(curSkel.sm.nep);
        allSkels{i}.swc=curSWC;
    end
end

tipArray=[10];
graphAll=1;
plot1=0;
randomTipOrder=0;
noPlot=1;
for tipIt=1:length(tipArray)
    tipDistCutoff=tipArray(tipIt);
    showTipType=1;
    if showTipType
        allTipHistDat={};
        if plot1
            f3=figure();
            fl2=tiledlayout('flow');
        end
        %title(tl2,"Order of Operations at tips");
        allDists=zeros(length(allSkels),1);
        for i=1:length(allSkels)
            %pause();
            quickHistDat=[];
            curCidHistDat={};
            curCid=curSkel.sm.cid;
            if plot1
                nexttile([4 1]);
                hold on
                title(num2str(curCid));
            end
            %nexttile
            colCodes=[[1 0 0];[0 1 0];[0 0 1];[1 0 1];[0 1 1]];
            curSkel=allSkels{i};
            curSM=curSkel.sm;
            curSWC=curSkel.swc;
            curPred=curSWC.pred(curSWC.arbor2swcID)+1;
            curPredInv=zeros(size(curPred));
            curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
            curSM.pred=curPred;
            allEdges=curSM.arbor.edges;
            uniqueNodes=unique(allEdges(:));
            nodeCounts=histcounts(allEdges,uniqueNodes);
            tipIDs=find(nodeCounts==1);
            forkIDs=find(nodeCounts>2);
            usedNodeList=[];
            sumVect=[0 0 0];
            %need to sort tipIDs
            distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
            [srtd srtIdx] = sort(distFromRoot,'descend');
            tipIDsSrtd=tipIDs(srtIdx');
            if graphAll
                graphTipsSrtd=tipIDsSrtd;
            else
                goodTips=TOI{i};
                goodTips=find(goodTips==1);
                graphTipsSrtd=tipIDsSrtd(goodTips);
            end
            if randomTipOrder
                graphTipsSrtd=graphTipsSrtd(randperm(length(graphTipsSrtd)));
            end
            totalDistCovered=0;
            for k=1:length(graphTipsSrtd)
                %curTipHistDat={};
                previousEdgeID=0;
                curTip=graphTipsSrtd(k);
                cidsTipsOrd={};
                if curTip~=1
                    if ~noPlot
                        drawnow
                        hold on
                    end
                    %xline(k);
                    curNode=curTip;
                    curEdgeLength=curSM.nep.props.nodeLength(curNode);
                    curTipDist=curEdgeLength;
                    curTipHistDat=[curNode,9,curEdgeLength,curTipDist];
                    
                    beginningPos=curSM.arbor.nodes.pos(curNode,:);
                    
                    abandon=0;
                    while ~abandon & curTipDist<tipDistCutoff
                        nextNode=curPredInv(curNode);
                        if nextNode==0 | ismember(nextNode,usedNodeList)
                            abandon=1;
                            break
                        end
                        if ismember(curNode,forkIDs)
                            if ~noPlot
                                scatter(k,curTipDist,25,'b_','LineWidth',1);
                            end
                            %                     scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                            %                         curSM.arbor.nodes.pos(curNode,2), ...
                            %                         curSM.arbor.nodes.pos(curNode,3),25,'bo','filled');
                            %                     quickHistDat=[quickHistDat;[curTip,6,curTipDist]];
                            endingPos=curSM.arbor.nodes.pos(curNode,:);
                            curTipHistDat=[curTipHistDat;[curNode,6,curEdgeLength,curTipDist]];
                            quickHistDat=[quickHistDat;[curTip,6,curEdgeLength,curTipDist]];
                            %break
                        end
                        curNodeRad=curSM.arbor.nodes.rad(curNode);
                        if curTipDist>20
                            graphTipDist=20;
                        else
                            graphTipDist=curTipDist;
                        end
                        %scatter(k,curTipDist,curNodeRad*50,cool8(round(graphTipDist*12.5)+1,:));
                        %                 scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                        %                     curSM.arbor.nodes.pos(curNode,2), ...
                        %                     curSM.arbor.nodes.pos(curNode,3), ...
                        %                     curNodeRad*10,cool8(round(graphTipDist*12.5)+1,:));
                        %curEdgeLength=curSM.skel2skel.linDist(nextNode,curNode);
                        
                        curTipDist=curTipDist+curEdgeLength;
                        %drawnow
                        closeSyns=find(curSM.syn2Skel.closest==curNode);
                        if ~isempty(closeSyns)
                            getSyns=1;
                            if getSyns
                                for m=1:length(closeSyns)
                                    curSyn=closeSyns(m);
                                    curSynType=getSMsynType(curSyn,curCid,curSM,curTis);
                                    curTipHistDat=[curTipHistDat;[curNode,curSynType,curEdgeLength,curTipDist]];
                                    quickHistDat=[quickHistDat;[curNode,curSynType,curEdgeLength,curTipDist]];
                                    if ~noPlot
                                        scatter(k,curTipDist,15,colCodes(curSynType,:),'filled');
                                    end
                                    %scatter3(curSM.arbor.nodes.pos(curNode,1), ...
                                    %    curSM.arbor.nodes.pos(curNode,2), ...
                                    %    curSM.arbor.nodes.pos(curNode,3), ...
                                    %    curNodeRad*12,cool8(round(graphTipDist*12.5)+1,:),'filled');
                                    
                                end
                            end
                            
                        end
                        usedNodeList=[usedNodeList curNode];
                        curNode=nextNode;
                        %pause();
                    end
                    
                    if curTipDist<tipDistCutoff | abandon
                        curTipHistDat=[];
                    else
                        totalDistCovered=totalDistCovered+curTipDist;
                    end
                    
                end
                curCidHistDat{k}=curTipHistDat;
            end
            curCid
            %sumVect
            
            if plot1
                title(num2str(curSM.cid));
                p1=scatter(0,0,1,colCodes(1,:),'filled');
                p2=scatter(0,0,1,colCodes(2,:),'filled');
                p3=scatter(0,0,1,colCodes(3,:),'filled');
                p4=scatter(0,0,1,colCodes(4,:),'filled');
                p5=scatter(0,0,1,colCodes(5,:),'filled');
                %legend([p1 p2 p3 p4 p5],{'bpc in','amc in','rgc out','amc out','unk'});
                ylim([0,11]);
            end
            winSize=20;
            %xlim([-winSize,winSize]);
            %ylim([-winSize,winSize]);
            %zlim([-winSize,winSize]);
            avgVect=sum(sumVect,1)/size(sumVect,1);
            
            %plot3([0 avgVect(1)*50],[0 avgVect(2)*50],[0 avgVect(3)*50],':','Color','b','LineWidth',5);
            %view(90,0);
            %xlabel('X');
            %ylabel('Y');
            %zlabel('Z');
            allTipHistDat{i}=curCidHistDat;
            allQuickHistDat{i}=quickHistDat;
            allDists(i)=totalDistCovered;
            if plot1
                a=gca;
                a.XGrid='on';
                a.XTick=[0:length(allTipHistDat{i})];
            end
        end
        
    end
end

%% This block is for drawing Skels (with node radiiii and everything)
f1=figure();
hold on
radMult=25;
skipMult=3;
i=2;
SknodeList=graphTipsSrtd(:);
SknodeList=tipIDsSrtd;
lengthList=zeros(length(SknodeList),1);
usedList=[];
allBranchNodes={};
for sknodeIt=1:length(SknodeList)
    curCol=rand([1 3]);
    curSknode=SknodeList(sknodeIt);
    curNodeList=curSknode;
    curTipDist=0;
    itr=1;
    while curSknode~=0
        if ~ismember(curSknode,usedList)
            curEdgeLength=curSM.nep.props.nodeLength(curSknode);
            curTipDist=curTipDist+curEdgeLength;
            
            usedList=[usedList curSknode];
            curSknode=curPredInv(curSknode);
            if curSknode==0
                break
            end
            curNodeList=[curNodeList curSknode];
            if mod(itr,skipMult)==0
                %curSknode
                curPos=curSM.nep.pos(curSknode,:);
                %curPos
                curRad=curSM.nep.nodeRad(curSknode);
                %curRad
                [a,b,curDepths]=getIPLdepth(curPos(:,3),curPos(:,1),curPos(:,2),[],[]);
                scatter3(curPos(1),curPos(2),curDepths,10,curCol,'filled');
                %scatter3(curPos(1),curPos(2),curDepths,curRad*radMult,curCol,'filled');
            end
            itr=itr+1;
        else
            %lengthList(sknodeIt)=curTipDist;
            break
        end
        %drawnow
    end
    drawnow
    allBranchNodes{sknodeIt}=curNodeList;
    lengthList(sknodeIt)=curTipDist;
end

%% This will plot out the branches


%%
branchLengthCutoff=5;
f2=figure();
hold on
i=1;
SknodeList=graphTipsSrtd(:);
%lengthList=zeros(length(SknodeList),1);
%usedList=[];
for sknodeIt=1:length(SknodeList)
    
    %curSknode=SknodeList(sknodeIt);
    %curNodeList=curSknode;
    %curTipDist=0;
    %while curSknode~=0
    %if ~ismember(curSknode,usedList)
    %curEdgeLength=curSM.nep.props.nodeLength(curSknode);
    %curTipDist=curTipDist+curEdgeLength;
    %curPos=curSM.nep.pos(curSknode,:);
    if lengthList(sknodeIt)>branchLengthCutoff
        curCol=rand([1 3]);
        nodeList=allBranchNodes{sknodeIt};
        nodeList(nodeList==0)=[];
        posList=curSM.nep.pos(nodeList,:);
        scatter3(posList(:,1),posList(:,2),posList(:,3),5,curCol,'filled');
        %for j=1:length(allBranchNodes{sknodeIt})
        %    curNode=allBranchNodes{sknodeIt}(j);
        %    curPos=curSM.nep.pos(curNode,:);
        %    scatter3(curPos(1),curPos(2),curPos(3),2,curCol,'filled');
        %drawnow
        %end
        %drawnow
        drawnow
    end
    %usedList=[usedList curSknode];
    %curSknode=curPredInv(curSknode);
    %curNodeList=[curNodeList curSknode];
    %drawnow
    %else
    %lengthList(sknodeIt)=curTipDist;
    %    break
    %end
end
%allBranchNodes{sknodeIt}=curNodeList;
%lengthList(sknodeIt)=curTipDist;
%end

%% physCorr
physDir='Y:\karlsRetina\EmilyImages\CorrectRegistratedFiles\';
analDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\';
datDir='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\';
physSaveDir='C:\work\phys\';
load([datDir 'ptDat.mat']);
loc1=[22,103];
loc2=[14,137];

testFile='Ai148_129SVG3_Translation_122618_1010.mat';
test=load([physDir testFile]);
if 0
    figure();
    hold on
    for i=400:1000
        imshow(test.I(:,:,i),[0 50], 'InitialMagnification', 400)
        hold on
        %image(test.I(:,:,i));
        plot(loc1(2),loc1(1), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        plot(loc2(2),loc2(1), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        drawnow
        pause(0.01)
    end
end

if 0
    fileNames={'Ai148_129SVG3_Translation_122618_1005.mat', ...
        'Ai148_129SVG3_Translation_122618_1006.mat', ...
        'Ai148_129SVG3_Translation_122618_1007.mat', ...
        'Ai148_129SVG3_Translation_122618_1008.mat', ...
        'Ai148_129SVG3_Translation_122618_1009.mat', ...
        'Ai148_129SVG3_Translation_122618_1010.mat', ...
        'Ai148_129SVG3_Translation_122618_2001.mat', ...
        'Ai148_129SVG3_Translation_122618_2002.mat', ...
        'Ai148_129SVG3_Translation_122618_2003.mat', ...
        'Ai148_129SVG3_Translation_122618_2004.mat', ...
        'Ai148_129SVG3_Translation_122618_2005.mat', ...
        'Ai148_129SVG3_Translation_122618_2006.mat'};
    saveNames={'1005.mat','1006.mat','1007.mat','1008.mat','1009.mat','1010.mat','2001.mat','2002.mat','2003.mat','2004.mat','2005.mat','2006.mat'};
    %uberDat=cell(length(saveNames),1); %zeros(32*256,32*256,length(saveNames));
    pixelcomps=nchoosek([1:32*256],2);
    for m=3:12
        fileName=fileNames{m};
        test=load([physDir fileName]);
        corrMat=zeros(32*256,32*256);
        sz=size(test.I);
        sz=sz(1:2);
        for i=1:size(pixelcomps,1)
            if mod(i,100000)==0
                i
            end
            curComp=pixelcomps(i,:);
            [ax,ay]=ind2sub(sz,curComp(1));
            [bx,by]=ind2sub(sz,curComp(2));
            atc=squeeze(test.I(ax,ay,:));
            btc=squeeze(test.I(bx,by,:));
            curCorr=corr(atc,btc);
            corrMat(curComp(1),curComp(2))=curCorr;
        end
        %uberDat{i}=corrMat;
        outputFilename=[physSaveDir saveNames{m}];
        save(outputFilename,'corrMat');
    end
    %load('C:\Users\karlf\Documents\work2022\corrMat.mat');
end


testPixel=squeeze(test.I(loc1(1),loc1(2),:));
testPixel2=squeeze(test.I(loc2(1),loc2(2),:));

ind1=sub2ind([32 256],loc1(1),loc1(2));
ind2=sub2ind([32 256],loc2(1),loc2(2));

figure(); plot(testPixel);
hold on
plot(testPixel2);

figure(); histogram(testPixel,[0:2:200]);
hold on
histogram(testPixel2,[0:2:200]);
legend();

corrCut=0.5;
f3=figure();
%tl3=tiledlayout(f3, 'flow');
subplot(3,1,1);
testPixCorr=corrMat(:,ind1)+corrMat(ind1,:)';
goodCorrInds=find(testPixCorr>corrCut);
[corrPxX,corrPxY]=ind2sub([32 256],goodCorrInds);
testPixCorrImg=reshape(testPixCorr,[32 256]);


testPixCorr2=corrMat(:,ind2)+corrMat(ind2,:)';
testPixCorrImg2=reshape(testPixCorr2,[32 256]);



%ax = nexttile(tl3);
imshow(testPixCorrImg,'InitialMagnification', 400 );
hold on
scatter(corrPxY,corrPxX);
subplot(3,1,2);
%ax = nexttile(tl3);
imshow(testPixCorrImg2);
subplot(3,1,3);
%ax = nexttile(tl3);
imshowpair(testPixCorrImg,testPixCorrImg2);

testPixTrace=testPixel(400:600);
testPix2Trace=testPixel2(400:600);

figure();
image(testPixCorrImg,'CDataMapping','scaled');
colorbar;

%get the best pixel correlation to the first one
testPixPalID=find(testPixCorr==max(testPixCorr));
[x1 y1]=ind2sub([32 256],testPixPalID);
testPixelPal=squeeze(test.I(x1,y1,:));

testPix2PalID=find(testPixCorr2==max(testPixCorr2));
[x2 y2]=ind2sub([32 256],testPix2PalID);
testPixel2Pal=squeeze(test.I(x2,y2,:));

figure();

plot(testPixel(400:600))
hold on
plot(testPixel2(400:600));
plot(testPixelPal(400:600));
plot(testPixel2Pal(400:600));
legend();

figure();
histogram(testPixCorr);

%% 07-05-22
% PCA test
% \
imgRshp=reshape(test.I,8192,[]);
%pcaTest=pca(imgRshp);
[coeff,score,latent,tsquared,explained]=pca(imgRshp');
figure();
histogram(coeff);
figure();
pcaComp1=reshape(coeff(:,1),[32 256]);
pcaComp2=reshape(coeff(:,2),[32 256]);
pcaComp3=reshape(coeff(:,3),[32 256]);
pcaComp4=reshape(coeff(:,4),[32 256]);
pcaComp5=reshape(coeff(:,5),[32 256]);
pcaComp6=reshape(coeff(:,6),[32 256]);
subplot(6,1,1)
image(pcaComp1(:,:),'CDataMapping','scaled');
colorbar
subplot(6,1,2)
image(pcaComp2(:,:),'CDataMapping','scaled');
colorbar
subplot(6,1,3)
image(pcaComp3(:,:),'CDataMapping','scaled');
colorbar
subplot(6,1,4)
image(pcaComp4(:,:),'CDataMapping','scaled');
colorbar
subplot(6,1,5)
image(pcaComp5(:,:),'CDataMapping','scaled');
colorbar
subplot(6,1,6)
image(pcaComp6(:,:),'CDataMapping','scaled');
colorbar

%% Skeleton stuff for after running the setup of the tip analysis
% requires graphTipsSrtd, which is just a list of the tips from the longest
% distance to root.
plotSkel=1;
for k=1:length(vgcCidList)
    %% copied from earlier to get the lists for the tips
    curSkel=allSkels{k};
    curSM=curSkel.sm;
    curSWC=curSkel.swc;
    curPred=curSWC.pred(curSWC.arbor2swcID)+1;
    curPredInv=zeros(size(curPred));
    curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
    curSM.pred=curPred;
    allEdges=curSM.arbor.edges;
    uniqueNodes=unique(allEdges(:));
    nodeCounts=histcounts(allEdges,uniqueNodes);
    tipIDs=find(nodeCounts==1);
    forkIDs=find(nodeCounts>2);
    usedNodeList=[];
    sumVect=[0 0 0];
    %need to sort tipIDs
    distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    graphTipsSrtd=tipIDsSrtd;
    
    if plotSkel
        figure();
        hold on
    end
    for l=1:length(graphTipsSrtd)
        previousEdgeID=0;
        curTip=graphTipsSrtd(k);
        cidsTipsOrd={};
        if curTip~=1
            curNode=curTip;
            curEdgeLength=curSM.nep.props.nodeLength(curNode);
            curTipDist=curEdgeLength;
            curTipHistDat=[curNode,9,curEdgeLength,curTipDist];
            beginningPos=curSM.arbor.nodes.pos(curNode,:);
            while ~abandon & curTipDist<tipDistCutoff
                nextNode=curPredInv(curNode);
                if nextNode==0 | ismember(nextNode,usedNodeList)
                    abandon=1;
                    break
                end
                
                
                
            end
        else
            endingPos=curSM.arbor.nodes.pos(curNode,:);
        end
    end
end

%% making the tables for the update


