%load skeletons if necessary
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

%% Main
%% quick check on the types of rgc targets and frequencies
targetTypeDat=cid2type(curTis.syn.edges(ismember(curTis.syn.edges(:,2),vgcCidList),1),curTis);
rgcTargSubs=targetTypeDat{3}(targetTypeDat{1}==1);
histDat=histcounts(rgcTargSubs,[1:length(curTis.cells.type.subTypeNames{1})]);
targSubIDs=find(histDat>10);
typesOfInterest={'4i','4ow','5ti','37','63','6sw'};
targSubIDs=find(ismember(curTis.cells.type.subTypeNames{1},typesOfInterest));
targSubsNames=curTis.cells.type.subTypeNames{1}(targSubIDs);
%targSub
% rgc with significant syn #
f0=figure();
hold on
bar(histDat(targSubIDs));
xticks([1:length(targSubIDs)]);
xticklabels(targSubsNames);
ylabel('# of synapses from VG3');
xlabel('RGC subtype');
title('VG3 to RGC outputs by subtype');
yline([10:10:120]);

%same check on the bipolar side
inputTypeDat=cid2type(curTis.syn.edges(ismember(curTis.syn.edges(:,1),vgcCidList),2),curTis);
bpcInputSubs=inputTypeDat{3}(inputTypeDat{1}==7);
bpcTypeBins=[-0.5:1:length(curTis.cells.type.subTypeNames{7})+0.5];
histDat=histcounts(bpcInputSubs,bpcTypeBins);
bpcOFFsubs=[1:5,15,18,19,24]+1;
bpcONsubs=[6:12,14,17,20,21,13,16]+1;

v2r=getClassConn('amc','vgc','rgc','all',curTis);
postVRcids=curTis.syn.edges(v2r,1);
postVRtypes=allPostTypes{3}(v2r);
uniquePostVRtypes=unique(postVRtypes);

rgcSubTypeList=[1:length(uniquePostVRtypes)];%[21 23 27 29 51];

bpcSubOrder=horzcat(bpcOFFsubs,bpcONsubs);
%targSubIDs=find(histDat>20);
targSubIDs=bpcSubOrder;
targSubsNames=curTis.cells.type.subTypeNames{7}(targSubIDs-1);

% BPC with significant number of inputs
f01=figure();
hold on
sigTarg=find(histDat>5);
smallTargSub=intersect(sigTarg,targSubIDs);
smallTargSubsNames=curTis.cells.type.subTypeNames{7}(smallTargSub-1);

bar(histDat(smallTargSub));
xticks([1:length(smallTargSub)]);
xticklabels(smallTargSubsNames);
ylabel('# of synapses to VG3');
xlabel('BPC input subtype');
title('BPC inputs to VG3 by subtypes');
%xline(9.5);
yline([10:50:250]);
yline([10:10:100]);

%% check the depths of the diferent classes of synapses and see which are the most overlappys
%whole dataset variables
allSynPos=curTis.syn.pos;
allSynEdges=curTis.syn.edges(:,1:2);
[n,m,allSynDepth]=getIPLdepth(allSynPos(:,3),allSynPos(:,1),allSynPos(:,2),[],[]);
bpcTypeSets={[[7 7 7 7];[ 3 4 5 19]],[[7 7 7 7 7 7];[6 7 8 19 20 21]]};
inputs=ismember(allSynEdges(:,1),vgcCidList);
outputs=ismember(allSynEdges(:,2),vgcCidList);
preTypeDat=cid2type(allSynEdges(:,2),curTis);
postTypeDat=cid2type(allSynEdges(:,1),curTis);

%These are all the bpc subtypes that have > 20 VGC synapses
bpcSubs=[3 4 5 6 7 8];
bpcSubNames=curTis.cells.type.subTypeNames{7}(bpcSubs);
%These are all the rgc subtypes that have > 20 VGC synapses
rgcSubs=[21 23 27 29 51 58];
rgcSubNames=curTis.cells.type.subTypeNames{1}(rgcSubs);

depthHistBins=[0:0.05:1];

%indices
%super sketch, but I'm just saying that input&~bpc = amc
bpcSynIDs=find(preTypeDat{1}==7 & inputs'==1);
amcSynIDs=find(inputs'==1 & preTypeDat{1}~=7);
onBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcONsubs));
offBpcSynIDs=find(preTypeDat{1}==7 & inputs'==1 & ismember(preTypeDat{3},bpcOFFsubs));
bpcIDs={};
bpcHistDat={};
for i=1:length(bpcSubs)
    curSub=bpcSubs(i);
    bpcIDs{i}=find(preTypeDat{1}==7 & inputs'==1 & preTypeDat{3}==curSub);
    curHist=histcounts(allSynDepth(bpcIDs{i}),depthHistBins);
    bpcHistDat{i}=curHist;
end
rgcIDs={};
rgcHistDat={};
for i=1:length(rgcSubs)
    curSub=rgcSubs(i);
    rgcIDs{i}=find(postTypeDat{1}==1 & outputs'==1 & postTypeDat{3}==curSub);
    curHist=histcounts(allSynDepth(rgcIDs{i}),depthHistBins);
    rgcHistDat{i}=curHist;
end

%%
% try to figure out the 'distances' between input and output synapse depth
% distributions
histDistMat=zeros(length(bpcHistDat),length(rgcHistDat));
for j=1:length(bpcHistDat)
    for k=1:length(rgcHistDat)
        histDistMat(j,k)=pdist2(bpcHistDat{j},rgcHistDat{k});
    end
end





%% depth dist plot
f02=figure();
hold on
subplot(2,1,1);
hold on
title('BPC input IPL depths');
for i=1:length(bpcIDs)
    curHistDat=bpcHistDat{i};
    a=area(curHistDat,'LineStyle','--');
    a.FaceAlpha = 0.2;
    %plot(1:length(curHistDat),curHistDat,'--');
end
xticks([0.5:2:19.5]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
set(gca, 'YDir','reverse');
legend(bpcSubNames);
xlim([4,15]);
%xline(9.5);
xline(10,'HandleVisibility','off');
xlabel('IPL depth');
ylabel('synapse #');

subplot(2,1,2);
hold on
title('RGC output IPL depths');
for i=1:length(rgcIDs)
    curHistDat=rgcHistDat{i};
    %plot(1:length(curHistDat),curHistDat,':');
    a=area(curHistDat,'LineStyle','--');
    a.FaceAlpha = 0.2;
end
legend(rgcSubNames);
xticks([0.5:2:19.5]);
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});
xline(10,'HandleVisibility','off');
xlim([4,15]);
xlabel('IPL depth');
ylabel('synapse #');
%legend({bpcSubNames{:},rgcSubNames{:}});


%% Go through skels
lc=20;
% 4ow=23; 51=24;
fig2=0;
fig3=0;
fig4=1;
sanityCheck=1;
sanityCheck2=0;
plotAll=1;
%bpcOFFsubs=[1:5,15,18,19,24];
%bpcONsubs=[6:12,14,17,20,21];

%f1=figure();
%hold on
%colMapA=hsv(8);

%f2=figure();
if fig3
    scf3=figure();
    title(['Input near VG3->RGC outputs by RGC subtype; lc=' num2str(lc)])
    hold on
end
if fig4
    if plotAll
        tp1=figure();
        tp2=figure();
    end
end

bpc2rgsColMap=hsv(16);
superBps2rgs=cell(length(rgcSubs),1);
bpcInputDat=struct();

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
    bpcSynTypeIDs=cell(length(curTis.cells.type.subTypeNames{7}),1);
    for o=1:length(bpcSynTypeIDs)
        bpcSynTypeIDs{o}=find(preTypeDat{1}==7 & inputs'==1 & preTypeDat{3}==o);
    end
    
    
    
    curDistMat=curSM.syn2Skel.syn2SynDist;
    vMat=zeros(size(curSM.syn2Skel.syn2SynDist));
    vMat(bpcSynIDs,:)=1;
    vMat(amcSynIDs,:)=-1;
    infMat=vMat.*exp(-curDistMat./lc);
    EIbalMat=sum(infMat,1);
    
    %figure(f1);
    %scatter(curSyns.depth(onBpcSynIDs),1+i/20,25,colMapA(i,:),'filled');
    %scatter(curSyns.depth(offBpcSynIDs),2+i/20,25,colMapA(i,:),'filled');
    %scatter(curSyns.depth(amcSynIDs),3+i/20,25,colMapA(i,:),'filled');
    
    %histogram(curSyns.depth);
    if 0
        scf1=figure();
        histogram(EIbalMat,50);
        title(['cid' num2str(curCid) ': E/I Bal; all syns; lc=' num2str(lc)]);
    end
    %curSM.syn2Skel.syn2SynDist;
    
    colMapB=hsv(256);
    
    
    searchRadius=50;
    for j=1:length(rgcSubs)
        curRGCsubType=rgcSubs(j);
        curSynIDs=find(postTypeDat{1}==1 & postTypeDat{3}==curRGCsubType & outputs'==1);
        curSynDists=curDistMat(curSynIDs,:);
        nearSynIDs=find(min(curSynDists,[],1)<searchRadius);
        nearBpcIDs=intersect(nearSynIDs,bpcSynIDs);
        nearAmcIDs=intersect(nearSynIDs,amcSynIDs);
        bps2rgs=[];
        for n=1:length(curSynIDs)
            curSID=curSynIDs(n);
            sD=curDistMat(curSID,:);
            nearSIDs=find(min(sD,[],1)<searchRadius);
            nearBIDs=intersect(nearSIDs,bpcSynIDs);
            nearBtypes=preTypeDat{3}(nearBIDs);
            bps2rgs=vertcat(bps2rgs,horzcat(nearBtypes',sD(nearBIDs)'));
        end
        bps2rgsrtd=sortrows(bps2rgs);
        if ~isempty(bps2rgsrtd)
        bps2rgsrtd=bps2rgsrtd(bps2rgsrtd(:,1)>0,:);
        superBps2rgs{j}=vertcat(superBps2rgs{j},bps2rgsrtd);
        typesPresent=unique(bps2rgsrtd(:,1));
        if plotAll
            if fig4
                figure(tp1);
                ars=2;
                scatter3(curSM.arbor.nodes.pos(1:ars:end,1)/10,curSM.arbor.nodes.pos(1:ars:end,2)/10, ...
                    curSM.arbor.nodes.pos(1:ars:end,3)/10,curSM.arbor.nodes.rad(1:ars:end)*4,'.k');
                hold on
                scatter3(curSyns.pos(curSynIDs,1),curSyns.pos(curSynIDs,2),curSyns.pos(curSynIDs,3),200,'mo');
                nearSizes=(20-min(curSynDists(:,nearBpcIDs),[],1))/0.1;
                scatter3(curSyns.pos(nearBpcIDs,1),curSyns.pos(nearBpcIDs,2),curSyns.pos(nearBpcIDs,3), ...
                    nearSizes,'co','filled');
                
                figure(tp2);
                %hold on
                
                pieSliceSize=2*pi/size(bps2rgsrtd,1);
                angV=0:pieSliceSize:2*pi-pieSliceSize;
                [sharedvals,idx] = ismember(bps2rgsrtd(:,1),typesPresent);
                %curHD=histcounts(idx,[0:length(typesPresent)]);
                colV=bpc2rgsColMap(idx,:);
                %tp2=figure();
                polarscatter(angV,20-bps2rgsrtd(:,2),(searchRadius-bps2rgsrtd(:,2))*10,colV,'o','filled');
                hold on
                ax = gca;
                ax.GridLineStyle='none';
                ax.RTickLabel={};
                ax.ThetaTickLabel={};
                pieCounter=0;
                lineDensity=16;
                for w=1:length(typesPresent)
                    curType=typesPresent(w);
                    lineMat=repmat([0;20],length(find(bps2rgsrtd(:,1)==curType))*lineDensity/2,1);
                    angMat=[pieCounter+pieSliceSize/lineDensity:pieSliceSize/lineDensity:pieSliceSize*length(find(bps2rgsrtd(:,1)==curType))+pieCounter];
                    bgcol=(bpc2rgsColMap(w,:)+[2 2 2])/3;
                    polarplot(angMat,lineMat,'LineWidth',2,'Color',bgcol);
                    %pieCounter=pieSliceSize*length(find(bps2rgsrtd(:,1),curType))+pieCounter;
                    pieCounter=angMat(end);
                end
                polarscatter(angV,20-bps2rgsrtd(:,2),(searchRadius-bps2rgsrtd(:,2))*10,colV,'o','filled');
                title(['vgc '+string(curCid) 'rgc '+string(curTis.cells.type.subTypeNames{1}(curRGCsubType)) ]);
                
                %polarscatter(bps2rgsrtd(:,2),bpc2rgsColMap(bps2rgsrtd(:,1),:),'*');
            end
        end
        end
        %for m=1:length(curSynIDs)
        
        
        
        
    end
end


%% super Pie
%figure();
circleDensity=180;
figure();
tL=tiledlayout(3,2,'Padding', 'none', 'TileSpacing', 'compact');
for j=1:length(rgcSubs)
    curRGCsubType=rgcSubs(j);
    nexttile
    %figure(sp1);
    %hold on
    bps2rgs=superBps2rgs{j};
    bps2rgsrtd=sortrows(bps2rgs);
    typesPresent=unique(bps2rgsrtd(:,1));
    pieSliceSize=2*pi/size(bps2rgsrtd,1);
    angV=0:pieSliceSize:2*pi-pieSliceSize;
    [sharedvals,idx] = ismember(bps2rgsrtd(:,1),typesPresent);
    %curHD=histcounts(idx,[0:length(typesPresent)]);
    colV=bpc2rgsColMap(idx,:);
    %tp2=figure();
    polarscatter(angV,20-bps2rgsrtd(:,2),(20-bps2rgsrtd(:,2))*10,colV,'.');
    hold on
    ax = gca;
    ax.GridLineStyle='none';
    %ax.RTick=([5 10 15 20]);
    ax.RTickLabel={};
    ax.ThetaTickLabel={};
    pieCounter=0;
    lineDensity=16;
    textDat=[];
    pieRack=[];
    typesPresent=typesPresent(typesPresent>0);
    for w=1:length(typesPresent)
        curType=typesPresent(w);
        if curType>0
            lineMat=repmat([0;20],length(find(bps2rgsrtd(:,1)==curType))*lineDensity/2,1);
            angMat=[pieCounter+pieSliceSize/lineDensity:pieSliceSize/lineDensity:pieSliceSize*length(find(bps2rgsrtd(:,1)==curType))+pieCounter];
            bgcol=(bpc2rgsColMap(w,:)+[2 2 2])/3;
            polarplot(angMat,lineMat,'LineWidth',2,'Color',bgcol);
            %pieCounter=pieSliceSize*length(find(bps2rgsrtd(:,1),curType))+pieCounter;
            pieCounter=angMat(end);
            textDat(w)=mean(angMat);
            pieRack(w)=pieCounter;
            %text(mean(angMat),10,string(curTis.cells.type.subTypeNames{7}(curType)));
        end
    end
    polarscatter(angV,bps2rgsrtd(:,2),(20-bps2rgsrtd(:,2))*20,colV,'.');
    %polarscatter(angV,20-bps2rgsrtd(:,2),(20-bps2rgsrtd(:,2))*10,colV,'o');
    title(['rgc '+string(curTis.cells.type.subTypeNames{1}(curRGCsubType)) ]);
    for w=1:length(typesPresent(typesPresent>0))
        curType=typesPresent(w);
        if curType>0
            curName=string(curTis.cells.type.subTypeNames{7}(curType));
            text(textDat(w),15,curName,'HorizontalAlignment','center','VerticalAlignment','middle');
        else
            curName='none';
        end
        text([2.8 2.8 2.8],[20 10 5],{'20um','10um','5um'},'HorizontalAlignment','center','VerticalAlignment','middle');
        polarplot([2*pi/circleDensity:2*pi/circleDensity:2*pi],repmat([10],circleDensity),'k:');
        polarplot([2*pi/circleDensity:2*pi/circleDensity:2*pi],repmat([5],circleDensity),'k:');
        polarplot([2*pi/circleDensity:2*pi/circleDensity:2*pi],repmat([20],circleDensity),'k:');
        polarplot([pieRack(w) pieRack(w)],[0 20],'k');
    end
    %polarscatter(bps2rgsrtd(:,2),bpc2rgsColMap(bps2rgsrtd(:,1),:),'*');
    %ax.GridLineStyle=(
    %ax.RTick=([5 10 15 20]);
    %ax.RTickLabel={'20','15','10','5'};
end
%for m=1:length(curSynIDs)



%% Swarm!!!
if fig5
    col1=[1 .5 0];
    col2=[0 .5 1];
    sp1=figure();
    hold on
    tiledlayout(2,3);
    %colorMaps=vertcat(repmat([1 0 1],3,1),repmat([0 1 1],3,1));
    %length(find(rgcTargSubs==21))
    rgcSynCounts=histcounts(rgcTargSubs,[0:59]);
        
for j=1:length(rgcSubs)
    curRGCsubType=rgcSubs(j);
    curRScount=length(find(rgcTargSubs==curRGCsubType));
    bps2rgs=superBps2rgs{j};
    infMat=exp(-bps2rgs(:,2)./lc);
    typesPresent=unique(bps2rgs(:,1));
    totInfMat=zeros(length(typesPresent),1);
    for u=1:length(typesPresent)
        curSType=typesPres
        
        ent(u);
        totInfMat(u)=sum(infMat(bps2rgs(:,1)==curSType));
    end
    infPercMat=totInfMat./sum(totInfMat);
    barCol=vertcat(repmat(col1,3,1),repmat(col2,4,1));
    colorMat=repmat(col1,size(bps2rgs,1),1);
    colorMat(find(ismember(bps2rgs(:,1),bpcONsubs)),:)= ...
        repmat(col2,length(find(ismember(bps2rgs(:,1),bpcONsubs))),1);
    nexttile
    hold on
    swarmchart(bps2rgs(:,1),infMat,infMat.^2*20+2,colorMat,'filled');
    curBar=bar([3 4 5 6 7 8 9],infPercMat(1:7));
    curBar.FaceColor='flat';
    curBar.CData([4:7],:)=repmat(col2,4,1);
    curBar.CData([1:3],:)=repmat(col1,3,1);
    curBar.FaceAlpha=0.2;
    curBar.EdgeAlpha=0;
    xlim([2 10]);
    title(['rgc '+string(curTis.cells.type.subTypeNames{1}(curRGCsubType)) ...
        '('+string(curRScount)+')']);
    xticks([3 4 5 6 7 8 9]);
    xticklabels(curTis.cells.type.subTypeNames{7}([3 4 5 6 7 8 9]));
    if j==1 | j==4
    ylabel('Influence');
    end
end
end

%% Everything after here can be ignored
if 0
    if 0
        %figure(f1);
        %scatter(curSyns.depth(curSynIDs),3+i/20+j,25,colMapA(i,:),'filled');
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
                    bpcInNear=bpcIn(bpcIn<20);
                    bpcOffInNear=bpcOffIn(bpcOffIn<20);
                    bpcOnInNear=bpcOnIn(bpcOnIn<20);
                    amcIn=abs(curDistMat(amcSynIDs,curSynIDs(k)));
                    amcInNear=amcIn(amcIn<20);
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
    if fig3
        figure(scf3);
        %xlabel('synapse #');
        ylabel('influence');
        set(gca, 'YDir','reverse');
        xticks([14:10:length(rgcSubs)*10+12]);
        xticklabels({'4i','4ow','51','63','37'});
        xlim([5,length(rgcSubs)*10+5]);
    end
    %     figure();
    %     hold on
    %     scatter3(curSM.arbor.nodes.pos(:,1),curSM.arbor.nodes.pos(:,2), ...
    %         curSM.arbor.nodes.pos(:,3),curSM.arbor.nodes.rad*5,'.k');
    %     scatter3(curSyns.pos(:,1)*10,curSyns.pos(:,2)*10,curSyns.pos(:,3)*10, ...
    %         25,synCol,'o','filled');
    %     title('Average ratio of topological to euclidean distance of synapse to all other synapses?');
    
    figure(f1);
    xline(0.4475);
    xlim([0.1,0.8]);
    yticks([1:length(rgcSubs)+3]+.2);
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
end
