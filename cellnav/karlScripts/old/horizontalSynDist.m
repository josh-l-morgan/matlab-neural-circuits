%% get the horizontal distributions of the synapse types to the center of the arbor in XY


global tis

vgcCidList=[2 3 4 5 13 14 20];
highcut=.65;
lowcut=.3;
fullCidList=tis.cells.cids(:);
fullTypeList=tis.cells.type.typeID(:);
allSynEdges=tis.syn.edges(:,[1 2]);
allSynPos=tis.syn.pos;

preType=zeros(length(allSynEdges(:,2)),1);
postType=zeros(length(allSynEdges(:,1)),1);

for k=1:length(preType)
    preCid=allSynEdges(k,2);
    if preCid==0
        preT=8;
    else
        preT=tis.cells.type.typeID(find(preCid==tis.cells.cids));
    end
    preType(k)=preT;
end

for k=1:length(postType)
    postCid=allSynEdges(k,1);
    if postCid==0
        postT=8;
    else
        postT=tis.cells.type.typeID(find(postCid==tis.cells.cids));
    end
    if isempty(postT)
        postT=8;
    end
    postType(k)=postT;
end
%% 
binVec=[0:2:70];
vgcDistDat=struct;
    figure();
    hold on
for i=[1 2 3 4 5]
    curVGC=vgcCidList(i);
    %get the center of the arbor
    sm=skelDat(i).SM.sm;
    meanXY = mean(sm.nep.pos(:,[1 2]));
    allSynDists=sqrt((allSynPos(:,1)-meanXY(1)).^2+(allSynPos(:,2)-meanXY(2)).^2);
    includeSynBool=synIPLdepths<highcut&synIPLdepths>lowcut;
    inputSynIDs=find(allSynEdges(:,1)==curVGC&includeSynBool==1);
    outputSynIDs=find(allSynEdges(:,2)==curVGC&includeSynBool==1);
    ACinIDs=find(allSynEdges(:,1)==curVGC & preType(:)==8);
    BPCinIDs=find(allSynEdges(:,1)==curVGC & preType(:)==7);
    ACoutIDs=find(allSynEdges(:,2)==curVGC & postType(:)==8);
    RGCoutIDs=find(allSynEdges(:,2)==curVGC & postType(:)==1);
    
    %get the distances of the syns to the center of the cell
    
    
    
    distMat={allSynDists(ACinIDs),allSynDists(BPCinIDs),allSynDists(ACoutIDs),allSynDists(RGCoutIDs)};
    %figure();
    %hold on
    subplot(5,2,i*2-1);
    hold on
    title([num2str(curVGC) 'input'])
    for q=1:2
        
        vgcDistDat(i).colNames=['ACin','BPCin','ACout','RGCout'];
        curHist=histogram(distMat{q},'BinEdges',binVec);
        vgcDistDat(i).counts{q}=curHist.BinCounts;
    end
    %figure();
    %hold on
    subplot(5,2,i*2);
    hold on
    title([num2str(curVGC) 'output'])
    for q=3:4
        
        vgcDistDat(i).colNames=['ACin','BPCin','ACout','RGCout'];
        curHist=histogram(distMat{q},'BinEdges',binVec);
        vgcDistDat(i).counts{q}=curHist.BinCounts;
    end
    
end

%% get vgc skel dat
inBin=zeros(length(binVec),1);
indInBin=inBin;
for i=[1 2 3 4 5]%length(vgcCidList)
    sm=skelDat(i).SM.sm;
meanPos = mean(cat(3,sm.nep.pos(sm.nep.edges(:,1),:),sm.nep.pos(sm.nep.edges(:,2),:)),3);
%%translate meanPos to new analysis frame
meanPosAdj = meanPos;
for i=1:length(meanPos(:,1))
 [g,x,meanPosAdj(i,3)] = getIPLdepth(meanPos(i,3),meanPos(i,1),meanPos(i,2),GCLplane,INLplane);
end
eLengths = sm.nep.props.edgeLength;
for i = 1:length(binVec)-1
   isBin = (meanPosAdj(:,3)>=binVec(i)) & (meanPosAdj(:,3)<binVec(i+1)) ;
   indInBin(i) = sum(eLengths(isBin));;
   inBin(i) = inBin(i)+sum(eLengths(isBin));
end
end

%% pooled results
allSynDists=sqrt((allSynPos(:,1)-meanXY(1)).^2+(allSynPos(:,2)-meanXY(2)).^2);
ACtoVGC=find(ismember(allSynEdges(:,1),vgcCidList)&preType(:)==8&includeSynBool==1);
BPCtoVGC=find(ismember(allSynEdges(:,1),vgcCidList)&preType(:)==7&includeSynBool==1);
VGCtoAC=find(ismember(allSynEdges(:,2),vgcCidList)&postType(:)==8&includeSynBool==1);
VGCtoRGC=find(ismember(allSynEdges(:,2),vgcCidList)&postType(:)==1&includeSynBool==1);

    
h1 = histogram(allSynDists(ACtoVGC),'BinEdges',binVec);
h1v = h1.Values;
h2 = histogram(allSynDists(BPCtoVGC),'BinEdges',binVec);
h2v = h2.Values;
h3 = histogram(allSynDists(VGCtoAC),'BinEdges',binVec);
h3v = h3.Values;
h4 = histogram(allSynDists(VGCtoRGC),'BinEdges',binVec);
h4v = h4.Values;
h5 = histogram(allSynDists(BPCtoVGC),'BinEdges',binVec);
h5v = h5.Values;
h6 = histogram(allSynDists(ACtoVGC),'BinEdges',binVec);
h6v = h6.Values;
h7 = histogram(allSynDists(VGCtoRGC),'BinEdges',binVec);
h7v = h7.Values;
h8 = histogram(allSynDists(VGCtoAC),'BinEdges',binVec);
h8v = h8.Values;

figure();
subplot(2,2,1);
hold on
title('INPUT Synapse distance from arbor center across all VGC');

bar(binVec(1:length(h1v)),h1v);
bar(binVec(1:length(h2v)),h2v);
%plot(binVec,inBin/10,'LineWidth',2)

subplot(2,2,2);
hold on
title('OUTPUT Synapse distance from arbor center across all VGC');

bar(binVec(1:length(h3v)),h3v);
bar(binVec(1:length(h4v)),h4v);

subplot(2,2,3);
hold on
title('INPUT E/I ratio for all VGCs');
curRatio = (h5v)./(h6v);
curRatio(~isfinite(curRatio))=0;
bar(binVec(1:length(h5v)),(h5v)./(h6v));
%bar(binVec,h2)
%plot(binVec,inBin/10,'LineWidth',2)

subplot(2,2,4);
hold on
title('OUTPUT E/I ratio for all VGCs');
curRatio = (h7v)./(h8v);
curRatio(~isfinite(curRatio))=0;
bar(binVec(1:length(h7v)),curRatio);
%bar(binVec,h2)
%plot(binVec,inBin/10,'LineWidth',2)

%% test plot to make sure that the plane fitting is working.
figure()
hold on
indBin=zeros(length(binVec),1);
for i=1:5 %length(vgcCidList)
    sm=skelDat(i).SM.sm;
    meanPos = mean(cat(3,sm.nep.pos(sm.nep.edges(:,1),:),sm.nep.pos(sm.nep.edges(:,2),:)),3);
    %%translate meanPos to new analysis frame
    meanPosAdj = meanPos;
    for i=1:length(meanPos(:,1))
        [g,x,meanPosAdj(i,3)] = getIPLdepth(meanPos(i,3),meanPos(i,1),meanPos(i,2),GCLplane,INLplane);
    end
    eLengths = sm.nep.props.edgeLength;
    for i = 1:length(binVec)-1
        isBin = (meanPosAdj(:,3)>=binVec(i)) & (meanPosAdj(:,3)<binVec(i+1)) ;
        indBin(i) = sum(eLengths(isBin));
    end
    plot(binVec,indBin,'LineWidth',2)
end
legend('2','3','4','5','13')
