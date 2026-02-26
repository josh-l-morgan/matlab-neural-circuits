load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\tis.mat');

fvDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\fvLibrary\';
plotDest=[fvDir '\dat\'];
if ~exist(plotDest,'dir')
    mkdir(plotDest)
end

if exist([fvDir 'ref_gcl nucEdge.mat'],'file')
    ipl_bord_GCL = load([fvDir 'ref_gcl nucEdge.mat']);
    ipl_bord_INL = load([fvDir 'ref_inl nucEdge.mat']);
    GCLbord=ipl_bord_GCL.fv.vertices(:,:);
    INLbord=ipl_bord_INL.fv.vertices(:,:);
else
    GCLbord = [0 0 0; 100 0 0; 0 100 0; 100 100 0];
    INLbord = [0 0 100; 100 0 100; 0 100 100; 100 100 100];
end

Locs={GCLbord;INLbord}; %These are in z,x,y I think.
for i=1:2
    P=Locs{i};
    B(:,i) = [P(:,3), P(:,2), ones(size(P,1),1)] \ P(:,1);
end
GCLplane=struct();
INLplane=struct();
GCLplane.Parameters=[-1 B(2,1) B(1,1) B(3,1)];
INLplane.Parameters=[-1 B(2,2) B(1,2) B(3,2)];

synIPLdepths=zeros(length(tis.syn.pos(:,1)),1);
for j=1:length(synIPLdepths)
    curPos=tis.syn.pos(j,:);
    [gcl,inl,curZ]=getIPLdepth(curPos(3),curPos(1),curPos(2),GCLplane,INLplane);
    synIPLdepths(j)=curZ;
end

vgcCids=[2 3 4 5 10 11 13 14 20];

%% making basic syn dist figures
if 0
     figure(); histogram(synIPLdepths,100); xlim([0 1]);
     
     %input vs output
     inputDepths=synIPLdepths(ismember(tis.syn.edges(:,1),vgcCids));
     outputDepths=synIPLdepths(ismember(tis.syn.edges(:,2),vgcCids));
     figure(); hold on; histogram(inputDepths,50); histogram(outputDepths,50);
     
     %excitatory vs inhibitory
     excInputDepths=synIPLdepths(ismember(tis.syn.edges(:,1),vgcCids)&tis.syn.preClass==8);
     InhInputDepths=synIPLdepths(ismember(tis.syn.edges(:,1),vgcCids)&tis.syn.preClass==7);
     figure(); hold on; histogram(InhInputDepths,50); histogram(excInputDepths,50);
     
     %VG3 outputs to different types
     downstreamRGCcids=tis.syn.edges(ismember(tis.syn.edges(:,2),vgcCids),1);
     downstreamTypeInfo=cid2type(downstreamRGCcids,tis);
     dsTypes=downstreamTypeInfo{1};
     dsSubTypes=downstreamTypeInfo{3};
     rgcSubTypes=dsSubTypes(dsTypes==1);
     rgcSubTypeNameList=tis.cells.type.subTypeNames{1};
     figure(); histo=histogram(rgcSubTypes,'BinEdges',[0:56]); xticks([0:56]), xticklabels(rgcSubTypeNameList);
end




global tis

vgcCidList=[2 3 4 5 13 14 20];

fullCidList=tis.cells.cids(:);
fullTypeList=tis.cells.type.typeID(:);
allSynEdges=tis.syn.edges(:,[1 2]);

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

binVec=[0:0.02:1];
vgcHistDat=struct;
for i=[1 2 3 4 5]
    curVGC=vgcCidList(i);
    inputSynIDs=find(allSynEdges(:,1)==curVGC);
    outputSynIDs=find(allSynEdges(:,2)==curVGC);
    ACinIDs=find(allSynEdges(:,1)==curVGC & preType(:)==8);
    BPCinIDs=find(allSynEdges(:,1)==curVGC & preType(:)==7);
    ACoutIDs=find(allSynEdges(:,2)==curVGC & postType(:)==8);
    RGCoutIDs=find(allSynEdges(:,2)==curVGC & postType(:)==1);
    depthMat={synIPLdepths(ACinIDs),synIPLdepths(BPCinIDs),synIPLdepths(ACoutIDs),synIPLdepths(RGCoutIDs)};
    figure();
    hold on
    title([num2str(curVGC) 'input'])
    for q=1:2
        
        vgcHistDat(i).colNames=['ACin','BPCin','ACout','RGCout'];
        curHist=histogram(depthMat{q},'BinEdges',binVec);
        vgcHistDat(i).counts{q}=curHist.BinCounts;
    end
    figure();
    hold on
    title([num2str(curVGC) 'output'])
    for q=3:4
        
        vgcHistDat(i).colNames=['ACin','BPCin','ACout','RGCout'];
        curHist=histogram(depthMat{q},'BinEdges',binVec);
        vgcHistDat(i).counts{q}=curHist.BinCounts;
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

%%
ACtoVGC=find(ismember(allSynEdges(:,1),vgcCidList)&preType(:)==8);
BPCtoVGC=find(ismember(allSynEdges(:,1),vgcCidList)&preType(:)==7);
VGCtoAC=find(ismember(allSynEdges(:,2),vgcCidList)&postType(:)==8);
VGCtoRGC=find(ismember(allSynEdges(:,2),vgcCidList)&postType(:)==1);

figure();
hold on
h1 = histogram(synIPLdepths(ACtoVGC),'BinEdges',binVec)
h2 = histogram(synIPLdepths(BPCtoVGC),'BinEdges',binVec)
bar(binVec,h1)
bar(binVec,h2)
plot(binVec,inBin/10,'LineWidth',2)

figure();
hold on
histogram(synIPLdepths(VGCtoAC),'BinEdges',binVec)
histogram(synIPLdepths(VGCtoRGC),'BinEdges',binVec)
bar(binVec,h1)
bar(binVec,h2)
plot(binVec,inBin/10,'LineWidth',2)

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


%%


%     
% %get the depth of a point in the IPL from the fitted planes
% function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCLplane,INLplane)
% zGCL=(-x*GCLplane.Parameters(2)-y*GCLplane.Parameters(3)-GCLplane.Parameters(4))/GCLplane.Parameters(1);
% zINL=(-x*INLplane.Parameters(2)-y*INLplane.Parameters(3)-INLplane.Parameters(4))/INLplane.Parameters(1);
% 
% IPLdepth=abs(z-zINL)/abs(zINL-zGCL);
% if eyewireFit==1
%     IPLdepth=IPLdepth*m+b;
% end
% end