


mpSize = matlabpool('size');
if ~mpSize
    'opening matlab pool'
    matlabpool close force
    matlabpool 
end

%% load cladeNetworks
load('D:\LGNs1\Analysis\cladistic\cladeNodeResults\cn_25_700.mat')


%% generic
MPN ='D:\LGNs1\Analysis\segmentations\vastMip3\collectVastSubs_mat\';

%{
  load([MPN 'cb2d.mat'])
    colorList = cb2d.IDs;
    colPropRaw = cb2d.areaUM;
        
    mtaCount = getList_mtaCount();
    colorList = mtaCount(:,1);
    colPropRaw = mtaCount(:,2);
    
    [colorList colPropRaw] = getList_giantBoutonsCounts(MPN);
    [colorList colPropRaw] = getList_biconality;
    [colorList2 colPropRaw2] = getList_giantBoutonsCounts(MPN);
[checkIDs checkCol]  = getList_axSeedCon();

%}
[checkIDs checkProp] 

    
    checkIDs = colorList;
checkProp = colPropRaw(:);

%% find clade differences





propRange = [min(checkProp(:)) max(checkProp(:))];
propSpread =  max(checkProp(:)) - min(checkProp(:));
histBin = [0:propSpread/10:propSpread];

%histBin = [0:.2:2];

for i = 1:length(realN)
    disp(sprintf('real %d of %d',i,length(realN)));
    realDifRes = cladePropDif(cn.real(1).cladeNodes,checkIDs,checkProp);
    useDifs = realDifRes.useDifs;
    realMeans(i) = mean(useDifs);
    realHist(i,:) = hist(useDifs,histBin);
end

meanRealHist = mean(realHist,1);

for i = 1:length(randN)
        disp(sprintf('rand %d of %d',i,length(randN)));
    randDifRes =cladePropDif(cn.rand(i).cladeNodes,checkIDs,checkProp);
    useDifs = randDifRes.useDifs;
    randMeans(i) = mean(useDifs);
    randHist(i,:) = hist(useDifs,histBin);
end

meanRandHist = mean(randHist,1);

%%
subplot(2,1,1)
bar(histBin,[meanRealHist/sum(meanRealHist); meanRandHist/sum(meanRandHist)]')


subplot(2,1,2)
histBin2 = [0 : .001 : 2];
histRandMeans = hist(randMeans,histBin2);
histRealMeans = hist(realMeans,histBin2);

bar(histBin2,[histRealMeans/sum(histRealMeans); histRandMeans/sum(histRandMeans)]')














% 
% 
% 
% 
% %% get nodeProp
% nodeProp = zeros(length(nodeOrder),3);
% nodeIDs = nodeOrder;
% 
% seedList = [ 108 201 109 907 903];
% useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
% synEdges = obI.nameProps.edges(:,[2 1]);
% 
% seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .0 1; 0 0 1];
% for s = 1:length(seedList)
%     useCol = seedCol(s,:);
%     for n = 1:length(nodeOrder)
%         nodeCol(n,:) = nodeCol(n,:) + useCol *  ...
%             double(sum((synEdges(:,1)==nodeOrder(n)) & (synEdges(:,2)) == seedList(s))>0);
%     end
% end
% 
% 
% for i = 1:length(nodeIDs)
%     for p = 1:length(nodeIDs)
%         propMat(i,p) = sum(nodeCol(i,:) & nodeCol(p,:));
%         compMat(i,p) = relMat(i,p);
%     end
% end
% 
% %% plot
% 
% rel1 = compMat(propMat>0);
% rel2 = compMat(propMat==0);
% 
% binRel = [0:.1:1];
% 
% histRel1 = hist(rel1,binRel)/length(rel1);
% histRel2 = hist(rel2,binRel)/length(rel2);
% 
% 
% bar(histRel1,'r')
% hold on
% bar(histRel2,'b')
% hold off
% 
% 
% bar([histRel1; histRel2]')
% 
% res.N1 = length(rel1);
% res.N2 = length(rel2);
% res.mean1 = mean(rel1);
% res.mean2 = mean(rel2);
% res.std1 = std(rel1);
% res.std2 = std(rel2);
% res.ranksumP = ranksum(rel1,rel2);
% 
% res





















