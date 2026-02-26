


%MPN = GetMyDir;
load('MPN.mat')
springRes = 'D:\LGNs1\Analysis\springDat\results\';

load([springRes 'res_snapLong01_randClade1.mat'])
randCladeResults = cladeResults;

load([springRes 'res_snapLong01_claded2.mat'])
realCladeResults = cladeResults;

nodeOrder = cladeResults.nodeOrder;
relMat = cladeResults.relMat;



%% node props
nodeIDs = nodeOrder;
load([MPN 'obI.mat'])
seedList = [ 108 201 109 907 903];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
synEdges = obI.nameProps.edges(:,[2 1]);

checkIDs = useList.preList;

seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .0 1; 0 0 1];
checkCol = zeros(length(checkIDs),3);
for s = 1:length(seedList)
    useCol = seedCol(s,:);
    for n = 1:length(checkIDs)
        checkCol(n,:) = checkCol(n,:) + useCol *  ...
            double(sum((synEdges(:,1)==checkIDs(n)) & (synEdges(:,2)) == seedList(s))>0);
    end
end

checkProp = checkCol;

%% find clade differences

realDifRes = cladePropDif(realCladeResults.cladedNodes,checkIDs,checkProp);
randDifRes =cladePropDif(randCladeResults.cladedNodes,checkIDs,checkProp);

histBin = [0:.2:2];
realHist = hist(realDifRes.useDifs,histBin);
randHist = hist(randDifRes.useDifs,histBin);

bar(histBin,[realHist/sum(realHist); randHist/sum(realHist)]')











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





















