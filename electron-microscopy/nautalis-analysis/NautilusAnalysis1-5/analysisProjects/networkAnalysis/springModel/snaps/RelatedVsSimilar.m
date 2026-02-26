


springRes = 'D:\LGNs1\Analysis\springDat\results\';
%load('D:\LGNs1\Analysis\springDat\results\res_snapLongNoSeedsAll_2.mat')
load([springRes 'res_snapLong01_baseClade1.mat'])

results  = cladeResults.results;
cellGroups = results.cellGroups;
nodeIDs = results.nodeIDs;



%% choose property list to test

    [colorList colPropRaw] = getList_pcLats;
    colProp = [colPropRaw(:,2)./colPropRaw(:,1) colPropRaw(:,3)./colPropRaw(:,2)];

   [colorList colPropRaw] = getList_sholl;
   colProp = colPropRaw(:,2);
   
%    
%    [colorList colPropRaw] = getList_meanDiam;
%    colProp = colPropRaw;
   
    %[colorList colProp] = getList_meanDiam();
    
   
   
    useList = colorList*0;
    for i = 1:length(colorList)
        if sum(nodeIDs==colorList(i))
            useList(i) = 1;
         end
    end
    percentWithConnections = sum(useList)/length(useList)
    
    colorList = colorList(useList>0);
    colProp = colProp(useList>0);
    
    
 %% Get similarity and relatedness matricies
 colormap(bluered(100))
 simMat = zeros( length(colorList));
    relMat = simMat;
 
   for y = 1:length(colProp)
       for x = 1:length(colProp);
           targY = find(nodeIDs == colorList(y));
           targX = find(nodeIDs == colorList(x));
           isSame = cellGroups(:,targY) == cellGroups(:,targX);
           if sum(isSame)
            relMat(y,x) = max(find(isSame));
           end
           
           isSame = cellGroups(:,y) == cellGroups(:,x);
           simMat(y,x) = abs(colProp(y) - colProp(x));
       end
   end
   simMat = abs(simMat);
   simMat = simMat./max(simMat(:));
   relMat = max(relMat(:))-relMat;
   relMat = relMat./max(relMat(:));
   
   subplot(3,1,1)
   image(relMat*100)
   subplot(3,1,2)
        image(simMat*50)    
        
%[shuffleGraph,newOrder] = reshuffleSimilarity(relMat);

%% correlation
RR = corrcoef(relMat(:),simMat(:));
realCC = RR(2,1)

reps = 1000;
for r = 1:reps
   pick = randperm(size(simMat,1)); 
   newRelMat = relMat(pick,:);
   newRelMat = newRelMat(:,pick);
   [R P] = corrcoef(newRelMat(:),simMat(:));
    randCC(r) = R(2,1);
end
subplot(3,1,3)
hist(randCC)
hold on
scatter(realCC,1,'r')
hold off
mean(randCC)
P = sum(randCC>=realCC)/length(randCC)


if 0
%% Plot vs
subplot(1,1,1)
scatter(relMat(:),simMat(:))

binRel = [1:50:300];
binSim = [0:.1:1];

clear histSim
for i = 1:length(binRel)-1
   useRel = (relMat>=binRel(i) ) & ( relMat<binRel(i+1));
   useSim = simMat(useRel>0); 
   histSim(i,:) = hist(useSim,binSim);
   subplot(length(binRel)-1,1,i)
   bar(histSim(i,:))
end
end






%%  Bar clades
load('.\data\clade\cladePick_six2.mat')
    
    
    
    colorTable = [ .3 .3 1; 1 1 0; 1 0 0; 1 1 1; 1 0 1 ; 0 1 0]
    col = [];
    allMembers = [];
    for i = 1:length(cladePick.members)
        members = cladePick.members{i};
        [members ia] = intersect(colorList,members);
        memProp{i} = colProp(ia,:);
%         allMembers = [allMembers members];
%         col = cat(1,col,repmat(colorTable(i,:),[length(members),1]));
        
    end

poolA = cat(1,memProp{[1, 2, 3]});
poolB = cat(1,memProp{[4 5 6]});

ranksum(poolA,poolB)
median(poolA)
rangeX(poolA)
median(poolB)
rangeX(poolB)












