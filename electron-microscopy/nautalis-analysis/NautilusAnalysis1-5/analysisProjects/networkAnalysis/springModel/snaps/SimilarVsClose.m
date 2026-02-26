


springRes = 'D:\LGNs1\Analysis\springDat\results\';
%load('D:\LGNs1\Analysis\springDat\results\res_snapLongNoSeedsAll_2.mat')
load([springRes 'res_snapLong01_baseClade1.mat'])

results  = cladeResults.results;
cellGroups = results.cellGroups;
nodeIDs = results.nodeIDs;



%% choose property list to test

%     [colorList colPropRaw] = getList_pcLats;
%     colProp = [colPropRaw(:,2)./colPropRaw(:,1) colPropRaw(:,3)./colPropRaw(:,2)];
%    % colProp = [colPropRaw(:,3)./colPropRaw(:,2) ];
%  
    
%     [colorList colPropRaw] =getList_cardShape
%     colProp = colPropRaw(:,2);
%     
%     
%     [colorList colPropRaw] =getList_cardSholl
%     colProp = colPropRaw(:,2);
    
   [shapeList colPropRaw] = getList_sholl;
   colProp = colPropRaw(:,2);
   shapeProp = colPropRaw(:,2);
   
    [posList cellPos] = getList_cellPositions;
    posProp = cellPos(:,1);
    
    
% %    
%    [colorList colPropRaw] = getList_meanDiam;
%    colProp = colPropRaw;
   
    %[colorList colProp] = getList_meanDiam();
    
   
   
    [colorList ia ib] = intersect(posList,shapeList);
    posProp = posProp(ia);
    shapeProp = shapeProp(ib);
    cellPos = cellPos(ia,:);
    scatter(posProp,shapeProp)
    
 %% Get similarity and relatedness matricies
 colormap(bluered(100))
 simMat = zeros( length(colorList));
    relMat = simMat;
 
   for y = 1:length(colProp)
       for x = 1:length(colProp);
           targY = y;
           targX = x;
           if sum(isSame)
            relMat(y,x) = sqrt((cellPos(targY,1)-cellPos(targX,1)).^2 +...
                (cellPos(targY,2)-cellPos(targX,2)).^2 +...
                (cellPos(targY,3)-cellPos(targX,3)).^2); 
           end
           
           isSame = cellGroups(:,y) == cellGroups(:,x);
           
           difSqr = 0;
           for s = 1:size(colProp,2)
               difSqr = difSqr+(colProp(y,s) - colProp(x,s)).^2;
           end
           difSqr = sqrt(difSqr);
           
           
           simMat(y,x) =difSqr;
       end
   end
   simMat = abs(simMat);
   simMat = 1-(simMat./max(simMat(:)));
   %relMat = max(relMat(:))-relMat;
   relMat = relMat./max(relMat(:));
   
   subplot(3,1,1)
   image(relMat*100)
   subplot(3,1,2)
        image(simMat*50)    
        
%[shuffleGraph,newOrder] = reshuffleSimilarity(relMat);

%% correlation
dmat = diag(1:size(relMat,1))==0;
[RR PP]= corrcoef(relMat(dmat),simMat(dmat));
realCC = RR(2,1)
ccP = PP(2,1)

reps = 10000;
for r = 1:reps
   pick = randperm(size(simMat,1)); 
   newRelMat = relMat(pick,:);
   newRelMat = newRelMat(:,pick);
   [R P] = corrcoef(newRelMat(dmat),simMat(dmat));
    randCC(r) = R(2,1);
end
subplot(3,1,3)
hist(randCC)
hold on
scatter(realCC,1,'r')
hold off
meanRand = mean(randCC)
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





















