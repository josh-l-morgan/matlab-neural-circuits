function[skelOri] = oriSkeletons(arbor,seedSub);

voxSize = [.512 .512 .480]; % in microns
seedSub = seedSub.*voxSize;
lookRange = [15 30]; % min and max value to look at
subplotY = 4;

%%
    %seedSub = arbor.seed .* voxSize;
    allEdgeMids = [];
    allEdgeLengths = [];
    for b = 1:length(arbor.branches)
       allEdgeMids = cat(1,allEdgeMids,arbor.branches(b).edgeMids);
       allEdgeLengths = cat(1,allEdgeLengths,arbor.branches(b).edgeLengths');
    end
    
    for d = 1:3
       allEdgeMids(:,d) = allEdgeMids(:,d) * voxSize(d); 
    end
    edgePos = cat(2,(allEdgeMids(:,1)-seedSub(:,1)), ...
        (allEdgeMids(:,2)-seedSub(:,2)),(allEdgeMids(:,3)-seedSub(:,3)) );
    
    edgeDist = sqrt(edgePos(:,1).^2 + edgePos(:,2).^2  + edgePos(:,3).^2);
    inRange = (edgeDist>lookRange(1)) & (edgeDist<=lookRange(2));
    
    [Coef,Score,latent,tsquare] = princomp(edgePos);
        skelOri.prinSTD = std(Score);
        skelOri.cartSTD = std(edgePos);
   %edgePos = Score;
    
    
    %% get rads
    binNum = 16;
    radDims = [1 2; 1 3; 2 3];
    radBins = (2*pi)/binNum*(0:1:binNum)-pi;

    edgeRads = allEdgeMids * 0;
    histRad = zeros(length(radBins)-1,3);
    for r = 1:3
        edgeRad = atan2(edgePos(:,radDims(r,1)),edgePos(:,radDims(r,2)));
        edgeRads(:,r) = edgeRad;
         for h = 1:length(radBins)-1
            isAngle = (edgeRad>radBins(h)) & (edgeRad<=radBins(h+1));
                        histRad(h,r) = sum(isAngle.*allEdgeLengths);
% 
%              subplot(1,2,1)
%             scatter(edgePos(isAngle,radDims(r,2)),edgePos(isAngle,radDims(r,1)),'.')
%                      xlim([-100 100]),ylim([-100 100])
%             	subplot(1,2,2)
%             %rose(edgeRad(isAngle),radBins+pi/8)
% 
%             myRose(histRad(:,r),radBins)
%             pause
         end
         
%          subplot(subplotY,3,r)
%          scatter(edgePos(:,radDims(r,2)),edgePos(:,radDims(r,1)),'.','b')
%          hold on
%                   scatter(edgePos(inRange,radDims(r,2)),edgePos(inRange,radDims(r,1)),'.','r')
%         hold off
%          xlim([-100 100]),ylim([-100 100])
%          subplot(subplotY,3,(r+3))
%          %rose(edgeRad(inRange)+pi/8,radBins)
%          myRose(histRad(:,r),radBins)
    end
    
    %%
     dirCount = zeros(3,2);
     useEdge = find(inRange);
        for e = 1:length(useEdge)
            edg = edgePos(useEdge(e),:);
           dimTarg = find(abs(edg) == max(abs(edg)));
           
           if edg(dimTarg)>0
               dirTarg = 1;
           else
               dirTarg = 2;               
           end
           dirCount(dimTarg,dirTarg) = dirCount(dimTarg,dirTarg) + allEdgeLengths(useEdge(e));
        end
    
%     subplot(subplotY,3,[7:9])
%     bar(dirCount(:,1)*-1)
%     hold on
%     bar(dirCount(:,2))
%     hold off
    
    skelOri.dirCount = dirCount;
    skelOri.oriCount = sum(dirCount,2);
    [sortedOri oriOrder] = sort(skelOri.oriCount,'descend');
    skelOri.oriOrder = oriOrder';
    skelOri.sortedOri = sortedOri;
    skelOri.oriRat = [(sortedOri(1)-sortedOri(2))/(sortedOri(1))...
        (sortedOri(1)-sortedOri(3))/(sortedOri(1))];
    skelOri.lengthUsed = sum(allEdgeLengths(useEdge));
    skelOri.totalLength = sum(allEdgeLengths);
    skelOri.lookRange = lookRange;
    
    
    
    %% Classify by model
    
    
    
    
    
    
    
    
    