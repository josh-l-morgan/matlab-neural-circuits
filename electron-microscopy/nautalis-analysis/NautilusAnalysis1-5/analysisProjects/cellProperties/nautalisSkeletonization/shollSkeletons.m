function[sholl] = shollSkeletons(cellStruc,seedSubRaw,lookRange);

%voxSize = [.512 .512 .480]; % in microns
%
% anchorScale = [.0184 0.016 0.030];
% voxelScale = [anchorScale(1) * 8  anchorScale(2) * 8  anchorScale(3) * 4];
% vastScale = [1 4/4.6 1];
% pixum = .0184 * 8;
%
% voxSize = voxelScale *2;
arbor = cellStruc.arbor;
voxSize = cellStruc.voxSize;
seedSub = seedSubRaw;

if ~exist('lookRange','var')
    lookRange = [5 10]; % min and max value to look at
end


%% Make bins

distBin = [lookRange(1):diff(lookRange)/4:lookRange(2)];
radBin = [0:pi/2:2*pi];
shiftAng = [0:pi/2/40:pi/2];
%shiftAng = [0];


%%
%seedSub = arbor.seed .* voxSize;
allEdgeMids = [];
allEdgeLengths = [];
branchIn1 = zeros(length(radBin)-1,length(distBin)-1,length(shiftAng));
branchIn2 = branchIn1;
minLength = 3;
clf
hold on
for b = 1:length(arbor.branches)
    allEdgeMids = arbor.branches(b).edgeMids;
    allEdgeLengths = arbor.branches(b).edgeLengths';
    allEdgeLengths = allEdgeLengths / cellStruc.voxSize(1);
    totEdge = sum(allEdgeLengths);
    
    if totEdge >=minLength
        
        for d = 1:3
            allEdgeMids(:,d) = allEdgeMids(:,d) * voxSize(d);
        end
        edgePos = cat(2,(allEdgeMids(:,1)-seedSub(:,1)), ...
            (allEdgeMids(:,2)-seedSub(:,2)),(allEdgeMids(:,3)-seedSub(:,3)) );
        
        edgeDist = sqrt(edgePos(:,1).^2 + edgePos(:,2).^2  + edgePos(:,3).^2);
        
        edgeDist1 = sqrt(edgePos(:,1).^2 + edgePos(:,3).^2 );
        edgeDist2 = sqrt(edgePos(:,2).^2 + edgePos(:,3).^2 );
        
        edgeRad1 = mod(atan2(edgePos(:,1),edgePos(:,3))+pi/4+pi,2*pi);
        edgeRad2 = mod(atan2(edgePos(:,2),edgePos(:,3))+pi/4+pi,2*pi);
        
        scatter(allEdgeMids(:,1),allEdgeMids(:,3),'.','k')
        
        funCol = {'r','g','b','m','y','c'};
        for r = 1:length(radBin)-1
            for d = 1:length(distBin)-1
                for s = 1:length(shiftAng)
                    test1 = mod(edgeRad1+shiftAng(s),2*pi);
                    test2 = mod(edgeRad2+shiftAng(s),2*pi);
                    ang1 = (test1>radBin(r)) & (test1<=radBin(r+1));
                    ang2 = (test2>radBin(r)) & (test2<=radBin(r+1));
                    indist1 = (edgeDist1>distBin(d)) & (edgeDist1<=distBin(d+1));
                    indist2 = (edgeDist2>distBin(d)) & (edgeDist2<=distBin(d+1));
                    branchIn1(r,d,s) = branchIn1(r,d,s) + sum(ang1 & indist1);
                    branchIn2(r,d,s) = branchIn2(r,d,s) + sum(ang2 & indist2);
                end
                if r == 1
                    scatter(allEdgeMids(ang1 & indist1,1),allEdgeMids(ang1 & indist1,3),'.',funCol{d})
                end
            end
        end
    end
end
hold off


sholl.distBin = distBin;
sholl.radBin = radBin;
sholl.branchIn1 = branchIn1;
sholl.branchIn2 = branchIn2;

%% find cardinal sholl

sum1 = squeeze(sum(branchIn1,2));
sum2 = squeeze(sum(branchIn2,2));
sholl.rat = [sum(sum1([2 4],1))/sum(sum1([1 3],1))    sum(sum2([2 4],1))/sum(sum2([1 3],1)) ]  ;



%% find best sholl

allRat = [(sum(sum1([2 4],:))./sum(sum1([1 3],:)))'    (sum(sum2([2 4],:))./sum(sum2([1 3],:)))' ];
min1 = find(allRat(:,1) == min(allRat(:,1)),1);
min2 = find(allRat(:,2) == min(allRat(:,2)),1);

sholl.bestAng = [shiftAng(min1) shiftAng(min2)];
sholl.bestRat = [allRat(min1,1) allRat(min2,2)];





