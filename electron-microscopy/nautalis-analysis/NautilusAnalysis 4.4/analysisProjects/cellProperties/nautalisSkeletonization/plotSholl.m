function[skelSholl] = plotSholl(arbor,seedSub,dim);

%{
Attempt to replicate Guido sholl analysis

1) compress into 2D
2) Draw 5 concentric circles whose width is defined by outermost dendrites
(15-27 um wide)
3) rings divided into four quadrants
4) quadrants positioned to maximize difference between orientations
5) count number of dendritic intersections passing through each ring
6) DOi = (a1 + a2) / (b1 + b2)

note: It is not clear what the z axis is in the guido data, but it is
probably rostral caudal.  His X cells may be the vertically oriented
rostral caudal cell.



%}



voxSize = [.512 .512 .480]; % in microns
voxSize = [.588 .512 .480]; % in microns adjusted for y compression

%dim = [1 2];
minLength = 3;


useSeed = seedSub(dim);
useSeed = useSeed.*voxSize(dim);


%%
%seedSub = arbor.seed .* voxSize;
allEdgeMids = [];
allEdgeLengths = [];
for b = 1:length(arbor.branches)
    allEdgeMids = cat(1,allEdgeMids,arbor.branches(b).edgeMids(:,dim));
    allEdgeLengths = cat(1,allEdgeLengths,arbor.branches(b).edgeLengths');
    
end
allEdgeMids(:,1) = allEdgeMids(:,1) * voxSize(dim(1)) - useSeed(1);
allEdgeMids(:,2) = allEdgeMids(:,2) * voxSize(dim(2)) - useSeed(2);

scatter(allEdgeMids(:,1),allEdgeMids(:,2),'.','b')
hold on
scatter(0,0,'.','r')
hold off

allDists = sqrt((allEdgeMids(:,1)).^2 ...
    + (allEdgeMids(:,2)).^2);

sortDists = sort(allDists,'ascend');
thresh95 = sortDists(round(length(sortDists)* 0.95));

outerRing = 100; % could use thresh95

showRad = outerRing* 1.2;
shollRings = outerRing .* [1:5]/5;




%%
crossPoint{5} = [];

for b = 1:length(arbor.branches)
    
    
    
    edges = arbor.branches(b).edges;
    
    startEdge = arbor.nodes.node2subs(edges(:,1),dim);
    startEdge(:,1) = startEdge(:,1) * voxSize(dim(1))-useSeed(1);
    startEdge(:,2) = startEdge(:,2) * voxSize(dim(2))-useSeed(2);
    
    stopEdge = arbor.nodes.node2subs(edges(:,2),dim);
    stopEdge(:,1) = stopEdge(:,1) * voxSize(dim(1))-useSeed(1);
    stopEdge(:,2) = stopEdge(:,2) * voxSize(dim(2))-useSeed(2);
    
    
    
    
    edgeLengths = sqrt((startEdge(:,1) - stopEdge(:,1)).^2 ...
        + (startEdge(:,2) - stopEdge(:,2)).^2);
    branchLength(b) = sum(edgeLengths);
    
    
        scatter(allEdgeMids(:,1),allEdgeMids(:,2),15,[.5 .5 .5],'.')
        hold on        
        scatter(0,0,'.','r')
        
        for s = 1:length(shollRings)
            circle(0,0,shollRings(s))
        end
        hold off
    
    
    if  branchLength(b) >= minLength
        
        meanEdge = (startEdge + stopEdge ) /2;
        startDist = sqrt((startEdge(:,1) ).^2 ...
            + (startEdge(:,2) ));
        stopDist = sqrt((stopEdge(:,1)).^2 ...
            + (stopEdge(:,2)));
        
%         scatter(startEdge(:,1),startEdge(:,2),'.','g')
%         scatter(allEdgeMids(:,1),allEdgeMids(:,2),'.','k')
%         hold on        
%         scatter(stopEdge(:,1),stopEdge(:,2),'.','c')
%         scatter(useSeed(1),useSeed(2),'.','r')
        
        for s = 1:length(shollRings)
%             circle(useSeed(1),useSeed(2),shollRings(s))
            
            whichSide = (startDist<=shollRings(s)) + (stopDist<=shollRings(s));
            crossed = find(whichSide == 1,1);
            if ~isempty(crossed)
                meanEdge(crossed,:);
                crossPoint{s} = cat(1,crossPoint{s},meanEdge(crossed,:));
            end
        end
        hold off
        pause(.01)
    end
end

crossPoints = cat(1,crossPoint{:});
if size(crossPoints,1) >0
crossPoints(:,1) = crossPoints(:,1);
crossPoints(:,2) = crossPoints(:,2);

crossRads = atan2(crossPoints(:,1),crossPoints(:,2))+ pi;

%%optimize
binRad = 2 * pi * [0:4]/4;
shiftR = 0 : .011 : (2*pi);
shiftHist = zeros(length(shiftR),length(binRad));
for r = 1:length(shiftR)
    shiftedCross = mod(crossRads + shiftR(r),2*pi);
    shiftHist(r,:) = histc(shiftedCross,binRad);
end

ab = [(shiftHist(:,1) + shiftHist(:,3))  (shiftHist(:,2) + shiftHist(:,4))];
minAB = min(ab,[],2);
maxAB = max(ab,[],2);
DOis = minAB./maxAB;
targDOis = find(DOis == min(DOis),1);

DOi = DOis(targDOis);
bestRot = shiftR(targDOis);



ylim([-showRad showRad])
xlim([-showRad showRad])
resString = sprintf('DOi = %01.3f \nrotation = %.1f rads \narbor = %.1fum',DOi,bestRot,sum(branchLength))
text(-showRad+showRad/10,showRad-showRad/3,resString,'FontWeight','Bold');



skelSholl.DOi = DOi;
skelSholl.arborSize = sum(branchLength);
skelSholl.bestRot = bestRot;
skelSholl.shiftR = shiftR;
skelSholl.bestShift = targDOis;
skelSholl.shiftHist = shiftHist(:,1:4);
skelSholl.ab = ab;
skelSholl.voxSize = voxSize;
skelSholl.useSeed = useSeed;
skelSholl.minLength = minLength;
skelSholl.dim =dim;
skelSholl.shollRings = shollRings;
skelSholl.crossPoint = crossPoint;
skelSholl.crossRads = crossRads;

else
   
    skelSholl.DOi = nan
    
end
