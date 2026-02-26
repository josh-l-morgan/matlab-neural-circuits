%function[] = springModel(uselist,modPar);



% %%
% F = -kX
% k = stiffness
% X = proportional to distance;

%% Get data
loadData = 0;
if loadData
    MPN = GetMyDir;
    load([MPN 'obI.mat']);
    seedList = [108 201];
    useList = obI2cellList_seedInput(obI,seedList);
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
end

springDir = 'D:\LGNs1\Analysis\springDat\'
%mkdir(springDir)

for rerun = 1:20

%% set variables
k = 1;
noise = .1;
damp = .5;
disperse = 5;
reps = 2000;
fsize = 200;
speedMax = 1;
centerSpring = 2;
imageFreq = 100;

startRepulse = reps/2;%
repulse = 5; %1
minDist = 5;

%% set fixed
wind = [0 fsize; 0 fsize];
center = [fsize/2 fsize/2];

%% configure nodes

postIDs = useList.postList;
preIDs = useList.preList;
nodeIDs = [preIDs postIDs];

nodeNum = length(nodeIDs);

postSynPref = seedPref.sharedSynNorm(1,:)./sum(seedPref.sharedSynNorm,1);
preSynPref = seedPref.ax2seed(1,:)./sum(seedPref.ax2seed,1);
synPref = [preSynPref postSynPref];
cellClass = [preIDs * 0 + 1   postIDs * 0 + 2];
isPre = cellClass ==1;
isPost = cellClass == 2;

rawCon = zeros(nodeNum,nodeNum);
useCon = rawCon * 0;
for i = 1:length(nodeIDs)
    for p = 1:length(nodeIDs)
        rawCon(i,p) = sum((allEdges(:,1) == nodeIDs(i)) & (allEdges(:,2) == nodeIDs(p)));
        useCon(i,p) = (p>i);
    end
end
con  = rawCon;

testDat = 0;
if  testDat
    con = [ 0 1; 0 0]
    useCon = con;
    nodeNum = size(con,1);
end


nodeX = rand(nodeNum,1)*fsize/2+fsize/4;
nodeY = rand(nodeNum,1)*fsize/2+fsize/4;
nodeXV = zeros(nodeNum,1);
nodeYV = zeros(nodeNum,1);

[ey ex] = find((con.*useCon)>0);
ev = con((con.*useCon)>0);

nodeCol(:,1) = synPref;
nodeCol(:,3) = 1-synPref;
set(0,'DefaultAxesColorOrder',nodeCol)


%%Zero out seeds
zeroSeeds = 1;
if zeroSeeds
    for s = 1:length(seedList);
        targ = find(nodeIDs == seedList(s))
        con(targ,:) = 0;
        con(:,targ) = 0;
    end
    ev = con((rawCon.*useCon)>0);
    
end

%con = con/max(con(:));

%%%%%
for r = 1:reps
    disp(sprintf('running %d of %d',r,reps))
    
    
    preXmat = repmat(nodeX,[1,nodeNum]);
    postXmat = repmat(nodeX',[nodeNum,1]);
    
    preYmat = repmat(nodeY,[1,nodeNum]);
    postYmat = repmat(nodeY',[nodeNum,1]);
    
    %%Independent dim calc
    difX = postXmat - preXmat;
    difY = postYmat - preYmat;
    
    dists = (sqrt(difX.^2 + difY.^2));
    ratX = difX./((abs(difX)+abs(difY)));
    ratY = difY./((abs(difX)+abs(difY)));
    ratX(~useCon) = 0;
    ratY(~useCon) = 0;
    
    springForce =  dists .* con * k;
    % springForce(springForce>springMax) = springMax;
    
    pullX = ratX .* springForce;
    pullY = ratY .* springForce;
    pullX(isnan(pullX)) = 0;
    pullY(isnan(pullY)) = 0;
    
    sumPullX = sum(pullX,2) - sum(pullX,1)';
    sumPullY = sum(pullY,2) - sum(pullY,1)';
    
    
    
    centDifX = center(1) - nodeX;
    centDifY = center(2) - nodeY;
    centRatX = centDifX./(abs(centDifX)+abs(centDifY));
    centRatY = centDifY./(abs(centDifX)+abs(centDifY));
    centerDist = sqrt(centDifX.^2 + centDifY.^2);
    centerForce = centerDist/fsize .* centerSpring;
    centerForce(centerDist>(fsize/2)) =  centerForce(centerDist>(fsize/2))*100;
    gravX = centRatX .* centerForce;
    gravY = centRatY .* centerForce;
    
    
    
    
    %repulsion = (1./dists .* repulse) + disperse/numNodes;
    %repulsion = (dists<minDist) * repulse + disperse/nodeNum;
    
    if r>=startRepulse
            repulsion = (dists<minDist) * repulse + (fsize./dists) * disperse/nodeNum;
    else
            (fsize./dists) * disperse/nodeNum;
    end
    
    %repulsion(repulsion>repulse) = repulse;
    repulsion(isnan(repulsion)) = 0;
    repulsion(~useCon) = 0;
    pushX = ratX .* repulsion;
    pushY = ratY .*repulsion;
    
    
    
    sumPushX = -sum(pushX,2) + sum(pushX,1)';
    sumPushY = -sum(pushY,2) + sum(pushY,1)';
    
    nodeXV = nodeXV + rand * noise - noise/2;
    nodeYV = nodeYV  + rand * noise - noise/2;
    nodeXV = (nodeXV + sumPullX + sumPushX + gravX);
    nodeYV = (nodeYV + sumPullY+ sumPushY + gravY);
    
    
    %%Regulate speed
    
    nodeXV = nodeXV * damp;
    nodeYV = nodeYV * damp;
    speed = sqrt(nodeXV.^2 + nodeYV.^2);
    speedLim = speed;
    speedLim(speed>speedMax) = speedMax;
    speedRat  = speedLim./speed;
    nodeXV = nodeXV.* speedRat;
    nodeYV = nodeYV.* speedRat;
    
    
    
    
    nodeX = nodeX + nodeXV;
    nodeY = nodeY + nodeYV;
    
    
    %%Position Seeds
    zeroSeeds = 1;
    if zeroSeeds
        for s = 1:length(seedList);
            targ = find(nodeIDs == seedList(s));
            if synPref(targ)>.5
                prefDif = synPref;
            else
                prefDif = 1- synPref;
            end
            meanX = sum(nodeX.*prefDif')/sum(prefDif);
            meanY = sum(nodeY.*prefDif')/sum(prefDif);
            rad  = atan2(meanY-center(1), meanX- center(2));
            newY = sin(rad) * (fsize/2)+center(1);
            newX = cos(rad) * (fsize/2) + center(2);
            
            nodeX(targ) = newX;
            nodeY(targ) = newY;
            
        end
    end
    
    
    %%Display
    showEdge = ~mod(r,imageFreq);
    if showEdge
        
        pl = plot([nodeX(ex) nodeX(ey)]',[nodeY(ex) nodeY(ey)]');
        hold on
        
        ei = sub2ind(size(con),ey,ex);
        ev = con(ei);
        edgeCol = nodeCol(ey,:);
        edgeCol(ev==0,:) = edgeCol(ev==0,:)*.2;
        for p = 1:length(pl)
            set(pl(p),'color',edgeCol(p,:));
        end
   
    sg = scatter(nodeX(isPre),nodeY(isPre),20,nodeCol(isPre,:),'^','filled');
    set(sg,'MarkerEdgeColor','w')
    
    sg = scatter(nodeX(isPost),nodeY(isPost),80,nodeCol(isPost,:),'o','filled');
    set(sg,'MarkerEdgeColor','w')
        
    
    ylim([wind(1,1) wind(1,2)])
    xlim([wind(2,1) wind(2,2)])
    set(gca,'color','k')
    
    hold off
    pause(.1)
     end
end


imageName = sprintf('%sspringRun_%03.0f.png',springDir,rerun);
saveas(gcf,imageName,'png')
%print(gcf, imageName, '-depsc2','-painters')
set(gcf, 'InvertHardCopy', 'off');
print(gcf, imageName, '-dpng')

%datName

end
















