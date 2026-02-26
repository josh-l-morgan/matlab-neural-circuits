function[] = makePOI()
%% Make points of interest that can be compared to vglut3 population polarities

%%Generate a list of POI that is similar to ROI lists, but whose position
%%matches known depths for polarity and whoes polarity matches that
%%observed in calcium imaging



global tis glob

if 1
    SPN =  [glob.datDir 'Analysis\Data\preproc\'];
    
    load([SPN 'ptDat.mat']);
    load([SPN 'ROI.mat']);
    load([SPN 'ROI2.mat']);
    load([SPN 'SOI.mat']);
    load([SPN 'NOI.mat']);
    load([SPN 'MOI.mat']);
    load('MPN.mat')
    swcDir = [WPN 'swc\'];

    datName = 'Circle_Polarity_Depth_MoveCircle072917.mat'
    popDat = load([SPN datName]); %load data from Kerschensteiner lab

end

testDepths = [0 .21 .29 .37 .44 .51 .60]; %Depths
popPol = [-.18 -.44 -.45 -.325 0.02 0.225 0.23];
poolDist = 5;

runCids = SOI.cids;%[2 3 4 ]; % !!!!!!!!!! add more
smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];

%% Parse pop data
D = popDat.Depth;
P = popDat.Polarity;
G = popDat.PolarityDepthGroup;

for i = 1:max(G)
    medD(i) = median(D(G==i));
    medP(i) = median(P(G==i));
end

cbP = median(P(G==1));
cbD = median(D(G==1));
isIPL = G>1;

ft = fittype( 'smoothingspline' );
[fitresult, gof] = fit( medD(2:end)', medP(2:end)', ft );

x = .2:.01:.8;
y = feval(fitresult,x);

scatter(D,P,'.','k')
hold on
scatter(x,y,'o','b')
scatter(medD,medP,'o','r','filled')
hold off
pause(1)

%% Make POI
clear POI
POI.cids = runCids; % create structure containing grouped rois
POI.poolDist = poolDist; 
POI.testDepths = testDepths;
POI.popPol = popPol;

for i = 1:length(runCids);
    
cid = runCids(i);
disp(sprintf('loading cell %d. %d of %d',cid,i,length(runCids)))

    
    
    %%Get distances between nodes
    fileName = sprintf('sm_cid%d.mat',cid);
    useSM(i) = 1;
    load([smDir fileName]);
    dAll = sm.skel2skel.linDist;
    pos = sm.nep.pos;
       
    %% find close
    
    %%Get depths for nodes
    GCLplane.Parameters = [-1 -0.0524864366532293 0.0261842280564985 10.3134566105665];
    INLplane.Parameters = [-1.0000 -0.0768 0.0327 43.0766];

    [zGCL,zINL,depthRaw] = getIPLdepth(pos(:,3),pos(:,1),pos(:,2),GCLplane,INLplane);
    
    %%Adjust Depth manually
    depth = depthRaw;
    lowMatch = [.3 .21];
    highMatch = [testDepths(end) testDepths(end)];
    oldRange = highMatch(1) - lowMatch(1);
    newRange = highMatch(2) - lowMatch(2);
    depth = depth-lowMatch(1);
    depth = depth* newRange/oldRange;
    depth = depth + lowMatch(2);


    %%Find depths relative to test depths
    depthDif = abs(depth - testDepths);
    minDepth = min(depthDif,[],2);

    checkN = minDepth * 0 + 1;
    useN = [];
    
    for n = 1:length(minDepth)
        
        minD = min(minDepth(checkN>0));
        nextN = find((minDepth==minD) & checkN,1);
        useN = [useN;nextN];
        checkN(nextN) = 0;
        tooClose = dAll(nextN,:) <= poolDist;
        checkN(tooClose) = 0;

        if 1;%isempty(nextN)
            scatter3(pos(:,1),pos(:,2),pos(:,3),'.','k');
            view([0 0])
            hold on
            scatter3(pos(useN,1),pos(useN,2),pos(useN,3),100,'r','filled')
            scatter3(pos(tooClose,1),pos(tooClose,2),pos(tooClose,3),100,'b','lineWidth',2)
            scatter3(pos(nextN,1),pos(nextN,2),pos(nextN,3),200,'r','lineWidth',4)
            hold off
            drawnow
            pause(.1)
            
        end

        if ~sum(checkN)
            break
        end
    end

    POI.c(i).closeNode = useN;
    POI.c(i).depth = depth(useN);

end
    
%%Combine cells

gNum = 0;
for c = 1:length(POI.c)
    for r = 1:length(POI.c(c).closeNode)
        gNum = gNum + 1;
        POI.closeNode(gNum,1) = POI.c(c).closeNode(r);
        POI.depth(gNum,1) = POI.c(c).depth(r);
        POI.roiCids(gNum,1) = POI.cids(c);
    end
end

save([SPN 'POI.mat'],'POI');

%% Get Physiology

%t = fittype( 'smoothingspline' );
%[fitresult, gof] = fit( testDepths(2:end)', popPol(2:end)', ft );
POI.Polarity = feval(fitresult,POI.depth);

lowPol = feval(fitresult,medD(2));
highPol = feval(fitresult,medD(end));

POI.Polarity(POI.Polarity<min(medP)) = min(medP);
POI.Polarity(POI.Polarity>max(medP)) = max(medP);
POI.fitmodel = fitresult;
 
nearP = abs(POI.depth - POI.testDepths);
% minP = min(nearP,[],2);
% isMin = nearP == repmat(minP,[1 size(nearP,2)]);
% tempInd = repmat([1:size(nearP,2)],[size(nearP,1) 1]);
% POI.plane = tempInd(isMin>0);

for i = 1:length(POI.depth)
    nearVec = nearP(i,:);
    POI.planeDist(i,1) = min(nearVec);
    POI.plane(i,1) = find(nearVec == min(nearVec),1);
end

x = .2:.01:.8;
y = feval(fitresult,x);
scatter(D,P,'.','k')
hold on
scatter(x,y,'o','b')
scatter(medD,medP,'o','r','filled')
scatter(POI.depth,POI.Polarity,'m','.')
hold off




%% Set first Polarity to cell body
POI.Polarity(POI.plane==1) = POI.popPol(1);

scatter(POI.depth,POI.plane)



save([SPN 'POI.mat'],'POI');









