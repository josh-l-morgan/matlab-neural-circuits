
%% load data
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'
TPN = [MPN 'skel\'];
load([TPN 'cellArbors'])
load([MPN 'obI.mat'])

%% Define variables
lookDist = 10;
gridSize = 10;

gridArbor.lookDist = lookDist;
gridArbor.gridSize = gridSize;

%% Get syn data and select synapses
targCell = 108;
synapses = obI.nameProps.edges;
preWithTarg = unique(synapses(synapses(:,1)==targCell,2));
preWithTarg = preWithTarg(preWithTarg>0);
synWithPre = [];
for i = 1:length(preWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==preWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);

synPos = ceil(obI.colStruc.anchors(synObj,:)/8);
synPos(:,1:2) = synPos(:,1:2); %%% Error in tif to sub assignment

gridArbor.targCell  = targCell;

%% find closest point and total length within 2x dist of closest point

preList = unique(synPre);
preList = preList(preList>0);
postList = unique(synPost);
postList = postList(postList>0);

goodCheck = synObj*0 + 1;
showPlots = 0;
synMat = zeros(length(preList),length(postList));
lengthMat = zeros(length(preList),length(postList),lookDist);

gridArbor.preList = preList;
gridArbor.postList = postList;


gridSyn = ceil(synPos/gridSize);
gridSyn(gridSyn==0) = 1;  %% error in anchor
maxGridSynPos = max(gridSyn,[],1)+gridSize;
synInd = sub2ind(maxGridSynPos,gridSyn(:,1),gridSyn(:,2),gridSyn(:,3));

for pre = 1 : length(preList)
    
    disp(sprintf('running pre %d of %d',pre,length(preList)))
    preCell = find(cellArbors.cellName == preList(pre));
    arbor1 = cellArbors.arbors(preCell);
    
    for post = 1: length(postList)
        
        postCell = find(cellArbors.cellName == postList(post));
        arbor2 = cellArbors.arbors(postCell);
        
        %%Compare arbors 
        if (length(arbor1)>0) & (length(arbor2)>0)

        mid1 = cat(1,arbor1.branches.edgeMids);
        mid2 = cat(1,arbor2.branches.edgeMids);
        
         mid1(:,1:2) = mid1(:,1:2) + 1024/8;
         mid2(:,1:2) = mid2(:,1:2) + 1024/8;
        
        else
            mid1 = [];
            mid2 = [];
        end
        sameGrid = countGrids(mid1,mid2,gridSize);
        trackArborOverlap(pre,post).subs = sameGrid.overlapSubs;
        touchMat(pre,post) = sameGrid.gridCount;
        
        %%Compare synapses
        rightPair = synPre==preList(pre) & synPost == postList(post);
        rightInd = synInd(rightPair);
        uniqueInd = unique(rightInd);
        
        if length(uniqueInd)>1
            histInd = hist(rightInd,uniqueInd);
        else
            histInd = length(rightInd);
        end
        
        [sY sX sZ] = ind2sub(maxGridSynPos,uniqueInd);
        
        synTrack(pre,post).subs = [sY sX sZ];
        synTrack(pre,post).synCount = histInd;
        
        if max(histInd)>1
           histInd
            pause(.1)
        end
        synMat(pre,post) = length(uniqueInd);  % how many unique synapse positions are there

    end
end

gridArbor.synMat = synMat;
gridArbor.touchMat = touchMat;
gridArbor.synTrack = synTrack;
gridArbor.arborTrack = trackArborOverlap;


%%

useVals = touchMat>0;
Xdat = touchMat(useVals);
Ydat = synMat(useVals);
[rho pc] = corr(Xdat,Ydat)

scatter(Xdat,Ydat,'.','r')
hold on

maxX = max(Xdat);
X = [0 maxX];
[fit1 gof] = fit(Xdat,Ydat,'poly1');
Y = X* fit1.p1 + fit1.p2;
line(X ,Y);
hold off

fitDat.Xdat = Xdat;
fitDat.Ydat = Ydat;
fitDat.fit1 = fit1;
fitDat.gof = gof;
fitDat.rho = rho;
fitDat.corrP = pc;

gridArbor.fitDat = fitDat;


%% Show syn vs length

colMat = synMat*256/max(synMat(:));
colMat(:,:,2) = touchMat * 256/max(touchMat(:));
colMat(:,:,3) = touchMat * 0;

image(uint8(colMat))

%%
save([TPN 'gridArbor.mat'],'gridArbor')


