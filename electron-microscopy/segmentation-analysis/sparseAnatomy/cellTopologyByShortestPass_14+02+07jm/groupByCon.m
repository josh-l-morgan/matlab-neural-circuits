function[] = downSampCon(allVox)

%% Set variables
dSamp = 10;


%% Get Data
subs = allVox.subs;
for m = 1:3
    subs(:,m) = subs(:,m) - allVox.minMax{1}(m);
end
conMat = allVox.conMat;


%% group data into new mat according to down sample factor.  Mat values = current voxel position in sub list.
groupSub = fix(subs/dSamp);
remainSub = subs-groupSub*dSamp+1;
groupSub = groupSub + 1;

maxGroup = max(groupSub,[],1);
groupInd = sub2ind(maxGroup,groupSub(:,1),groupSub(:,2),groupSub(:,3)); %index of voxel group (down samp space)
voxInd = sub2ind([dSamp dSamp dSamp],remainSub(:,1),remainSub(:,2),remainSub(:,3)); %index of voxel within group

maxGrouped = [max(groupInd) max(voxInd)];
uGroupInd = unique(groupInd);

numGroup = length(uGroupInd);
lookupInd(uGroupInd) = 1:numGroup;
groupPos = lookupInd(groupInd);  %get new position of each voxel within new group ordering

groupVoxMat = zeros(numGroup,dSamp^3);
gvMatInd = sub2ind(size(groupVoxMat),groupPos,voxInd');
groupVoxMat(gvMatInd) = 1:length(gvMatInd);



%% Check for multiple objects
lookupGroupVox = zeros(length(groupInd)+1,1);
tic
groupIDs = zeros(numGroup,1);
extraIDs = zeros(numGroup,1);
extraIDnum = 0;
extraGroupVoxMat = groupVoxMat * 0;
for i = 1:size(groupVoxMat,1);
    if ~mod(i,1000)
        disp(sprintf('%d of %d',i,numGroup))
        toc
    end
    groupVox = groupVoxMat(i,:);
    groupVox = groupVox(groupVox>0);
    
    
    lookupGroupVox(groupVox+1) = 1:length(groupVox);
    
    
    groupCon = conMat(groupVox,:);
    groupCon = lookupGroupVox(groupCon+1);
    lookupGroupVox(groupVox+1) = 0;
    
    
    if length(groupCon(:)) > 26
        vProp = rand(size(groupVox));
        reps = 10;
        
        [passProp] = quickPassMax(groupCon,vProp);
        uSeg = unique(passProp);
        
        if length(uSeg)==1
            groupIDs(i) = i;
            
        elseif length(uSeg)>1
            for us = 1:length(uSeg)
                extraIDnum = extraIDnum + 1;
                goodVox = groupVox(passProp == uSeg(us));
                extraGroupVoxMat(extraIDnum,1:length(goodVox)) = goodVox;
                extraGroupIDs(extraIDnum) = i + us/(length(uSeg)+1);
            end
        end
        
        %
        %         testVox.conMat = groupCon;
        %         testVox.conDat = allVox.conDat;
        %         testVox.name = 'hi'
        
        %         tic
        %         [surfClose2node] = passMaxAndDist(testVox,vProp,reps);
        %         toc
    end
end
toc

allGroupVoxMat = cat(1,groupVoxMat(groupIDs>0,:),extraGroupVoxMat(1:extraIDnum,:));
allGroupIDs = [groupIDs(groupIDs>0); extraGroupIDs(1:extraIDnum)'];


%% Get subs of nodes

lookupSubs = cat(1,[0 0 0],subs);
groupVoxSubs = zeros(size(allGroupVoxMat,1),size(allGroupVoxMat,2),3);
grabSub = allGroupVoxMat * 0;

for i = 1:3
    grabSub(:)= lookupSubs(allGroupVoxMat+1,i) ;
    groupVoxSubs(:,:,i) = grabSub;
end

numVoxInGroup = sum(allGroupVoxMat>0,2);
meanSub = squeeze(sum(groupVoxSubs,2) ./ repmat(numVoxInGroup,[1 1 3]));

%% Get connectivity

conMat = cat(1, zeros(1,26), conMat);

for i = 1: 26
   
    grabCon =  conMat(allGroupVoxMat+1,:);
    
end








