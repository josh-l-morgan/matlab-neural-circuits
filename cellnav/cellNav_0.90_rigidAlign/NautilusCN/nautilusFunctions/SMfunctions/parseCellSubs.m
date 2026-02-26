function[allVox] = parseCellSubs(allVox)




%% Clean object subs
allVox.name = 'All object voxels';
numVox = size(allVox.subs,1);

%% Find Dimensions
midObj = median(allVox.subs,1);
allVox.mid = midObj;

minAV = min(allVox.subs,[],1);
offsetSubs = minAV-1;
maxAV = max(allVox.subs,[],1);
imageDims = maxAV-offsetSubs + 1;

minMax{1} = offsetSubs;
minMax{2} = imageDims;
allVox.minMax = minMax;



allVox.voxFV = renderCon(allVox.subs,[],[1 .5 0],.3);








