%%Read in segments, segment the segments, turn segments into points. 

%%Load data
global glob tis
load('MPN.mat');
load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

%%Find segments
segName = 'pre cid1007'

oldNames = obI.nameProps.oldNames;
isSeg = zeros(length(oldNames),1);
for i = 1:length(oldNames)
    isReg = regexp(oldNames{i},segName);
    if ~isempty(isReg)
        isSeg(i) = isReg(1);
    end
end

obTarg = find(isSeg>0);
subs = cat(1,dsObj(obTarg).subs);

%%Find centroids
cutSubs = min(subs,[],1)-1;
shortSubs = subs - repmat(cutSubs,[size(subs,1) 1]);
maxSubs = max(shortSubs,[],1);
vInd = sub2ind(maxSubs,shortSubs(:,1),shortSubs(:,2),shortSubs(:,3));

vol = zeros(maxSubs,'logical');
vol(vInd) = 1;
vLab = bwlabeln(vol,6);
props = regionprops(vLab,'centroid');
cents = cat(1,props.Centroid);

cents = cents + double(repmat(cutSubs([2 1 3]),[size(cents,1) 1]));

vCent = round(cents * obI.em.dsRes * 1000 ./ repmat(obI.em.res,[size(cents,1) 1]));
vCent = vCent(:,[1 2 3]);




















