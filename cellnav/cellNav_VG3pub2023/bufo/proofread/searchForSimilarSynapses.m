

tooClose = .2; % distance for synapses too close to one another

%
% for i  = 1:length(sms)
%
%     sm = sms(i).sm;
%     nep = sm.nep;
%     closest = sm.syn2Skel.closest;
%
%
syn = tis.syn;

pos = syn.pos;
numSyn = size(pos,1);
dif1 = pos(:,1)-pos(:,1)';
dif2 = pos(:,2)-pos(:,2)';
dif3 = pos(:,3)-pos(:,3)';
dif = sqrt(dif1.^2 + dif2.^2 + dif3.^2);
difMask = repmat([1:numSyn]',[ 1 numSyn]) > repmat([1:numSyn],[numSyn 1]);
isClose = dif<=tooClose;
isId = (syn.pre == syn.pre') & (syn.post == syn.post');% & (syn.obID == syn.obID');
[sameY sameX] = find(difMask & isClose & isId);
isSame = [sameY sameX];
samePos = syn.synPosRaw(sameY,[2 1 3]);
samePos1 = syn.synPosRaw(sameY,[2 1 3]);
samePos2 = syn.synPosRaw(sameX,[2 1 3]);
samePos = cat(3,samePos1,samePos2);
samePos = permute(samePos,[3 2 1])

tis.obI.colStruc.names{syn.obID(sameY(1))}

clear sp
for n = 1:length(sameY)

    obId1 = syn.obID(sameY(n));
    obId2 = syn.obID(sameX(n));
    sp(n).name1 = tis.obI.colStruc.names{obId1};
    sp(n).name2 = tis.obI.colStruc.names{obId2};
    sp(n).source1 = tis.obI.fuse.exportDir{tis.obI.fuse.obSource(obId1)};
    sp(n).source2 = tis.obI.fuse.exportDir{tis.obI.fuse.obSource(obId2)};
    sp(n).obId1 = obId1;
    sp(n).obId2 = obId2;
    sp(n).pos1 = samePos1(n,:);
    sp(n).pos2 = samePos2(n,:);
end


for n = 1:length(sp)

    sp(n)
    pause
end




