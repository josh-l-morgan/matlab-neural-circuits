function[vSubs] = smSub2um(sm)




vSubs = sm.arbor.subs .* repmat(sm.arbor.voxSize,[size(sm.arbor.subs,1) 1]);
vSubs = vSubs + repmat(sm.arbor.offset .* sm.arbor.voxSize,...
    [size(sm.arbor.subs,1) 1]);
