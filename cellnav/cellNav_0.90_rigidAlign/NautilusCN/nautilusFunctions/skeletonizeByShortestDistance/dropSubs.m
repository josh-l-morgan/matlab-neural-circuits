function[subs] = dropSubs(subs);

subs = subs-repmat(min(subs,[],1),[size(subs,1) 1]) + 1;
