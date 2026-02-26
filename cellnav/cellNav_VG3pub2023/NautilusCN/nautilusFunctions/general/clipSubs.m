function[inSubs outSubs] = clipSubs(subs,clip)

if ~isempty(clip);
    subsG = subs >= repmat(clip(1,:),[size(subs,1) 1]);
    subsL = subs <= repmat(clip(2,:),[size(subs,1) 1]);
    subsOK = sum(subsG & subsL,2)==3;
    inSubs = subs(subsOK,:);
    outSubs = subs(~subsOK,:);
else
    inSubs = subs;
    outSubs = [];
end