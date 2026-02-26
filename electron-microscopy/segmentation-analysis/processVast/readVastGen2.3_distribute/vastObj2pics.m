function[] = vastObj2pics(TPN)

load([TPN 'vastOb.mat']);

IPN = [TPN 'sumImages\']
if ~exist(IPN,'dir'),mkdir(IPN); end

I = zeros(vastOb.size(1:2),'uint8');
for i = 1:length(vastOb)
    
    I = I *0;
    ind = sub2ind(vastOb.size(1:2),vastOb.subs{i}(:,1),...
        vastOb.subs{i}(:,2));
    uind = unique(ind);
    hind = hist(ind,uind);
    I(uind) = hind;
    
    nam = sprintf('sum%d_o%04.0f.png',3,vastOb.uniqueIds(i));
    
end

