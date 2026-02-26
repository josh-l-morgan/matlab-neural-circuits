function[pred] = edges2pred(edges,seed)
%%generate list of upstream nodes relative to a seed point using an edge
%%list

oldN = seed; % choose starting point

fresh = ones(size(edges,1),1); % define edges as unused
uN = unique(edges(:));
pred = uN * 0;

while ~isempty(oldN)
    newN = [];
    hits = [];
    for n = 1:length(oldN)
        is1 = (edges(:,1)== oldN(n)) & (fresh);
        is2 = (edges(:,2) == oldN(n)) & (fresh);
        fresh(is1) = 0;
        fresh(is2) = 0;
        hit = cat(1,edges(is1,2),edges(is2,1));
        hits = cat(1,hits,hit);
        newN = cat(1,newN,hit);
        pred(hit) = oldN(n);        
    end
    oldN = newN;
end

