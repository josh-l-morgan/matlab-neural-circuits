function[cids]  = checkAlias(cids,aliases)


for i = 1:length(aliases)
    alias = aliases{i};
    for a = 2:length(alias)
        cids(cids==alias(a)) = alias(1);
    end
end
