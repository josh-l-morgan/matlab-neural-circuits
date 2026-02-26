function[subs] = indToSub(siz,ind)

if length(siz) == 3
    [y x z] = ind2sub(siz,ind);
    subs = [y x z];
else
    
    [y x ] = ind2sub(siz,ind);
    subs = [y x ];
end