function[cols] = val2cmap(vals,cmap);

if ~exist('cmap','var')
    cmap = jet(100);
end

wall = size(cmap,1);
cVal = round(vals);
cVal(cVal>wall) = wall;
cVal(cVal<1) = 1;
cols = cmap(cVal,:);