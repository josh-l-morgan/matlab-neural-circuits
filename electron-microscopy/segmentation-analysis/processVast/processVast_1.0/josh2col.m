function[cmap] = josh2col(typeLab)

cmap = zeros(length(typeLab),3);

for i = 1:length(typeLab)
    if typeLab(i) == 1;
        cmap(i,1) = 1;
       % cmap(i,3) = rand;
    elseif typeLab(i) == 2
        cmap(i,2) = 1;
        cmap(i,3) = rand;
    end
    
end
cmap = cat(1, [ 0 0 0], cmap);