function[subs] = scaleNodes(subs,scaleVec);

%%Rescale nx3 sub list accorting to 1x3 scale vector

%{
Reference

Mip8
anchorScale = [-.0184 0.016 0.030];
voxelScale = [anchorScale(1) * 8 * 4 anchorScale(2) * 8 * 4 anchorScale(3)* 4 * 4];

%}

if ~isempty(subs)
    for d = 1:size(subs,2)
        subs(:,d) = subs(:,d) * scaleVec(d);
    end   
end