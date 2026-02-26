function[pathLength voxLength] = path2seed(propPath)


pathLength = zeros(length(propPath.tips),1);
bases = pathLength;
owner = zeros(length(propPath.owner),1);
for t = 1:length(propPath.tips)
    t
    sourceTip = propPath.tips(t);
    tip = sourceTip;
    pathLength(t) = 0;
    bases(t) = tip;
    for r = 1:length(propPath.pred)*2
        bases(t) = tip;
        if tip<1, break, end
        owner(tip) = t;
        pathLength(t) = pathLength(t) + propPath.predLength(tip);
        tip =  propPath.pred(tip);
    end
end


voxLength = owner*0;
voxLength(owner>0) = pathLength(owner(owner>0)); %flawed list







