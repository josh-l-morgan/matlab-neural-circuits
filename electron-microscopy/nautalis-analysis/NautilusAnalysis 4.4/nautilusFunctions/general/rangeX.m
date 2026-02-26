function[rangeXX] = rangeX(props,cI)


if ~exist('cI','var')
    cI = 0.95;
end
L = length(props);
sortProps = sort(props);
shiftPick = ceil(L*((1-cI)/2));
shiftPick = max(1,shiftPick);
shiftPick = min(length(sortProps),shiftPick);
rangeXX = [sortProps(shiftPick) sortProps(L-shiftPick+1)];
