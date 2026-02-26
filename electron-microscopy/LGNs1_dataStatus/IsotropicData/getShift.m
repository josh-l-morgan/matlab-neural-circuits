function[shiftInd] = getShift(sizeI,indList)


[y x] = ind2sub(sizeI, indList);

yshift = [-1 0 1 -1  1 -1 0 1];
xshift = [-1 -1 -1 0  0 1 1 1];

yshifted = repmat(y,[1 8]) + repmat(yshift,[length(y) 1]);
xshifted = repmat(x,[1 8]) + repmat(xshift,[length(x) 1]);

ywall = (yshifted<1) | (yshifted>sizeI(1));
xwall = (xshifted<1) | (xshifted>sizeI(2));

yshifted(ywall) = y(1);
xshifted(xwall) = x(1);

shiftInd = sub2ind(sizeI, yshifted,xshifted);