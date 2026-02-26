function[maxObj objVol] = showMaxProp(moveObj,prop,minMax)


if exist('minMax','var')
    maxSubs  = minMax{2};
else
    maxSubs = max(moveObj,[],1);
end

objVol = zeros(maxSubs,'uint8');
objVol(sub2ind(maxSubs,moveObj(:,1),moveObj(:,2),moveObj(:,3))) = prop;
maxObj = max(objVol,[],3);
image(maxObj), pause(.01)