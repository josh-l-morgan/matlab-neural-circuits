function[maxObj objVol] = showMaxProp(moveObj,prop)

maxSubs = max(moveObj,[],1);
objVol = zeros(maxSubs,'uint8');
objVol(sub2ind(maxSubs,moveObj(:,1),moveObj(:,2),moveObj(:,3))) = prop;
maxObj = max(objVol,[],3);
image(maxObj), pause(.01)