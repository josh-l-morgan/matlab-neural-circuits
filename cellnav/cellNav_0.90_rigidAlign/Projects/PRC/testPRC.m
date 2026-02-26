
dilateRad  = 2; %radius to dilate PRC



global tis

load('MPN.mat')
load([MPN 'dsObj.mat'])


obIsPrc = find(tis.obI.nameProps.tag.prc);
prcNum = length(obIsPrc);
obIsGlia = find(tis.obI.nameProps.tag.glia);
ob2cid = tis.obI.nameProps.cellNum;

prcCids = ob2cid(obIsPrc)';
typeID = tis.cells.type.typeID;

%% Get Segmentation dimensions
minSeg = [inf inf inf];
maxSeg = [0 0 0];
for i = 1:length(dsObj)
    subs = dsObj(i).subs;
    cMin = min(subs,[],1);
    cMax = max(subs,[],1);
    minSeg = min([minSeg;cMin],[],1);
    maxSeg = max([maxSeg;cMax],[],1);
end

%% Make Segmentation Volume
vBuf = 10;
vSize = maxSeg-minSeg + vBuf*2;
numOb = length(dsObj);

segV = zeros(vSize);
for i = 1:numOb
    subs = dsObj(i).subs;
    subs = subs - repmat(minSeg,[size(subs,1) 1]) +vBuf;
    subInd = sub2ind(vSize,subs(:,1),subs(:,2),subs(:,3));
    segV(subInd) = i;
end

cmap = hsv(numOb);
cmap = [0 0 0;cmap(randperm(size(cmap,1)),:)];
clf
colormap(cmap)
for i = 1:size(segV,3)
    image(segV(:,:,i)+1);
    drawnow
end

sumSegV = sum(segV>0,3);



%% Dilate PRC
se = strel('sphere',dilateRad);
showDilation = 1;
clear dVals
for i = 1:prcNum

    [y x z] = find(segV==obIsPrc(i));
!!!!!!!!!!!!!! Make cut out for dilation
    prcV1 = (segV==obIsPrc(i));
    dPrcV = imdilate(prcV1,se);
    cutVals = find((dPrcV>0) &(prcV1==0));

    dVals{i} = segV(cutVals);
    if showDilation
        if sum(dVals{i})
            cutV = segV* 0;
            cutV(cutVals)  = segV(cutVals);
            sumCutV = sum(cutV>0,3);
            sumDPrcV = sum(dPrcV,3);

            colSum = zeros(size(sumPrcV,1),size(prcV,2),3);
            colSum(:,:,1) = sumPrcV;
            colSum(:,:,2) = sumCutV;
            colSum(:,:,3) = sumSegV;

            maxCutV = max(cutV,[],3);
            [y x] = find(maxCutV>0);
            minCut = min([y x],[],1);
            maxCut = max([y x],[],1);
            colormap(cmap)

            minCut = max(minCut-10,[1 1]);
            maxCut = min(maxCut+10,size(maxCutV));

            subplot(1,2,1)

            image(uint8(colSum(minCut(1):maxCut(1),minCut(2):maxCut(2),:)*50))
            subplot(1,2,2)
            image(maxCutV(minCut(1):maxCut(1),minCut(2):maxCut(2)))

            drawnow
        end
    end
end


%% Parse dVals
clear n
for i = 1:length(dVals);
    
    v = dVals{i};
    n(i).v = v;
    n(i).cutVoxNum = length(v);
    n(i).emptyVox = sum(v==0);
    cutCidsV  = ob2cid(v(v>0));
    n(i).cutCids = unique(cutCidsV);
    n(i).cutCidsCount = histc(cutCidsV,n(i).cutCids);
    n(i).types = typeID(n(i).cutCids);

end

%% Check cell types
checkTypeNames = {'rgc' 'bpc' 'amc' 'mul'}; % N + 1 will be other
clear checkTypeId
for i = 1:length(checkTypeNames)
    checkTypeId(i) = find(strcmp(tis.cells.type.typeNames,checkTypeNames{i}));
end

for i = 1:length(n)
    types = n(i).types;
    n(i).checkTypeNames = cat(2,checkTypeNames,'other');
    n(i).checkTypeIds = checkTypeId;
    for t = 1:length(checkTypeId)

        isT = find(types==checkTypeId(t));
        n(i).typeCount(t) = length(isT);
        n(i).typeVCount(t) = sum(n(i).cutCidsCount(isT));
        types(isT) = 0;
    end
    isOther = find(isT>0);
    n(i).typeCount(length(checkTypeId)+1) = length(isOther);
    n(i).typeVCount(length(checkTypeId)+1) = sum(n(i).cutCidsCount(isOther));


end


%% Combine cells

typeCellNum = zeros(prcNum,length(checkTypeNames)+1);
typeVoxNum = typeCellNum;
totVox = zeros(prcNum,1);
for i = 1:length(n)
   typeCellNum(i,:) = n(i).typeCount;
   typeVoxNum(i,:) = n(i).typeVCount;
   totVox(i) = n(i).cutVoxNum;
end











