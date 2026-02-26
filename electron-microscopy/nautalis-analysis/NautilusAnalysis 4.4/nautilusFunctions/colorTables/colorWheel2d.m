function[colMat] = colorWheel2d(colNum,satVal);

%%
if ~exist('colNum','var')
    colNum = 100;
end
if ~exist('satVal','var')
    satVal = colNum/2;
end


field = ones(colNum);
fInd = find(field);
[y x] = find(field);

ydif = y-mean(y);
xdif = x-mean(x);

dists = sqrt(ydif.^2+xdif.^2);
rads = atan2(xdif,ydif);
rads = rads - min(rads);
rads = round(rads * 99/max(rads))+1;



idxField = field;
idxField(fInd) = rads;
image(idxField);
colormap hsv(100)



hmap = hsv(100);
satmap = hmap(rads,:);
unsatmap = satmap;
for i = 1:100;
    isCol = find(rads == i);
    getDists = dists(isCol).^1.5;
    %maxDist = max(getDists);
    maxDist = (satVal).^1.5;
    getDists = min(getDists/maxDist,1);
    getCol = satmap(isCol,:);
    getGrey = getCol*0+.5;
    newCol = getCol .* repmat(getDists,[1 3]) + ...
        getGrey .* repmat(1-getDists,[1 3]) ;
    unsatmap(isCol,:) = newCol;
end

for c = 1 : 3
   
    tempI = field;
    tempI(:) = unsatmap(:,c);
    colMat(:,:,c) = tempI;
    
end
image(uint8(colMat*256))

