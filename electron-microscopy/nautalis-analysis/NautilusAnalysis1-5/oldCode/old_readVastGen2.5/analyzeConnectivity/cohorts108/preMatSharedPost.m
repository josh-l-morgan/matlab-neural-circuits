function[sharedMat] = getMinSharedPre(inMat);

showShared = 0;
sharedMat = zeros(size(inMat,1),size(inMat,1),size(inMat,2));
zeroX = [];
[sizeY sizeX] = size(inMat);

for x = 1:sizeX
    for y = 1:sizeY
            sharedMat(y,:,x) = min(inMat(:,x),inMat(y,x));
    end
    sharedMat(sub2ind(size(sharedMat),1:sizeY,1:sizeY,ones(1,sizeY)*x))=0;
end

%%

if showShared
    for i = 1:size(sharedMat,3)
       image(sharedMat(:,:,i)*50)
       pause(.01)
    end
end
