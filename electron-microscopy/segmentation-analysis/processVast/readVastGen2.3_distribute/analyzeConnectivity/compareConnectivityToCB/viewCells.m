
%function[Ic,maxHeight] = showAllObjects(showObj, col,dim,fsize);

dim = 3;
fsize = showObj.ImageSize;
for i = 1:length(showObj.PixelIdxList);
    i
    [y x z] = ind2sub(fsize,showObj.PixelIdxList{i});
    obj(i).subs = [y x z];
end

objNum = length(obj);
cellId = 1:objNum;
col = rand(length(cellId),3);
colMap = hsv(256);
col = colMap(ceil((1:objNum)*256/objNum),:);
col = col(randperm(objNum),:);


minInt = 50;
contrastFactor = 10;


%%

if ~exist('dim','var')
    dim = 1;
end

if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end


%%
%fsize = [1700 1700 1300];
I = zeros(fsize(dims));
newHeight = I;
maxHeight = I;
Ic = cat(3,I,I,I);


uCellId = unique(cellId(cellId>0));

for i = 1:length(uCellId)
    i
    sub = obj(i).subs;
    if ~isempty(sub)
        inds = sub2ind(fsize(dims),sub(:,dims(1)),sub(:,dims(2)));
        uinds = unique(inds);
        
        if length(uinds)>1
            hinds = hist(inds,uinds);
        else
            hinds = length(inds);
        end
        
        I(uinds) = I(uinds) + hinds';
        newHeight(inds) = max(newHeight(inds),sub(:,dim));
    end
    
    I = I * contrastFactor;
    I((I>0)) = I(I>0) + minInt;
    I(I>255) = 255;
    I((I>0) & (I<minInt)) = minInt;
    
    
    useI = newHeight>maxHeight;
    maxHeight(useI) = newHeight(useI);
    
    for c = 1:3
        Itemp = Ic(:,:,c);
        Itemp(useI) = I(useI)*col(i,c);
        Ic(:,:,c) = Itemp;
    end
    
    image(uint8(Ic))
    pause(.001)
    I = I * 0;
    
end

Ic = uint8(Ic);



