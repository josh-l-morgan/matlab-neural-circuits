function[croppedI] = cropI(I)


[y x z] = ind2sub(size(I),find(I>0));
croppedI = I(min(y):max(y), min(x):max(x),1:end);