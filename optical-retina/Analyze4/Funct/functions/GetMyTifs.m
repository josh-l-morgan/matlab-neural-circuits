function[I] = GetMyTifs

TPN = GetMyDir;
Idir=dir(TPN);
In={};
for i = 1: length(Idir)
    name=Idir(i).name;
    LN=length(name);
    if LN>=3
        if length(name)>=4
        if sum(name(LN-3:LN)=='.tif')==4;
            In{length(In)+1}=name;
        end
        end
    end
end

I = imread([TPNi '\' In{1}]);
I = max(I,[],3);
siz = [size(I,1) size(I,2) length(In)];

if siz(3)>1
    I = zeros(siz,'uint8');
for i = 1:siz(3)
    I(:,:,i)=imread([TPNi '\' In{i}]);
end
end