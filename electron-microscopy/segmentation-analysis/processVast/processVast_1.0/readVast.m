function[S] = readVast(SPN)

%% find images
dSPN = dir(SPN);
iNam = {};
for i =1:length(dSPN)
    if sum(regexp(dSPN(i).name,'.tif')) |...
            sum(regexp(dSPN(i).name,'.png'));
        iNam{length(iNam)+1,1} = dSPN(i).name;
    end
end

%% Read Data
try imageInfo = imfinfo([SPN iNam{1}]);
    ys = imageInfo.Height;
    xs = imageInfo.Width;
catch err
    testI = imread([SPN iNam{1}]);
    [ys xs] = size(testI);
end

zs = length(iNam);

S = zeros(ys,xs,zs);
for i = 1:length(iNam)
    I = imread([SPN iNam{i}]);
    if size(I,3) == 1;
        S(:,:,i) = I; 
    else
            dI = double(I);
     indI = dI(:,:,1)* 256^2 + dI(:,:,2) * 256 + dI(:,:,3);
     S(:,:,i) = indI; 
    end
   
end
