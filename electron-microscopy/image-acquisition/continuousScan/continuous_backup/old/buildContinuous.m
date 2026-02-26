

%TPN = GetMyDir;

picNames = getPics(TPN);

for i = 1 : length(picNames)
    I = imread([TPN picNames{i}]);
    
    image(I)
end