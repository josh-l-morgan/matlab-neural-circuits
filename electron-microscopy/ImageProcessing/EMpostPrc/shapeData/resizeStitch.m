

TPN = GetMyDir;


dTPN = dir(TPN);
dTPN = dTPN(3:end);

iNam = [];
for i = 1:length(dTPN)
    nam = dTPN(i).name;
    if length(nam)>4
        if strcmp(nam(end-3:end), '.tif' )
            iNam{length(iNam)+1} = nam;
        end
    end
end

newTPN = [TPN(1:end-1) '_resized\'];
if ~exist(newTPN)
    mkdir(newTPN)
end

%% Get sizes
for i = 1:length(iNam)
   Iinfo = imfinfo([TPN iNam{i}]);
   iWidth(i) = Iinfo.Width;
   iHeight(i) = Iinfo.Height;
end

maxW = max(iWidth);
maxH = max(iHeight);

%% Resize images
newI = zeros(maxH,maxW,'uint8');
for i = 1:length(iNam)
    newI = newI * 0;
    I = imread([TPN iNam{i}]);
    [ys xs] = size(I);
    shiftY = fix((maxH - ys)/2);
    shiftX = fix((maxW - xs)/2);
    newI(shiftY + 1:shiftY + size(I,1),...
        shiftX + 1:shiftX + size(I,2)) = I;
    imwrite(newI,[newTPN iNam{i}],'compression','none');
end







