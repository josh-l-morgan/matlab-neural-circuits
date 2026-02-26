


SPN{1}  = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_plusBips2b\';
SPN{2} = 'X:\Active\morganLab\ANALYSIS\IxQ\Analysis\movies\IxQ\cid4_plusBips2a\';
SPN{3} = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_plusBips2\';
TPN = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\collect_plusBips2\';

usePlanes{1} = 1:359;
usePlanes{2} = 1:359;
usePlanes{3} = 1:359;

xWindow = [300 1650];
yWindow = [260 1500];

fileType = '.png';
c = 0;
for s = 1:length(SPN)
    dSPN = dir([SPN{s} '*' fileType]);
    nams = {dSPN.name};
    for n = 1:length(nams)
        c = c+1;
        I = imread([SPN{s} nams{n}]);
        I = I(yWindow(1):yWindow(2),xWindow(1):xWindow(2));
        fileName = sprintf('%s%d.png',TPN,c);
        imwrite(I,fileName);
    end
end
    