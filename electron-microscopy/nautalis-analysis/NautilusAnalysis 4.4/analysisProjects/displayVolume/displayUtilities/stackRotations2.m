


SPN{1}  = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_B_VGnoMarkers\';
SPN{2} = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_B_VGmarkers\';
SPN{3} = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_B_VGBip\';
SPN{4} = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_B_VGfunc\';

TPN = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\collect_VG4_Bfunc2\';
if ~exist(TPN,'dir'), mkdir(TPN),end

usePlanes{1} = [1 359];
usePlanes{2} = [1 359];
usePlanes{3} = [1 359];
usePlanes{4} = [1 359];


xWindow = [300 1650];
yWindow = [260 1500];

fileType = '.png';
step = 1;

c = 0;
for s = 1:length(SPN)
    dSPN = dir([SPN{s} '*' fileType]);
    nams = {dSPN.name};
    for n = usePlanes{s}(1):step:usePlanes{s}(2)
        if n<=length(nams)
        if exist([SPN{s} nams{n}],'file')
            c = c+1;
            I = imread([SPN{s} nams{n}]);
            I = I(yWindow(1):yWindow(2),xWindow(1):xWindow(2),:);
            fileName = sprintf('%s%d.png',TPN,c);
            imwrite(uint8(I),fileName);
        end
        end
    end
end
