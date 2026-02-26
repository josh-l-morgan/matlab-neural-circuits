

SPN = 'J:\Chas P11 ixD-LGN\W002_retakes_aligned_tif\'
TPN = [SPN(1:end-1) '_invert\'];
if ~exist(TPN,'dir'),mkdir(TPN), end

clear Idir
Idir = dir([SPN '*.tif']);
Inam = {Idir.name};

for i = 1:length(Inam)
    i
    if ~exist([TPN Inam{i}],'file')
        I = imread([SPN Inam{i}]);
        I = 255-I;
        imwrite(I,[TPN Inam{i}]);
    end
end









