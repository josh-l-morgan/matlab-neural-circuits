


% SPN1 = 'D:\LGNs1\Analysis\movies\subLin125_Region1_SynOFF_scaleBar2\'
% SPN2 = 'D:\LGNs1\Analysis\movies\subLin125_Region1_SynON_scaleBar2\'
% TPN = 'D:\LGNs1\Analysis\movies\subLin125_Region1_SynONOFF_scaleBar2\'

SPN1 = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_A_VGonly\'
SPN2 = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_A_VGmarkers\'
TPN = 'E:\IxQ_KarlsRetinaVG3_2019\Analysis\movies\IxQ\cid4_VG4_A_withBips2\'


SPN1 = 'D:\LGNs1\Analysis\movies\subLin125_Glom1c_SynOFF_scaleBar2\'
SPN2 = 'D:\LGNs1\Analysis\movies\subLin125_Glom1c_SynON_scaleBar2\'
TPN = 'D:\LGNs1\Analysis\movies\subLin125_Glom1c_SynONFF_scaleBar2\'

if ~exist(TPN,'dir'),mkdir(TPN); end
tag = 'region1_'

c = 0;

dir1 = dir([SPN1 '*.png']);
inam1 = {dir1.name};


dir2 = dir([SPN2 '*.png']);
inam2 = {dir2.name};

for i = 1:length(inam1)
   c = c+1; 
    fileName = sprintf('%s_%05.0f.png',tag,c)
    copyfile([SPN1 inam1{i}],[TPN fileName]);
end


for i = 1:length(inam2)
   c = c+1; 
    fileName = sprintf('%s_%05.0f.png',tag,c)
    copyfile([SPN2 inam2{i}],[TPN fileName]);
end























