SPN = GetMyDir
%SPN = 'F:\gxM\gxM_2p2p_vect_num1\CF\bigstack\bigstack_60_1600.oif.files\'
TPN = [SPN(1:end-1) '_reverse\'];
if ~exist(TPN,'dir'), mkdir(TPN);end


inams = dir([SPN '*.tif']);
maxI = length(inams)*2;

for i = 1:length(inams)
    i
   newNam = sprintf('withReverse_%04.0f.tif',i);
    I = imread([SPN inams(i).name]);
    imwrite(I,[TPN newNam]);
    newNam2 = sprintf('withReverse_%4.0f.tif',maxI-i+1);
    imwrite(I,[TPN newNam2]);
    
    
end
