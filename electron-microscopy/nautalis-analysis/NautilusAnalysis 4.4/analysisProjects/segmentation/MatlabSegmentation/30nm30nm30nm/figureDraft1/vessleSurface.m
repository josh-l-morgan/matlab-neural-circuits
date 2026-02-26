
SPN = 'D:\LGNs1\PSC_alignments\joshmProcess\S32intermediateSingleList_ds_medZ_filtBlood4_filtCBV3d_colorCBV\'

iNams = dir([SPN '*.tif']);
Iraw = imread([SPN iNams(i).name]);
[ys xs zs] = size(Iraw);
I = zeros(ys xs length(iNams));

for i = 1:length(iNams)
    
    Iraw = imread([SPN iNams(i).name]);
    
    
end



































