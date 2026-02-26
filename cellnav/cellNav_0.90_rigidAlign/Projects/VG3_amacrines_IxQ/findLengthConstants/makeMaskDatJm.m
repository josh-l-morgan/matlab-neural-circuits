global tis glob
SPN = [glob.datDir 'Analysis\Data\preproc\'];


%%  Ai148_129SVG3_AutoRregMask_122618_1005

iSize = [32 256]

load([SPN 'ptDat.mat']);

ps = ptDat(:,3);
cs = ptDat(:,4);
rs = ptDat(:,5);

maskDat = zeros(iSize(1),iSize(2),length(ps));
for i = 1:length(ps)
    maskDat(rs(i),cs(i),i) = 1;
end

clf
colormap gray(255)
sumMask = sum(maskDat,3);
imshow(sumMask * 255/max(sumMask(:)));

save([SPN 'maskDatPts.mat'],'maskDat')














