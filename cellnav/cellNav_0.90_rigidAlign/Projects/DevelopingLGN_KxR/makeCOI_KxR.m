function[] = makeCOI_LGN()

global tis glob
datFold = [glob.datDir 'Analysis\Data\'];
if ~exist(datFold,'dir'),mkdir(datFold);end
clear COI

COI.targetTCs = [3004 3107 3210 3211 3213 3018 3097 3032 3098 3200 3202 3203 3204 3206 3209 3219];

save([datFold 'COI.mat'],'COI')









