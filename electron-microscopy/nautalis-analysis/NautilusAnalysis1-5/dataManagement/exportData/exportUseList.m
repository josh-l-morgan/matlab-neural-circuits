

TPN = 'D:\LGNs1\Analysis\exports\';
fileName = 'useList_seed_RGC_TCR.mat';

seedList = [108 201 109 903 907];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);

save([TPN fileName],'useList')