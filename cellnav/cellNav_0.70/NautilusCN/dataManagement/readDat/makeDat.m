function[] = makeDat

load('MPN.mat')

%% Data for IxQ
dat.googleLinks.datSheet = '14LMWg7-75tov0vBcVP9OaYf8_6EkEFOpE3hGJm5YVcY';
dat.googleLinks.datGID = '1560450437';
dat.googleLinks.aliasGID = '70625524';

dat = parseGoogleDat(dat.googleLinks)

save([MPN 'dat.mat'],'dat')