function[] = makeTisDat()

load('MPN.mat')

load([WPN 'tis.mat'])

fvDir = [WPN 'fvLibrary\']
tisDat.created = datestr(clock);
tisDat.depthHist = getHisto(tis,fvDir);



save([fvDir 'tisDat.mat'])


























