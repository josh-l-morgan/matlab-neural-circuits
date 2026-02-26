function[tisDat] = makeTisDat(tis)

global glob

if ~exist('tis','var')
    MPN = glob.NA.MPN;
    WPN = glob.NA.WPN;
    load([WPN 'tis.mat'])
end

fvDir = glob.fvDir;
%load('MPN.mat')
%fvDir = [WPN 'fvLibrary\']

tisDat.created = datestr(clock);
tisDat.depthHist = getHisto(tis,fvDir);



save([fvDir 'tisDat.mat'],'tisDat')


























