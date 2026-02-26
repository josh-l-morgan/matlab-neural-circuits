function[tisDat] = makeTisDat(tis)

global glob
MPN = glob.NA.MPN;
WPN = glob.NA.WPN;
if ~exist('tis','var')
    load([WPN 'tis.mat'])
    saveDat = 1;
else
    saveDat = 0;
end

fvDir = [WPN 'fvLibrary\'];
%load('MPN.mat')
%fvDir = [WPN 'fvLibrary\']

tisDat.created = datestr(clock);
try
    tisDat.depthHist = getHisto(tis,fvDir);
catch
    errStr = 'Did not successfully run getHisto';
    set(glob.NA.export.handles.textOut,'String',errStr)
    for c = 1:10
        set(glob.NA.export.handles.textOut,'backgroundcolor',[1 .5 .5])
        pause(.2)
        set(glob.NA.export.handles.textOut,'backgroundcolor',[1 1 1])
        pause(.2)
    end
    tisDat.getHistoFailed = 1;
    disp('getHisto failed')
end

if saveDat
    save([fvDir 'tisDat.mat'],'tisDat')
end


























