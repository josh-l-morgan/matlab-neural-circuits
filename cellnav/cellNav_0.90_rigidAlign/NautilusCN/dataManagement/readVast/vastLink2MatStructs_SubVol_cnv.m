function[] = vastLink2MatStructs_SubVol_cnv

global glob

mipLevel = glob.NA.export.mipLevel;

clear global vdata
if ~exist('vdata','var')
    vasttools
end

global vdata

if ~vdata.state.isconnected
    vdata.vast.connect('127.0.0.1',22081,1000)
end


try
    segNum = vdata.vast.getnumberofsegments
catch
    segNum = 0;
    disp('No segments viewable from VAST.  Make sure VAST Remote Control API Server is enabled.')
end

glob.NA.export.segNum = segNum;

if segNum
    MPN = glob.NA.MPN;
    WPN = glob.NA.WPN;
    
    %{
slash = regexp(MPN,'\');
EPN = [MPN(1:slash(end-1)) 'Export\'];
SPN=uigetdir(EPN)
info=vdata.vast.getinfo();
[selectedlayernr, selectedemlayernr, selectedsegmentlayernr, res] = vdata.vast.getselectedlayernr()
layerNum = vdata.vast.getnroflayers;
[linfo res] = vdata.vast.getlayerinfo(selectedsegmentlayernr);
    %}
    
    startTime = clock;
    
    %     n = vdata.vast.getnroflayers()
    %     [a b] = vdata.vast.getlayerinfo(n)
    %     vdata.vast.getselectedsegmentnr()
    %     vdata.vast.getallsegmentnames
    TPN = [glob.NA.MPN glob.NA.export.exportName '\'];
    
    useSaved = 1;
    
    
    % mpSize = parpool('size');
    % if ~mpSize
    %     'opening matlab pool'
    %     parpool close force
    %     parpool
    % end
    %}
    
    
    %% Read in vast Colors
    %Please connect to VAST with vasttools first!
    
    getRLEtoSubs(TPN,mipLevel);
    
    rleTemp2Subs(TPN);
    
    obI = makeOBI(TPN)
    
    
    %% Down sample vastSubs
    
    
    res = glob.NA.export.volRes;
    vRes = [1 1 1] .* [2^mipLevel 2^mipLevel 1].* res;
    dsRes = repmat(glob.NA.export.dsRes,[1 3]);
    dsDim = dsRes./vRes;
    
    useSaved = 0;
    dsObj = downSampObj(TPN, dsDim);
    
    
    obI.em.res = res; %or [4.6 4 30]?
    obI.em.vRes =vRes;
    obI.em.dsRes =dsRes/1000;
    obI.em.mipLevel = mipLevel;
    
    
    %% Get seed
    
    
    save([TPN 'obI.mat'],'obI')
    
end
   
    close(vdata.fh)
    
end














