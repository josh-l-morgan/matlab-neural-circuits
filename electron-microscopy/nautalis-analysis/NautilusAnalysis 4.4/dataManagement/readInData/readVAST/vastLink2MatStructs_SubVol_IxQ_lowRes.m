

mipLevel = 3;


if ~exist('vdata','var')
    vasttools
end

global vdata

if ~vdata.state.isconnected
    vdata.vast.connect('127.0.0.1',22081,1000)
end

load('MPN.mat')




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

TPN = GetMyDir;

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

vdata.vast.disconnect
close(vdata.fh)


%% Down sample vastSubs


res = [20 20 40]; 
vRes = [1 1 1] .* [2^mipLevel 2^mipLevel 1].* res;
dsRes = [100 100 100];
dsDim = dsRes./vRes;

useSaved = 0;
dsObj = downSampObj(TPN, dsDim);


obI.em.res = res; %or [4.6 4 30]?
obI.em.vRes =vRes;
obI.em.dsRes =dsRes/1000;
obI.em.mipLevel = mipLevel;


%% Get seed


save([TPN 'obI.mat'],'obI')
















